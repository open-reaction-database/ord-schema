# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Scoring natural-language translations against expected interpretations.

The eval harness has two halves. Deciding whether a translation is *right* needs only
the cases and the model's output, so it lives here alongside the prompt it grades: you
change a prompt and a case together, and a scoring rule drifting from the prompt is the
bug this placement prevents. Actually *running* the resulting query needs a database, so
that half stays with the server that has one.

Identifiers are compared by structure rather than by string, so a prompt change that
emits an equivalent-but-differently-spelled pattern does not read as a regression.
"""

from importlib import resources

import yaml
from pydantic import BaseModel
from rdkit import Chem

from ord_schema.agent.nl_query import NLQuery

# Numeric filters the model must populate only when the question asks for them; an
# unexpected value here is over-extraction and counts as a miss.
_NUMERIC_FIELDS = (
    "min_yield",
    "max_yield",
    "min_conversion",
    "max_conversion",
    "similarity_threshold",
    "limit",
)


class ComponentExpectation(BaseModel):
    """Expected component constraint.

    Identifiers are compared by structure, not by exact string: SMARTS and SMILES are
    canonicalized with RDKit so equivalent patterns match (e.g. ``cB(O)O`` vs
    ``[c]B(O)O``), while names that RDKit cannot parse fall back to a case- and
    whitespace-insensitive string compare (so "4-aminophenol" still differs from
    "aminophenol").
    """

    identifier: str
    target: str
    mode: str


def _canonical_identifier(identifier: str, mode: str) -> str:
    """Canonicalizes a component identifier for structure-aware comparison.

    Args:
        identifier: The component identifier (a name, SMILES, or SMARTS).
        mode: The match mode; ``SMARTS`` is parsed as SMARTS, everything else as SMILES.

    Returns:
        The RDKit-canonical SMARTS or SMILES, or the lowercased, stripped identifier
        RDKit cannot parse it (e.g. a compound name awaiting resolution).
    """
    if mode == "SMARTS":
        mol = Chem.MolFromSmarts(identifier)
        if mol is not None:
            return Chem.MolToSmarts(mol)
    else:
        mol = Chem.MolFromSmiles(identifier)
        if mol is not None:
            return Chem.MolToSmiles(mol)
    return identifier.strip().lower()


class CaseExpectation(BaseModel):
    """Per-case expectations checked against the model's interpretation."""

    components: list[ComponentExpectation] = []
    min_yield: float | None = None
    max_yield: float | None = None
    min_conversion: float | None = None
    max_conversion: float | None = None
    similarity_threshold: float | None = None
    limit: int | None = None


class EvalCase(BaseModel):
    """A single evaluation example: a question and its expected interpretation."""

    question: str
    expect: CaseExpectation


def load_cases() -> list[EvalCase]:
    """Loads the evaluation cases bundled alongside this module."""
    raw = (resources.files("ord_schema.agent") / "nl_query_eval_cases.yaml").read_text(
        encoding="utf-8"
    )
    return [EvalCase.model_validate(case) for case in yaml.safe_load(raw)]


def check_interpretation(expect: CaseExpectation, interpretation: NLQuery) -> list[str]:
    """Returns a list of mismatch messages between expectation and interpretation.

    An empty list means the interpretation matched all expectations.

    Args:
        expect: The case's expected interpretation.
        interpretation: The structured query the model produced.

    Returns:
        Human-readable mismatch descriptions; empty if the case passed.
    """
    mismatches = []
    remaining = list(interpretation.components)
    for wanted in expect.components:
        for candidate in remaining:
            if (
                candidate.target == wanted.target
                and candidate.mode == wanted.mode
                and _canonical_identifier(candidate.identifier, candidate.mode)
                == _canonical_identifier(wanted.identifier, wanted.mode)
            ):
                remaining.remove(candidate)
                break
        else:
            mismatches.append(
                f"missing component {wanted.identifier!r} "
                f"({wanted.target}/{wanted.mode})"
            )
    mismatches.extend(
        f"unexpected component {extra.identifier!r} ({extra.target}/{extra.mode})"
        for extra in remaining
    )
    for field in _NUMERIC_FIELDS:
        wanted_value = getattr(expect, field)
        actual_value = getattr(interpretation, field)
        if wanted_value is None and actual_value is not None:
            mismatches.append(f"unexpected {field}={actual_value}")
        elif wanted_value is not None and actual_value != wanted_value:
            mismatches.append(f"{field}: expected {wanted_value}, got {actual_value}")
    return mismatches
