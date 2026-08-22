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

"""Builds the eval cases, taking reaction IDs from a corpus rather than inventing them.

Each case names a reference query -- what a correct translation looks like -- and, where
one exists, the near-miss a model plausibly writes instead. Reactions the reference
returns become ``must_return``. Reactions the near-miss returns *and the reference does
not* become ``must_not_return``, and those are what make a wrong translation fail rather
than merely differ.

The counterexamples are found by asking for them -- ``near_miss AND NOT reference`` --
rather than by differencing two capped result sets, which compares samples instead of
sets and produced two cases that failed correct translations.

Rerun this when the corpus changes; the IDs are only meaningful against the corpus they
were drawn from.
"""

import argparse
import json
import pathlib
from collections.abc import Sequence
from typing import Any

from ord_schema.search import execute, query

PYRIDINE = {"op": "eq", "path": "smiles", "value": {"compound": "pyridine"}}
SOLVENT = {"op": "eq", "path": "reaction_role", "value": {"literal": "SOLVENT"}}
DESIRED = {"op": "eq", "path": "is_desired_product", "value": {"literal": True}}
YIELD_OVER_50 = {
    "op": "and",
    "clauses": [
        {"op": "eq", "path": "type", "value": {"literal": "YIELD"}},
        {"op": "gt", "path": "percentage.value", "value": {"literal": 50}},
    ],
}

CASES: list[dict[str, Any]] = [
    {
        "question": "which reactions use pyridine as a solvent?",
        "why": "two conditions on one element, which a wrong translation splits in two",
        "reference": {
            "op": "exists",
            "path": "inputs.components",
            "where": {"op": "and", "clauses": [SOLVENT, PYRIDINE]},
        },
        # Pyridine somewhere and a solvent somewhere, which need not be the same
        # component: the reaction using pyridine as a reactant in toluene matches.
        "near_miss": {
            "op": "and",
            "clauses": [
                {"op": "exists", "path": "inputs.components", "where": PYRIDINE},
                {"op": "exists", "path": "inputs.components", "where": SOLVENT},
            ],
        },
    },
    {
        "question": "reactions run above 350 K",
        "why": "a scalar comparison needing no quantifier at all",
        "reference": {
            "op": "gt",
            "path": "conditions.temperature.setpoint_kelvin",
            "value": {"literal": 350},
        },
        "near_miss": None,
    },
    {
        "question": "reactions where a desired product has a yield above 50%",
        "why": (
            "correlation: the yield has to belong to the desired product, not to "
            "whichever product happens to carry one"
        ),
        "reference": {
            "op": "exists",
            "path": "outcomes.products",
            "where": {
                "op": "and",
                "clauses": [
                    DESIRED,
                    {"op": "exists", "path": "measurements", "where": YIELD_OVER_50},
                ],
            },
        },
        "near_miss": {
            "op": "and",
            "clauses": [
                {"op": "exists", "path": "outcomes.products", "where": DESIRED},
                {
                    "op": "exists",
                    "path": "outcomes.products.measurements",
                    "where": YIELD_OVER_50,
                },
            ],
        },
    },
    {
        "question": "reactions with no solvent at all",
        "why": "a forall, which the occurrence index cannot answer",
        "reference": {
            "op": "forall",
            "path": "inputs.components",
            "where": {
                "op": "ne",
                "path": "reaction_role",
                "value": {"literal": "SOLVENT"},
            },
        },
        "near_miss": None,
    },
]

# A question the grammar cannot express: a value is a literal or a compound, never
# another column, so no comparison between two columns can be written.
INEXPRESSIBLE = {
    "question": "reactions that ran longer than their workup took",
    "why": (
        "comparing two columns, which the grammar cannot express: a value is a literal "
        "or a compound, never another column"
    ),
    "compiles": False,
    "must_return": [],
    "must_not_return": [],
}


def build(corpus: execute.Corpus, *, examples: int = 3) -> list[dict[str, Any]]:
    """Returns the cases, with reaction IDs drawn from ``corpus``.

    Args:
        corpus: The corpus to draw IDs from.
        examples: How many reactions to require of a correct translation.

    Returns:
        Cases ready to serialize.

    Raises:
        ValueError: If a near-miss returns nothing its reference excludes, which makes
            it no near-miss and would leave the case with no counterexample.
    """
    built = []
    for case in CASES:
        rows = corpus.search(
            query.Query.model_validate({"where": case["reference"], "limit": 200})
        )
        right = sorted(rows.column("reaction_id").to_pylist())
        wrong: list[str] = []
        if case["near_miss"] is not None:
            counterexamples = corpus.search(
                query.Query.model_validate(
                    {
                        "where": {
                            "op": "and",
                            "clauses": [
                                case["near_miss"],
                                {"op": "not", "clause": case["reference"]},
                            ],
                        },
                        "limit": 2,
                    }
                )
            )
            wrong = sorted(counterexamples.column("reaction_id").to_pylist())
            if not wrong:
                raise ValueError(
                    f"{case['question']}: the near-miss returns nothing the reference "
                    f"excludes, so it is not a near-miss"
                )
        built.append(
            {
                "question": case["question"],
                "why": case["why"],
                "must_return": right[:examples],
                "must_not_return": wrong,
            }
        )
    built.append(INEXPRESSIBLE)
    return built


def main(argv: Sequence[str] | None = None) -> None:
    """Writes the cases file.

    Args:
        argv: Command-line arguments; ``sys.argv`` when omitted.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projections", required=True, help="Glob for projections")
    parser.add_argument("--structures", required=True, help="Glob for structures")
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path(__file__).parent / "nl_cases.json",
        help="Where to write the cases",
    )
    parser.add_argument(
        "--require-current",
        action="store_true",
        help="Refuse artifacts not written by this version of the library",
    )
    args = parser.parse_args(argv)
    with execute.Corpus(
        args.projections,
        args.structures,
        require_current=args.require_current,
        pivot_budget_bytes=0,
        # The only name any reference query resolves, kept local so building the cases
        # needs no external service.
        resolver={"pyridine": "c1ccncc1"}.__getitem__,
    ) as corpus:
        cases = build(corpus)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(cases, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
