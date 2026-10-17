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

"""Structures a person drew or pasted into a question, named by placeholders.

A question can carry a structure that has no usable name -- a drawn product, a scaffold,
a pattern with an atom list, a transformation. The person's interface names each one
with a placeholder, ``@A``, written into the question, and sends the structure beside
it. The model sees each placeholder's name and format and never its text: it writes the
name where a structure goes, and :func:`substitute` puts the text there before the
query is validated.

Substituting rather than resolving keeps the grammar unaware of placeholders. The query
that runs names its structures outright, so it can be run again as it stands; the
grammar's own checks apply to the text; and a placeholder never reaches a name
resolver, where ``@A`` would be looked up as a compound.
"""

import enum
import re
from collections.abc import Mapping
from typing import Any

from pydantic import BaseModel, ConfigDict, ValidationError, model_validator
from rdkit import Chem, rdBase

from ord_schema.search import query

# Each structure predicate is a pass over the corpus's structures, so a question's cost
# grows with the placeholders it carries.
MOST_PLACEHOLDERS = 10
MOST_CHARACTERS = 2000
PLACEHOLDERS_HEADING = "Placeholders in this question:"

_NAME = re.compile(r"@[A-Z]+")
# A marker is a name standing alone in prose. The lookbehind keeps a SMILES chirality
# mark, [C@@H] or C[C@H](O), and an address, user@Example.com, from reading as one.
_MARKER = re.compile(r"(?<![\w@])@[A-Z]+(?!\w)")
# The keys a structure predicate names its structure under.
_SLOTS = ("compound", "smiles", "smarts")


class StructureFormat(enum.StrEnum):
    """What a placeholder's text is."""

    COMPOUND_SMILES = "compound_smiles"
    COMPOUND_SMARTS = "compound_smarts"
    REACTION_SMARTS = "reaction_smarts"

    @property
    def label(self) -> str:
        """Returns the format as a person reads it, such as "compound SMILES"."""
        return _LABELS[self]


_LABELS = {
    StructureFormat.COMPOUND_SMILES: "compound SMILES",
    StructureFormat.COMPOUND_SMARTS: "compound SMARTS",
    StructureFormat.REACTION_SMARTS: "reaction SMARTS",
}

# The field each structure predicate takes its structure in, and the formats it takes.
# A compound SMILES read as a substructure is read as SMARTS, as the search form reads
# a drawn substructure.
_FIELDS = {
    "substructure": "smarts",
    "similarity": "smiles",
    "same_compound": "smiles",
    "same_parent": "smiles",
    "reaction_smarts": "smarts",
}
_FORMATS = {
    "substructure": (StructureFormat.COMPOUND_SMILES, StructureFormat.COMPOUND_SMARTS),
    "similarity": (StructureFormat.COMPOUND_SMILES,),
    "same_compound": (StructureFormat.COMPOUND_SMILES,),
    "same_parent": (StructureFormat.COMPOUND_SMILES,),
    "reaction_smarts": (StructureFormat.REACTION_SMARTS,),
}


class PlaceholderError(ValueError):
    """A placeholder, or a question's use of one, does not fit."""


def _grammar_check(model: type[BaseModel], payload: dict[str, Any]) -> None:
    """Raises the grammar's own reason when a structure would not validate.

    Args:
        model: The predicate that would hold the structure.
        payload: The predicate, holding it.

    Raises:
        ValueError: With the grammar's messages, unwrapped from pydantic's report.
    """
    try:
        model.model_validate(payload)
    except ValidationError as error:
        raise ValueError(
            "; ".join(
                str(detail["msg"]).removeprefix("Value error, ")
                for detail in error.errors()
            )
        ) from None


def _has_wildcard(molecule: Chem.Mol) -> bool:
    """Returns whether a molecule holds a dummy atom, ``*``."""
    return any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms())


class Placeholder(BaseModel):
    """One structure a question names by placeholder.

    The text is held to what the grammar accepts for its format, so a placeholder that
    exists is one a query can use.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    format: StructureFormat
    text: str

    @model_validator(mode="after")
    def _check(self) -> "Placeholder":
        if len(self.text) > MOST_CHARACTERS:
            raise ValueError(
                f"a placeholder's text is at most {MOST_CHARACTERS} characters"
            )
        match self.format:
            case StructureFormat.COMPOUND_SMILES:
                _grammar_check(
                    query.SameCompound,
                    {"op": "same_compound", "path": "smiles", "smiles": self.text},
                )
                # A dummy atom matches nothing as a compound's identity, but any atom
                # as a pattern, so a structure holding one is a pattern.
                if _has_wildcard(Chem.MolFromSmiles(self.text)):
                    raise ValueError(
                        "a compound SMILES holds no wildcard atom; "
                        f"{self.text!r} is a compound SMARTS"
                    )
            case StructureFormat.COMPOUND_SMARTS:
                _grammar_check(
                    query.Substructure,
                    {"op": "substructure", "path": "smiles", "smarts": self.text},
                )
            case StructureFormat.REACTION_SMARTS:
                _grammar_check(
                    query.ReactionSmarts,
                    {"op": "reaction_smarts", "smarts": self.text},
                )
        return self


def is_name(value: object) -> bool:
    """Returns whether a value is a placeholder name, such as ``@A``."""
    return isinstance(value, str) and _NAME.fullmatch(value) is not None


def markers(question: str) -> set[str]:
    """Returns the placeholder names a question writes.

    Args:
        question: The question, in English.

    Returns:
        Each name standing alone in the text, once.
    """
    return set(_MARKER.findall(question))


def classify(text: str) -> Placeholder:
    """Returns a structure as a placeholder, in the format its text is written in.

    A reaction is a reaction SMARTS. A molecule RDKit reads as SMILES is a compound
    SMILES unless it holds a wildcard, and anything else RDKit reads as SMARTS is a
    compound SMARTS.

    Args:
        text: The structure.

    Returns:
        The placeholder.

    Raises:
        PlaceholderError: If the text is none of the three.
        ValueError: If it reads as a reaction but the grammar refuses it, such as one
            with too many templates on a side.
    """
    if ">" in text:
        return Placeholder(format=StructureFormat.REACTION_SMARTS, text=text)
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(text)
        pattern = Chem.MolFromSmarts(text)
    if molecule is not None and molecule.GetNumAtoms() and not _has_wildcard(molecule):
        return Placeholder(format=StructureFormat.COMPOUND_SMILES, text=text)
    if pattern is not None and pattern.GetNumAtoms():
        return Placeholder(format=StructureFormat.COMPOUND_SMARTS, text=text)
    raise PlaceholderError(f"{text!r} is not a SMILES, SMARTS, or reaction SMARTS")


def check(question: str, structures: Mapping[str, Placeholder]) -> None:
    """Raises unless a question and its placeholders fit each other.

    Args:
        question: The question, in English.
        structures: Each placeholder's name, mapped to its structure.

    Raises:
        PlaceholderError: If a name is not a placeholder name, there are more than
            ``MOST_PLACEHOLDERS``, the question writes a name with no structure, or a
            structure's name is not in the question.
    """
    for name in structures:
        if not is_name(name):
            raise PlaceholderError(
                f"{name!r} is not a placeholder name: write @ and capital letters, "
                "such as @A"
            )
    if len(structures) > MOST_PLACEHOLDERS:
        raise PlaceholderError(
            f"a question holds at most {MOST_PLACEHOLDERS} placeholders, and this one "
            f"has {len(structures)}"
        )
    written = markers(question)
    unbound = _ordered(written - set(structures))
    if unbound:
        raise PlaceholderError(f"{unbound[0]} is in the question but has no structure")
    unused = _ordered(set(structures) - written)
    if unused:
        raise PlaceholderError(
            f"{unused[0]} has a structure but is not in the question"
        )


def _ordered(names: set[str]) -> list[str]:
    """Returns names in placeholder order: ``@Z`` before ``@AA``."""
    return sorted(names, key=lambda name: (len(name), name))


def describe(structures: Mapping[str, Placeholder]) -> str:
    """Returns what the model is told about a question's placeholders.

    Args:
        structures: Each placeholder's name, mapped to its structure.

    Returns:
        A heading and one line per placeholder naming its format, never its text; empty
        for a question with none.
    """
    if not structures:
        return ""
    lines = [
        f"{name} is a {structures[name].format.label}."
        for name in _ordered(set(structures))
    ]
    return "\n".join([PLACEHOLDERS_HEADING, *lines])


def dump(structures: Mapping[str, Placeholder]) -> dict[str, dict[str, str]] | None:
    """Returns a question's placeholders as JSON-ready values, for a record of it.

    Args:
        structures: Each placeholder's name, mapped to its structure.

    Returns:
        Each name mapped to its ``format`` and ``text``, or None for a question with no
        placeholders.
    """
    if not structures:
        return None
    return {name: value.model_dump(mode="json") for name, value in structures.items()}


def substitute(value: Any, structures: Mapping[str, Placeholder]) -> Any:
    """Returns a query with each placeholder replaced by its structure.

    A placeholder is the whole value of a structure predicate's ``compound``,
    ``smiles``, or ``smarts``, whichever the model wrote. Its text lands in the field
    the predicate takes -- ``smarts`` for ``substructure`` and ``reaction_smarts``,
    ``smiles`` for the rest -- and the key the model used is removed.

    Args:
        value: A query as the model wrote it, already decoded.
        structures: Each placeholder's name, mapped to its structure.

    Returns:
        A new query. One that names no placeholder comes back equal to ``value``.

    Raises:
        PlaceholderError: If a placeholder has no structure, its format is not one the
            predicate takes, or a name appears anywhere but as a structure slot's whole
            value.
    """
    return _substitute(value, structures, context=None)


def _substitute(
    value: Any, structures: Mapping[str, Placeholder], *, context: str | None
) -> Any:
    """Returns ``value`` substituted, naming ``context`` in any error.

    Args:
        value: Part of a query.
        structures: Each placeholder's name, mapped to its structure.
        context: The nearest enclosing ``op``, or the key above ``value`` where there is
            none, for an error to say where a misplaced name was.

    Returns:
        The part, substituted.

    Raises:
        PlaceholderError: As ``substitute`` does.
    """
    if isinstance(value, list):
        return [_substitute(item, structures, context=context) for item in value]
    if isinstance(value, str):
        written = markers(value)
        if written:
            raise PlaceholderError(
                f"{_ordered(written)[0]} appears in {context}, but a placeholder is "
                "only the whole value of compound, smiles, or smarts in substructure, "
                "similarity, same_compound, same_parent, or reaction_smarts"
            )
        return value
    if not isinstance(value, dict):
        return value
    op = value.get("op")
    result: dict[str, Any] = {}
    for key, item in value.items():
        if op in _FIELDS and key in _SLOTS and is_name(item):
            field, text = _structure(op, item, structures)
            if field in result or (field != key and field in value):
                raise PlaceholderError(f"{op} names {item} beside another structure")
            result[field] = text
        else:
            nested = op if isinstance(op, str) else (context or key)
            result[key] = _substitute(item, structures, context=nested)
    return result


def _structure(
    op: str, name: str, structures: Mapping[str, Placeholder]
) -> tuple[str, str]:
    """Returns the field a predicate takes a placeholder's structure in, and the text.

    Args:
        op: The structure predicate.
        name: The placeholder it names.
        structures: Each placeholder's name, mapped to its structure.

    Returns:
        The field and the text.

    Raises:
        PlaceholderError: If the placeholder has no structure, or its format is not one
            the predicate takes.
    """
    if name not in structures:
        listed = ", ".join(_ordered(set(structures))) or "none"
        raise PlaceholderError(
            f"{name} is not a placeholder in this question; the placeholders are "
            f"{listed}"
        )
    placeholder = structures[name]
    taken = _FORMATS[op]
    if placeholder.format not in taken:
        labels = " or a ".join(fmt.label for fmt in taken)
        raise PlaceholderError(
            f"{name} is a {placeholder.format.label}, which {op} does not take; {op} "
            f"takes a {labels}"
        )
    return _FIELDS[op], placeholder.text
