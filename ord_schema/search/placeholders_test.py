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

"""Tests for ord_schema.search.placeholders."""

import re

import pytest
from pydantic import ValidationError

from ord_schema.search import placeholders

_COMPOUND_SMILES = placeholders.StructureFormat.COMPOUND_SMILES
_COMPOUND_SMARTS = placeholders.StructureFormat.COMPOUND_SMARTS
_REACTION_SMARTS = placeholders.StructureFormat.REACTION_SMARTS
_SMILES = placeholders.Placeholder(format=_COMPOUND_SMILES, text="CCO")


def test_a_placeholder_is_validated_by_its_format():
    assert placeholders.Placeholder(format=_COMPOUND_SMILES, text="c1ccncc1")
    with pytest.raises(ValidationError, match="SMILES does not parse"):
        placeholders.Placeholder(format=_COMPOUND_SMILES, text="c1ccnc")
    with pytest.raises(ValidationError, match="wildcard"):
        placeholders.Placeholder(format=_COMPOUND_SMILES, text="c1ccccc1*")
    with pytest.raises(ValidationError, match="SMARTS has no atoms"):
        placeholders.Placeholder(format=_COMPOUND_SMARTS, text="")
    with pytest.raises(ValidationError, match="at most 5 are supported"):
        placeholders.Placeholder(format=_REACTION_SMARTS, text="C.C.C.C.C.C>>C")
    with pytest.raises(ValidationError, match="at most 2000 characters"):
        placeholders.Placeholder(format=_COMPOUND_SMILES, text="C" * 2001)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("c1ccccc1Br", _COMPOUND_SMILES),
        ("c1ccccc1*", _COMPOUND_SMARTS),
        ("[c;R1]-[Br,I]", _COMPOUND_SMARTS),
        ("C(=O)Cl.N>>C(=O)N", _REACTION_SMARTS),
    ],
)
def test_classify(text, expected):
    assert placeholders.classify(text) == placeholders.Placeholder(
        format=expected, text=text
    )


@pytest.mark.parametrize("text", ["2-pyridone", "Pd(OAc)2", "CO2", "70%"])
def test_classify_refuses_what_is_not_a_structure(text):
    with pytest.raises(
        placeholders.PlaceholderError,
        match="not a SMILES, SMARTS, or reaction SMARTS",
    ):
        placeholders.classify(text)


def test_markers():
    assert placeholders.markers("make (@A) from @B, and @AA's yield") == {
        "@A",
        "@B",
        "@AA",
    }
    assert placeholders.markers("[C@@H] and C[C@H](O) from user@Example.com") == set()


def test_check():
    placeholders.check("reactions using @A", {"@A": _SMILES})
    placeholders.check("reactions using ethanol", {})
    with pytest.raises(
        placeholders.PlaceholderError,
        match="@B is in the question but has no structure",
    ):
        placeholders.check("reactions using @A and @B", {"@A": _SMILES})
    with pytest.raises(
        placeholders.PlaceholderError,
        match="@B has a structure but is not in the question",
    ):
        placeholders.check("reactions using @A", {"@A": _SMILES, "@B": _SMILES})
    with pytest.raises(
        placeholders.PlaceholderError, match="'@a' is not a placeholder name"
    ):
        placeholders.check("reactions using @a", {"@a": _SMILES})
    many: dict[str, placeholders.Placeholder] = {
        f"@{letter}": _SMILES for letter in "ABCDEFGHIJK"
    }
    with pytest.raises(placeholders.PlaceholderError, match="at most 10 placeholders"):
        placeholders.check(" ".join(many), many)


def test_describe():
    assert placeholders.describe({}) == ""
    assert placeholders.describe(
        {
            "@AA": placeholders.Placeholder(format=_REACTION_SMARTS, text="C>>C"),
            "@B": placeholders.Placeholder(format=_COMPOUND_SMARTS, text="[Br,I]"),
        }
    ) == (
        "Placeholders in this question:\n"
        "@B is a compound SMARTS.\n"
        "@AA is a reaction SMARTS."
    )


def test_substitute_moves_text_into_the_field_each_predicate_takes():
    structures = {
        "@A": placeholders.Placeholder(format=_COMPOUND_SMILES, text="c1ccncc1"),
        "@B": placeholders.Placeholder(format=_COMPOUND_SMARTS, text="c[Br,I]"),
        "@E": placeholders.Placeholder(
            format=_REACTION_SMARTS, text="cB(O)O.c[Br,I]>>c-c"
        ),
    }
    query = {
        "where": {
            "op": "and",
            "clauses": [
                {"op": "same_compound", "path": "smiles", "compound": "@A"},
                {"op": "substructure", "path": "smiles", "smiles": "@A"},
                {"op": "substructure", "path": "smiles", "compound": "@B"},
                {"op": "reaction_smarts", "smarts": "@E"},
            ],
        }
    }
    assert placeholders.substitute(query, structures)["where"]["clauses"] == [
        {"op": "same_compound", "path": "smiles", "smiles": "c1ccncc1"},
        {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
        {"op": "substructure", "path": "smiles", "smarts": "c[Br,I]"},
        {"op": "reaction_smarts", "smarts": "cB(O)O.c[Br,I]>>c-c"},
    ]


def test_substitute_leaves_a_query_without_placeholders_alone():
    query = {
        "where": {"op": "same_compound", "path": "smiles", "compound": "pyridine"},
        "limit": 3,
    }
    assert placeholders.substitute(query, {}) == query


@pytest.mark.parametrize(
    ("clause", "message"),
    [
        (
            {
                "op": "similarity",
                "path": "smiles",
                "compound": "@B",
                "threshold": 0.6,
            },
            "@B is a compound SMARTS, which similarity does not take; similarity "
            "takes a compound SMILES",
        ),
        (
            {"op": "same_compound", "path": "smiles", "compound": "@D"},
            "@D is not a placeholder in this question",
        ),
        (
            {"op": "reaction_smarts", "smarts": "@A>>@B"},
            "@A appears in reaction_smarts",
        ),
        (
            {"op": "eq", "path": "smiles", "value": {"literal": "@A"}},
            "@A appears in eq",
        ),
    ],
)
def test_substitute_refuses(clause, message):
    structures = {
        "@A": placeholders.Placeholder(format=_COMPOUND_SMILES, text="CCO"),
        "@B": placeholders.Placeholder(format=_COMPOUND_SMARTS, text="c[Br,I]"),
    }
    with pytest.raises(placeholders.PlaceholderError, match=re.escape(message)):
        placeholders.substitute({"where": clause}, structures)
