# Copyright 2022 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for ord_schema.orm.reaction_class."""

import pytest

# Skip the whole module unless the optional 'reaction-class' extra is installed;
# importing reaction_class pulls in rxn_insight at module load.
pytest.importorskip("rxn_insight")

from rxnmapper import RXNMapper  # ty: ignore[unresolved-import]

from ord_schema.orm.reaction_class import classify_reaction_smiles

_SUZUKI = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
_DEOXYFLUORINATION = "CC(O)CCc1ccccc1.O=S(=O)(F)c1ccc(Cl)cc1>>CC(F)CCc1ccccc1"


@pytest.fixture(scope="module")
def rxn_mapper():
    """Shared RXNMapper so the transformer model loads once for the module."""
    return RXNMapper()


def test_classify_named_reaction(rxn_mapper):
    reaction_class, reaction_name = classify_reaction_smiles(_SUZUKI, rxn_mapper)
    assert reaction_class == "C-C Coupling"
    assert reaction_name is not None
    assert "Suzuki" in reaction_name


def test_unnamed_reaction_name_is_none(rxn_mapper):
    # Rxn-INSIGHT returns the sentinel name "OtherReaction" for this reaction;
    # classify_reaction_smiles normalizes that to None while keeping the class.
    reaction_class, reaction_name = classify_reaction_smiles(
        _DEOXYFLUORINATION, rxn_mapper
    )
    assert reaction_class is not None
    assert reaction_name is None


def test_unclassifiable_returns_none(rxn_mapper):
    assert classify_reaction_smiles("foo>>bar", rxn_mapper) == (None, None)
