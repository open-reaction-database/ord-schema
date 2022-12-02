# Copyright 2020 Open Reaction Database Project Authors
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
"""Tests for ord_schema.units."""

import pytest
from rdkit import Chem

from ord_schema import resolvers
from ord_schema.proto import reaction_pb2


class TestNameResolvers:
    def test_resolve_names(self):
        roundtrip_smi = lambda smi: Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        message = reaction_pb2.Reaction()
        message.inputs["test"].components.add().identifiers.add(type="NAME", value="aspirin")
        assert resolvers.resolve_names(message)
        resolved_smi = roundtrip_smi(message.inputs["test"].components[0].identifiers[1].value)
        assert resolved_smi == roundtrip_smi("CC(=O)Oc1ccccc1C(O)=O")
        assert (
            message.inputs["test"].components[0].identifiers[1].type
            == reaction_pb2.CompoundIdentifier.IdentifierType.SMILES
        )
        assert "NAME resolved" in message.inputs["test"].components[0].identifiers[1].details


class TestInputResolvers:
    def test_input_resolve(self):
        roundtrip_smi = lambda smi: Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        string = "10 g of THF"
        reaction_input = resolvers.resolve_input(string)
        assert len(reaction_input.components) == 1
        assert reaction_input.components[0].amount.mass == reaction_pb2.Mass(value=10, units="GRAM")
        assert reaction_input.components[0].identifiers[0] == reaction_pb2.CompoundIdentifier(type="NAME", value="THF")
        assert reaction_pb2.CompoundIdentifier.SMILES == reaction_input.components[0].identifiers[1].type
        assert roundtrip_smi(reaction_input.components[0].identifiers[1].value) == roundtrip_smi("C1COCC1")

        string = "100 mL of 5.0uM sodium hydroxide in water"
        reaction_input = resolvers.resolve_input(string)
        assert len(reaction_input.components) == 2
        assert reaction_input.components[0].amount.moles == reaction_pb2.Moles(value=500, units="NANOMOLE")
        assert reaction_input.components[0].identifiers[0] == reaction_pb2.CompoundIdentifier(
            type="NAME", value="sodium hydroxide"
        )
        assert reaction_pb2.CompoundIdentifier.SMILES == reaction_input.components[0].identifiers[1].type
        assert roundtrip_smi(reaction_input.components[0].identifiers[1].value) == roundtrip_smi("[Na+].[OH-]")
        assert reaction_input.components[1].amount.volume == reaction_pb2.Volume(value=100, units="MILLILITER")
        assert reaction_input.components[1].amount.volume_includes_solutes == True
        assert reaction_input.components[1].identifiers[0] == reaction_pb2.CompoundIdentifier(
            type="NAME", value="water"
        )
        assert reaction_pb2.CompoundIdentifier.SMILES == reaction_input.components[1].identifiers[1].type
        assert roundtrip_smi(reaction_input.components[1].identifiers[1].value) == roundtrip_smi("O")

    @pytest.mark.parametrize(
        "string,expected",
        (
            ("100 g of 5.0uM sodium hydroxide in water", "amount of solution must be a volume"),
            ("100 L of 5 grapes in water", "String did not match template"),
        ),
    )
    def test_input_resolve_should_fail(self, string, expected):
        with pytest.raises((KeyError, ValueError), match=expected):
            resolvers.resolve_input(string)
