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

from absl.testing import absltest
from absl.testing import parameterized
from rdkit import Chem

from ord_schema import resolvers
from ord_schema.proto import reaction_pb2


class NameResolversTest(absltest.TestCase):

    def test_resolve_names(self):
        roundtrip_smi = lambda smi: Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        message = reaction_pb2.Reaction()
        message.inputs['test'].components.add().identifiers.add(type='NAME',
                                                                value='aspirin')
        self.assertTrue(resolvers.resolve_names(message))
        resolved_smi = roundtrip_smi(
            message.inputs['test'].components[0].identifiers[1].value)
        self.assertEqual(resolved_smi, roundtrip_smi('CC(=O)Oc1ccccc1C(O)=O'))
        self.assertEqual(
            message.inputs['test'].components[0].identifiers[1].type,
            reaction_pb2.CompoundIdentifier.IdentifierType.SMILES)
        self.assertRegex(
            message.inputs['test'].components[0].identifiers[1].details,
            'NAME resolved')


class InputResolversTest(parameterized.TestCase, absltest.TestCase):

    def test_input_resolve(self):
        string = '10 g of THF'
        reaction_input = resolvers.resolve_input(string)
        self.assertEqual(reaction_input.components[0].mass,
                         reaction_pb2.Mass(value=10, units='GRAM'))
        self.assertEqual(
            reaction_input.components[0].identifiers[0],
            reaction_pb2.CompoundIdentifier(type='NAME', value='THF'))
        self.assertEqual(len(reaction_input.components), 1)

        string = '100 mL of 5.0uM sodium hydroxide  in water'
        reaction_input = resolvers.resolve_input(string)
        self.assertEqual(len(reaction_input.components), 2)
        self.assertEqual(reaction_input.components[0].moles,
                         reaction_pb2.Moles(value=500, units='NANOMOLE'))
        self.assertEqual(
            reaction_input.components[0].identifiers[0],
            reaction_pb2.CompoundIdentifier(type='NAME',
                                            value='sodium hydroxide'))
        self.assertEqual(reaction_input.components[1].volume,
                         reaction_pb2.Volume(value=100, units='MILLILITER'))
        self.assertEqual(reaction_input.components[1].volume_includes_solutes,
                         True)
        self.assertEqual(
            reaction_input.components[1].identifiers[0],
            reaction_pb2.CompoundIdentifier(type='NAME', value='water'))

    @parameterized.named_parameters(
        ('bad amount', '100 g of 5.0uM sodium hydroxide  in water',
         'amount of solution must be a volume'),
        ('missing concentration', '100 L of 5 grapes in water',
         'String did not match template'),
    )
    def test_input_resolve_should_fail(self, string, expected_error):
        with self.assertRaisesRegex((KeyError, ValueError), expected_error):
            resolvers.resolve_input(string)


if __name__ == '__main__':
    absltest.main()
