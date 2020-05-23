# Copyright 2020 The Open Reaction Database Authors
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

"""Tests for ord_schema.proto.reaction_pb2."""

from absl.testing import absltest
try:
    from rdkit import Chem
except ImportError:
    Chem = None

from ord_schema.proto import reaction_pb2


class ReactionPb2Test(absltest.TestCase):

    def test_simple(self):
        reaction = reaction_pb2.Reaction()
        reaction.identifiers.add(value='C(C)Cl.Br>>C(C)Br.Cl',
                                 type='REACTION_SMILES')
        self.assertTrue(reaction.IsInitialized())
        self.assertLen(reaction.identifiers, 1)
        self.assertFalse(reaction.HasField('setup'))
        with self.assertRaisesRegex(ValueError,
                                    'Reaction has no field not_a_field'):
            reaction.HasField('not_a_field')

    @absltest.skipIf(Chem is None, 'no rdkit')
    def test_rdkit_binary_compound_identifier(self):
        mol = Chem.MolFromSmiles('COO')
        identifier = reaction_pb2.CompoundIdentifier(type='RDKIT_BINARY',
                                                     bytes_value=mol.ToBinary())
        self.assertEqual(Chem.MolToSmiles(mol),
                         Chem.MolToSmiles(Chem.Mol(identifier.bytes_value)))


if __name__ == '__main__':
    absltest.main()
