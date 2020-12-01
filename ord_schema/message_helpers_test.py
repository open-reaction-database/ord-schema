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
"""Tests for ord_schema.message_helpers."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import parameterized
from google.protobuf import json_format
from google.protobuf import text_format
from rdkit import Chem

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.proto import test_pb2

_BENZENE_MOLBLOCK = """241
  -OEChem-07232015262D

 12 12  0     0  0  0  0  0  0999 V2000
    2.8660    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7320    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7320   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660    1.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4631    0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690    0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4631   -0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690   -0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -1.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  1  7  1  0  0  0  0
  2  4  1  0  0  0  0
  2  8  1  0  0  0  0
  3  5  2  0  0  0  0
  3  9  1  0  0  0  0
  4  6  2  0  0  0  0
  4 10  1  0  0  0  0
  5  6  1  0  0  0  0
  5 11  1  0  0  0  0
  6 12  1  0  0  0  0
M  END"""


class MessageHelpersTest(parameterized.TestCase, absltest.TestCase):

    @parameterized.parameters(
        ('ord-1234567890', 'data/12/ord-1234567890'),
        ('test/ord-foo.pbtxt', 'data/fo/ord-foo.pbtxt'),
        ('ord_dataset-f00.pbtxt', 'data/f0/ord_dataset-f00.pbtxt'),
        ('ord_data-123456foo7.jpg', 'data/12/ord_data-123456foo7.jpg'))
    def test_id_filename(self, filename, expected):
        self.assertEqual(message_helpers.id_filename(filename), expected)

    @parameterized.named_parameters(
        ('SMILES', 'c1ccccc1', 'SMILES', 'c1ccccc1'),
        ('INCHI', 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H', 'INCHI', 'c1ccccc1'),
        ('MOLBLOCK', _BENZENE_MOLBLOCK, 'MOLBLOCK', 'c1ccccc1'))
    def test_smiles_from_compound(self, value, identifier_type, expected):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value=value, type=identifier_type)
        self.assertEqual(message_helpers.smiles_from_compound(compound),
                         expected)

    @parameterized.named_parameters(
        ('SMILES', 'c1ccccc1', 'SMILES', 'c1ccccc1'),
        ('INCHI', 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H', 'INCHI', 'c1ccccc1'),
        ('MOLBLOCK', _BENZENE_MOLBLOCK, 'MOLBLOCK', 'c1ccccc1'))
    def test_mol_from_compound(self, value, identifier_type, expected):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value=value, type=identifier_type)
        mol = message_helpers.mol_from_compound(compound)
        self.assertEqual(Chem.MolToSmiles(mol), expected)
        mol, identifier = message_helpers.mol_from_compound(
            compound, return_identifier=True)
        self.assertEqual(Chem.MolToSmiles(mol), expected)
        self.assertEqual(identifier, compound.identifiers[0])

    @parameterized.named_parameters(('SMILES', 'invalid_smiles', 'SMILES'),
                                    ('INCHI', 'invalid_inchi', 'INCHI'))
    def test_mol_from_compound_failures(self, value, identifier_type):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value=value, type=identifier_type)
        with self.assertRaisesRegex(ValueError,
                                    'invalid structural identifier'):
            message_helpers.mol_from_compound(compound)

    def test_get_product_yield(self):
        product = reaction_pb2.ProductCompound()
        self.assertIsNone(message_helpers.get_product_yield(product))
        product.measurements.add(type='YIELD', percentage=dict(value=23))
        self.assertEqual(23, message_helpers.get_product_yield(product))

    def test_check_compound_identifiers(self):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value='c1ccccc1', type='SMILES')
        message_helpers.check_compound_identifiers(compound)
        compound.identifiers.add(value='InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H',
                                 type='INCHI')
        message_helpers.check_compound_identifiers(compound)
        compound.identifiers.add(value='c1ccc(O)cc1', type='SMILES')
        with self.assertRaisesRegex(ValueError, 'inconsistent'):
            message_helpers.check_compound_identifiers(compound)

    def test_get_reaction_smiles(self):
        reaction = reaction_pb2.Reaction()
        reactant1 = reaction.inputs['reactant1']
        reactant1.components.add(reaction_role='REACTANT').identifiers.add(
            value='c1ccccc1', type='SMILES')
        reactant1.components.add(reaction_role='SOLVENT').identifiers.add(
            value='N', type='SMILES')
        self.assertEqual(message_helpers.get_reaction_smiles(reaction),
                         'c1ccccc1>N>')
        reactant2 = reaction.inputs['reactant2']
        reactant2.components.add(reaction_role='REACTANT').identifiers.add(
            value='Cc1ccccc1', type='SMILES')
        reactant2.components.add(reaction_role='SOLVENT').identifiers.add(
            value='N', type='SMILES')
        reaction.outcomes.add().products.add(
            reaction_role=reaction_pb2.ReactionRole.PRODUCT).identifiers.add(
                value='O=C=O', type='SMILES')
        self.assertEqual(message_helpers.get_reaction_smiles(reaction),
                         'Cc1ccccc1.c1ccccc1>N>O=C=O')

    def test_get_reaction_smiles_failure(self):
        reaction = reaction_pb2.Reaction()
        reactant = reaction.inputs['reactant']
        component = reactant.components.add(reaction_role='REACTANT')
        component.identifiers.add(value='benzene', type='NAME')
        with self.assertRaisesRegex(ValueError,
                                    'no valid reactants or products'):
            message_helpers.get_reaction_smiles(reaction)
        component.identifiers.add(value='c1ccccc1', type='SMILES')
        with self.assertRaisesRegex(ValueError, 'must contain at least one'):
            message_helpers.get_reaction_smiles(reaction,
                                                allow_incomplete=False)
        reaction.outcomes.add().products.add(
            reaction_role=reaction_pb2.ReactionRole.PRODUCT).identifiers.add(
                value='invalid', type='SMILES')
        with self.assertRaisesRegex(ValueError, 'reaction contains errors'):
            message_helpers.get_reaction_smiles(reaction,
                                                allow_incomplete=False)


class FindSubmessagesTest(absltest.TestCase):

    def test_scalar(self):
        message = test_pb2.Scalar(int32_value=5, float_value=6.7)
        self.assertEmpty(
            message_helpers.find_submessages(message, test_pb2.Scalar))
        with self.assertRaisesRegex(TypeError, 'must be a Protocol Buffer'):
            message_helpers.find_submessages(message, float)

    def test_nested(self):
        message = test_pb2.Nested()
        self.assertEmpty(
            message_helpers.find_submessages(message, test_pb2.Nested.Child))
        message.child.value = 5.6
        submessages = message_helpers.find_submessages(message,
                                                       test_pb2.Nested.Child)
        self.assertLen(submessages, 1)
        # Show that the returned submessages work as references.
        submessages[0].value = 7.8
        self.assertAlmostEqual(message.child.value, 7.8, places=4)

    def test_repeated_nested(self):
        message = test_pb2.RepeatedNested()
        message.children.add().value = 1.2
        message.children.add().value = 3.4
        self.assertLen(
            message_helpers.find_submessages(message,
                                             test_pb2.RepeatedNested.Child), 2)

    def test_map_nested(self):
        message = test_pb2.MapNested()
        message.children['one'].value = 1.2
        message.children['two'].value = 3.4
        self.assertLen(
            message_helpers.find_submessages(message, test_pb2.MapNested.Child),
            2)

    def test_compounds(self):
        message = reaction_pb2.Reaction()
        message.inputs['test'].components.add().identifiers.add(type='NAME',
                                                                value='aspirin')
        self.assertLen(
            message_helpers.find_submessages(message, reaction_pb2.Compound), 1)


class BuildDataTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        self.data = b'test data'
        self.filename = os.path.join(self.test_subdirectory, 'test.data')
        with open(self.filename, 'wb') as f:
            f.write(self.data)

    def test_build_data(self):
        message = message_helpers.build_data(self.filename,
                                             description='binary data')
        self.assertEqual(message.bytes_value, self.data)
        self.assertEqual(message.description, 'binary data')
        self.assertEqual(message.format, 'data')

    def test_bad_filename(self):
        with self.assertRaisesRegex(ValueError,
                                    'cannot deduce the file format'):
            message_helpers.build_data('testdata', 'no description')


class BuildCompoundTest(parameterized.TestCase, absltest.TestCase):

    def test_smiles_and_name(self):
        compound = message_helpers.build_compound(smiles='c1ccccc1',
                                                  name='benzene')
        expected = reaction_pb2.Compound(identifiers=[
            reaction_pb2.CompoundIdentifier(value='c1ccccc1', type='SMILES'),
            reaction_pb2.CompoundIdentifier(value='benzene', type='NAME')
        ])
        self.assertEqual(compound, expected)

    @parameterized.named_parameters(
        ('mass', '1.2 g', reaction_pb2.Mass(value=1.2, units='GRAM')),
        ('moles', '3.4 mol', reaction_pb2.Moles(value=3.4, units='MOLE')),
        ('volume', '5.6 mL', reaction_pb2.Volume(value=5.6, units='MILLILITER'))
    )
    def test_amount(self, amount, expected):
        compound = message_helpers.build_compound(amount=amount)
        self.assertEqual(
            getattr(compound.amount, compound.amount.WhichOneof('kind')),
            expected)

    @parameterized.named_parameters(('missing_units', '1.2'),
                                    ('negative_mass', '-3.4 g'))
    def test_bad_amount(self, amount):
        with self.assertRaises((KeyError, ValueError)):
            message_helpers.build_compound(amount=amount)

    def test_role(self):
        compound = message_helpers.build_compound(role='solvent')
        self.assertEqual(compound.reaction_role,
                         reaction_pb2.ReactionRole.SOLVENT)

    def test_bad_role(self):
        with self.assertRaisesRegex(KeyError, 'not a supported type'):
            message_helpers.build_compound(role='flavorant')

    def test_is_limiting(self):
        self.assertTrue(
            message_helpers.build_compound(is_limiting=True).is_limiting)
        self.assertFalse(
            message_helpers.build_compound(is_limiting=False).is_limiting)
        self.assertFalse(
            message_helpers.build_compound().HasField('is_limiting'))

    @parameterized.named_parameters(
        ('prep_without_details', 'dried', None,
         reaction_pb2.CompoundPreparation(type='DRIED')),
        ('prep_with_details', 'dried', 'in the fire of the sun',
         reaction_pb2.CompoundPreparation(type='DRIED',
                                          details='in the fire of the sun')),
        ('custom_prep_with_details', 'custom', 'threw it on the ground',
         reaction_pb2.CompoundPreparation(type='CUSTOM',
                                          details='threw it on the ground')))
    def test_prep(self, prep, details, expected):
        compound = message_helpers.build_compound(prep=prep,
                                                  prep_details=details)
        self.assertEqual(compound.preparations[0], expected)

    def test_bad_prep(self):
        with self.assertRaisesRegex(KeyError, 'not a supported type'):
            message_helpers.build_compound(prep='shaken')

    def test_prep_details_without_prep(self):
        with self.assertRaisesRegex(ValueError, 'prep must be provided'):
            message_helpers.build_compound(prep_details='rinsed gently')

    def test_custom_prep_without_details(self):
        with self.assertRaisesRegex(ValueError,
                                    'prep_details must be provided'):
            message_helpers.build_compound(prep='custom')

    def test_vendor(self):
        self.assertEqual(
            message_helpers.build_compound(vendor='Sally').source.vendor,
            'Sally')


class SetSoluteMolesTest(parameterized.TestCase, absltest.TestCase):

    def test_set_solute_moles_should_fail(self):
        solute = message_helpers.build_compound(name='Solute')
        solvent = message_helpers.build_compound(name='Solvent')
        with self.assertRaisesRegex(ValueError, 'defined volume'):
            message_helpers.set_solute_moles(solute, [solvent], '10 mM')

        solute = message_helpers.build_compound(name='Solute', amount='1 mol')
        solvent = message_helpers.build_compound(name='Solvent', amount='1 L')
        with self.assertRaisesRegex(ValueError, 'overwrite'):
            message_helpers.set_solute_moles(solute, [solvent], '10 mM')

    def test_set_solute_moles(self):
        solute = message_helpers.build_compound(name='Solute')
        solvent2 = message_helpers.build_compound(name='Solvent',
                                                  amount='100 mL')
        message_helpers.set_solute_moles(solute, [solvent2], '1 molar')
        self.assertEqual(solute.amount.moles,
                         reaction_pb2.Moles(units='MILLIMOLE', value=100))
        solvent3 = message_helpers.build_compound(name='Solvent',
                                                  amount='75 uL')
        message_helpers.set_solute_moles(solute, [solvent3],
                                         '3 mM',
                                         overwrite=True)
        self.assertEqual(solute.amount.moles,
                         reaction_pb2.Moles(units='NANOMOLE', value=225))
        solvent4 = message_helpers.build_compound(name='Solvent',
                                                  amount='0.2 uL')
        message_helpers.set_solute_moles(solute, [solvent4],
                                         '30 mM',
                                         overwrite=True)
        self.assertEqual(solute.amount.moles,
                         reaction_pb2.Moles(units='NANOMOLE', value=6))
        solvent5 = message_helpers.build_compound(name='Solvent',
                                                  amount='0.8 uL')
        message_helpers.set_solute_moles(solute, [solvent4, solvent5],
                                         '30 mM',
                                         overwrite=True)
        self.assertEqual(solute.amount.moles,
                         reaction_pb2.Moles(units='NANOMOLE', value=30))


class CompoundIdentifiersTest(absltest.TestCase):

    def test_identifier_setters(self):
        compound = reaction_pb2.Compound()
        identifier = message_helpers.set_compound_name(compound, 'water')
        self.assertEqual(
            identifier,
            reaction_pb2.CompoundIdentifier(type='NAME', value='water'))
        self.assertEqual(
            compound.identifiers[0],
            reaction_pb2.CompoundIdentifier(type='NAME', value='water'))
        message_helpers.set_compound_smiles(compound, 'O')
        self.assertEqual(
            compound.identifiers[1],
            reaction_pb2.CompoundIdentifier(type='SMILES', value='O'))
        identifier = message_helpers.set_compound_name(compound, 'ice')
        self.assertEqual(
            identifier, reaction_pb2.CompoundIdentifier(type='NAME',
                                                        value='ice'))
        self.assertEqual(
            compound.identifiers[0],
            reaction_pb2.CompoundIdentifier(type='NAME', value='ice'))
        compound = reaction_pb2.Compound()
        identifier = message_helpers.set_compound_molblock(
            compound, _BENZENE_MOLBLOCK)
        self.assertEqual(_BENZENE_MOLBLOCK, compound.identifiers[0].value)

    def test_identifier_getters(self):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(type='NAME', value='water')
        self.assertEqual(message_helpers.get_compound_name(compound), 'water')
        self.assertIsNone(message_helpers.get_compound_smiles(compound))
        compound.identifiers.add(type='SMILES', value='O')
        self.assertEqual(message_helpers.get_compound_smiles(compound), 'O')
        self.assertEqual(message_helpers.smiles_from_compound(compound), 'O')
        compound = reaction_pb2.Compound()
        compound.identifiers.add(type='MOLBLOCK', value=_BENZENE_MOLBLOCK)
        self.assertEqual(message_helpers.get_compound_molblock(compound),
                         _BENZENE_MOLBLOCK)
        self.assertEqual(message_helpers.molblock_from_compound(compound),
                         _BENZENE_MOLBLOCK)


class SetDativeBondsTest(parameterized.TestCase, absltest.TestCase):

    def test_has_transition_metal(self):
        self.assertFalse(
            message_helpers.has_transition_metal(Chem.MolFromSmiles('P')))
        self.assertTrue(
            message_helpers.has_transition_metal(
                Chem.MolFromSmiles('Cl[Pd]Cl')))

    @parameterized.named_parameters(
        ('Pd(PH3)(NH3)Cl2', '[PH3][Pd](Cl)(Cl)[NH3]', ('N', 'P'),
         'N->[Pd](<-P)(Cl)Cl'))
    def test_set_dative_bonds(self, smiles, from_atoms, expected):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        dative_mol = message_helpers.set_dative_bonds(mol,
                                                      from_atoms=from_atoms)
        self.assertEqual(Chem.MolToSmiles(dative_mol), expected)


class LoadAndWriteMessageTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.messages = [
            test_pb2.Scalar(int32_value=3, float_value=4.5),
            test_pb2.RepeatedScalar(values=[1.2, 3.4]),
            test_pb2.Enum(value='FIRST'),
            test_pb2.RepeatedEnum(values=['FIRST', 'SECOND']),
            test_pb2.Nested(child=test_pb2.Nested.Child(value=1.2)),
        ]

    @parameterized.parameters(message_helpers.MessageFormat)
    def test_round_trip(self, message_format):
        for message in self.messages:
            with tempfile.NamedTemporaryFile(suffix=message_format.value) as f:
                message_helpers.write_message(message, f.name)
                f.flush()
                self.assertEqual(
                    message,
                    message_helpers.load_message(f.name, type(message)))

    def test_bad_binary(self):
        with tempfile.NamedTemporaryFile(suffix='.pb') as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(message.SerializeToString())
            f.flush()
            # NOTE(kearnes): The decoder is not perfect; for example, it will
            # not be able to distinguish from a message with the same tags and
            # types (e.g. test_pb2.Scalar and test_pb2.RepeatedScalar).
            with self.assertRaisesRegex(ValueError, 'Error parsing message'):
                message_helpers.load_message(f.name, test_pb2.Nested)

    def test_bad_json(self):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.json') as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(json_format.MessageToJson(message))
            f.flush()
            with self.assertRaisesRegex(ValueError, 'no field named "values"'):
                message_helpers.load_message(f.name, test_pb2.Nested)

    def test_bad_pbtxt(self):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.pbtxt') as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(text_format.MessageToString(message))
            f.flush()
            with self.assertRaisesRegex(ValueError, 'no field named "values"'):
                message_helpers.load_message(f.name, test_pb2.Nested)

    def test_bad_suffix(self):
        message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
        with self.assertRaisesRegex(ValueError, 'not a valid MessageFormat'):
            message_helpers.write_message(message, 'test.proto')


class CreateMessageTest(parameterized.TestCase, absltest.TestCase):

    @parameterized.named_parameters(
        ('reaction', 'Reaction', reaction_pb2.Reaction),
        ('temperature', 'Temperature', reaction_pb2.Temperature),
        ('temperature measurement', 'TemperatureConditions.Measurement',
         reaction_pb2.TemperatureConditions.Measurement))
    def test_valid_messages(self, message_name, expected_class):
        message = message_helpers.create_message(message_name)
        self.assertIsInstance(message, expected_class)

    @parameterized.named_parameters(('bad case', 'reaction'),
                                    ('gibberish', 'aosdifjasdf'))
    def test_invalid_messages(self, message_name):
        with self.assertRaisesRegex(ValueError, 'Cannot resolve'):
            message_helpers.create_message(message_name)


if __name__ == '__main__':
    absltest.main()
