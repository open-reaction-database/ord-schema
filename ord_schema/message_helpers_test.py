"""Tests for ord_schema.message_helpers."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import parameterized
from google import protobuf
from google.protobuf import json_format
from google.protobuf import text_format

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.proto import test_pb2

try:
    from rdkit import Chem
except ImportError:
    Chem = None


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
        submessages = message_helpers.find_submessages(
            message, test_pb2.Nested.Child)
        self.assertLen(submessages, 1)
        # Show that the returned submessages work as references.
        submessages[0].value = 7.8
        self.assertAlmostEqual(message.child.value, 7.8, places=4)

    def test_repeated_nested(self):
        message = test_pb2.RepeatedNested()
        message.children.add().value = 1.2
        message.children.add().value = 3.4
        self.assertLen(
            message_helpers.find_submessages(
                message, test_pb2.RepeatedNested.Child),
            2)

    def test_map_nested(self):
        message = test_pb2.MapNested()
        message.children['one'].value = 1.2
        message.children['two'].value = 3.4
        self.assertLen(
            message_helpers.find_submessages(
                message, test_pb2.MapNested.Child),
            2)


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


class WriteDataTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)

    def test_string_value(self):
        message = reaction_pb2.Data(value='test value')
        filename = message_helpers.write_data(message, self.test_subdirectory)
        expected = os.path.join(
            self.test_subdirectory,
            'ord_data-'
            '47d1d8273710fd6f6a5995fac1a0983fe0e8828c288e35e80450ddc5c4412def'
            '.txt')
        self.assertEqual(filename, expected)
        # NOTE(kearnes): Open with 'r' to get the decoded string.
        with open(filename, 'r') as f:
            self.assertEqual(message.value, f.read())

    def test_bytes_value(self):
        message = reaction_pb2.Data(bytes_value=b'test value')
        filename = message_helpers.write_data(message, self.test_subdirectory)
        expected = os.path.join(
            self.test_subdirectory,
            'ord_data-'
            '47d1d8273710fd6f6a5995fac1a0983fe0e8828c288e35e80450ddc5c4412def'
            '.txt')
        self.assertEqual(filename, expected)
        with open(filename, 'rb') as f:
            self.assertEqual(message.bytes_value, f.read())

    def test_url_value(self):
        message = reaction_pb2.Data(url='test value')
        self.assertIsNone(
            message_helpers.write_data(message, self.test_subdirectory))

    def test_missing_value(self):
        message = reaction_pb2.Data()
        with self.assertRaisesRegex(ValueError, 'no value to write'):
            message_helpers.write_data(message, self.test_subdirectory)

    def test_min_size(self):
        message = reaction_pb2.Data(value='test_value')
        self.assertIsNone(
            message_helpers.write_data(
                message, self.test_subdirectory, min_size=1.0))

    def test_max_size(self):
        message = reaction_pb2.Data(value='test value')
        with self.assertRaisesRegex(ValueError, 'larger than max_size'):
            message_helpers.write_data(
                message, self.test_subdirectory, max_size=1e-6)

    def test_find_data_messages(self):
        message = reaction_pb2.Reaction()
        self.assertEmpty(
            message_helpers.find_submessages(message, reaction_pb2.Data))
        message = reaction_pb2.ReactionObservation()
        message.image.value = 'not an image'
        self.assertLen(
            message_helpers.find_submessages(message, reaction_pb2.Data), 1)
        message = reaction_pb2.ReactionSetup()
        message.automation_code['test1'].value = 'test data 1'
        message.automation_code['test2'].bytes_value = b'test data 2'
        self.assertLen(
            message_helpers.find_submessages(message, reaction_pb2.Data), 2)
        message = reaction_pb2.Reaction()
        message.observations.add().image.value = 'not an image'
        message.setup.automation_code['test1'].value = 'test data 1'
        message.setup.automation_code['test2'].bytes_value = b'test data 2'
        self.assertLen(
            message_helpers.find_submessages(message, reaction_pb2.Data), 3)


class BuildCompoundTest(parameterized.TestCase, absltest.TestCase):

    def test_smiles_and_name(self):
        compound = message_helpers.build_compound(
            smiles='c1ccccc1', name='benzene')
        expected = reaction_pb2.Compound(
            identifiers=[
                reaction_pb2.CompoundIdentifier(value='c1ccccc1',
                                                type='SMILES'),
                reaction_pb2.CompoundIdentifier(value='benzene', type='NAME')])
        self.assertEqual(compound, expected)

    @parameterized.named_parameters(
        ('mass', '1.2 g', reaction_pb2.Mass(value=1.2, units='GRAM')),
        ('moles', '3.4 mol', reaction_pb2.Moles(value=3.4, units='MOLES')),
        ('volume', '5.6 mL',
         reaction_pb2.Volume(value=5.6, units='MILLILITER')))
    def test_amount(self, amount, expected):
        compound = message_helpers.build_compound(amount=amount)
        self.assertEqual(getattr(compound, compound.WhichOneof('amount')),
                         expected)

    @parameterized.named_parameters(('missing_units', '1.2'),
                                    ('negative_mass', '-3.4 g'))
    def test_bad_amount(self, amount):
        with self.assertRaises((KeyError, ValueError)):
            message_helpers.build_compound(amount=amount)

    def test_role(self):
        compound = message_helpers.build_compound(role='solvent')
        self.assertEqual(compound.reaction_role,
                         reaction_pb2.Compound.ReactionRole.SOLVENT)

    def test_bad_role(self):
        with self.assertRaisesRegex(KeyError, 'not a supported type'):
            message_helpers.build_compound(role='flavorant')

    def test_is_limiting(self):
        self.assertTrue(
            message_helpers.build_compound(is_limiting=True).is_limiting)
        self.assertFalse(
            message_helpers.build_compound(is_limiting=False).is_limiting)
        self.assertFalse(message_helpers.build_compound().is_limiting)

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
        self.assertEqual(compound.preparation, expected)

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
            message_helpers.build_compound(vendor='Sally').vendor_source,
            'Sally')


class GetCompoundSmilesTest(absltest.TestCase):

    def test_get_compound_smiles(self):
        compound = message_helpers.build_compound(
            smiles='c1ccccc1', name='benzene')
        self.assertEqual(message_helpers.get_compound_smiles(compound),
                         'c1ccccc1')


class GetCompoundMolTest(absltest.TestCase):

    @absltest.skipIf(Chem is None, 'no rdkit')
    def test_get_compound_mol(self):
        mol = Chem.MolFromSmiles('c1ccccc1')
        compound = message_helpers.build_compound(
            smiles='c1ccccc1', name='benzene')
        identifier = compound.identifiers.add()
        identifier.type = identifier.RDKIT_BINARY
        identifier.bytes_value = mol.ToBinary()
        self.assertEqual(
            Chem.MolToSmiles(mol),
            Chem.MolToSmiles(message_helpers.get_compound_mol(compound)))


class LoadMessageTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.messages = [
            test_pb2.Scalar(int32_value=3, float_value=4.5),
            test_pb2.RepeatedScalar(values=[1.2, 3.4]),
            test_pb2.Enum(value='FIRST'),
            test_pb2.RepeatedEnum(values=['FIRST', 'SECOND']),
            test_pb2.Nested(child=test_pb2.Nested.Child(value=1.2)),
        ]

    def test_binary(self):
        for message in self.messages:
            with tempfile.NamedTemporaryFile() as f:
                f.write(message.SerializeToString())
                f.flush()
                self.assertEqual(
                    message,
                    message_helpers.load_message(
                        f.name, type(message), 'binary'))

    def test_json(self):
        for message in self.messages:
            with tempfile.NamedTemporaryFile(mode='w+') as f:
                f.write(json_format.MessageToJson(message))
                f.flush()
                self.assertEqual(
                    message,
                    message_helpers.load_message(
                        f.name, type(message), 'json'))

    def test_pbtxt(self):
        for message in self.messages:
            with tempfile.NamedTemporaryFile(mode='w+') as f:
                f.write(text_format.MessageToString(message))
                f.flush()
                self.assertEqual(
                    message,
                    message_helpers.load_message(
                        f.name, type(message), 'pbtxt'))

    def test_bad_binary(self):
        with tempfile.NamedTemporaryFile() as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(message.SerializeToString())
            f.flush()
            # NOTE(kearnes): The decoder is not perfect; for example, it will
            # not be able to distinguish from a message with the same tags and
            # types (e.g. test_pb2.Scalar and test_pb2.RepeatedScalar).
            with self.assertRaisesRegex(protobuf.message.DecodeError,
                                        'Error parsing message'):
                message_helpers.load_message(f.name, test_pb2.Nested, 'binary')

    def test_bad_json(self):
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(json_format.MessageToJson(message))
            f.flush()
            # NOTE(kearnes): The decoder is not perfect; for example, it will
            # not be able to distinguish from a message with the same tags and
            # types (e.g. test_pb2.Scalar and test_pb2.RepeatedScalar).
            with self.assertRaisesRegex(json_format.ParseError,
                                        'no field named "values"'):
                message_helpers.load_message(f.name, test_pb2.Nested, 'json')

    def test_bad_pbtxt(self):
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(text_format.MessageToString(message))
            f.flush()
            # NOTE(kearnes): The decoder is not perfect; for example, it will
            # not be able to distinguish from a message with the same tags and
            # types (e.g. test_pb2.Scalar and test_pb2.RepeatedScalar).
            with self.assertRaisesRegex(text_format.ParseError,
                                        'no field named "values"'):
                message_helpers.load_message(f.name, test_pb2.Nested, 'pbtxt')


if __name__ == '__main__':
    absltest.main()
