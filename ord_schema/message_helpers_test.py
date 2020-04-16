"""Tests for ord_schema.message_helpers."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import parameterized

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2


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
            message_helpers.build_compound(role='product')

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


if __name__ == '__main__':
    absltest.main()
