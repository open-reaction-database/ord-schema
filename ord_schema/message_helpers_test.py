"""Tests for ord_schema.message_helpers."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema import message_helpers
from ord_schema import testing
from ord_schema.proto import ord_schema_pb2 as schema


class BuildCompoundTest(parameterized.TestCase, testing.TestCase):

    def test_smiles_and_name(self):
        compound = message_helpers.build_compound(
            smiles='c1ccccc1', name='benzene')
        expected = schema.Compound(
            identifiers=[
                schema.CompoundIdentifier(value='c1ccccc1', type='SMILES'),
                schema.CompoundIdentifier(value='benzene', type='NAME')])
        self.assertEqual(compound, expected)

    @parameterized.named_parameters(
        ('mass', '1.2 g', schema.Mass(value=1.2, units='GRAM')),
        ('moles', '3.4 mol', schema.Moles(value=3.4, units='MOLES')),
        ('volume', '5.6 mL', schema.Volume(value=5.6, units='MILLILITER')))
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
                         schema.Compound.ReactionRole.SOLVENT)

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
         schema.CompoundPreparation(type='DRIED')),
        ('prep_with_details', 'dried', 'in the fire of the sun',
         schema.CompoundPreparation(type='DRIED',
                                    details='in the fire of the sun')),
        ('custom_prep_with_details', 'custom', 'threw it on the ground',
         schema.CompoundPreparation(type='CUSTOM',
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
