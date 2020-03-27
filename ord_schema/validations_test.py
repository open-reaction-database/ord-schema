"""Tests for ord_schema.units."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema import validations
from ord_schema.proto import ord_schema_pb2 as schema


class UnitsTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.validations = validations

    @parameterized.named_parameters(
        ('volume', schema.Volume(value=15.0, units=schema.Volume.MILLILITER)),
        ('time', schema.Time(value=24, units=schema.Time.HOUR)),
        ('mass', schema.Mass(value=32.1, units=schema.Mass.GRAM)),
    )
    def test_units(self, message):
        self.assertEqual(self.validations.ValidateMessage(message), message)

    @parameterized.named_parameters(
        ('neg volume', 
            schema.Volume(value=-15.0, units=schema.Volume.MILLILITER),
            'non-negative'),
        ('neg time', schema.Time(value=-24, units=schema.Time.HOUR), 
            'non-negative'),
        ('neg mass', schema.Mass(value=-32.1, units=schema.Mass.GRAM),
            'non-negative'),
        ('no units', schema.FlowRate(value=5), 'units'),
        ('percentage out of range', schema.Percentage(value=200), 'between'),
    )
    def test_units_should_fail(self, message, expected_error):
        with self.assertRaisesRegex((ValueError), expected_error):
            self.validations.ValidateMessage(message)

    def test_orcid(self):
        message = schema.Person(orcid='0000-0001-2345-678X')
        self.assertEqual(self.validations.ValidateMessage(message), message)

    def test_orcid_should_fail(self):
        message = schema.Person(orcid='abcd-0001-2345-678X')
        with self.assertRaisesRegex((ValueError), 'Invalid'):
            self.validations.ValidateMessage(message)

    def test_reaction(self):
        message = schema.Reaction()
        with self.assertRaisesRegex((ValueError), 'reaction input'):
            self.validations.ValidateReaction(message)
        

    def test_reaction_recursive(self):
        message = schema.Reaction()
        with self.assertRaisesRegex((ValueError), 'reaction input'):
            self.validations.ValidateMessage(message, recurse=False)
        dummy_input = message.inputs['dummy_input']
        with self.assertRaisesRegex((ValueError), 'component'):
            self.validations.ValidateMessage(message)
        dummy_component = dummy_input.components.add()
        with self.assertRaisesRegex((ValueError), 'identifier'):
            self.validations.ValidateMessage(message)
        dummy_component.identifiers.add(type='CUSTOM')
        with self.assertRaisesRegex((ValueError), 'details'):
            self.validations.ValidateMessage(message)
        dummy_component.identifiers[0].details = 'custom_identifier'
        self.assertEqual(self.validations.ValidateMessage(message), message)
        outcome = message.outcomes.add()
        analysis = outcome.analyses['dummy_analysis']
        self.assertEqual(self.validations.ValidateMessage(message), message)

if __name__ == '__main__':
    absltest.main()
