"""Tests for ord_schema.units."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema import validations
from ord_schema.proto import ord_schema_pb2 as schema


class ValidationsTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()

    @parameterized.named_parameters(
        ('volume', schema.Volume(value=15.0, units=schema.Volume.MILLILITER)),
        ('time', schema.Time(value=24, units=schema.Time.HOUR)),
        ('mass', schema.Mass(value=32.1, units=schema.Mass.GRAM)),
    )
    def test_units(self, message):
        self.assertEqual(validations.ValidateMessage(message), message)

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
        ('low temperature', schema.Temperature(value=-5, units='KELVIN'),
            'between'),
        ('low temperature 2', schema.Temperature(value=-500, units='CELSIUS'),
            'between'),
    )
    def test_units_should_fail(self, message, expected_error):
        with self.assertRaisesRegex((ValueError), expected_error):
            validations.ValidateMessage(message)

    def test_orcid(self):
        message = schema.Person(orcid='0000-0001-2345-678X')
        self.assertEqual(validations.ValidateMessage(message), message)

    def test_orcid_should_fail(self):
        message = schema.Person(orcid='abcd-0001-2345-678X')
        with self.assertRaisesRegex((ValueError), 'Invalid'):
            validations.ValidateMessage(message)

    def test_reaction(self):
        message = schema.Reaction()
        with self.assertRaisesRegex((ValueError), 'reaction input'):
            validations.ValidateReaction(message)

    def test_reaction_recursive(self):
        message = schema.Reaction()
        with self.assertRaisesRegex((ValueError), 'reaction input'):
            validations.ValidateMessage(message)
        dummy_input = message.inputs['dummy_input']
        self.assertEqual(validations.ValidateMessage(message, recurse=False), 
            message)
        with self.assertRaisesRegex((ValueError), 'component'):
            validations.ValidateMessage(message)
        dummy_component = dummy_input.components.add()
        with self.assertRaisesRegex((ValueError), 'identifier'):
            validations.ValidateMessage(message)
        dummy_component.identifiers.add(type='CUSTOM')
        with self.assertRaisesRegex((ValueError), 'details'):
            validations.ValidateMessage(message)
        dummy_component.identifiers[0].details = 'custom_identifier'
        self.assertEqual(validations.ValidateMessage(message), message)
        outcome = message.outcomes.add()
        analysis = outcome.analyses['dummy_analysis']
        self.assertEqual(validations.ValidateMessage(message), message)

    def test_datetimes(self):
        message = schema.ReactionProvenance()
        message.experiment_start.value = '11 am'
        message.record_created.value = '10 am'
        with self.assertRaisesRegex((ValueError), 'after'):
            validations.ValidateMessage(message)
        message.record_created.value = '11:15 am'
        self.assertEqual(validations.ValidateMessage(message), message)

if __name__ == '__main__':
    absltest.main()
