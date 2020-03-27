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
        ('orcid', schema.Person(orcid='0000-0001-2345-678X'))
    )
    def test_resolve(self, message):
        self.assertEqual(self.validations.ValidateMessage(message), message)

    @parameterized.named_parameters(
        ('neg volume', 
            schema.Volume(value=-15.0, units=schema.Volume.MILLILITER),
            'nonnegative'),
        ('neg time', schema.Time(value=-24, units=schema.Time.HOUR), 
            'nonnegative'),
        ('neg mass', schema.Mass(value=-32.1, units=schema.Mass.GRAM),
            'nonnegative'),
        ('invalid ORCID', schema.Person(orcid='abcd-0001-2345-678X'), 
            'invalid'),
    )
    def test_resolve_should_fail(self, message, expected_error):
        with self.assertRaisesRegex((ValueError), expected_error):
            self.validations.ValidateMessage(message)


if __name__ == '__main__':
    absltest.main()
