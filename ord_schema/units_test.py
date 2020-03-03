"""Tests for ord_schema.units."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema import units
from ord_schema.proto import ord_schema_pb2 as schema


class UnitsTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self._resolver = units.UnitResolver()

    @parameterized.named_parameters(
        ('capitalized', '15.0 ML',
         schema.Volume(value=15.0, units=schema.Volume.VOLUME_UNIT_MILLILITER)),
        ('integer', '24 H',
         schema.Time(value=24, units=schema.Time.TIME_UNIT_HOUR)),
        ('no space', '32.1g',
         schema.Mass(value=32.1, units=schema.Mass.MASS_UNIT_GRAM)),
        ('extra space', '   32.1      \t   g  ',
         schema.Mass(value=32.1, units=schema.Mass.MASS_UNIT_GRAM)),
    )
    def test_resolve(self, string, expected):
        self.assertEqual(self._resolver.resolve(string), expected)

    @parameterized.named_parameters(
        ('bad units', '1.21 GW', 'unrecognized units'),
        ('multiple matches', '15.0 ML 20.0 L',
         'string does not contain a value with units'),
        ('extra period', '15.0. ML',
         'string does not contain a value with units'),
    )
    def test_resolve_should_fail(self, string, expected_error):
        with self.assertRaisesRegex((KeyError, ValueError), expected_error):
            self._resolver.resolve(string)


if __name__ == '__main__':
    absltest.main()
