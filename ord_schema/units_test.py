"""Tests for ord_schema.units."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema import units
from ord_schema.proto import reaction_pb2


class UnitsTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self._resolver = units.UnitResolver()

    @parameterized.named_parameters(
        ('capitalized', '15.0 ML',
         reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER)),
        ('integer', '24 H',
         reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR)),
        ('no space', '32.1g',
         reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM)),
        ('extra space', '   32.1      \t   g  ',
         reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM)),
        ('lengths', ' 10 millimeter',
         reaction_pb2.Length(value=10, units=reaction_pb2.Length.MILLIMETER)),
    )
    def test_resolve(self, string, expected):
        self.assertEqual(self._resolver.resolve(string), expected)

    @parameterized.named_parameters(
        ('bad units', '1.21 GW', 'unrecognized units'),
        ('multiple matches', '15.0 ML 20.0 L',
         'string does not contain a value with units'),
        ('extra period', '15.0. ML',
         'string does not contain a value with units'),
        ('ambiguous units', '5.2 m', 'ambiguous'),
    )
    def test_resolve_should_fail(self, string, expected_error):
        with self.assertRaisesRegex((KeyError, ValueError), expected_error):
            self._resolver.resolve(string)


if __name__ == '__main__':
    absltest.main()
