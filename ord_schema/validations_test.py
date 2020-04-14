"""Tests for ord_schema.units."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema import validations
from ord_schema.proto import reaction_pb2

try:
    from rdkit import Chem
except ImportError:
    Chem = None


class ValidationsTest(parameterized.TestCase, absltest.TestCase):

    @parameterized.named_parameters(
        ('volume',
         reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER)),
        ('time', reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR)),
        ('mass', reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM)),
    )
    def test_units(self, message):
        self.assertEqual(validations.validate_message(message), message)

    @parameterized.named_parameters(
        ('neg volume',
         reaction_pb2.Volume(
             value=-15.0, units=reaction_pb2.Volume.MILLILITER),
         'non-negative'),
        ('neg time', reaction_pb2.Time(value=-24, units=reaction_pb2.Time.HOUR),
         'non-negative'),
        ('neg mass',
         reaction_pb2.Mass(value=-32.1, units=reaction_pb2.Mass.GRAM),
         'non-negative'),
        ('no units', reaction_pb2.FlowRate(value=5), 'units'),
        ('percentage out of range', reaction_pb2.Percentage(value=200),
         'between'),
        ('low temperature', reaction_pb2.Temperature(value=-5, units='KELVIN'),
         'between'),
        ('low temperature 2',
         reaction_pb2.Temperature(value=-500, units='CELSIUS'),
         'between'),
    )
    def test_units_should_fail(self, message, expected_error):
        with self.assertRaisesRegex(ValueError, expected_error):
            validations.validate_message(message)

    def test_orcid(self):
        message = reaction_pb2.Person(orcid='0000-0001-2345-678X')
        self.assertEqual(validations.validate_message(message), message)

    def test_orcid_should_fail(self):
        message = reaction_pb2.Person(orcid='abcd-0001-2345-678X')
        with self.assertRaisesRegex(ValueError, 'Invalid'):
            validations.validate_message(message)

    def test_reaction(self):
        message = reaction_pb2.Reaction()
        with self.assertRaisesRegex(ValueError, 'reaction input'):
            validations.validate_reaction(message)

    def test_reaction_recursive(self):
        message = reaction_pb2.Reaction()
        with self.assertRaisesRegex(ValueError, 'reaction input'):
            validations.validate_message(message)
        dummy_input = message.inputs['dummy_input']
        self.assertEqual(validations.validate_message(message, recurse=False),
                         message)
        with self.assertRaisesRegex(ValueError, 'component'):
            validations.validate_message(message)
        dummy_component = dummy_input.components.add()

        with self.assertRaisesRegex(ValueError, 'identifier'):
            validations.validate_message(message)
        dummy_component.identifiers.add(type='CUSTOM')
        with self.assertRaisesRegex(ValueError, 'details'):
            validations.validate_message(message)
        dummy_component.identifiers[0].details = 'custom_identifier'
        dummy_component.identifiers[0].value = 'custom_value'
        with self.assertRaisesRegex(ValueError, 'require an amount'):
            validations.validate_message(message)
        dummy_component.mass.value = 1
        dummy_component.mass.units = reaction_pb2.Mass.GRAM
        self.assertEqual(validations.validate_message(message), message)
        outcome = message.outcomes.add()
        _ = outcome.analyses['dummy_analysis']
        self.assertEqual(validations.validate_message(message), message)

    def test_datetimes(self):
        message = reaction_pb2.ReactionProvenance()
        message.experiment_start.value = '11 am'
        message.record_created.time.value = '10 am'
        with self.assertRaisesRegex(ValueError, 'after'):
            validations.validate_message(message)
        message.record_created.time.value = '11:15 am'
        self.assertEqual(validations.validate_message(message), message)

    def test_compound_name_resolver(self):
        message = reaction_pb2.Compound()
        identifier = message.identifiers.add()
        identifier.type = identifier.NAME
        identifier.value = 'aspirin'
        self.assertEqual(
            validations.validate_message(message).identifiers[1],
            reaction_pb2.CompoundIdentifier(type='SMILES',
                                            value='CC(=O)OC1=CC=CC=C1C(=O)O',
                                            details='NAME resolved by PubChem'))

    @absltest.skipIf(Chem is None, 'no rdkit')
    def test_compound_rdkit_binary(self):
        mol = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
        message = reaction_pb2.Compound()
        identifier = message.identifiers.add()
        identifier.type = identifier.SMILES
        identifier.value = Chem.MolToSmiles(mol)
        self.assertEqual(
            validations.validate_message(message).identifiers[1],
            reaction_pb2.CompoundIdentifier(type='RDKIT_BINARY',
                                            bytes_value=mol.ToBinary()))

    def test_data(self):
        message = reaction_pb2.Data()
        with self.assertRaisesRegex(ValueError, 'requires one of'):
            validations.validate_message(message)
        message.bytes_value = b'test data'
        with self.assertRaisesRegex(ValueError, 'requires a description'):
            validations.validate_message(message)
        message.description = 'test data'
        with self.assertRaisesRegex(ValueError, 'format is required'):
            validations.validate_message(message)
        message.value = 'test data'
        validations.validate_message(message)


if __name__ == '__main__':
    absltest.main()
