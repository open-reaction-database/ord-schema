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
        self.assertEmpty(validations.validate_message(message))

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
        with self.assertRaisesRegex(
                validations.ValidationError, expected_error):
            validations.validate_message(message)

    def test_orcid(self):
        message = reaction_pb2.Person(orcid='0000-0001-2345-678X')
        self.assertEmpty(validations.validate_message(message))

    def test_orcid_should_fail(self):
        message = reaction_pb2.Person(orcid='abcd-0001-2345-678X')
        with self.assertRaisesRegex(validations.ValidationError, 'Invalid'):
            validations.validate_message(message)

    def test_reaction(self):
        message = reaction_pb2.Reaction()
        with self.assertRaisesRegex(
                validations.ValidationError, 'reaction input'):
            validations.validate_message(message)

    def test_reaction_recursive(self):
        message = reaction_pb2.Reaction()
        # Reactions must have at least one input
        with self.assertRaisesRegex(
                validations.ValidationError, 'reaction input'):
            validations.validate_message(message, recurse=False)
        dummy_input = message.inputs['dummy_input']
        # Reactions must have at least one outcome
        with self.assertRaisesRegex(
                validations.ValidationError, 'reaction outcome'):
            validations.validate_message(message, recurse=False)
        outcome = message.outcomes.add()
        self.assertEmpty(validations.validate_message(message, recurse=False))
        # Inputs must have at least one component
        with self.assertRaisesRegex(validations.ValidationError, 'component'):
            validations.validate_message(message)
        dummy_component = dummy_input.components.add()
        # Components must have at least one identifier
        with self.assertRaisesRegex(validations.ValidationError, 'identifier'):
            validations.validate_message(message)
        dummy_component.identifiers.add(type='CUSTOM')
        # Custom identifiers must have details specified
        with self.assertRaisesRegex(validations.ValidationError, 'details'):
            validations.validate_message(message)
        dummy_component.identifiers[0].details = 'custom_identifier'
        dummy_component.identifiers[0].value = 'custom_value'
        # Components of reaction inputs must have a defined amount
        with self.assertRaisesRegex(
                validations.ValidationError, 'require an amount'):
            validations.validate_message(message)
        dummy_component.mass.value = 1
        dummy_component.mass.units = reaction_pb2.Mass.GRAM
        # Reactions must have defined products or conversion
        with self.assertRaisesRegex(
                validations.ValidationError, 'products or conversion'):
            validations.validate_message(message)
        outcome.conversion.value = 75
        # If converseions are defined, must have limiting reagent flag
        with self.assertRaisesRegex(
                validations.ValidationError, 'is_limiting'):
            validations.validate_message(message)
        dummy_component.is_limiting = True
        self.assertEmpty(validations.validate_message(message))

        # If an analysis uses an internal standard, a component must have
        # an INTERNAL_STANDARD reaction role
        outcome.analyses['dummy_analysis'].uses_internal_standard = True
        with self.assertRaisesRegex(
                validations.ValidationError, 'INTERNAL_STANDARD'):
            validations.validate_message(message)
        # Assigning internal standard role to input should resolve the error
        message_input_istd = reaction_pb2.Reaction()
        message_input_istd.CopyFrom(message)
        message_input_istd.inputs['dummy_input'].components[0].reaction_role = (
            reaction_pb2.Compound.ReactionRole.INTERNAL_STANDARD)
        self.assertEmpty(validations.validate_message(message_input_istd))
        # Assigning internal standard role to workup should resolve the error
        message_workup_istd = reaction_pb2.Reaction()
        message_workup_istd.CopyFrom(message)
        workup = message_workup_istd.workup.add()
        istd = workup.components.add()
        istd.identifiers.add(type='SMILES', value='CCO')
        istd.mass.value = 1
        istd.mass.units = reaction_pb2.Mass.GRAM
        istd.reaction_role = istd.ReactionRole.INTERNAL_STANDARD
        self.assertEmpty(validations.validate_message(message_workup_istd))

    def test_reaction_recursive_noraise_on_error(self):
        message = reaction_pb2.Reaction()
        message.inputs['dummy_input'].components.add()
        errors = validations.validate_message(message, raise_on_error=False)
        expected = [
            'Compounds must have at least one identifier',
            "Reaction input's components require an amount",
            'Reactions should have at least 1 reaction outcome',
        ]
        self.assertEqual(errors, expected)

    def test_datetimes(self):
        message = reaction_pb2.ReactionProvenance()
        message.experiment_start.value = '11 am'
        message.record_created.time.value = '10 am'
        with self.assertRaisesRegex(validations.ValidationError, 'after'):
            validations.validate_message(message)
        message.record_created.time.value = '11:15 am'
        self.assertEmpty(validations.validate_message(message))

    def test_record_id(self):
        message = reaction_pb2.ReactionProvenance()
        message.record_created.time.value = '10 am'
        message.record_id = 'ord-c0bbd41f095a44a78b6221135961d809'
        self.assertEmpty(validations.validate_message(message))

    @parameterized.named_parameters(
        ('too short', 'ord-c0bbd41f095a4'),
        ('too long', 'ord-c0bbd41f095a4c0bbd41f095a4c0bbd41f095a4'),
        ('bad prefix', 'foo-c0bbd41f095a44a78b6221135961d809'),
        ('bad capitalization', 'ord-C0BBD41F095A44A78B6221135961D809'),
        ('bad characters', 'ord-h0bbd41f095a44a78b6221135961d809'),
        ('bad characters 2', 'ord-notARealId'),
    )
    def test_bad_record_id(self, record_id):
        message = reaction_pb2.ReactionProvenance()
        message.record_created.time.value = '10 am'
        message.record_id = record_id
        with self.assertRaisesRegex(validations.ValidationError, 'malformed'):
            validations.validate_message(message)

    def test_compound_name_resolver(self):
        message = reaction_pb2.Compound()
        identifier = message.identifiers.add()
        identifier.type = identifier.NAME
        identifier.value = 'aspirin'
        validations.validate_message(message)  # Message is modified in place.
        self.assertEqual(
            message.identifiers[1],
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
        validations.validate_message(message)  # Message is modified in place.
        self.assertEqual(
            message.identifiers[1],
            reaction_pb2.CompoundIdentifier(type='RDKIT_BINARY',
                                            bytes_value=mol.ToBinary()))

    def test_data(self):
        message = reaction_pb2.Data()
        with self.assertRaisesRegex(
                validations.ValidationError, 'requires one of'):
            validations.validate_message(message)
        message.bytes_value = b'test data'
        with self.assertRaisesRegex(
                validations.ValidationError, 'format is required'):
            validations.validate_message(message)
        message.value = 'test data'
        self.assertEmpty(validations.validate_message(message))


if __name__ == '__main__':
    absltest.main()
