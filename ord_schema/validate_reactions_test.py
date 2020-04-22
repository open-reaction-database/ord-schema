"""Tests for ord_schema.validation_reactions."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from google.protobuf import text_format

from ord_schema import validate_reactions
from ord_schema.proto import reaction_pb2


class ValidateReactionsTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        reaction1 = reaction_pb2.Reaction()
        dummy_input = reaction1.inputs['dummy_input']
        dummy_component = dummy_input.components.add()
        dummy_component.identifiers.add(type='CUSTOM')
        dummy_component.identifiers[0].details = 'custom_identifier'
        dummy_component.identifiers[0].value = 'custom_value'
        dummy_component.is_limiting = True
        dummy_component.mass.value = 1
        dummy_component.mass.units = reaction_pb2.Mass.GRAM
        reaction1.outcomes.add().conversion.value = 75
        with open(os.path.join(self.test_subdirectory, 'reaction1.pbtxt'),
                  'w') as f:
            f.write(text_format.MessageToString(reaction1))
        reaction2 = reaction_pb2.Reaction()
        with open(os.path.join(self.test_subdirectory, 'reaction2.pbtxt'),
                  'w') as f:
            f.write(text_format.MessageToString(reaction2))
        with open(os.path.join(self.test_subdirectory, 'reaction3.pbtxt'),
                  'w') as f:
            f.write('garbage that is not a reaction proto')

    def test_main(self):
        input_pattern = os.path.join(self.test_subdirectory, 'reaction1.pbtxt')
        output = os.path.join(self.test_subdirectory, 'output.txt')
        with flagsaver.flagsaver(input_pattern=input_pattern, output=output):
            validate_reactions.main(())
        self.assertFalse(os.path.exists(output))

    def test_main_with_input_file(self):
        input_file = os.path.join(self.test_subdirectory, 'input_file.txt')
        output = os.path.join(self.test_subdirectory, 'output.txt')
        with open(input_file, 'w') as f:
            filename = os.path.join(self.test_subdirectory, 'reaction1.pbtxt')
            f.write(f'{filename}\n')
        with flagsaver.flagsaver(input_file=input_file, output=output):
            validate_reactions.main(())
        self.assertFalse(os.path.exists(output))

    def test_main_with_errors(self):
        input_pattern = os.path.join(self.test_subdirectory, 'reaction*.pbtxt')
        output = os.path.join(self.test_subdirectory, 'output.txt')
        with flagsaver.flagsaver(input_pattern=input_pattern, output=output):
            with self.assertRaisesRegex(SystemExit,
                                        'validation encountered errors'):
                validate_reactions.main(())
        self.assertTrue(os.path.exists(output))
        expected_output = [
            'reaction2.pbtxt: Reactions should have '
            'at least 1 reaction input\n',
            'reaction2.pbtxt: Reactions should have '
            'at least 1 reaction outcome\n',
            'reaction3.pbtxt: 1:1 : Message type "ord.Reaction" '
            'has no field named "garbage".\n',
        ]
        with open(output) as f:
            self.assertEqual(f.readlines(), expected_output)

    def test_main_no_input(self):
        with flagsaver.flagsaver():
            with self.assertRaisesRegex(ValueError, 'is required'):
                validate_reactions.main(())

    def test_main_ambiguous_input(self):
        with flagsaver.flagsaver(input_pattern='foo', input_file='bar'):
            with self.assertRaisesRegex(ValueError, 'not both'):
                validate_reactions.main(())


if __name__ == '__main__':
    absltest.main()
