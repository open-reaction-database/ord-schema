"""Tests for ord_schema.scripts.build_dataset."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import build_dataset


class BuildDatasetTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        reaction1 = reaction_pb2.Reaction()
        dummy_input = reaction1.inputs['dummy_input']
        dummy_component = dummy_input.components.add()
        dummy_component.identifiers.add(type='CUSTOM')
        dummy_component.identifiers[0].details = 'custom_identifier'
        dummy_component.identifiers[0].value = 'custom_value'
        dummy_component.is_limiting = reaction_pb2.Boolean.TRUE
        dummy_component.mass.value = 1
        dummy_component.mass.units = reaction_pb2.Mass.GRAM
        reaction1.outcomes.add().conversion.value = 75
        message_helpers.write_message(
            reaction1,
            os.path.join(self.test_subdirectory, 'reaction-1.pbtxt'))
        # reaction2 is empty.
        reaction2 = reaction_pb2.Reaction()
        message_helpers.write_message(
            reaction2,
            os.path.join(self.test_subdirectory, 'reaction-2.pbtxt'))

    def test_simple(self):
        input_pattern = os.path.join(self.test_subdirectory, 'reaction-1.pbtxt')
        output_filename = os.path.join(self.test_subdirectory, 'dataset.pbtxt')
        with flagsaver.flagsaver(
                input=input_pattern,
                name='test dataset',
                description='this is a test dataset',
                output=output_filename):
            build_dataset.main(())
        self.assertTrue(os.path.exists(output_filename))
        dataset = message_helpers.load_message(
            output_filename, dataset_pb2.Dataset)
        self.assertEqual(dataset.name, 'test dataset')
        self.assertEqual(dataset.description, 'this is a test dataset')
        self.assertLen(dataset.reactions, 1)

    def test_validation(self):
        input_pattern = os.path.join(self.test_subdirectory, 'reaction-?.pbtxt')
        output_filename = os.path.join(self.test_subdirectory, 'dataset.pbtxt')
        with flagsaver.flagsaver(
                input=input_pattern,
                name='test dataset',
                description='this is a test dataset',
                output=output_filename):
            with self.assertRaisesRegex(
                    validations.ValidationError,
                    'Reactions should have at least 1 reaction input'):
                build_dataset.main(())
        # Make sure disabling validation works.
        with flagsaver.flagsaver(
                input=input_pattern,
                name='test dataset',
                description='this is a test dataset',
                output=output_filename,
                validate=False):
            build_dataset.main(())


if __name__ == '__main__':
    absltest.main()
