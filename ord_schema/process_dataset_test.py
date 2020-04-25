"""Tests for ord_schema.process_dataset."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver

from ord_schema import process_dataset
from ord_schema.proto import dataset_pb2
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
        dataset1 = dataset_pb2.Dataset(reactions=[reaction1])
        self.dataset1_filename = os.path.join(self.test_subdirectory,
                                              'dataset1.pb')
        with open(self.dataset1_filename, 'wb') as f:
            f.write(dataset1.SerializeToString())
        # reaction2 is empty.
        reaction2 = reaction_pb2.Reaction()
        dataset2 = dataset_pb2.Dataset(reactions=[reaction1, reaction2])
        self.dataset2_filename = os.path.join(self.test_subdirectory,
                                              'dataset2.pb')
        with open(self.dataset2_filename, 'wb') as f:
            f.write(dataset2.SerializeToString())

    def test_main(self):
        output = os.path.join(self.test_subdirectory, 'output.pb')
        with flagsaver.flagsaver(input_pattern=self.dataset1_filename,
                                 output=output):
            process_dataset.main(())
        self.assertTrue(os.path.exists(output))

    def test_main_with_errors(self):
        output = os.path.join(self.test_subdirectory, 'output.pb')
        with flagsaver.flagsaver(input_pattern=self.dataset2_filename,
                                 write_errors=True,
                                 output=output):
            with self.assertRaisesRegex(ValueError,
                                        'validation encountered errors'):
                process_dataset.main(())
        error_filename = f'{self.dataset2_filename}.error'
        self.assertTrue(os.path.exists(error_filename))
        expected_output = [
            'Reactions should have at least 1 reaction input\n',
            'Reactions should have at least 1 reaction outcome\n',
        ]
        with open(error_filename) as f:
            self.assertEqual(f.readlines(), expected_output)


if __name__ == '__main__':
    absltest.main()
