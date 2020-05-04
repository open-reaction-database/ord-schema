"""Tests for ord_schema.process_dataset."""

import glob
import os
import subprocess
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver

from ord_schema import message_helpers
from ord_schema import process_dataset
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


class ProcessDatasetTest(absltest.TestCase):

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
                                              'dataset1.pbtxt')
        message_helpers.write_message(
            dataset1, self.dataset1_filename, output_format='pbtxt')
        # reaction2 is empty.
        reaction2 = reaction_pb2.Reaction()
        dataset2 = dataset_pb2.Dataset(reactions=[reaction1, reaction2])
        self.dataset2_filename = os.path.join(self.test_subdirectory,
                                              'dataset2.pbtxt')
        message_helpers.write_message(
            dataset2, self.dataset2_filename, output_format='pbtxt')

    def test_main_with_input_pattern(self):
        with flagsaver.flagsaver(input_pattern=self.dataset1_filename):
            process_dataset.main(())

    def test_main_with_input_file(self):
        input_file = os.path.join(self.test_subdirectory, 'input_file.txt')
        with open(input_file, 'w') as f:
            f.write(f'{self.dataset1_filename}\n')
        with flagsaver.flagsaver(input_file=input_file):
            process_dataset.main(())

    def test_main_with_validation_errors(self):
        with flagsaver.flagsaver(input_pattern=self.dataset2_filename,
                                 write_errors=True):
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

    def test_main_with_updates(self):
        output = os.path.join(self.test_subdirectory, 'output.pbtxt')
        with flagsaver.flagsaver(input_pattern=self.dataset1_filename,
                                 update=True,
                                 output=output):
            process_dataset.main(())
        self.assertTrue(os.path.exists(output))
        dataset = message_helpers.load_message(
            output, dataset_pb2.Dataset, input_format='pbtxt')
        self.assertLen(dataset.reactions, 1)
        self.assertStartsWith(dataset.reactions[0].provenance.record_id, 'ord-')

    def test_main_with_too_many_flags(self):
        with flagsaver.flagsaver(input_pattern=self.dataset1_filename,
                                 input_file=self.dataset2_filename):
            with self.assertRaisesRegex(ValueError, 'not both'):
                process_dataset.main(())

    def test_preserve_existing_dataset_id(self):
        dataset = dataset_pb2.Dataset(
            reactions=[reaction_pb2.Reaction()],
            dataset_id='64b14868c5cd46dd8e75560fd3589a6b')
        filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(
            dataset, filename, output_format='pbtxt')
        expected_filename = os.path.join(
            self.test_subdirectory, 'data', '64',
            '64b14868c5cd46dd8e75560fd3589a6b.pbtxt')
        self.assertFalse(os.path.exists(expected_filename))
        with flagsaver.flagsaver(root=self.test_subdirectory,
                                 input_pattern=filename,
                                 validate=False,
                                 update=True):
            process_dataset.main(())
        self.assertTrue(os.path.exists(expected_filename))

    def test_bad_dataset_id(self):
        dataset = dataset_pb2.Dataset(
            reactions=[reaction_pb2.Reaction()],
            dataset_id='not-a-real-dataset-id')
        filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(
            dataset, filename, output_format='pbtxt')
        with flagsaver.flagsaver(root=self.test_subdirectory,
                                 input_pattern=filename,
                                 validate=False,
                                 update=True):
            with self.assertRaisesRegex(ValueError, 'malformed dataset ID'):
                process_dataset.main(())


class SubmissionWorkflowTest(absltest.TestCase):
    """Test suite for the ORD submission workflow.

    setUp() starts each test with a clean git environment containing some
    data. To create a new test, you should made a modification to the git
    repo (e.g. adding a new dataset or editing an existing one) and call
    self._run_main() to commit the changes and run process_datasets.py.
    """

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        os.chdir(self.test_subdirectory)
        subprocess.run(['git', 'init'], check=True)
        # Add some initial data.
        reaction = reaction_pb2.Reaction()
        methylamine = reaction.inputs['methylamine']
        component = methylamine.components.add()
        component.identifiers.add(type='SMILES', value='CN')
        component.is_limiting = True
        component.moles.value = 1
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 75
        dataset_id = '64b14868c5cd46dd8e75560fd3589a6b'
        dataset = dataset_pb2.Dataset(reactions=[reaction],
                                      dataset_id=dataset_id)
        # Make sure the initial dataset is valid.
        validations.validate_message(dataset)
        os.makedirs(os.path.join('data', '64'))
        self.dataset_filename = os.path.join(
            'data', f'{dataset_id[:2]}', f'{dataset_id}.pbtxt')
        message_helpers.write_message(
            dataset, self.dataset_filename, output_format='pbtxt')
        subprocess.run(['git', 'add', 'data'], check=True)
        subprocess.run(['git', 'commit', '-m', 'Initial commit'], check=True)

    def _run_main(self):
        subprocess.run(['git', 'add', '*.pbtxt', 'data/*/*.pbtxt'], check=True)
        changed = subprocess.run(['git', 'diff', '--name-only', '--staged'],
                                 check=True,
                                 capture_output=True)
        with open('changed.txt', 'wb') as f:
            f.write(changed.stdout)
        subprocess.run(['git', 'commit', '-m', 'Submission'], check=True)
        with flagsaver.flagsaver(
                input_file='changed.txt', update=True, cleanup=True):
            process_dataset.main(())
        return glob.glob(
            os.path.join(self.test_subdirectory, '**/*.pbtxt'),
            recursive=True)

    def test_new_dataset(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(
            dataset, dataset_filename, output_format='pbtxt')
        filenames = self._run_main()
        self.assertLen(filenames, 2)
        self.assertFalse(os.path.exists(dataset_filename))


if __name__ == '__main__':
    absltest.main()
