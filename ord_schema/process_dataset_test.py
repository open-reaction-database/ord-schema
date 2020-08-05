# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for ord_schema.process_dataset."""

import glob
import os
import subprocess
import tempfile

from absl import flags
from absl import logging
from absl.testing import absltest
from absl.testing import flagsaver
from rdkit import RDLogger

from ord_schema import message_helpers
from ord_schema import process_dataset
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


class ProcessDatasetTest(absltest.TestCase):
    def setUp(self):
        super().setUp()
        # Suppress RDKit warnings to clean up the test output.
        RDLogger.logger().setLevel(RDLogger.CRITICAL)
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
        reaction1.provenance.record_created.time.value = '11 am'
        reaction1.provenance.record_created.person.username = 'test'
        reaction1.provenance.record_created.person.email = 'test@example.com'
        dataset1 = dataset_pb2.Dataset(reactions=[reaction1])
        self.dataset1_filename = os.path.join(self.test_subdirectory,
                                              'dataset1.pbtxt')
        message_helpers.write_message(dataset1, self.dataset1_filename)
        # reaction2 is empty.
        reaction2 = reaction_pb2.Reaction()
        dataset2 = dataset_pb2.Dataset(reactions=[reaction1, reaction2])
        self.dataset2_filename = os.path.join(self.test_subdirectory,
                                              'dataset2.pbtxt')
        message_helpers.write_message(dataset2, self.dataset2_filename)

    def test_main_with_input_pattern(self):
        with flagsaver.flagsaver(input_pattern=self.dataset1_filename):
            process_dataset.main(())

    def test_main_with_input_file(self):
        input_file = os.path.join(self.test_subdirectory, 'input_file.txt')
        with open(input_file, 'w') as f:
            f.write(f'A\t{self.dataset1_filename}\n')
        with flagsaver.flagsaver(input_file=input_file):
            process_dataset.main(())

    def test_main_with_validation_errors(self):
        with flagsaver.flagsaver(input_pattern=self.dataset2_filename,
                                 write_errors=True):
            with self.assertRaisesRegex(validations.ValidationError,
                                        'validation encountered errors'):
                process_dataset.main(())
        error_filename = f'{self.dataset2_filename}.error'
        self.assertTrue(os.path.exists(error_filename))
        expected_output = [
            'Reaction: Reactions should have at least 1 reaction input\n',
            'Reaction: Reactions should have at least 1 reaction outcome\n',
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
        dataset = message_helpers.load_message(output, dataset_pb2.Dataset)
        self.assertLen(dataset.reactions, 1)
        self.assertStartsWith(dataset.reactions[0].reaction_id, 'ord-')

    def test_main_with_too_many_flags(self):
        with flagsaver.flagsaver(input_pattern=self.dataset1_filename,
                                 input_file=self.dataset2_filename):
            with self.assertRaisesRegex(ValueError, 'not both'):
                process_dataset.main(())


class SubmissionWorkflowTest(absltest.TestCase):
    """Test suite for the ORD submission workflow.

    setUp() starts each test with a clean git environment containing some
    data. To create a new test, you should made a modification to the git
    repo (e.g. adding a new dataset or editing an existing one) and call
    self._run() to commit the changes and run process_datasets.py.
    """
    _DEFAULT_BRANCH = 'main'

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        os.chdir(self.test_subdirectory)
        subprocess.run(['git', 'init', '-b', self._DEFAULT_BRANCH], check=True)
        subprocess.run(
            ['git', 'config', '--local', 'user.email', 'test@ord-schema'],
            check=True)
        subprocess.run(
            ['git', 'config', '--local', 'user.name', 'Test Runner'],
            check=True)
        # Add some initial data.
        reaction = reaction_pb2.Reaction()
        methylamine = reaction.inputs['methylamine']
        component = methylamine.components.add()
        component.identifiers.add(type='SMILES', value='CN')
        component.is_limiting = True
        component.moles.value = 1
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 75
        reaction.provenance.record_created.time.value = '2020-01-01'
        reaction.provenance.record_created.person.username = 'test'
        reaction.provenance.record_created.person.email = 'test@example.com'
        reaction.reaction_id = 'ord-10aed8b5dffe41fab09f5b2cc9c58ad9'
        dataset_id = 'ord_dataset-64b14868c5cd46dd8e75560fd3589a6b'
        dataset = dataset_pb2.Dataset(reactions=[reaction],
                                      dataset_id=dataset_id)
        # Make sure the initial dataset is valid.
        validations.validate_message(dataset)
        os.makedirs(os.path.join('data', '64'))
        self.dataset_filename = os.path.join(self.test_subdirectory, 'data',
                                             '64', f'{dataset_id}.pbtxt')
        message_helpers.write_message(dataset, self.dataset_filename)
        subprocess.run(['git', 'add', 'data'], check=True)
        subprocess.run(['git', 'commit', '-m', 'Initial commit'], check=True)
        # Use a new branch for tests.
        subprocess.run(['git', 'checkout', '-b', 'test'], check=True)

    def _run(self, **kwargs):
        """Runs process_dataset.main().

        Args:
            **kwargs: Keyword arguments to set flags.

        Returns:
            added: Set of added reaction IDs.
            removed: Set of deleted reaction IDs.
            filenames: List of .pbtxt filenames in the updated database.
        """
        subprocess.run(['git', 'add', '*.pbtxt', 'data/*/*.pbtxt'], check=True)
        changed = subprocess.run(
            ['git', 'diff', '--name-status', self._DEFAULT_BRANCH],
            check=True,
            capture_output=True,
            text=True)
        with open('changed.txt', 'w') as f:
            f.write(changed.stdout)
        logging.info(f'Changed files:\n{changed.stdout}')
        subprocess.run(['git', 'commit', '-m', 'Submission'], check=True)
        run_flags = {
            'input_file': 'changed.txt',
            'update': True,
            'cleanup': True
        }
        run_flags.update(kwargs)
        with flagsaver.flagsaver(**run_flags):
            added, removed = process_dataset.run()
        filenames = glob.glob(os.path.join(self.test_subdirectory,
                                           '**/*.pbtxt'),
                              recursive=True)
        return added, removed, filenames

    def test_add_dataset(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = '2020-01-01'
        reaction.provenance.record_created.person.username = 'test'
        reaction.provenance.record_created.person.email = 'test@example.com'
        reaction.reaction_id = 'test'
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(dataset, dataset_filename)
        added, removed, filenames = self._run()
        self.assertEqual(added, {'test'})
        self.assertEmpty(removed)
        self.assertLen(filenames, 2)
        self.assertFalse(os.path.exists(dataset_filename))
        # Check for assignment of dataset and reaction IDs.
        filenames.pop(filenames.index(self.dataset_filename))
        self.assertLen(filenames, 1)
        dataset = message_helpers.load_message(filenames[0],
                                               dataset_pb2.Dataset)
        self.assertNotEmpty(dataset.dataset_id)
        self.assertLen(dataset.reactions, 1)
        self.assertNotEmpty(dataset.reactions[0].reaction_id)

    def test_add_sharded_dataset(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = '2020-01-02'
        reaction.provenance.record_created.person.username = 'test2'
        reaction.provenance.record_created.person.email = 'test2@example.com'
        reaction.reaction_id = 'test1'
        dataset1 = dataset_pb2.Dataset(reactions=[reaction])
        dataset1_filename = os.path.join(self.test_subdirectory, 'test1.pbtxt')
        message_helpers.write_message(dataset1, dataset1_filename)
        reaction.provenance.record_created.time.value = '2020-01-03'
        reaction.provenance.record_created.person.username = 'test3'
        reaction.provenance.record_created.person.email = 'test3@example.com'
        reaction.reaction_id = 'test2'
        dataset2 = dataset_pb2.Dataset(reactions=[reaction])
        dataset2_filename = os.path.join(self.test_subdirectory, 'test2.pbtxt')
        message_helpers.write_message(dataset2, dataset2_filename)
        added, removed, filenames = self._run()
        self.assertEqual(added, {'test1', 'test2'})
        self.assertEmpty(removed)
        self.assertLen(filenames, 2)
        filenames.pop(filenames.index(self.dataset_filename))
        self.assertLen(filenames, 1)
        dataset = message_helpers.load_message(filenames[0],
                                               dataset_pb2.Dataset)
        self.assertLen(dataset.reactions, 2)
        self.assertFalse(os.path.exists(dataset1_filename))
        self.assertFalse(os.path.exists(dataset2_filename))

    def test_add_dataset_with_existing_reaction_ids(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction_id = 'ord-10aed8b5dffe41fab09f5b2cc9c58ad9'
        reaction.reaction_id = reaction_id
        reaction.provenance.record_created.time.value = '2020-01-01 11 am'
        reaction.provenance.record_created.person.username = 'test'
        reaction.provenance.record_created.person.email = 'test@example.com'
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(dataset, dataset_filename)
        added, removed, filenames = self._run()
        self.assertEqual(added, {'ord-10aed8b5dffe41fab09f5b2cc9c58ad9'})
        self.assertEmpty(removed)
        self.assertLen(filenames, 2)
        self.assertFalse(os.path.exists(dataset_filename))
        filenames.pop(filenames.index(self.dataset_filename))
        self.assertLen(filenames, 1)
        dataset = message_helpers.load_message(filenames[0],
                                               dataset_pb2.Dataset)
        # Check that existing record IDs for added datasets are not overridden.
        self.assertEqual(dataset.reactions[0].reaction_id, reaction_id)
        self.assertLen(dataset.reactions[0].provenance.record_modified, 0)

    def test_modify_dataset(self):
        dataset = message_helpers.load_message(self.dataset_filename,
                                               dataset_pb2.Dataset)
        # Modify the existing reaction...
        dataset.reactions[0].inputs['methylamine'].components[
            0].moles.value = 2
        # ...and add a new reaction.
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = '2020-01-01'
        reaction.provenance.record_created.person.username = 'test'
        reaction.provenance.record_created.person.email = 'test@example.com'
        reaction.reaction_id = 'test'
        dataset.reactions.add().CopyFrom(reaction)
        message_helpers.write_message(dataset, self.dataset_filename)
        added, removed, filenames = self._run()
        self.assertEqual(added, {'test'})
        self.assertEmpty(removed)
        self.assertCountEqual([self.dataset_filename], filenames)
        # Check for preservation of dataset and record IDs.
        updated_dataset = message_helpers.load_message(self.dataset_filename,
                                                       dataset_pb2.Dataset)
        self.assertLen(updated_dataset.reactions, 2)
        self.assertEqual(dataset.dataset_id, updated_dataset.dataset_id)
        self.assertEqual(dataset.reactions[0].reaction_id,
                         updated_dataset.reactions[0].reaction_id)
        self.assertNotEmpty(updated_dataset.reactions[1].reaction_id)

    def test_modify_reaction_id(self):
        dataset = message_helpers.load_message(self.dataset_filename,
                                               dataset_pb2.Dataset)
        dataset.reactions[0].reaction_id = 'test_rename'
        message_helpers.write_message(dataset, self.dataset_filename)
        added, removed, filenames = self._run()
        self.assertEqual(added, {'test_rename'})
        self.assertEqual(removed, {'ord-10aed8b5dffe41fab09f5b2cc9c58ad9'})
        self.assertCountEqual([self.dataset_filename], filenames)

    def test_add_dataset_with_validation_errors(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='C#O')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(dataset, dataset_filename)
        with self.assertRaisesRegex(validations.ValidationError,
                                    'could not validate SMILES'):
            self._run()

    def test_add_sharded_dataset_with_validation_errors(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        dataset1 = dataset_pb2.Dataset(reactions=[reaction])
        dataset1_filename = os.path.join(self.test_subdirectory, 'test1.pbtxt')
        message_helpers.write_message(dataset1, dataset1_filename)
        reaction.inputs['ethylamine'].components[0].identifiers[
            0].value = 'C#O'
        dataset2 = dataset_pb2.Dataset(reactions=[reaction])
        dataset2_filename = os.path.join(self.test_subdirectory, 'test2.pbtxt')
        message_helpers.write_message(dataset2, dataset2_filename)
        with self.assertRaisesRegex(validations.ValidationError,
                                    'could not validate SMILES'):
            self._run()

    def test_modify_dataset_with_validation_errors(self):
        dataset = message_helpers.load_message(self.dataset_filename,
                                               dataset_pb2.Dataset)
        dataset.reactions[0].inputs['methylamine'].components[
            0].moles.value = (-2)
        message_helpers.write_message(dataset, self.dataset_filename)
        with self.assertRaisesRegex(validations.ValidationError,
                                    'must be non-negative'):
            self._run()

    def test_add_dataset_with_large_data(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = '2020-01-01'
        reaction.provenance.record_created.person.username = 'test'
        reaction.provenance.record_created.person.email = 'test@example.com'
        reaction.reaction_id = 'test'
        image = reaction.observations.add().image
        image.bytes_value = b'test data value'
        image.format = 'png'
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(dataset, dataset_filename)
        added, removed, filenames = self._run(min_size=0.0)
        self.assertEqual(added, {'test'})
        self.assertEmpty(removed)
        self.assertLen(filenames, 2)
        filenames.pop(filenames.index(self.dataset_filename))
        dataset = message_helpers.load_message(filenames[0],
                                               dataset_pb2.Dataset)
        relative_path = (
            'data/36/ord_data-'
            '36443a1839bf1160087422b7468a93c7b97dac7eea423bfac189208a15823139'
            '.png')
        expected = ('https://github.com/Open-Reaction-Database/'
                    'ord-submissions-test/tree/' + relative_path)
        self.assertEqual(dataset.reactions[0].observations[0].image.url,
                         expected)
        with open(os.path.join(self.test_subdirectory, relative_path),
                  'rb') as f:
            self.assertEqual(b'test data value', f.read())

    def test_add_dataset_with_too_large_data(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='SMILES', value='CCN')
        component.is_limiting = True
        component.moles.value = 2
        component.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        image = reaction.observations.add().image
        image.bytes_value = b'test data value'
        image.format = 'png'
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(self.test_subdirectory, 'test.pbtxt')
        message_helpers.write_message(dataset, dataset_filename)
        with self.assertRaisesRegex(ValueError, 'larger than max_size'):
            self._run(min_size=0.0, max_size=0.0)


if __name__ == '__main__':
    absltest.main()
