# Copyright 2020 The Open Reaction Database Authors
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

"""Tests for ord_schema.scripts.validate_dataset."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import validate_dataset


class ValidateDatasetTest(absltest.TestCase):

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
        dataset1 = dataset_pb2.Dataset(reactions=[reaction1])
        dataset1_filename = os.path.join(
            self.test_subdirectory, 'dataset1.pbtxt')
        message_helpers.write_message(dataset1, dataset1_filename)
        # reaction2 is empty.
        reaction2 = reaction_pb2.Reaction()
        dataset2 = dataset_pb2.Dataset(reactions=[reaction1, reaction2])
        dataset2_filename = os.path.join(
            self.test_subdirectory, 'dataset2.pbtxt')
        message_helpers.write_message(dataset2, dataset2_filename)

    def test_simple(self):
        input_pattern = os.path.join(self.test_subdirectory, 'dataset1.pbtxt')
        with flagsaver.flagsaver(input=input_pattern):
            validate_dataset.main(())

    def test_validation_errors(self):
        input_pattern = os.path.join(self.test_subdirectory, 'dataset*.pbtxt')
        with flagsaver.flagsaver(input=input_pattern):
            with self.assertRaisesRegex(
                    validations.ValidationError,
                    'Reactions should have at least 1 reaction input'):
                validate_dataset.main(())


if __name__ == '__main__':
    absltest.main()
