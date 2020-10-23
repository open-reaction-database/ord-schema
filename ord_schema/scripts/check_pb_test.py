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
"""Tests for ord_schema.scripts.check_pb."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.scripts import check_pb


class CheckPbTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        self.pb_filename = os.path.join(self.test_subdirectory, 'test.pb')
        self.pbtxt_filename = os.path.join(self.test_subdirectory, 'test.pbtxt')

    def _run(self):
        with flagsaver.flagsaver(pb=self.pb_filename,
                                 pbtxt=self.pbtxt_filename):
            check_pb.main(())

    def test_main_pass(self):
        dataset = dataset_pb2.Dataset()
        reaction = dataset.reactions.add()
        component = reaction.inputs['test'].components.add()
        component.identifiers.add(value='c1ccccc1', type='SMILES')
        message_helpers.write_message(dataset, self.pb_filename)
        message_helpers.write_message(dataset, self.pbtxt_filename)
        self._run()

    def test_main_fail(self):
        dataset = dataset_pb2.Dataset()
        reaction = dataset.reactions.add()
        component = reaction.inputs['test'].components.add()
        component.identifiers.add(value='c1ccccc1', type='SMILES')
        message_helpers.write_message(dataset, self.pb_filename)
        component.identifiers.add(value='benzene', type='NAME')
        message_helpers.write_message(dataset, self.pbtxt_filename)
        with self.assertRaisesRegex(ValueError, 'Datasets differ'):
            self._run()


if __name__ == '__main__':
    absltest.main()
