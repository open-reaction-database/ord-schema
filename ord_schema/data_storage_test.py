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
"""Tests for ord_schema.data_storage."""

import os
import tempfile

from absl import flags
from absl.testing import absltest

from ord_schema import data_storage
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2


class WriteDataTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)

    def test_string_value(self):
        message = reaction_pb2.Data(string_value='test value')
        filename, data_size = data_storage.write_data(message,
                                                      self.test_subdirectory)
        expected = os.path.join(
            self.test_subdirectory, 'ord_data-'
            '47d1d8273710fd6f6a5995fac1a0983fe0e8828c288e35e80450ddc5c4412def'
            '.txt')
        self.assertEqual(filename, expected)
        self.assertAlmostEqual(data_size, 0.000043)
        # NOTE(kearnes): Open with 'r' to get the decoded string.
        with open(filename, 'r') as f:
            self.assertEqual(message.string_value, f.read())

    def test_bytes_value(self):
        message = reaction_pb2.Data(bytes_value=b'test value')
        filename, _ = data_storage.write_data(message, self.test_subdirectory)
        expected = os.path.join(
            self.test_subdirectory, 'ord_data-'
            '47d1d8273710fd6f6a5995fac1a0983fe0e8828c288e35e80450ddc5c4412def'
            '.txt')
        self.assertEqual(filename, expected)
        with open(filename, 'rb') as f:
            self.assertEqual(message.bytes_value, f.read())

    def test_url_value(self):
        message = reaction_pb2.Data(url='test value')
        filename, _ = data_storage.write_data(message, self.test_subdirectory)
        self.assertIsNone(filename)

    def test_missing_value(self):
        message = reaction_pb2.Data()
        with self.assertRaisesRegex(ValueError, 'no value to write'):
            data_storage.write_data(message, self.test_subdirectory)

    def test_min_max_size(self):
        message = reaction_pb2.Data(string_value='test_value')
        with self.assertRaisesRegex(ValueError, 'must be less than or equal'):
            data_storage.write_data(message,
                                    self.test_subdirectory,
                                    min_size=2.0,
                                    max_size=1.0)

    def test_min_size(self):
        message = reaction_pb2.Data(string_value='test_value')
        filename, _ = data_storage.write_data(message,
                                              self.test_subdirectory,
                                              min_size=1.0)
        self.assertIsNone(filename)

    def test_max_size(self):
        message = reaction_pb2.Data(string_value='test value')
        with self.assertRaisesRegex(ValueError, 'larger than max_size'):
            data_storage.write_data(message,
                                    self.test_subdirectory,
                                    max_size=1e-6)


class ExtractDataTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)

    def test_find_data_messages(self):
        message = reaction_pb2.Reaction()
        self.assertEmpty(
            message_helpers.find_submessages(message, reaction_pb2.Data))
        message = reaction_pb2.ReactionObservation()
        message.image.string_value = 'not an image'
        self.assertLen(
            message_helpers.find_submessages(message, reaction_pb2.Data), 1)
        message = reaction_pb2.ReactionSetup()
        message.automation_code['test1'].string_value = 'test data 1'
        message.automation_code['test2'].bytes_value = b'test data 2'
        self.assertLen(
            message_helpers.find_submessages(message, reaction_pb2.Data), 2)
        message = reaction_pb2.Reaction()
        message.observations.add().image.string_value = 'not an image'
        message.setup.automation_code['test1'].string_value = 'test data 1'
        message.setup.automation_code['test2'].bytes_value = b'test data 2'
        self.assertLen(
            message_helpers.find_submessages(message, reaction_pb2.Data), 3)

    def test_extract_data(self):
        message = reaction_pb2.ReactionObservation()
        message.image.string_value = 'not an image'
        data_storage.extract_data(message, root=self.test_subdirectory)
        relative_path = (
            'data/54/ord_data-'
            '5464533c9647b67eb320c40ccc5959537c09102ae75388f6a7675b433e745c9d'
            '.txt')
        expected = ('https://github.com/Open-Reaction-Database/'
                    'ord-submissions-test/tree/' + relative_path)
        self.assertEqual(message.image.url, expected)
        with open(os.path.join(self.test_subdirectory, relative_path)) as f:
            self.assertEqual(f.read(), 'not an image')


if __name__ == '__main__':
    absltest.main()
