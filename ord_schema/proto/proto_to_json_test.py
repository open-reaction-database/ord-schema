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
"""Tests for ord_schema.proto.proto_to_json."""

import base64

from absl.testing import absltest

from ord_schema.proto import proto_to_json
from ord_schema.proto import test_pb2


class GetDatabaseJsonTest(absltest.TestCase):
    def test_scalar(self):
        message = test_pb2.Scalar(int32_value=1,
                                  int64_value=2,
                                  float_value=3.4,
                                  string_value='five',
                                  bytes_value=b'six')
        record = proto_to_json.get_database_json(message)
        # NOTE(kearnes): Using individual tests here to avoid introducing a
        # dependency on pandas for assert_series_equal.
        self.assertLen(record, 5)
        self.assertEqual(record['int32_value'], 1)
        self.assertEqual(record['int64_value'], 2)
        self.assertAlmostEqual(record['float_value'], 3.4, places=3)
        self.assertEqual(record['string_value'], 'five')
        # Note that bytes values are converted to base64 and then decoded.
        self.assertEqual(record['bytes_value'], 'c2l4')
        self.assertEqual(
            base64.b64decode(record['bytes_value'].encode('utf-8')), b'six')

    def test_repeated_scalar(self):
        message = test_pb2.RepeatedScalar(values=[1.1, 2.2, 3.3])
        record = proto_to_json.get_database_json(message)
        # NOTE(kearnes): Using individual tests here to avoid introducing a
        # dependency on numpy for assert_almost_equal.
        self.assertLen(record, 1)
        self.assertLen(record['values'], 3)
        self.assertAlmostEqual(record['values'][0], 1.1, places=3)
        self.assertAlmostEqual(record['values'][1], 2.2, places=3)
        self.assertAlmostEqual(record['values'][2], 3.3, places=3)

    def test_enum(self):
        message = test_pb2.Enum(value='FIRST')
        record = proto_to_json.get_database_json(message)
        expected = {'value': 'FIRST'}
        self.assertEqual(record, expected)

    def test_repeated_enum(self):
        message = test_pb2.RepeatedEnum(values=['FIRST', 'FIRST', 'SECOND'])
        record = proto_to_json.get_database_json(message)
        expected = {'values': ['FIRST', 'FIRST', 'SECOND']}
        self.assertEqual(record, expected)

    def test_nested(self):
        message = test_pb2.Nested()
        message.child.value = 1.2
        record = proto_to_json.get_database_json(message)
        self.assertLen(record, 1)
        self.assertLen(record['child'], 1)
        self.assertAlmostEqual(record['child']['value'], 1.2, places=3)

    def test_repeated_nested(self):
        message = test_pb2.RepeatedNested()
        message.children.add().value = 1.2
        message.children.add().value = 3.4
        record = proto_to_json.get_database_json(message)
        self.assertLen(record, 1)
        self.assertLen(record['children'], 2)
        self.assertAlmostEqual(record['children'][0]['value'], 1.2, places=3)
        self.assertAlmostEqual(record['children'][1]['value'], 3.4, places=3)

    def test_map(self):
        message = test_pb2.Map()
        message.values['one'] = 1.2
        message.values['two'] = 3.4
        record = proto_to_json.get_database_json(message)
        self.assertLen(record, 1)
        self.assertLen(record['values'], 2)
        self.assertEqual(record['values'][0]['key'], 'one')
        self.assertAlmostEqual(record['values'][0]['value'], 1.2, places=3)
        self.assertEqual(record['values'][1]['key'], 'two')
        self.assertAlmostEqual(record['values'][1]['value'], 3.4, places=3)

    def test_nested_map(self):
        message = test_pb2.MapNested()
        message.children['one'].value = 1.2
        message.children['two'].value = 3.4
        record = proto_to_json.get_database_json(message)
        self.assertLen(record, 1)
        self.assertLen(record['children'], 2)
        self.assertEqual(record['children'][0]['key'], 'one')
        self.assertLen(record['children'][0]['value'], 1)
        self.assertAlmostEqual(record['children'][0]['value']['value'],
                               1.2,
                               places=3)
        self.assertEqual(record['children'][1]['key'], 'two')
        self.assertLen(record['children'][1]['value'], 1)
        self.assertAlmostEqual(record['children'][1]['value']['value'],
                               3.4,
                               places=3)


if __name__ == '__main__':
    absltest.main()
