"""Tests for ord_schema.proto.json_schema."""

from absl.testing import absltest

from ord_schema.proto import bq_schema
from ord_schema.proto import test_pb2


class BqSchemaTest(absltest.TestCase):

    def test_scalar(self):
        schema = bq_schema.get_schema(test_pb2.Scalar.DESCRIPTOR)
        expected = [
            {'name': 'int32_value', 'type': 'INT64', 'mode': 'NULLABLE'},
            {'name': 'int64_value', 'type': 'INT64', 'mode': 'NULLABLE'},
            {'name': 'float_value', 'type': 'FLOAT64', 'mode': 'NULLABLE'},
            {'name': 'string_value', 'type': 'STRING', 'mode': 'NULLABLE'},
            {'name': 'bytes_value', 'type': 'BYTES', 'mode': 'NULLABLE'},
        ]
        self.assertEqual(schema, expected)

    def test_repeated_scalar(self):
        schema = bq_schema.get_schema(test_pb2.RepeatedScalar.DESCRIPTOR)
        expected = [
            {'name': 'values', 'type': 'FLOAT64', 'mode': 'REPEATED'}
        ]
        self.assertEqual(schema, expected)

    def test_enum(self):
        schema = bq_schema.get_schema(test_pb2.Enum.DESCRIPTOR)
        expected = [
            {'name': 'value', 'type': 'STRING', 'mode': 'NULLABLE'}
        ]
        self.assertEqual(schema, expected)

    def test_repeated_enum(self):
        schema = bq_schema.get_schema(test_pb2.RepeatedEnum.DESCRIPTOR)
        expected = [
            {'name': 'values', 'type': 'STRING', 'mode': 'REPEATED'}
        ]
        self.assertEqual(schema, expected)

    def test_nested(self):
        schema = bq_schema.get_schema(test_pb2.Nested.DESCRIPTOR)
        expected = [
            {
                'name': 'child',
                'type': 'RECORD',
                'mode': 'NULLABLE',
                'fields': [
                    {
                        'name': 'value',
                        'type': 'FLOAT64',
                        'mode': 'NULLABLE'
                    }
                ]
            }
        ]
        self.assertEqual(schema, expected)

    def test_repeated_nested(self):
        schema = bq_schema.get_schema(test_pb2.RepeatedNested.DESCRIPTOR)
        expected = [
            {
                'name': 'children',
                'type': 'RECORD',
                'mode': 'REPEATED',
                'fields': [
                    {
                        'name': 'value',
                        'type': 'FLOAT64',
                        'mode': 'NULLABLE'
                    }
                ]
            }
        ]
        self.assertEqual(schema, expected)

    def test_map(self):
        schema = bq_schema.get_schema(test_pb2.Map.DESCRIPTOR)
        expected = [
            {
                'name': 'values',
                'type': 'RECORD',
                'mode': 'REPEATED',
                'fields': [
                    {
                        'name': 'key',
                        'type': 'STRING',
                        'mode': 'NULLABLE'
                    },
                    {
                        'name': 'value',
                        'type': 'FLOAT64',
                        'mode': 'NULLABLE'
                    }
                ]
            }
        ]
        self.assertEqual(schema, expected)

    def test_map_nested(self):
        schema = bq_schema.get_schema(test_pb2.MapNested.DESCRIPTOR)
        expected = [
            {
                'name': 'children',
                'type': 'RECORD',
                'mode': 'REPEATED',
                'fields': [
                    {
                        'name': 'key',
                        'type': 'STRING',
                        'mode': 'NULLABLE'
                    },
                    {
                        'name': 'value',
                        'type': 'RECORD',
                        'mode': 'NULLABLE',
                        'fields': [
                            {
                                'name': 'value',
                                'type': 'FLOAT64',
                                'mode': 'NULLABLE'
                            }
                        ]
                    }
                ]
            }
        ]
        self.assertEqual(schema, expected)


if __name__ == '__main__':
    absltest.main()
