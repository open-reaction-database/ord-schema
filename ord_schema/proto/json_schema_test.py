"""Tests for ord_schema.proto.json_schema."""

from absl.testing import absltest

from ord_schema.proto import json_schema
from ord_schema.proto import test_pb2


class JsonSchemaTest(absltest.TestCase):

    def test_scalar(self):
        schema = json_schema.get_schema(test_pb2.Scalar.DESCRIPTOR)
        expected = {
            'title': 'Scalar',
            'type': 'object',
            'properties': {
                'int32_value': {'type': 'number', 'multipleOf': 1.0},
                'int64_value': {'type': 'number', 'multipleOf': 1.0},
                'float_value': {'type': 'number'},
                'string_value': {'type': 'string'},
                'bytes_value': {'type': 'string'}
            },
            'additionalProperties': False,
        }
        self.assertEqual(schema, expected)

    def test_repeated_scalar(self):
        schema = json_schema.get_schema(test_pb2.RepeatedScalar.DESCRIPTOR)
        expected = {
            'title': 'RepeatedScalar',
            'type': 'object',
            'properties': {
                'value': {
                    'type': 'array',
                    'items': {'type': 'number'}
                }
            },
            'additionalProperties': False,
        }
        self.assertEqual(schema, expected)

    def test_enum(self):
        schema = json_schema.get_schema(test_pb2.Enum.DESCRIPTOR)
        expected = {
            'title': 'Enum',
            'type': 'object',
            'properties': {
                'value': {
                    'type': 'string',
                    'enum': ['UNSPECIFIED', 'FIRST', 'SECOND']
                }
            },
            'additionalProperties': False,
        }
        self.assertEqual(schema, expected)

    def test_repeated_enum(self):
        schema = json_schema.get_schema(test_pb2.RepeatedEnum.DESCRIPTOR)
        expected = {
            'title': 'RepeatedEnum',
            'type': 'object',
            'properties': {
                'value': {
                    'type': 'array',
                    'items': {
                        'type': 'string',
                        'enum': ['UNSPECIFIED', 'FIRST', 'SECOND']
                    },
                }
            },
            'additionalProperties': False,
        }
        self.assertEqual(schema, expected)

    def test_nested(self):
        schema = json_schema.get_schema(test_pb2.Nested.DESCRIPTOR)
        expected = {
            'title': 'Nested',
            'type': 'object',
            'properties': {
                'child': {
                    'title': 'Child',
                    'type': 'object',
                    'properties': {
                        'value': {'type': 'number'}
                    },
                    'additionalProperties': False,
                }
            },
            'additionalProperties': False,
        }
        self.assertEqual(schema, expected)

    def test_repeated_nested(self):
        schema = json_schema.get_schema(test_pb2.RepeatedNested.DESCRIPTOR)
        expected = {
            'title': 'RepeatedNested',
            'type': 'object',
            'properties': {
                'children': {
                    'type': 'array',
                    'items': {
                        'title': 'Child',
                        'type': 'object',
                        'properties': {
                            'value': {'type': 'number'}
                        },
                        'additionalProperties': False,
                    }
                }
            },
            'additionalProperties': False,
        }
        self.assertEqual(schema, expected)

    def test_map(self):
        # NOTE(kearnes): Proto maps are represented in JSON as arrays of
        # <map_name>Entry objects, which have 'key' and 'value' fields.
        # In other words:
        #
        # map<string, float> values = 1;
        #
        # is equivalent to
        #
        # message ValuesEntry {
        #   string key = 1;
        #   float value = 2;
        # }
        # repeated ValuesEntry = 1;
        schema = json_schema.get_schema(test_pb2.Map.DESCRIPTOR)
        expected = {
            'title': 'Map',
            'type': 'object',
            'properties': {
                'values': {
                    'type': 'array',
                    'items': {
                        'title': 'ValuesEntry',
                        'type': 'object',
                        'properties': {
                            'key': {'type': 'string'},
                            'value': {'type': 'number'}
                        },
                        'additionalProperties': False,
                    }
                }
            },
            'additionalProperties': False,
        }
        self.assertEqual(schema, expected)


if __name__ == '__main__':
    absltest.main()
