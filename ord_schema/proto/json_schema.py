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
"""Generates a JSON schema to handle Reaction protos.

See:
* https://json-schema.org/
* https://json-schema.org/understanding-json-schema/index.html

Inspired by https://github.com/GoogleCloudPlatform/protoc-gen-bq-schema.
"""

import json

from absl import app
from absl import flags

from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('output', None, 'Output filename (*.json).')


def get_schema(message_descriptor):
    """Builds the JSON schema for a Message.

    Args:
        message_descriptor: google.protobuf.descriptor.Descriptor.

    Returns:
        Dict of dicts; a JSON schema.

    Raises:
        NotImplementedError: If a FieldDescriptor type is not supported.
    """
    schema = {
        # NOTE(kearnes): Setting additionalProperties to False here ensures
        # that only the listed properties are allowed. See
        # https://json-schema.org/understanding-json-schema/reference/object.html#properties.
        'additionalProperties': False,
        'type': 'object',
        'properties': {},
        'title': message_descriptor.name,
    }
    for field in message_descriptor.fields:
        field_properties = {}
        if field.type == field.TYPE_BOOL:
            field_properties.update({'type': 'boolean'})
        elif field.type == field.TYPE_FLOAT:
            field_properties.update({'type': 'number'})
        elif field.type in [field.TYPE_INT32, field.TYPE_INT64]:
            # NOTE(kearnes): We avoid 'integer' since its implementation differs
            # between validators; see
            # https://json-schema.org/understanding-json-schema/reference/numeric.html#integer.
            field_properties.update({'type': 'number', 'multipleOf': 1.0})
        elif field.type in [field.TYPE_STRING, field.TYPE_BYTES]:
            field_properties.update({'type': 'string'})
        elif field.type == field.TYPE_ENUM:
            enum_descriptor = field.enum_type
            enum_names = enum_descriptor.values_by_name.keys()
            field_properties.update({'type': 'string', 'enum': enum_names})
        elif field.type == field.TYPE_MESSAGE:
            field_properties.update(get_schema(field.message_type))
        else:
            raise NotImplementedError(f'unsupported type: {field.type}')
        if field.label == field.LABEL_REPEATED:
            schema['properties'][field.name] = {
                'type': 'array',
                'items': field_properties,
            }
        else:
            schema['properties'][field.name] = field_properties
    return schema


def main(argv):
    del argv  # Only used by app.run().
    schema = get_schema(reaction_pb2.Reaction.DESCRIPTOR)
    with open(FLAGS.output, 'w') as f:
        json.dump(schema, f, indent=2)


if __name__ == '__main__':
    flags.mark_flag_as_required('output')
    app.run(main)
