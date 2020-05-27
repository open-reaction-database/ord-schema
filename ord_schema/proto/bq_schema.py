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
"""Generates a BigQuery table schema to handle Reaction protos.

See:
* https://cloud.google.com/bigquery/docs/schemas
* https://cloud.google.com/bigquery/docs/nested-repeated

Inspired by https://github.com/GoogleCloudPlatform/protoc-gen-bq-schema.
"""

import json

from absl import app
from absl import flags

from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('output', None, 'Output filename (*.json).')


def get_schema(message_descriptor):
    """Builds the table schema for a Message.

    Args:
        message_descriptor: google.protobuf.descriptor.Descriptor.

    Returns:
        List of dictionaries containing schema information.

    Raises:
        NotImplementedError: If a FieldDescriptor type is not supported.
    """
    schema = []
    for field in message_descriptor.fields:
        if field.label == field.LABEL_REPEATED:
            mode = 'REPEATED'
        else:
            mode = 'NULLABLE'  # All fields are optional.
        if field.type == field.TYPE_BOOL:
            schema.append({'name': field.name, 'type': 'BOOL', 'mode': mode})
        elif field.type == field.TYPE_FLOAT:
            schema.append({
                'name': field.name,
                'type': 'FLOAT64',
                'mode': mode
            })
        elif field.type in [field.TYPE_INT32, field.TYPE_INT64]:
            schema.append({'name': field.name, 'type': 'INT64', 'mode': mode})
        elif field.type == field.TYPE_STRING:
            schema.append({'name': field.name, 'type': 'STRING', 'mode': mode})
        elif field.type == field.TYPE_BYTES:
            schema.append({'name': field.name, 'type': 'BYTES', 'mode': mode})
        elif field.type == field.TYPE_ENUM:
            schema.append({'name': field.name, 'type': 'STRING', 'mode': mode})
        elif field.type == field.TYPE_MESSAGE:
            schema.append({
                'name': field.name,
                'type': 'RECORD',
                'mode': mode,
                'fields': get_schema(field.message_type)
            })
        else:
            raise NotImplementedError(f'unsupported type: {field.type}')
    return schema


def main(argv):
    del argv  # Only used by app.run().
    schema = get_schema(reaction_pb2.Reaction.DESCRIPTOR)
    # Add fields for the dataset ID.
    schema.append({
        'name': '_dataset_id',
        'type': 'STRING',
        'mode': 'NULLABLE'
    })
    # Add a field for the full serialized message.
    schema.append({'name': '_serialized', 'type': 'BYTES', 'mode': 'NULLABLE'})
    with open(FLAGS.output, 'w') as f:
        json.dump(schema, f, indent=2)


if __name__ == '__main__':
    flags.mark_flag_as_required('output')
    app.run(main)
