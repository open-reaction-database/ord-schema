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
"""Converts serialized protocol buffers to JSON for use in BigQuery.

The output is a JSON-lines file, where newlines are used to separate records.
See http://jsonlines.org/ for more details.
"""

import base64
import glob
import json

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern (glob).')
flags.DEFINE_string('output', None, 'Output filename (*.jsonl).')


def encode_bytes(value):
    """Encodes bytes values for BigQuery."""
    # Decode to UTF8 since JSON does not support bytes.
    return base64.b64encode(value).decode('utf-8')


def get_processed_value(field, value):
    """Returns a processed field value.

    Args:
        field: FieldDescriptor for the current field.
        value: The field value.

    Returns:
        A processed version of the value. For instance, enum values are returned
        as strings.
    """
    if field.type == field.TYPE_MESSAGE:
        return get_database_json(value)
    if field.type == field.TYPE_ENUM:
        return field.enum_type.values_by_number[value].name
    if field.type == field.TYPE_BYTES:
        return encode_bytes(value)
    return value


def get_database_json(message):
    """Generates JSON that conforms to JSON schema.

    In particular, this means that proto maps are treated as arrays of objects
    with (key, value) fields. This is critical for inclusion in structured
    databases like BigQuery.

    Args:
        message: Protocol buffer.

    Returns:
        Dict matching the JSON schema specification.
    """
    record = {}
    for field, value in message.ListFields():
        if (field.type == field.TYPE_MESSAGE
                and field.message_type.GetOptions().map_entry):
            # Convert proto maps to lists of (key, value) pairs.
            field_key = field.message_type.fields_by_name['key']
            field_value = field.message_type.fields_by_name['value']
            entries = []
            # Order is not preserved in the map, so we attempt to sort here.
            try:
                keys = sorted(value)
            except TypeError:
                keys = value.keys()
            for map_key in keys:
                entry = {
                    'key': get_processed_value(field_key, map_key),
                    'value': get_processed_value(field_value, value[map_key]),
                }
                entries.append(entry)
            record[field.name] = entries
        elif field.label == field.LABEL_REPEATED:
            entries = []
            for entry in value:
                entries.append(get_processed_value(field, entry))
            record[field.name] = entries
        else:
            record[field.name] = get_processed_value(field, value)
    return record


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input)
    logging.info('Found %d datasets', len(filenames))
    records = []
    for filename in filenames:
        dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
        for reaction in dataset.reactions:
            record_dict = get_database_json(reaction)
            record_dict['_dataset_id'] = dataset.dataset_id
            record_dict['_serialized'] = encode_bytes(
                reaction.SerializeToString())
            records.append(json.dumps(record_dict))
    with open(FLAGS.output, 'w') as f:
        for record in records:
            f.write(f'{record}\n')


if __name__ == '__main__':
    flags.mark_flags_as_required(['input', 'output'])
    app.run(main)
