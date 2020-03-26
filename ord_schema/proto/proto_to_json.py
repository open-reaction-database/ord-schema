"""Converts serialized protocol buffers to JSON.

The output is a JSON-lines file, where newlines are used to separate records.
See http://jsonlines.org/ for more details.
"""

import glob
import json

from absl import app
from absl import flags
from absl import logging

from google.protobuf import json_format
from google.protobuf.pyext import _message
from ord_schema.proto import ord_schema_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern (glob).')
flags.DEFINE_string('output', None, 'Output filename (*.jsonl).')
flags.DEFINE_boolean('database', False,
                     'If True, maps will be treated as repeated values. '
                     'This allows the JSON to be used as input to BigQuery.')


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
    elif field.type == field.TYPE_ENUM:
        return field.enum_type.values_by_number[value].name
    else:
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
        if isinstance(value, (
               _message.ScalarMapContainer, _message.MessageMapContainer)):
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
            for v in value:
                entries.append(get_processed_value(field, v))
            record[field.name] = entries
        else:
            record[field.name] = get_processed_value(field, value)
    return record


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input)
    logging.info('Found %d records', len(filenames))
    records = []
    for filename in filenames:
        with open(filename, 'rb') as f:
            reaction = ord_schema_pb2.Reaction.FromString(f.read())
        if FLAGS.database:
            record_dict = get_database_json(reaction)
            record = json.dumps(record_dict)
        else:
            record = json_format.MessageToJson(
                reaction,
                preserving_proto_field_name=True,
                indent=None)
        records.append(record)
    with open(FLAGS.output, 'w') as f:
        for record in records:
            f.write(f'{record}\n')


if __name__ == '__main__':
    flags.mark_flags_as_required(['input', 'output'])
    app.run(main)
