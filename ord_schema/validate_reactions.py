"""Runs validations on Reaction messages.

Reaction protos are expected to be stored as one message per file.
"""

import glob
import os
import sys

from absl import app
from absl import flags
from absl import logging
from google import protobuf
from google.protobuf import json_format
from google.protobuf import text_format

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input_pattern', None,
                    'Pattern (glob) matching input records.')
flags.DEFINE_enum('input_format', 'pbtxt', ['pbtxt', 'json'],
                  'Input record format.')
flags.DEFINE_string('output', None, 'Output file for validation errors.')


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input_pattern)
    errors = []
    for filename in sorted(filenames):
        basename = os.path.basename(filename)
        try:
            reaction = message_helpers.load_message(
                filename, reaction_pb2.Reaction, FLAGS.input_format)
        except (json_format.ParseError,
                protobuf.message.DecodeError,
                text_format.ParseError) as error:
            errors.append((basename, error))
            logging.warning('Parsing error for %s: %s', basename, error)
            continue
        for error in validations.validate_message(
                reaction, raise_on_error=False):
            errors.append((basename, error))
            logging.warning('Validation error for %s: %s', basename, error)
    logging.info('Validation summary: %d/%d successful (%d failures)',
                 len(filenames) - len(errors), len(filenames), len(errors))
    if errors and FLAGS.output:
        with open(FLAGS.output, 'w') as f:
            for filename, error in errors:
                f.write(f'{filename}: {error}\n')
    if errors:
        # Tell the shell that validation failed; see
        # https://docs.python.org/3/library/sys.html#sys.exit.
        sys.exit('validation encountered errors')


if __name__ == '__main__':
    flags.mark_flags_as_required(['input_pattern'])
    app.run(main)
