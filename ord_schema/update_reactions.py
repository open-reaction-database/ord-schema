"""Updates Reaction protos prior to submission.

Current updates:
  * Add record IDs.

Planned updates:
  * Canonicalize structure identifiers.
  * Replace large bytes data with URLs.
"""

import os
import subprocess
import uuid

from absl import app
from absl import flags
from google.protobuf import text_format

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string(
    'input_file', None,
    'Input file containing the git status of changed files. '
    'This should be the output of `git diff --name-status <head_ref>`.')
flags.DEFINE_string('root_dir', None, 'Root directory for output records.')
flags.DEFINE_boolean('cleanup', True, 'If True, remove the original protos.')


def add_record_id(message):
    """Adds a record ID to the given Reaction message.

    The message is updated in-place with the record ID.

    Args:
        message: reaction_pb2.Reaction proto.

    Returns:
        Text record ID.

    Raises:
        ValueError: if record_id is already set.
    """
    if message.provenance.record_id:
        raise ValueError(
            f'record_id is already set: {message.provenance.record_id}')
    record_id = f'ord-{uuid.uuid4().hex}'
    message.provenance.record_id = record_id
    return record_id


def process(filename, status):
    """Processes a single Reaction proto.

    Args:
        filename: Text filename containing a serialized Reaction proto.
        status: Text git status, e.g. 'A' for 'added' or 'M' for 'modified'.

    Returns:
        Updated Reaction proto.
    """
    message = message_helpers.load_message(
        filename, reaction_pb2.Reaction, input_format='pbtxt')
    #########################
    # BEGIN MESSAGE UPDATES #
    #########################
    if status == 'A':
        add_record_id(message)
    #######################
    # END MESSAGE UPDATES #
    #######################
    # TODO(kearnes): Run validations on the updated message?
    return message


def main(argv):
    del argv  # Only used by app.run().
    statuses = {}
    with open(FLAGS.input_file) as f:
        for line in f:
            # NOTE(kearnes): This doesn't work for renames; expand?
            status, filename = line.strip().split()
            if status == 'D':
                continue
            if not filename.endswith('.pbtxt'):
                continue
            statuses[filename] = status
    records = {}
    for filename, status in statuses.items():
        records[filename] = process(filename, status)
    # Write updated protos.
    for original_filename, message in records.items():
        record_id = message.provenance.record_id
        if not record_id.startswith('ord-'):
            raise ValueError(f'malformed record_id: {record_id}')
        # Paths are sharded into subdirectories as "${ROOT}/ab/abcd.pbtxt".
        filename = os.path.join(FLAGS.root_dir,
                                record_id[len('ord-'):len('ord-') + 2],
                                f'{record_id}.pbtxt')
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        if FLAGS.cleanup:
            try:
                # Use git to branch/move the original file.
                subprocess.check_output(
                    ['git', 'mv', original_filename, filename],
                    stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as error:
                # The only acceptable error is that we are not in a git repo.
                if b'outside repository' not in error.output:
                    raise error
                os.rename(original_filename, filename)
        with open(filename, 'w') as f:
            f.write(text_format.MessageToString(message))


if __name__ == '__main__':
    flags.mark_flag_as_required(['input_file', 'output_dir'])
    app.run(main)
