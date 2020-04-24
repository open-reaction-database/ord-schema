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

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string(
    'input_file', None,
    'Input file containing the git status of changed files. '
    'This should be the output of `git diff --name-status <base_ref>`.')
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
        filename, reaction_pb2.Reaction, input_format='binary')
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


def get_diff():
    """Fetches a list of added or changed files.

    Returns:
        Dict mapping text filenames to status indicators; e.g. 'A' for 'added'
        or 'M' for 'modified'.
    """
    with open(FLAGS.input_file) as f:
        input_file = f.readlines()
    statuses = {}
    for line in input_file:
        # NOTE(kearnes): This doesn't work for renames; expand?
        status, filename = line.strip().split()
        if status == 'D':
            continue
        if not filename.endswith('.pb'):
            continue
        statuses[filename] = status
    return statuses


def main(argv):
    del argv  # Only used by app.run().
    statuses = get_diff()
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
                                f'{record_id}.pb')
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        if FLAGS.cleanup and original_filename != filename:
            try:
                # Use git to branch/move the original file.
                subprocess.check_output(
                    ['git', 'mv', original_filename, filename],
                    stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as error:
                # The only acceptable error is that we are not in a git repo.
                if b'outside repository' not in error.output:
                    raise error
                os.unlink(original_filename)
        with open(filename, 'wb') as f:
            f.write(message.SerializeToString())


if __name__ == '__main__':
    flags.mark_flag_as_required('output_dir')
    app.run(main)
