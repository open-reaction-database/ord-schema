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
"""Dataset processing script for database submissions.

This script is meant to be a one-stop shop for preparing submissions to the
Open Reaction Database.

By default, the script only validates the input Dataset messages. Validation
may introduce changes to the Reaction messages, such as the addition of SMILES
for compounds identified only by NAME. Users should frequently run these checks
as they are preparing a dataset for submission.

With the optional --update flag, the script also performs database-specific
updates (such as adding record IDs). These operations are meant to be run as
part of the submission process and not as part of the pre-submission validation
cycle.

Example usage:

* For normal validation-only operation:
  $ python process_dataset.py --input_pattern=my_dataset.pb
* To write the validated protos to disk:
  $ python process_dataset.py --input_pattern=my_dataset.pb --output=out.pb
* To write errors to disk:
  $ python process_dataset.py --input_pattern=my_dataset.pb --write_errors
* To process multiple Dataset protos:
  $ python process_dataset.py --input_pattern="my_dataset-*.pb"
"""

import dataclasses
import glob
import os
import subprocess
import uuid

from absl import app
from absl import flags
from absl import logging
import github
from google.protobuf import text_format

from ord_schema import data_storage
from ord_schema import message_helpers
from ord_schema import updates
from ord_schema import validations
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('root', '', 'Root of the repository.')
flags.DEFINE_string('input_pattern', None,
                    'Pattern (glob) matching input Dataset protos.')
flags.DEFINE_string('input_file', None,
                    'Filename containing Dataset proto filenames.')
flags.DEFINE_boolean('write_errors', False,
                     'If True, errors will be written to <filename>.error.')
flags.DEFINE_string('output', None, 'Filename for output Dataset.')
flags.DEFINE_boolean('validate', True, 'If True, validate Reaction protos.')
flags.DEFINE_boolean('update', False, 'If True, update Reaction protos.')
flags.DEFINE_boolean('cleanup', False, 'If True, use git to clean up.')
flags.DEFINE_float('min_size', 0.1,
                   'Minimum size (in MB) for offloading Data values.')
flags.DEFINE_float('max_size', 100.0,
                   'Maximum size (in MB) for offloading Data values.')
flags.DEFINE_string('base', None, 'Git branch to diff against.')
flags.DEFINE_integer(
    'issue', None,
    'GitHub pull request number. If provided, a comment will be added.')
flags.DEFINE_string('token', None, 'GitHub authentication token.')


@dataclasses.dataclass(eq=True, frozen=True, order=True)
class FileStatus:
    """A filename and its status in Git."""
    filename: str
    status: str
    original_filename: str

    def __post_init__(self):
        if self.status[0] not in ['A', 'D', 'M', 'R']:
            raise ValueError(f'unsupported file status: {self.status}')


def _get_inputs():
    """Gets a list of Dataset proto filenames to process.

    Returns:
        List of FileStatus objects.

    Raises:
        ValueError: If a git-diff status is not one of {'A', 'D', 'M', 'R'}.
    """
    if FLAGS.input_pattern and FLAGS.input_file:
        raise ValueError(
            'one of --input_pattern or --input_file is required, not both')
    if FLAGS.input_pattern:
        # Setting recursive=True allows recursive matching with '**'.
        filenames = glob.glob(FLAGS.input_pattern, recursive=True)
        return [FileStatus(filename, 'A', '') for filename in filenames]
    if FLAGS.input_file:
        inputs = []
        with open(FLAGS.input_file) as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) == 3:
                    status, original_filename, filename = fields
                    if not status.startswith('R'):
                        raise ValueError(
                            f'malformed status line: {line.strip()}')
                else:
                    status, filename = fields
                    if not status.startswith(('A', 'D', 'M')):
                        raise ValueError(
                            f'unsupported git-diff statue: {status}')
                    original_filename = ''
                inputs.append(FileStatus(filename, status, original_filename))
        return inputs
    raise ValueError('one of --input_pattern or --input_file is required')


def _combine_datasets(datasets):
    """Combines multiple Dataset messages into a single Dataset.

    Args:
        datasets: Dict mapping FileStatus objects to Dataset messages.

    Returns:
        Combined Dataset message.

    Raises:
        ValueError: if the combined dataset ID is malformed.
    """
    combined = None
    for file_status in sorted(datasets):
        dataset = datasets[file_status]
        if combined is None:
            combined = dataset
        else:
            combined.reactions.extend(dataset.reactions)
    if not combined.dataset_id or len(datasets) > 1:
        combined.dataset_id = f'ord_dataset-{uuid.uuid4().hex}'
    return combined


def cleanup(inputs, output_filename):
    """Removes and/or renames the input Dataset files.

    Args:
        inputs: List of FileStatus objects; the input Datasets.
        output_filename: Text filename for the output Dataset.
    """
    if len(inputs) == 1 and inputs[0].filename == output_filename:
        logging.info('editing an existing dataset; no cleanup needed')
        return  # Reuse the existing dataset ID.
    # Branch the first input file...
    args = ['git', 'mv', inputs[0].filename, output_filename]
    logging.info('Running command: %s', ' '.join(args))
    subprocess.run(args, check=True)
    # ...and remove the others.
    for file_status in inputs[1:]:
        args = ['git', 'rm', file_status.filename]
        logging.info('Running command: %s', ' '.join(args))
        subprocess.run(args, check=True)


def _get_reaction_ids(dataset):
    """Returns a list of reaction IDs in a Dataset."""
    reaction_ids = []
    for reaction in dataset.reactions:
        reaction_ids.append(reaction.reaction_id)
    return reaction_ids


def _load_base_dataset(file_status, base):
    """Loads a Dataset message from another branch."""
    if file_status.status.startswith('A'):
        return None  # Dataset only exists in the submission.
    # NOTE(kearnes): Use --no-pager to avoid a non-zero exit code.
    args = ['git', '--no-pager', 'show']
    if file_status.status.startswith('R'):
        args.append(f'{base}:{file_status.original_filename}')
    else:
        args.append(f'{base}:{file_status.filename}')
    logging.info('Running command: %s', ' '.join(args))
    dataset_pbtxt = subprocess.run(args,
                                   capture_output=True,
                                   check=True,
                                   text=True)
    return text_format.Parse(dataset_pbtxt.stdout, dataset_pb2.Dataset())


def get_change_stats(datasets, inputs, base):
    """Computes diff statistics for the submission.

    Args:
        datasets: Dict mapping filenames to Dataset messages.
        inputs: List of FileStatus objects.
        base: Git branch to diff against.

    Returns:
        added: Set of added reaction IDs.
        removed: Set of deleted reaction IDs.
    """
    old, new = set(), set()
    for file_status in inputs:
        if not file_status.status.startswith('D'):
            new.update(_get_reaction_ids(datasets[file_status.filename]))
        dataset = _load_base_dataset(file_status, base)
        if dataset is not None:
            old.update(_get_reaction_ids(dataset))
    return new - old, old - new


def _run_updates(inputs, datasets):
    """Updates the submission files.

    Args:
        inputs: List of FileStatus objects.
        datasets: Dict mapping filenames to Dataset messages.
    """
    for dataset in datasets.values():
        # Set reaction_ids, resolve names, fix cross-references, etc.
        updates.update_dataset(dataset)
        # Offload large Data values.
        data_filenames = data_storage.extract_data(dataset,
                                                   FLAGS.root,
                                                   min_size=FLAGS.min_size,
                                                   max_size=FLAGS.max_size)
        if data_filenames:
            args = ['git', 'add'] + data_filenames
            logging.info('Running command: %s', ' '.join(args))
            subprocess.run(args, check=True)
    combined = _combine_datasets(datasets)
    # Final validation (incl. IDs) to make sure we didn't break anything.
    validations.validate_datasets({'_COMBINED': combined},
                                  FLAGS.write_errors,
                                  validate_ids=True)
    if FLAGS.output:
        output_filename = FLAGS.output
    else:
        _, suffix = os.path.splitext(inputs[0].filename)
        output_filename = os.path.join(
            FLAGS.root,
            message_helpers.id_filename(f'{combined.dataset_id}{suffix}'))
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    if FLAGS.cleanup:
        cleanup(inputs, output_filename)
    logging.info('writing combined Dataset to %s', output_filename)
    message_helpers.write_message(combined, output_filename)


def run():
    """Main function that returns added/removed reaction ID sets.

    This function should be called directly by tests to get access to the
    return values. If main() returns something other than None it will break
    shell error code logic downstream.

    Returns:
        added: Set of added reaction IDs.
        removed: Set of deleted reaction IDs.
    """
    inputs = sorted(_get_inputs())
    if not inputs:
        logging.info('nothing to do')
        return set(), set()  # Nothing to do.
    datasets = {}
    for file_status in inputs:
        if file_status.status == 'D':
            continue  # Nothing to do for deleted files.
        datasets[file_status.filename] = message_helpers.load_message(
            file_status.filename, dataset_pb2.Dataset)
    if FLAGS.validate:
        # Note: this does not check if IDs are malformed.
        validations.validate_datasets(datasets, FLAGS.write_errors)
    if FLAGS.base:
        added, removed = get_change_stats(datasets, inputs, base=FLAGS.base)
        logging.info('Summary: +%d -%d reaction IDs', len(added), len(removed))
        if FLAGS.issue and FLAGS.token:
            client = github.Github(FLAGS.token)
            repo = client.get_repo(os.environ['GITHUB_REPOSITORY'])
            issue = repo.get_issue(FLAGS.issue)
            issue.create_comment(
                f'Summary: +{len(added)} -{len(removed)} reaction IDs')
    else:
        added, removed = None, None
    if FLAGS.update:
        _run_updates(inputs, datasets)
    else:
        logging.info('nothing else to do; use --update for more')
    return added, removed


def main(argv):
    del argv  # Only used by app.run().
    run()


if __name__ == '__main__':
    app.run(main)
