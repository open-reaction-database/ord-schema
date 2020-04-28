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

import glob
import os
import subprocess
import uuid

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input_pattern', None,
                    'Pattern (glob) matching input Dataset protos.')
flags.DEFINE_string('input_file', None,
                    'Filename containing Dataset proto filenames.')
flags.DEFINE_enum('input_format', 'binary',
                  [f.value for f in message_helpers.MessageFormats],
                  'Input message format.')
flags.DEFINE_boolean('write_errors', False,
                     'If True, errors will be written to <filename>.error.')
flags.DEFINE_string('output', None, 'Filename for output Dataset.')
flags.DEFINE_boolean('validate', True, 'If True, validate Reaction protos.')
flags.DEFINE_boolean('update', False, 'If True, update Reaction protos.')
flags.DEFINE_boolean('cleanup', False, 'If True, use git to clean up.')


def validate(datasets):
    """Runs validation for a set of datasets.

    Args:
        datasets: Dict mapping filenames to Dataset protos.
    """
    validation_errors = False
    for filename, dataset in datasets.items():
        errors = _validate_dataset(filename, dataset)
        if errors:
            validation_errors = True
            if FLAGS.write_errors:
                with open(f'{filename}.error', 'w') as f:
                    for error in errors:
                        f.write(f'{error}\n')
    # NOTE(kearnes): Run validation for all datasets before
    # exiting if there are errors.
    if validation_errors:
        raise ValueError('validation encountered errors; '
                         're-run with `--write_errors` to save errors to disk.')


def _validate_dataset(filename, dataset):
    """Validates Reaction messages in a Dataset.

    Note that validation may change the message. For example, NAME
    identifiers will be resolved to structures.

    Args:
        filename: Text filename; the dataset source.
        dataset: dataset_pb2.Dataset message.

    Returns:
        List of validation error messages.
    """
    basename = os.path.basename(filename)
    errors = []
    num_bad_reactions = 0
    for i, reaction in enumerate(dataset.reactions):
        reaction_errors = validations.validate_message(
            reaction, raise_on_error=False)
        if reaction_errors:
            num_bad_reactions += 1
        for error in reaction_errors:
            errors.append(error)
            logging.warning('Validation error for %s[%d]: %s',
                            basename, i, error)
    logging.info('Validation summary for %s: %d/%d successful (%d failures)',
                 basename,
                 len(dataset.reactions) - num_bad_reactions,
                 len(dataset.reactions),
                 num_bad_reactions)
    return errors


def _update_reaction(reaction):
    """Updates a Reaction message.

    Current updates:
      * Sets record_id if not already set.

    Args:
        reaction: reaction_pb2.Reaction message.
    """
    # TODO(kearnes): There's more complexity here regarding record_id values
    # that are set outside of the submission pipeline. It's not as simple as
    # requiring them to be empty, since a submission may include edits to
    # existing reactions.
    if not reaction.provenance.record_id:
        record_id = f'ord-{uuid.uuid4().hex}'
        reaction.provenance.record_id = record_id


def _get_filenames():
    """Gets a list of Dataset proto filenames to process.

    Returns:
        List of text filenames.
    """
    if FLAGS.input_pattern and FLAGS.input_file:
        raise ValueError(
            'one of --input_pattern or --input_file is required, not both')
    if FLAGS.input_pattern:
        # Setting recursive=True allows recursive matching with '**'.
        return glob.glob(FLAGS.input_pattern, recursive=True)
    if FLAGS.input_file:
        with open(FLAGS.input_file) as f:
            return [line.strip() for line in f]
    raise ValueError('one of --input_pattern or --input_file is required')


def main(argv):
    del argv  # Only used by app.run().
    filenames = sorted(_get_filenames())
    if not filenames:
        logging.info('nothing to do')
        return  # Nothing to do.
    datasets = {}
    for filename in filenames:
        datasets[filename] = message_helpers.load_message(
            filename, dataset_pb2.Dataset, FLAGS.input_format)
    if FLAGS.validate:
        validate(datasets)
    if not FLAGS.update:
        logging.info('nothing else to do; use --update for more')
        return  # Nothing else to do.
    for dataset in datasets.values():
        for reaction in dataset.reactions:
            _update_reaction(reaction)
    # Combine all datasets into a single object.
    combined = None
    for key in sorted(datasets):
        dataset = datasets[key]
        if combined is None:
            combined = dataset
        else:
            combined.reactions.extend(dataset.reactions)
    if FLAGS.output:
        output_filename = FLAGS.output
    else:
        suffix = message_helpers.get_suffix(FLAGS.input_format)
        dataset_id = uuid.uuid4().hex
        output_filename = os.path.join(
            'data', dataset_id[:2], f'{dataset_id}{suffix}')
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    if FLAGS.cleanup:
        # Branch the first input file...
        subprocess.run(['git', 'mv', filenames[0], output_filename], check=True)
        # ...and remove the others.
        for filename in filenames[1:]:
            subprocess.run(['git', 'rm', filename], check=True)
    logging.info('writing combined Dataset to %s', output_filename)
    message_helpers.write_message(
        combined, output_filename, FLAGS.input_format)


if __name__ == '__main__':
    app.run(main)
