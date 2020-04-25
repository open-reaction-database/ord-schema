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
flags.DEFINE_enum('input_format', 'binary',
                  [f.value for f in message_helpers.MessageFormats],
                  'Input message format.')
flags.DEFINE_boolean('write_errors', False,
                     'If True, errors will be written to <filename>.error.')
flags.DEFINE_string('output', None, 'Filename for output Dataset.')
flags.DEFINE_boolean('validate', True, 'If True, validate Reaction protos.')
flags.DEFINE_boolean('update', False, 'If True, update Reaction protos.')


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
        raise ValueError('validation encountered errors')


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
    for reaction in dataset.reactions:
        for error in validations.validate_message(
                reaction, raise_on_error=False):
            errors.append(error)
            logging.warning('Validation error for %s: %s', basename, error)
    logging.info('Validation summary for %s: %d/%d successful (%d failures)',
                 basename,
                 len(dataset.reactions) - len(errors),
                 len(dataset.reactions),
                 len(errors))
    return errors


def _update_reaction(reaction):
    """Updates a Reaction message.

    Current updates:
      * Sets record_id if not already set.

    Args:
        reaction: reaction_pb2.Reaction message.
    """
    if not reaction.provenance.record_id:
        record_id = f'ord-{uuid.uuid4().hex}'
        reaction.provenance.record_id = record_id


def main(argv):
    del argv  # Only used by app.run().
    # Setting recursive=True allows recursive matching with '**'.
    filenames = glob.glob(FLAGS.input_pattern, recursive=True)
    datasets = {}
    for filename in filenames:
        datasets[filename] = message_helpers.load_message(
            filename, dataset_pb2.Dataset, FLAGS.input_format)
    if FLAGS.validate:
        validate(datasets)
    if FLAGS.update:
        for dataset in datasets.values():
            for reaction in dataset.reactions:
                _update_reaction(reaction)
    # Combine all datasets into a single object.
    combined = None
    for dataset in datasets.values():
        if combined is None:
            combined = dataset
        else:
            combined.reactions.extend(dataset.reactions)
    if FLAGS.output:
        message_helpers.write_message(
            combined, FLAGS.output, FLAGS.input_format)


if __name__ == '__main__':
    flags.mark_flag_as_required('input_pattern')
    app.run(main)
