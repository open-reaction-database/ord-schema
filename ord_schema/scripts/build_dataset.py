"""Builds a Dataset from a set of Reaction protos.

TODO(kearnes): Add support for automatic sharding?

Example usage:
$ python build_dataset.py \
  --input="reaction-*.pbtxt" \
  --name="My dataset" \
  --description="Experiments from DOI 10.1000/xyz123" \
  --output=my_dataset.pbtxt
"""

import glob

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern for Reaction protos.')
flags.DEFINE_string('output', None, 'Output Dataset filename (*.pbtxt).')
flags.DEFINE_string('name', None, 'Name for this dataset.')
flags.DEFINE_string('description', None, 'Description for this dataset.')
flags.DEFINE_boolean('validate', True, 'If True, run validations on Reactions.')


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input, recursive=True)
    logging.info('Found %d Reaction protos', len(filenames))
    reactions = []
    for filename in filenames:
        reactions.append(
            message_helpers.load_message(filename, reaction_pb2.Reaction))
    if not FLAGS.name:
        logging.warning('Consider setting the dataset name with --name')
    if not FLAGS.description:
        logging.warning(
            'Consider setting the dataset description with --description')
    dataset = dataset_pb2.Dataset(
        name=FLAGS.name,
        description=FLAGS.description,
        reactions=reactions)
    if FLAGS.validate:
        validations.validate_datasets({'_COMBINED': dataset})
    message_helpers.write_message(dataset, FLAGS.output)


if __name__ == '__main__':
    flags.mark_flags_as_required(['input', 'output'])
    app.run(main)
