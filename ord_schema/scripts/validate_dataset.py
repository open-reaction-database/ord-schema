"""Validates a set of Dataset protocol buffers.

Example usage:
$ python validate.py --input="my_dataset.pbtxt"
"""

import glob

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern for Dataset protos.')


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input, recursive=True)
    logging.info('Found %d datasets', len(filenames))
    datasets = {}
    for filename in filenames:
        logging.info('Validating %s', filename)
        datasets[filename] = message_helpers.load_message(
            filename, dataset_pb2.Dataset)
    validations.validate_datasets(datasets)


if __name__ == '__main__':
    flags.mark_flag_as_required('input')
    app.run(main)
