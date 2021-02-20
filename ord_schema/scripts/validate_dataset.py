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
"""Validates a set of Dataset protocol buffers.

Example usage:
$ python validate.py --input="my_dataset.pbtxt"
"""

import glob

from absl import app
from absl import flags
from absl import logging
from rdkit import RDLogger

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern for Dataset protos.')


def main(argv):
    del argv  # Only used by app.run().
    filenames = sorted(glob.glob(FLAGS.input, recursive=True))
    logging.info('Found %d datasets', len(filenames))
    for filename in filenames:
        logging.info('Validating %s', filename)
        dataset = message_helpers.load_message(
            filename, dataset_pb2.Dataset)
        validations.validate_datasets({filename: dataset})


if __name__ == '__main__':
    flags.mark_flag_as_required('input')
    RDLogger.DisableLog('rdApp.*')  # Disable RDKit logging.
    app.run(main)
