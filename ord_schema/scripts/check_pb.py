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
"""Compares pbtxt and pb Datasets.

Specifically, checks that a pb Dataset:
    - Can be read
    - Matches the contents of a pbtxt ground truth
"""

import difflib
import pprint

from absl import app
from absl import flags
from google.protobuf import text_format

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('pb', None, 'Path to *.pb Dataset.')
flags.DEFINE_string('pbtxt', None, 'Path to *.pbtxt Dataset.')


def main(argv):
    del argv  # Only used by app.run().
    dataset = message_helpers.load_message(FLAGS.pb, dataset_pb2.Dataset)
    pb_data = text_format.MessageToString(dataset)
    with open(FLAGS.pbtxt) as f:
        pbtxt_data = f.read()
    if pb_data != pbtxt_data:
        diff = difflib.context_diff(pb_data.splitlines(),
                                    pbtxt_data.splitlines())
        raise ValueError(f'Datasets differ:\n{pprint.pformat(list(diff))}')


if __name__ == '__main__':
    flags.mark_flags_as_required(['pb', 'pbtxt'])
    app.run(main)
