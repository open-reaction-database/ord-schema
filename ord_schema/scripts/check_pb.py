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

import argparse
import difflib
import pprint

from google.protobuf import text_format

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Compare pbtxt and pb Datasets")
    parser.add_argument("--pb", required=True, help="Path to *.pb Dataset")
    parser.add_argument("--pbtxt", required=True, help="Path to *.pbtxt Dataset")
    return parser.parse_args(argv)


def main(args):
    dataset = message_helpers.load_message(args.pb, dataset_pb2.Dataset)
    pb_data = text_format.MessageToString(dataset)
    with open(args.pbtxt) as f:
        pbtxt_data = f.read()
    if pb_data != pbtxt_data:
        diff = difflib.context_diff(pb_data.splitlines(), pbtxt_data.splitlines())
        raise ValueError(f"Datasets differ:\n{pprint.pformat(list(diff))}")


if __name__ == "__main__":
    main(parse_args())
