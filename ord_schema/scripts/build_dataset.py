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
"""Builds a Dataset from a set of Reaction protos.

TODO(kearnes): Add support for automatic sharding?
"""

import argparse
import glob

from ord_schema import message_helpers, validations
from ord_schema.logging import get_logger
from ord_schema.proto import dataset_pb2, reaction_pb2

logger = get_logger(__name__)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Build a Dataset from Reaction protos")
    parser.add_argument("--input", required=True, help="Input pattern for Reaction protos")
    parser.add_argument("--output", required=True, help="Output Dataset filename (*.pbtxt)")
    parser.add_argument("--name", required=True, help="Name for this dataset")
    parser.add_argument("--description", required=True, help="Description for this dataset")
    parser.add_argument("--no-validate", action="store_true", help="If set, do not run validations on reactions")
    return parser.parse_args(argv)


def main(args):
    filenames = glob.glob(args.input, recursive=True)
    logger.info("Found %d Reaction protos", len(filenames))
    reactions = []
    for filename in filenames:
        reactions.append(message_helpers.load_message(filename, reaction_pb2.Reaction))
    dataset = dataset_pb2.Dataset(name=args.name, description=args.description, reactions=reactions)
    if not args.no_validate:
        validations.validate_datasets({"_COMBINED": dataset})
    message_helpers.write_message(dataset, args.output)


if __name__ == "__main__":
    main(parse_args())
