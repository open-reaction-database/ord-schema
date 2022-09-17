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

Usage:
    build_dataset.py --input=<str> --output=<str> --name=<str> --description=<str> [--no-validate]

Options:
    --input=<str>           Input pattern for Reaction protos
    --output=<str>          Output Dataset filename (*.pbtxt)
    --name=<str>            Name for this dataset
    --description=<str>     Description for this dataset
    --no-validate           If set, do run validations on reactions
"""
import glob

import docopt

from ord_schema.logging import get_logger
from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)


def main(kwargs):
    filenames = glob.glob(kwargs["--input"], recursive=True)
    logger.info("Found %d Reaction protos", len(filenames))
    reactions = []
    for filename in filenames:
        reactions.append(message_helpers.load_message(filename, reaction_pb2.Reaction))
    dataset = dataset_pb2.Dataset(name=kwargs["--name"], description=kwargs["--description"], reactions=reactions)
    if not kwargs["--no-validate"]:
        validations.validate_datasets({"_COMBINED": dataset})
    message_helpers.write_message(dataset, kwargs["--output"])


if __name__ == "__main__":
    main(docopt.docopt(__doc__))
