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

Usage:
    validate.py --input=<str> [--filter=<str>]

Options:
    --input=<str>       Input pattern for Dataset protos
    --filter=<str>      Regex filename filter
"""
import glob
import re
from collections.abc import Iterable

import docopt
from rdkit import RDLogger

from ord_schema.logging import get_logger
from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def filter_filenames(filenames: Iterable[str], pattern: str) -> list[str]:
    """Filters filenames according to a regex pattern."""
    filtered_filenames = []
    for filename in filenames:
        if re.search(pattern, filename):
            filtered_filenames.append(filename)
    return filtered_filenames


def main(kwargs):
    filenames = sorted(glob.glob(kwargs["--input"], recursive=True))
    logger.info("Found %d datasets", len(filenames))
    if kwargs["--filter"]:
        filenames = filter_filenames(filenames, kwargs["--filter"])
        logger.info("Filtered to %d datasets", len(filenames))
    for filename in filenames:
        logger.info("Validating %s", filename)
        dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
        validations.validate_datasets({filename: dataset})


if __name__ == "__main__":
    RDLogger.DisableLog("rdApp.*")  # Disable RDKit logging.
    main(docopt.docopt(__doc__))
