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
"""Creates reaction-related badges for ord-data.

Usage:
    reactions.py --root=<str> --output=<str>

Options:
    --root=<str>        ORD root
    --output=<str>      Output SVG filename
"""
import glob
import os
import requests

import docopt

from ord_schema import message_helpers
from ord_schema.logging import get_logger
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def main(kwargs):
    num_reactions = 0
    for filename in glob.glob(os.path.join(kwargs["--root"], "*", "*.pb*")):
        dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
        logger.info("%s:\t%d", filename, len(dataset.reactions))
        num_reactions += len(dataset.reactions)
    args = {
        "label": "Reactions",
        "message": num_reactions,
        "color": "informational",
    }
    response = requests.get("https://img.shields.io/static/v1", params=args)
    with open(kwargs["--output"], "w") as f:
        f.write(response.text)


if __name__ == "__main__":
    main(docopt.docopt(__doc__))
