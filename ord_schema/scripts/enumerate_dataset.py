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
"""Creates a Dataset by enumerating a template with a spreadsheet.

Usage:
    dataset_templating.py --template=<str> --spreadsheet=<str> --output=<str> [--no-validate]

Options:
    --template=<str>        Path to a Reaction pbtxt file defining a template
    --spreadsheet=<str>     Path to a spreadsheet file with a header row matching template placeholders
    --output=<str>          Filename for output Dataset
    --no-validate           If set, do not validate Reaction protos
"""
import docopt

from ord_schema.logging import get_logger
from ord_schema import message_helpers
from ord_schema import templating

logger = get_logger(__name__)


def main(kwargs):
    with open(kwargs["--template"]) as f:
        template_string = f.read()
    df = templating.read_spreadsheet(kwargs["--spreadsheet"])
    logger.info(
        "generating new Dataset from %s and %s",
        kwargs["--template"],
        kwargs["--spreadsheet"],
    )
    dataset = templating.generate_dataset(template_string, df, validate=(not kwargs["--no-validate"]))
    logger.info("writing new Dataset to %s", kwargs["--output"])
    message_helpers.write_message(dataset, kwargs["--output"])


if __name__ == "__main__":
    main(docopt.docopt(__doc__))
