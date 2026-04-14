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
"""Creates a Dataset by enumerating a template with a spreadsheet."""

import argparse

from ord_schema import message_helpers, templating
from ord_schema.logging import get_logger

logger = get_logger(__name__)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Enumerate a template with a spreadsheet")
    parser.add_argument("--name", required=True, help="Dataset name")
    parser.add_argument("--description", required=True, help="Dataset description")
    parser.add_argument("--template", required=True, help="Path to a Reaction pbtxt file defining a template")
    parser.add_argument(
        "--spreadsheet",
        required=True,
        help="Path to a spreadsheet file with a header row matching template placeholders",
    )
    parser.add_argument("--output", required=True, help="Filename for output Dataset")
    parser.add_argument("--no-validate", action="store_true", help="If set, do not validate Reaction protos")
    return parser.parse_args(argv)


def main(args):
    with open(args.template) as f:
        template_string = f.read()
    df = templating.read_spreadsheet(args.spreadsheet)
    logger.info(
        "generating new Dataset from %s and %s",
        args.template,
        args.spreadsheet,
    )
    dataset = templating.generate_dataset(
        name=args.name,
        description=args.description,
        template_string=template_string,
        df=df,
        validate=not args.no_validate,
    )
    logger.info("writing new Dataset to %s", args.output)
    message_helpers.write_message(dataset, args.output)


if __name__ == "__main__":
    main(parse_args())
