# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Walks ord-data's history into the index of each reaction's first commit.

The index is what ``derive_arrivals --index`` reads; see
``ord_schema.artifacts.history`` for what a first commit is::

    index_arrivals.py --ord_data=ord-data --output=arrivals_index.parquet

An existing index at ``--output`` is extended from its head, so only the commits since
are read; ``--rebuild`` walks from the root instead, which an index walked under an
older definition, or along a chain the branch no longer holds, requires. Every LFS
object the walk reads has to be fetched first.
"""

import argparse
import pathlib

from ord_schema.artifacts import history
from ord_schema.logging import get_logger

logger = get_logger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ord_data", required=True, help="The ord-data checkout")
    parser.add_argument("--output", required=True, help="The index file to write")
    parser.add_argument(
        "--revision", default="main", help="The branch or commit to walk to"
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Walk from the root rather than extending the index at --output",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Walks the history and writes the index.

    Args:
        args: Parsed arguments, read for ``ord_data``, ``output``, ``revision``, and
            ``rebuild``.
    """
    previous = None
    if pathlib.Path(args.output).exists() and not args.rebuild:
        previous = history.read_index(args.output)
    index = history.walk(args.ord_data, revision=args.revision, previous=previous)
    history.write_index(index, args.output)
    logger.info(
        "%d reaction IDs across %d commits, walked to %s",
        index.table.num_rows,
        len(index.chain),
        index.head,
    )


if __name__ == "__main__":
    main(parse_args())
