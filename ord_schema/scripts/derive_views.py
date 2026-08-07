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

"""Derives tabular views from Parquet datasets.

Each input dataset yields one view (see ``ord_schema.views`` for what a view contains,
and what it deliberately does not). Outputs mirror the inputs' directory layout beneath
``--output_dir``, rooted at the leading components of ``--input_pattern`` that hold no
wildcard::

    derive_views.py --input_pattern='data/*/*.parquet' --output_dir=views

    data/aa/ord_dataset-<id>.parquet  ->  views/aa/ord_dataset-<id>.parquet

Views whose footer already records the current source content, library version, and
artifact version are skipped, so re-running is cheap; ``--force`` rewrites them anyway.
A match that is itself a derived artifact is ignored rather than derived from, so
``--output_dir`` may sit inside a recursive pattern's reach; an ``--output_dir`` that
would write over any source dataset is an error, as is a run that derives nothing.
"""

import argparse

from ord_schema import artifacts, views
from ord_schema.logging import silence_rdkit_logs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_pattern", required=True, help="Input pattern for Parquet datasets"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to write views beneath"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite views that are already current",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Derives a view for every dataset matching the input pattern.

    Raises:
        ValueError: If the pattern matched only derived artifacts, which means it was
            aimed at an output tree rather than a source tree. Silence there would let a
            pipeline step downstream proceed as though the views had been built.
    """
    silence_rdkit_logs()
    written, skipped, ignored = artifacts.derive_tree(
        args.input_pattern,
        args.output_dir,
        artifact=views.ARTIFACT,
        write=views.write_view,
        force=args.force,
    )
    if not written and not skipped:
        raise ValueError(
            f"no sources derived: all {ignored} matches for {args.input_pattern!r} are "
            "already derived artifacts"
        )


if __name__ == "__main__":
    main(parse_args())
