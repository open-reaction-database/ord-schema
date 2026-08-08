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

"""Derives tabular views from Parquet projections.

A view is a narrowing of the projection, so this reads projections rather than source
datasets: run ``derive_projection.py`` first (see ``ord_schema.views`` for what a view
contains, and what it deliberately does not). Each input projection yields one view, and
outputs mirror the inputs' directory layout beneath ``--output_dir``, rooted at the
leading components of ``--input_pattern`` that hold no wildcard::

    derive_views.py --input_pattern='projections/*/*.parquet' --output_dir=views

    projections/aa/ord_dataset-<id>.parquet  ->  views/aa/ord_dataset-<id>.parquet

A view records the hash of the *source dataset* its projection was built from, not of
the projection, so one comparison answers "is this current for that dataset?". Views
already recording the current source content, library version, and artifact version are
skipped, so re-running is cheap; --force rewrites them anyway.

A match that is not a projection -- a source dataset, or a view from an earlier run --
is ignored rather than derived from, so --output_dir may sit inside a recursive
pattern's reach. These are errors: an --output_dir that would write over any input, a
match that cannot be read as Parquet at all, and a run that derives nothing.
"""

import argparse

from ord_schema import artifacts, projection, views
from ord_schema.logging import silence_rdkit_logs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input_pattern", required=True, help="Input pattern for Parquet projections"
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
    """Derives a view for every projection matching the input pattern.

    Args:
        args: Parsed command-line arguments.

    Raises:
        ValueError: If the pattern matched no projections, which usually means it was
            aimed at the source tree, or at an output tree, rather than at the
            projections. Silence there would let a pipeline step downstream proceed as
            though the views had been built.
    """
    silence_rdkit_logs()
    written, skipped, ignored = artifacts.derive_tree(
        args.input_pattern,
        args.output_dir,
        artifact=views.ARTIFACT,
        write=views.write_view,
        force=args.force,
        parent_artifact=projection.ARTIFACT,
    )
    if not written and not skipped:
        raise ValueError(
            f"no views derived: none of the {ignored} matches for "
            f"{args.input_pattern!r} are projections"
        )


if __name__ == "__main__":
    main(parse_args())
