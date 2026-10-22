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

r"""Derives each projection's arrivals artifact from the index of first commits.

Each projection yields one artifact: its reactions' first commits to ord-data, looked up
in the index ``index_arrivals`` writes; see ``ord_schema.artifacts.arrivals``.
Outputs mirror the inputs' directory layout beneath ``--output_dir``::

    derive_arrivals.py --input_pattern='projections/*/*.parquet' \
        --output_dir=arrivals --index=arrivals_index.parquet

    projections/aa/ord_dataset-<id>.parquet -> arrivals/aa/ord_dataset-<id>.parquet

An artifact whose stamps record the current source content, versions, and a head on the
index's chain is skipped, so re-running after the index is extended rewrites nothing;
--force rewrites them anyway. A projection holding a reaction the index lacks is an
error, as is a run that finds no projections.
"""

import argparse
import functools

from ord_schema.artifacts import arrivals, base, history, projection


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input_pattern", required=True, help="Input pattern for Parquet projections"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to write artifacts beneath"
    )
    parser.add_argument(
        "--index",
        required=True,
        help="The index of first commits, from index_arrivals",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite artifacts that are already current",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Derives an arrivals artifact for every projection matching the pattern.

    Args:
        args: Parsed arguments, read for ``input_pattern``, ``output_dir``, ``index``,
            and ``force``.

    Raises:
        ValueError: If the pattern matched no projections, which usually means it was
            aimed at the source tree or an output tree, or if a projection holds a
            reaction the index lacks.
    """
    index = history.read_index(args.index)
    written, skipped, ignored = base.derive_tree(
        args.input_pattern,
        args.output_dir,
        artifact=arrivals.ARTIFACT,
        write=functools.partial(arrivals.write_arrivals, index=index),
        schema=arrivals.SCHEMA,
        force=args.force,
        parent_artifact=projection.ARTIFACT,
        current=index.covers,
    )
    if not written and not skipped:
        raise ValueError(
            f"no arrivals artifacts derived: none of the {ignored} matches for "
            f"{args.input_pattern!r} are projections"
        )


if __name__ == "__main__":
    main(parse_args())
