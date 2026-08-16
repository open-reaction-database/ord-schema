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

r"""Derives pivot artifacts from projections, one per repeated level.

A quantifier over a repeated level is answered by a semi-join against a pivot -- one
row per element, carrying the ordinals that say which element it was and the element's
own fields with the repeated ones removed; see ``ord_schema.agent.pivot`` for the shape
and why. Building one unnests the projection, which is minutes per level over ORD, so
it belongs here rather than in whichever query asks first.

Each level gets its own subdirectory, and each projection its own file within it::

    derive_pivots.py --input_pattern='projections/*/*.parquet' --output_dir=pivots

    projections/aa/ord_dataset-<id>.parquet
        -> pivots/outcomes.products/aa/ord_dataset-<id>.parquet
        -> pivots/inputs.components/aa/ord_dataset-<id>.parquet

That is the layout ``Corpus(..., pivots_dir=...)`` reads. Levels are chosen with
--levels; the default is every level the schema has, which is far more than a corpus is
usually asked about, so naming the ones a deployment answers questions over is the
cheaper run.

Artifacts whose footer already records the current source content, library version, and
artifact version are skipped, so re-running is cheap; --force rewrites them anyway.

A match that is not a projection -- a source dataset, or an artifact from an earlier
run -- is ignored rather than derived from, so --output_dir may sit inside a recursive
pattern's reach. These are errors: an --output_dir that would write over any input, a
match that cannot be read as Parquet at all, a level the schema does not have, and a
run that finds no projections at all.
"""

import argparse
import functools
import pathlib

from ord_schema import artifacts, projection
from ord_schema.agent import pivot
from ord_schema.logging import get_logger

logger = get_logger(__name__)


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
        "--levels",
        nargs="+",
        default=sorted(pivot.LEVELS),
        help="Repeated levels to pivot; every level in the schema by default",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite artifacts that are already current",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Derives a pivot artifact per level for every projection matching the pattern.

    Args:
        args: Parsed arguments, read for ``input_pattern``, ``output_dir``, ``levels``,
            and ``force``.

    Raises:
        ValueError: If a requested level is not one the projection schema has, or if
            the pattern matched no projections -- which usually means it was aimed at
            the source tree, or at an output tree, rather than at the projections.
            Silence there would let a pipeline step downstream proceed as though the
            artifacts had been built.
    """
    unknown = sorted(set(args.levels) - set(pivot.LEVELS))
    if unknown:
        raise ValueError(
            f"not repeated levels of the projection schema: {unknown}; "
            f"known levels are {sorted(pivot.LEVELS)}"
        )
    for level_path in args.levels:
        written, skipped, ignored = artifacts.derive_tree(
            args.input_pattern,
            str(pathlib.Path(args.output_dir) / level_path),
            artifact=pivot.ARTIFACT,
            write=functools.partial(pivot.write_pivot, level_path=level_path),
            force=args.force,
            parent_artifact=projection.ARTIFACT,
        )
        if not written and not skipped:
            raise ValueError(
                f"no pivots derived for {level_path}: none of the {ignored} matches "
                f"for {args.input_pattern!r} are projections"
            )
        logger.info("%s: %d written, %d already current", level_path, written, skipped)


if __name__ == "__main__":
    main(parse_args())
