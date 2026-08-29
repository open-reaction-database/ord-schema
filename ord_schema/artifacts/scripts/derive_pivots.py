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

A quantifier over a repeated level is answered by a semi-join against a pivot -- one row
per element, carrying the ordinals that say which element it was and the element's own
fields with the repeated ones removed; see ``ord_schema.artifacts.pivot`` for the shape
and why. Building one unnests the projection, which is minutes per level over ORD, so it
belongs here rather than in whichever query asks first.

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
pattern's reach. These are errors: an --output_dir that would write over any input, an
input that cannot be read as Parquet at all, a level the schema does not have or that
is named twice, a projection that changed while the run was deriving from it or whose
count and unnest disagreed, and a run that finds no projections at all.
"""

import argparse
import functools
import pathlib

from ord_schema.artifacts import base, pivot, projection
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


def _write_counted(
    source: pathlib.Path,
    output: pathlib.Path,
    /,
    *,
    level_path: str,
    level_paths: tuple[str, ...],
    source_md5: str,
    source_dataset_id: str | None,
) -> int:
    """Writes one level's artifact, told how many elements it holds.

    Args:
        source: The projection to pivot.
        output: Where the artifact goes.
        level_path: The level to pivot.
        level_paths: Every level this run derives, counted together on first need.
        source_md5: Hash of the source dataset to stamp.
        source_dataset_id: Source dataset ID to stamp.

    Returns:
        The number of rows written.
    """
    return pivot.write_pivot(
        source,
        output,
        level_path=level_path,
        counts=_counts(source, pivot.file_identity(source), level_paths),
        source_md5=source_md5,
        source_dataset_id=source_dataset_id,
    )


@functools.cache
def _counts(
    source: pathlib.Path,
    identity: pivot.FileIdentity,
    level_paths: tuple[str, ...],
) -> pivot.Counts:
    """Returns the element count per level for one projection, read once.

    Cached because derive_tree drives one level at a time: the levels share ancestors,
    so counting them together reads the projection once where a count per level reads
    it once per level -- measured at 15x deriving every level of a two-shard tree of
    projections that record them all. The margin grows with the bytes each level
    carries, since that is what the count reads and what the extra reads repeat.

    Asked for every level being derived rather than only the one at hand, which is the
    trade: a run rebuilding one stale level of a projection reads the columns of all of
    them, against a run rebuilding many reading those columns once instead of once per
    level. Only projections that need at least one level written are counted at all,
    since nothing calls this for a source derive_tree skips.

    Args:
        source: The projection to count.
        identity: ``pivot.file_identity(source)``, which is what makes the key name the
            file rather than the path. The dataset hash will not do: it is stamped
            through from the projection's parent, so a projection rebuilt from an
            unchanged dataset carries the same one and would hit an entry counted over
            the previous bytes.
        level_paths: Every level being derived, so one read answers for all of them.

    Returns:
        The counts, carrying the identity of the file they were read from.
    """
    return pivot.count_levels(source, level_paths)


def main(args: argparse.Namespace) -> None:
    """Derives a pivot artifact per level for every projection matching the pattern.

    Args:
        args: Parsed arguments, read for ``input_pattern``, ``output_dir``, ``levels``,
            and ``force``.

    Raises:
        ValueError: If a requested level is not one the projection schema has or is
            named twice; if a projection changed while it was being derived, or its
            count and its unnest disagreed; or if the pattern matched no projections --
            which usually means it was aimed at the source tree, or at an output tree,
            rather than at the projections. Silence there would let a pipeline step
            downstream proceed as though the artifacts had been built.
    """
    levels = pivot.check_levels(args.levels)
    for level_path in levels:
        directory = pathlib.Path(args.output_dir) / level_path
        written, skipped, ignored = base.derive_tree(
            args.input_pattern,
            str(directory),
            artifact=pivot.ARTIFACT,
            write=functools.partial(
                _write_counted,
                level_path=level_path,
                level_paths=tuple(levels),
            ),
            # Reaches the ordinals, which are top-level columns: a repeated level added
            # to the schema above this one gives every level beneath it another one,
            # and the artifacts carrying the old set are derived again rather than
            # stamping as current. What it does not reach is the element struct, which
            # is one column however its fields change; that is what a version is for.
            schema=pivot.schema(pivot.LEVELS[level_path]),
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
