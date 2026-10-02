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
count and unnest disagreed, a run that finds no projections at all, a level whose
artifacts cannot be found once it has derived them, and a level the tree already holds
that this run did not name and whose artifacts the current schema has outgrown -- a
subset run is the cheaper run, but the levels it leaves behind still have to match the
projections a corpus reads them beside.

Every run reports, per level, how many of the artifacts on disk hold no rows, and warns
when all of them do. That is ordinary for a level nothing in the corpus records, and it
is also what a wrong count looks like, which nothing downstream can tell apart. Counted
from the artifacts rather than from what the run wrote, so a re-run that skips
everything reports the same thing as the run that built them, over exactly the files a
corpus reads. A file found there that is not this level's pivot -- unreadable,
unstamped, or another level's -- is reported and left out of the tally rather than
ending the run, since the output tree is not this script's to police.
"""

import argparse
import functools
import pathlib

import pyarrow as pa
import pyarrow.parquet as pq

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


def _empty_artifacts(
    pivots_dir: str, level_path: str
) -> tuple[int, list[str], list[str]]:
    """Counts a level's artifacts, and says which hold no rows and which were skipped.

    Taken over exactly the files ``pivot.artifact_paths`` lists, which is what a corpus
    reads: a report over any other set describes a tree nobody queries. Their rows are
    read from the Parquet footers rather than tallied as the run writes, so the answer
    covers the artifacts this run skipped as already current, and a re-run says what the
    run that built them said.

    A file the reader would pick up that is not this level's pivot is skipped rather
    than counted -- an artifact of another level, or something left in the tree by hand
    -- and so is one that cannot be read at all, a truncated copy or a bad block. A run
    killed mid-write is not among them: atomic_io names its temp for the destination
    with a .tmp suffix, which neither this nor a reader ever lists.

    Args:
        pivots_dir: The directory the levels sit under.
        level_path: The level whose artifacts to look for.

    Returns:
        How many pivots of ``level_path`` there are, the paths of those holding no rows,
        and the paths skipped -- each named rather than tallied, since a count alone
        sends whoever reads it looking through the whole tree.
    """
    found, empty, skipped = 0, [], []
    for path in pivot.artifact_paths(pivots_dir, level_path):
        try:
            metadata = pq.read_metadata(path)
        except (OSError, pa.ArrowInvalid, pa.ArrowNotImplementedError):
            skipped.append(f"{path} cannot be read as Parquet")
            continue
        stamped = (metadata.metadata or {}).get(pivot.META_PIVOT_PATH.encode())
        if stamped is None:
            skipped.append(f"{path} carries no pivot level")
            continue
        if stamped.decode() != level_path:
            skipped.append(f"{path} holds {stamped.decode()}")
            continue
        found += 1
        if metadata.num_rows == 0:
            empty.append(str(path))
    return found, empty, skipped


def main(args: argparse.Namespace) -> None:
    """Derives a pivot artifact per level for every projection matching the pattern.

    Args:
        args: Parsed arguments, read for ``input_pattern``, ``output_dir``, ``levels``,
            and ``force``.

    Raises:
        ValueError: If a requested level is not one the projection schema has or is
            named twice; if a projection changed while it was being derived, or its
            count and its unnest disagreed; if the pattern matched no projections --
            which usually means it was aimed at the source tree, or at an output tree,
            rather than at the projections; or if a level's artifacts cannot be found
            after the run derived them, which means they were written somewhere no
            reader looks; or if the tree holds a level this run did not derive whose
            artifacts no longer match the schema. Silence in any of those would let a
            pipeline step downstream proceed as though the artifacts had been built.
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
            # Reaches the ordinals, which are top-level columns, and the fields of the
            # element struct, which are compared down through it: a level added above
            # this one gives every level beneath it another ordinal, and a field added
            # to the element gives it another leaf. Artifacts carrying the old set of
            # either are derived again rather than stamping as current.
            schema=pivot.schema(pivot.LEVELS[level_path]),
            force=args.force,
            parent_artifact=projection.ARTIFACT,
        )
        if not written and not skipped:
            raise ValueError(
                f"no pivots derived for {level_path}: none of the {ignored} matches "
                f"for {args.input_pattern!r} are projections"
            )
        artifacts, empty, unread = _empty_artifacts(args.output_dir, level_path)
        logger.info(
            "%s: %d written, %d already current, %d of %d artifacts empty",
            level_path,
            written,
            skipped,
            len(empty),
            artifacts,
        )
        if len(empty) < artifacts:
            # Named one by one where only some are empty, which is the shape that needs
            # looking at. Where all of them are the warning below says so once, and a
            # line per artifact would bury it on a corpus whose levels are mostly empty.
            for path in empty:
                logger.info("%s: no elements at %s", path, level_path)
        for reason in unread:
            logger.warning("%s: not counted at %s", reason, level_path)
        if artifacts < written + skipped:
            # Every artifact this run wrote or found current, and a reader cannot list
            # them all -- so some of them went somewhere nothing looks, and the tree is
            # not the one a query will read. Orphans left by a larger earlier run push
            # this the other way and cannot raise it.
            #
            # Raised rather than reported, because the denominator below shrinks to
            # match whatever was found: a level half of which is invisible reports the
            # visible half as the whole, and a level entirely invisible reports zero of
            # zero, which reads as every artifact being empty. That is the one thing
            # this report exists to say, said about the search rather than the corpus.
            raise ValueError(
                f"{directory} holds {artifacts} pivot artifacts for {level_path}, "
                f"fewer than the {written + skipped} this run wrote or found current; "
                "the rest are where no reader looks"
                + (f" ({len(unread)} not counted: {unread})" if unread else "")
            )
        if len(empty) == artifacts:
            # Ordinary for a level nothing in this corpus records, and identical on
            # disk to a level whose count came back wrong: an empty artifact is stamped
            # current like any other, so no later run revisits it and no reader tells
            # the two apart. Nothing else in the chain aggregates rows across a tree,
            # or says anything about them above INFO.
            logger.warning("%s: every artifact is empty at this level", level_path)
    behind = base.left_behind(
        args.output_dir,
        pivot.ARTIFACT,
        {level: pivot.schema(pivot.LEVELS[level]) for level in pivot.LEVELS},
        pivot.artifact_paths,
        levels,
    )
    if behind:
        # Raised rather than reported: the run succeeded at what it was asked to do and
        # still left a tree a corpus refuses, and a pipeline that reads the exit status
        # would go on to publish it. Naming the levels says which to add to --levels.
        raise ValueError(
            f"{args.output_dir} holds {len(behind)} levels this run did not derive "
            "whose artifacts no longer match the schema: "
            + "; ".join(f"{level} ({reason})" for level, reason in behind.items())
        )


if __name__ == "__main__":
    main(parse_args())
