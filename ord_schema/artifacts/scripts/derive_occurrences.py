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

r"""Derives occurrence artifacts from pivots, one per structure-bearing path.

A structure query finds the reactions holding a match by semi-joining a relation of one
row per structure occurrence; see ``ord_schema.artifacts.occurrences`` for the shape and
why it holds dataset-local IDs. Built in process that relation is 1.19 GiB and a minute
over ORD, and wants several gigabytes of DuckDB memory to assemble, so it belongs here.

Each indexed path gets its own subdirectory, mirroring the pivots it descends from::

    derive_occurrences.py --pivots_dir=pivots --output_dir=occurrences

    pivots/inputs.components/aa/ord_dataset-<id>.parquet
        -> occurrences/inputs.components/aa/ord_dataset-<id>.parquet
    pivots/outcomes.products.measurements/aa/ord_dataset-<id>.parquet
        -> occurrences/outcomes.products.measurements.authentic_standard/aa/...

Which paths those are is not a choice: a corpus checks that its index reaches every
structure it holds, so a tree missing one leaves structures no query can find. Deriving
a subset with --paths is for rebuilding one, not for choosing what to carry.

The pivots for the levels those paths range within have to exist first --
``derive_pivots.py --levels`` names them, and ``--print_levels`` here prints exactly
that list. Deriving from the pivot rather than from the projection is what keeps the two
from disagreeing about which elements exist: the occurrences are that artifact projected
down to three columns.

Artifacts whose footer already records the current source content, library version, and
artifact version are skipped, so re-running is cheap; --force rewrites them anyway.

These are errors: an --output_dir that would write over any input, a path this schema
does not carry a structure at, a pivot that is stale or holds another level, a run that
finds no pivots for a path at all, and a path the tree already holds that this run did
not name whose artifacts are short a column the schema declares.
"""

import argparse
import functools
import glob
import pathlib

import pyarrow.parquet as pq

from ord_schema.artifacts import base, occurrences, pivot
from ord_schema.logging import get_logger

logger = get_logger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--pivots_dir", help="Directory holding the pivot artifacts to derive from"
    )
    parser.add_argument("--output_dir", help="Directory to write artifacts beneath")
    parser.add_argument(
        "--paths",
        nargs="+",
        default=sorted(occurrences.PATHS),
        help="Indexed paths to derive; every one the schema carries by default",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite artifacts that are already current",
    )
    parser.add_argument(
        "--print_levels",
        action="store_true",
        help="Print the pivot levels these paths need, for derive_pivots --levels",
    )
    return parser.parse_args(argv)


def check_paths(paths: list[str]) -> list[str]:
    """Returns the requested paths, refusing any this artifact does not cover.

    Args:
        paths: The paths asked for.

    Returns:
        Them, deduplicated and sorted.

    Raises:
        ValueError: If a path is not one the schema carries a structure at, or is named
            twice -- which would derive it twice and report the second as already
            current.
    """
    unknown = sorted(set(paths) - set(occurrences.PATHS))
    if unknown:
        raise ValueError(
            f"no structure sits at {unknown}, so there are no occurrences to derive "
            f"there. Known paths: {sorted(occurrences.PATHS)}"
        )
    repeated = sorted({path for path in paths if paths.count(path) > 1})
    if repeated:
        raise ValueError(f"named more than once: {repeated}")
    return sorted(set(paths))


def levels_for(paths: list[str]) -> list[str]:
    """Returns the pivot levels the given paths are read from.

    Args:
        paths: Indexed paths, already checked.

    Returns:
        The levels, deduplicated and sorted, as ``derive_pivots --levels`` takes them.
    """
    return sorted({occurrences.PATHS[path][0].path for path in paths})


def _audit(
    output_dir: str, path: str
) -> tuple[int, list[str], list[tuple[str, str | None]]]:
    """Reads back what a reader would find for one path, and says what is wrong with it.

    Taken over exactly the files ``occurrences.artifact_paths`` lists, which is what a
    corpus reads: a report over any other set describes a tree nobody queries. Rows come
    from the Parquet footers rather than being tallied as the run writes, so the answer
    covers artifacts this run skipped as already current, and a re-run says what the run
    that built them said.

    Args:
        output_dir: The directory the paths sit under.
        path: The indexed path to audit.

    Returns:
        How many artifacts a reader finds, which of them hold no rows, and which are
        stamped with a path other than this one -- the last carrying what each holds,
        which is None for a file that is not an occurrences artifact at all.
    """
    empty: list[str] = []
    mismatched: list[tuple[str, str | None]] = []
    found = occurrences.artifact_paths(output_dir, path)
    for artifact in found:
        held = occurrences.occurrence_path(artifact)
        if held != path:
            mismatched.append((str(artifact), held))
            continue
        if not pq.read_metadata(artifact).num_rows:
            empty.append(str(artifact))
    return len(found), empty, mismatched


def main(args: argparse.Namespace) -> None:
    """Derives the occurrences at each requested path.

    Args:
        args: Parsed command-line arguments.

    Raises:
        ValueError: If a requested path is not one the schema carries a structure at or
            is named twice; if a pivot is stale or holds another level; or if a path's
            level has no pivots at all -- which usually means they have not been derived
            yet, and would otherwise leave a tree a corpus refuses for reaching too few
            structures; or if the tree holds a path this run did not derive whose
            artifacts are short a column the schema declares.
    """
    paths = check_paths(args.paths)
    if args.print_levels:
        print(" ".join(levels_for(paths)))
        return
    if not args.pivots_dir or not args.output_dir:
        raise ValueError("--pivots_dir and --output_dir are both required")
    for path in paths:
        level = occurrences.PATHS[path][0].path
        # --pivots_dir is a directory someone named, not a pattern they wrote, so a
        # bracket or a star in it is part of the name: escaped, it matches that tree
        # rather than a sibling fitting the character class. The root goes separately
        # because glob_root reads the escape as the wildcard it is spelled with, and
        # would stop the mirror short of it.
        root = pathlib.Path(args.pivots_dir) / level
        pattern = f"{glob.escape(str(root))}/**/*.parquet"
        directory = pathlib.Path(args.output_dir) / path
        # Globbed exactly as derive_tree will glob it, rather than through
        # pivot.artifact_paths: a check that finds pivots the derivation then cannot
        # would turn this friendly message into derive_tree's report of an empty
        # pattern.
        if not glob.glob(pattern, recursive=True):
            raise ValueError(
                f"no pivots over {level} under {args.pivots_dir}, so there are no "
                f"occurrences to derive for {path}. Derive them first with "
                f"derive_pivots.py --levels {level}"
            )
        written, skipped, ignored = base.derive_tree(
            pattern,
            str(directory),
            artifact=occurrences.ARTIFACT,
            write=functools.partial(occurrences.write_occurrences, path=path),
            schema=occurrences.SCHEMA,
            force=args.force,
            parent_artifact=pivot.ARTIFACT,
            root=root,
        )
        if not written and not skipped:
            raise ValueError(
                f"no occurrences derived for {path}: none of the {ignored} matches for "
                f"{pattern!r} are pivots over {level}. Derive them first with "
                f"derive_pivots.py --levels {level}"
            )
        found, empty, mismatched = _audit(args.output_dir, path)
        for name, held in mismatched:
            logger.error("%s: holds %s, not %s", name, held or "no path", path)
        if mismatched:
            # Every indexed path writes the same three columns, so a file filed under
            # the wrong one passes the currency check: same artifact name, same source
            # hash, same versions, same columns. It is then skipped by every later run
            # and read as this path's occurrences ever after, answering a quantifier
            # with another level's elements. The stamped path is the only thing that
            # tells them apart, and this is the only place that reads it.
            raise ValueError(
                f"{directory} holds {len(mismatched)} artifacts stamped with another "
                f"path; they are current for {path} by every other measure and would "
                "answer its quantifiers with the wrong elements"
            )
        logger.info(
            "%s: %d written, %d already current, %d of %d artifacts empty",
            path,
            written,
            skipped,
            len(empty),
            found,
        )
        if len(empty) < found:
            # Named one by one where only some are empty, which is the shape worth
            # looking at. Where all of them are, the warning below says so once.
            for name in empty:
                logger.info("%s: no occurrences at %s", name, path)
        elif found:
            # Ordinary for a path nothing in this corpus records -- most corpora have no
            # authentic standards -- and identical on disk to a derivation that filtered
            # wrongly: an empty artifact is stamped current like any other, so no later
            # run revisits it and no reader tells the two apart.
            logger.warning("%s: every artifact is empty at this path", path)
        if found < written + skipped:
            # Every artifact this run wrote or found current, and a reader cannot list
            # them all, so some went where nothing looks and the tree is not the one a
            # corpus will read.
            raise ValueError(
                f"{directory} holds {found} occurrence artifacts for {path}, fewer "
                f"than the {written + skipped} this run wrote or found current; the "
                "rest are where no reader looks"
            )
    behind = base.left_behind(
        args.output_dir,
        occurrences.ARTIFACT,
        # One schema for every path: the three columns are the same wherever a
        # structure sits, so nothing on disk is short one until this schema grows.
        dict.fromkeys(occurrences.PATHS, occurrences.SCHEMA),
        occurrences.artifact_paths,
        paths,
    )
    if behind:
        # Raised rather than reported: the run succeeded at what it was asked to do and
        # still left a tree a corpus refuses. Naming the paths says which to add to
        # --paths.
        raise ValueError(
            f"{args.output_dir} holds {len(behind)} paths this run did not derive "
            "whose artifacts no longer match the schema: "
            + "; ".join(f"{path} ({reason})" for path, reason in behind.items())
        )


if __name__ == "__main__":
    main(parse_args())
