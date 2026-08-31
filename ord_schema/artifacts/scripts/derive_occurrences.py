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
does not carry a structure at, a pivot that is stale or holds another level, and a run
that finds no pivots for a path at all.
"""

import argparse
import functools
import pathlib

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


def main(args: argparse.Namespace) -> None:
    """Derives the occurrences at each requested path.

    Args:
        args: Parsed command-line arguments.

    Raises:
        ValueError: If a requested path is not one the schema carries a structure at or
            is named twice; if a pivot is stale or holds another level; or if a path's
            level has no pivots at all -- which usually means they have not been derived
            yet, and would otherwise leave a tree a corpus refuses for reaching too few
            structures.
    """
    paths = check_paths(args.paths)
    if args.print_levels:
        print(" ".join(levels_for(paths)))
        return
    if not args.pivots_dir or not args.output_dir:
        raise ValueError("--pivots_dir and --output_dir are both required")
    for path in paths:
        level = occurrences.PATHS[path][0].path
        # Escaped for the same reason artifact_paths escapes: only the last segments
        # are a pattern, and a shard directory holding a bracket is a name, not a class.
        pattern = f"{pathlib.Path(args.pivots_dir) / level}/**/*.parquet"
        directory = pathlib.Path(args.output_dir) / path
        if not pivot.artifact_paths(args.pivots_dir, level):
            # Said here rather than left to derive_tree, which reports an empty glob
            # against the pattern it was handed and cannot know what would fill it.
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
        )
        if not written and not skipped:
            raise ValueError(
                f"no occurrences derived for {path}: none of the {ignored} matches for "
                f"{pattern!r} are pivots over {level}. Derive them first with "
                f"derive_pivots.py --levels {level}"
            )
        found = len(occurrences.artifact_paths(args.output_dir, path))
        logger.info(
            "%s: %d written, %d already current, %d artifacts",
            path,
            written,
            skipped,
            found,
        )
        if found < written + skipped:
            # Every artifact this run wrote or found current, and a reader cannot list
            # them all, so some went where nothing looks and the tree is not the one a
            # corpus will read.
            raise ValueError(
                f"{directory} holds {found} occurrence artifacts for {path}, fewer "
                f"than the {written + skipped} this run wrote or found current; the "
                "rest are where no reader looks"
            )


if __name__ == "__main__":
    main(parse_args())
