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

"""One row per structure occurrence, for the paths a structure query can be routed to.

A quantifier whose body is a structure predicate on the element's own ``smiles`` is
answered by a semi-join rather than by reading the projection: which reactions hold a
match is a separate cost from the chemistry, and the larger one. The relation that
answers it is one row per occurrence, carrying the reaction, the structure, and the
element's own ``reaction_role`` -- keeping the role beside the structure is what makes
"pyridine as the solvent" a condition on one row rather than an intersection of two
reaction sets, which over-returns.

**Dataset-local IDs, always.** An occurrence names its structure by the ``structure_id``
its own dataset assigned, never by a corpus-wide one. A corpus-wide ID is that plus an
offset that is a running total over the corpus's datasets, so it belongs to no artifact:
the same rows read beside a different set of datasets need a different offset, and an ID
baked in here would stay in range while naming another dataset's molecule -- passing on
the whole corpus and failing only on a subset of it. The corpus joins the offset on at
open, keyed by the source dataset every artifact names; see
:mod:`ord_schema.search.execute`.

That rule is what lets this be an artifact at all rather than a relation assembled at
open, which is what :mod:`ord_schema.search.execute` does: over ORD it builds
18,847,978 rows holding 1.19 GiB, wanting 5-6.5 GB of DuckDB memory and 16-25 GB of
temporary files to get there, before the first query can be answered. The same rows
derived here are 242 MB of Parquet. **Nothing reads them yet** -- the corpus that does
is a separate change, and until it lands a tree derived from this is inert.

Derived from the **pivot** over the level an indexed path ranges within, which already
holds one row per element: this artifact is that projected down to three columns, so the
two cannot disagree about which elements exist. A path descending from its level through
singular struct fields -- an authentic standard is one compound per measurement -- is
read from that level's pivot through the remainder; see ``ord_schema.artifacts.pivot``.
"""

import glob
import os
import pathlib

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from ord_schema import atomic_io
from ord_schema.artifacts import base, pivot, projection, structures
from ord_schema.logging import get_logger

logger = get_logger(__name__)

ARTIFACT = "occurrences"

# The indexed path a file holds, stamped so a reader can tell one from another without
# trusting where it sits; see ``occurrence_path``.
META_OCCURRENCE_PATH = "ord.occurrence_path"

# The element field an occurrence carries beside its structure. One rather than all of
# them: a row is what a semi-join reads, and every column on it is paid for by every
# occurrence in the corpus.
INDEXED_FIELD = "reaction_role"

SCHEMA = pa.schema(
    [
        pa.field("reaction_id", pa.string(), nullable=False),
        # Local to the dataset this was derived from, and meaningless beside another's
        # structures artifact; see the module docstring.
        pa.field("structure_id", pa.uint32(), nullable=False),
        pa.field(INDEXED_FIELD, pa.string()),
    ]
)


def indexed_paths(
    schema: pa.Schema = projection.SCHEMA,
) -> dict[str, tuple[pivot.RepeatedLevel, tuple[str, ...]]]:
    """Returns the paths this artifact covers, and the pivot each is read from.

    Taken from the same schema walk the structures artifact is built from, so the
    occurrences reach the elements the projection does by construction rather than by a
    list here agreeing with one there. Every path the schema can carry a structure at
    has to come out covered: one dropped would leave the structures sitting only there
    reachable by no query, which reads as a corpus holding none of them.

    Args:
        schema: The projection schema.

    Returns:
        Each structure-bearing path, mapped to the repeated level whose pivot holds it
        and the segments from that level's element down to the structure.

    Raises:
        ValueError: If a structure-bearing path is reached by no repeated level, or
            carries no ``INDEXED_FIELD``, so no pivot can answer for it; or if the
            schema carries no structure at all. Raised at import, naming the path,
            because the answer is to change this module rather than to retry.
    """
    paths: dict[str, tuple[pivot.RepeatedLevel, tuple[str, ...]]] = {}
    # Levels from the schema the caller passed, so asking what a written file covers is
    # answered about that file. Reaching against the module's default levels would
    # answer about the schema this library was imported with -- the disagreement between
    # two walks that this function exists to prevent.
    levels = pivot.repeated_levels(schema)
    for _, path, dtype in structures.structure_levels(schema):
        if INDEXED_FIELD not in [field.name for field in dtype]:
            raise ValueError(
                f"{path} carries a structure but no {INDEXED_FIELD}, so an occurrence "
                "row cannot bind a query to the element it matched, and a structure "
                "sitting only there would look like one this artifact lost"
            )
        reached = pivot.reach(path, levels)
        if reached is None:
            raise ValueError(
                f"{path} carries a structure but no repeated level reaches it, so no "
                "pivot holds its elements and they can be read only by unnesting the "
                "projection"
            )
        level, remainder, _ = reached
        paths[path] = (level, remainder)
    if not paths:
        raise ValueError(
            "the schema carries no structure at any path, so this artifact has nothing "
            "to hold and a corpus reading it would index nothing"
        )
    return paths


PATHS = indexed_paths()


def _select(level: pivot.RepeatedLevel, remainder: tuple[str, ...]) -> str:
    """Returns the SELECT projecting a pivot's rows down to occurrences.

    Args:
        level: The repeated level the pivot holds.
        remainder: Segments from that level's element down to the structure.

    Returns:
        The SELECT, over a view named ``elements``.
    """
    element = ".".join([pivot.ELEMENT, *remainder])
    # An element carrying no structure is no occurrence: the pivot is complete, holding
    # a row for every element including one whose fields are all NULL, and this is not.
    # S608: the level and its remainder come from this module's walk of the schema.
    return f"""
        SELECT {pivot.REACTION_ID}, {element}.structure_id AS structure_id,
               {element}.{INDEXED_FIELD} AS {INDEXED_FIELD}
        FROM elements
        WHERE {element}.structure_id IS NOT NULL
        """  # noqa: S608


def _check_pivot(source: str | os.PathLike[str], level: str) -> base.Stamps:
    """Returns the stamps of a current pivot over ``level``, refusing anything else.

    Args:
        source: Path to the artifact to derive from.
        level: The level it has to hold.

    Returns:
        Its stamps, which this artifact passes through.

    Raises:
        ValueError: If it is not a pivot, holds another level, or is stale.
    """
    stamps = base.load_stamps(source)
    if stamps.artifact != pivot.ARTIFACT:
        raise ValueError(f"{source} is a {stamps.artifact}, not a {pivot.ARTIFACT}")
    held = pivot.pivot_path(source)
    if held != level:
        raise ValueError(
            f"{source} holds the pivot over {held}, but the occurrences at this path "
            f"are read from {level}"
        )
    if not base.stamps_are_current(stamps, pivot.ARTIFACT):
        raise ValueError(
            f"{source} is stale, and occurrences derived from it would claim a "
            "provenance they do not have; derive the pivot again first"
        )
    return stamps


def write_occurrences(
    source: str | os.PathLike[str],
    output: str | os.PathLike[str],
    /,
    *,
    path: str,
    source_md5: str | None = None,
    source_dataset_id: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 100_000,
) -> int:
    """Derives the occurrences at one indexed path from a pivot artifact.

    Args:
        source: The pivot over the level ``path`` ranges within.
        output: Where to write.
        path: The indexed path, which must be a key of ``PATHS``.
        source_md5: The source dataset's hash; taken from ``source`` when unset.
        source_dataset_id: The source dataset's ID; taken from ``source`` when unset.
        compression: Parquet codec.
        row_group_size: Rows per group.

    Returns:
        How many occurrences were written. Zero is ordinary -- most corpora record no
        authentic standards -- and produces a valid empty file rather than none.

    Raises:
        ValueError: If ``path`` is not an indexed path, or if ``source`` is not a
            current pivot over the level that path ranges within.
    """
    reached = PATHS.get(path)
    if reached is None:
        raise ValueError(
            f"{path} is not a path this artifact covers, so it has nothing to hold. "
            f"Known paths: {sorted(PATHS)}"
        )
    level, remainder = reached
    parent = _check_pivot(source, level.path)
    if source_md5 is None:
        source_md5 = parent.source_md5
    if source_dataset_id is None:
        source_dataset_id = parent.source_dataset_id
    stamps = base.current_stamps(ARTIFACT, source_dataset_id, source_md5)
    metadata = base.to_metadata(stamps) | {META_OCCURRENCE_PATH: path}
    target = SCHEMA.with_metadata(metadata)
    connection = duckdb.connect()
    written = 0
    try:
        # No row order to keep: a reader semi-joins this on reaction_id, and which
        # occurrence came back first says nothing. Same trade the pivot makes.
        connection.execute("SET preserve_insertion_order=false")
        # Through the relational API rather than a path interpolated into SQL, which a
        # filename holding a quote would otherwise close.
        connection.read_parquet(str(pathlib.Path(source))).create_view("elements")
        with (
            atomic_io.atomic_path(output) as temp_path,
            pq.ParquetWriter(temp_path, target, compression=compression) as writer,
        ):
            # A path with no occurrences writes no batches, and the writer's close still
            # produces a valid zero-row file carrying the schema and stamps -- what a
            # reader globs and finds nothing in, rather than a file that is missing.
            reader = connection.execute(_select(level, remainder)).to_arrow_reader(
                row_group_size
            )
            for batch in reader:
                # Cast rather than trust: DuckDB's own widths are its business, and the
                # artifact's schema is a promise to whoever reads it later.
                writer.write_table(
                    pa.Table.from_batches([batch])
                    .cast(SCHEMA)
                    .replace_schema_metadata(metadata)
                )
                written += batch.num_rows
    finally:
        connection.close()
    logger.info("%s: %d occurrences at %s", source, written, path)
    return written


def occurrence_path(path: str | os.PathLike[str]) -> str | None:
    """Returns the indexed path an artifact holds, or None if it holds none.

    Args:
        path: Path to a Parquet file.

    Returns:
        The dotted path stamped in the footer, or None if the file is not an
        occurrences artifact.
    """
    metadata = pq.read_schema(path).metadata or {}
    value = metadata.get(META_OCCURRENCE_PATH.encode())
    return value.decode() if value is not None else None


def artifact_paths(
    occurrences_dir: str | os.PathLike[str], path: str
) -> list[pathlib.Path]:
    """Returns the artifacts a reader finds for one indexed path, in a stable order.

    The one definition of what a path's artifacts are, for the same reason
    ``pivot.artifact_paths`` is: a build judging its own output by a wider rule than a
    reader applies calls a tree healthy that no corpus can read.

    Recursive, because the tree mirrors the projections it descends from: a sharded
    corpus puts its files under ``<path>/<shard>/`` rather than directly under the path.

    Args:
        occurrences_dir: The directory the indexed paths sit under.
        path: The indexed path, as the query grammar names it.

    Returns:
        The paths, sorted. Empty if the path has no directory at all.
    """
    directory = pathlib.Path(occurrences_dir) / path
    # Escaped, because only the last two segments are a pattern: a directory is a name
    # someone chose, and one holding a bracket or a star would otherwise be read as a
    # character class and match somewhere else, or nowhere.
    return sorted(
        pathlib.Path(found)
        for found in glob.glob(
            f"{glob.escape(str(directory))}/**/*.parquet", recursive=True
        )
    )
