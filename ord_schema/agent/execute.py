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

"""Running a compiled query against a corpus of projections and structures artifacts.

The compiler (:mod:`ord_schema.agent.query`) leaves two kinds of work undone on
purpose, and this module is where both happen:

* **Names become structures.** Compound names bind as resolved SMILES, through
  :mod:`ord_schema.resolvers` by default or any injected resolver.
* **Structure predicates become bitmaps.** Each compiles to a bitmap test on the
  element's ``structure_id``, and the chemistry runs here. Substructure goes through an
  RDKit ``SubstructLibrary`` built over the corpus: it screens on the pattern
  fingerprints the artifact already stores and matches the survivors exactly, across
  its own threads with the GIL released, so a server can share one corpus between
  concurrent requests without forking. Similarity stays in SQL and has no verification
  step at all -- Tanimoto is defined on the Morgan fingerprint, so the screen is the
  whole answer.

Structure IDs are dataset-local, so a corpus-wide bitmap needs an offset per file.
``Corpus`` pairs every projection with its structures artifact, refuses a pair whose
stamps disagree about the source dataset -- the two files are two statements of one
derivation, and a mismatched pair would join IDs to the wrong molecules silently --
and publishes a ``reactions`` relation carrying each row's offset in the column the
compiled SQL expects.

A ``Corpus`` is safe to share between concurrent searches. A ``DuckDBPyConnection``
holds the pending result of its last ``execute``, so two threads stepping on one read
each other's rows rather than raising; each search therefore takes its own cursor. A
cursor in turn sees only the catalog, not the objects ``register`` attaches to the
connection that registered them, so every relation here is a catalog table or view.

Finding the reactions that hold a match is a scan of every reaction, which costs more
than the chemistry for all but the narrowest queries. An **occurrence index** -- one row
per structure occurrence, carrying the corpus-wide ID, the path it sat at, and the
element's own ``reaction_role`` -- turns that scan into a semi-join against a narrow
table. It is spent through the compiler's ``ElementIndex`` hook, one quantifier at a
time: an ``exists`` whose body asks one occurrence row's worth -- a structure predicate
on the element's own ``smiles``, at most a ``reaction_role`` equality beside it --
becomes a condition on the row, and everything around it compiles exactly as it would
have. An aggregate over the quantifier, a scalar beside it, an ordering, a limit, a
negation, a disjunction, or a second quantifier each compose with the condition rather
than disqualifying it. Keeping the role beside the structure in the index is what
preserves element binding: "pyridine as the solvent" stays a condition on a single row
rather than an intersection of two reaction sets, which over-returns.

What still compiles over the elements, and so reads the projection: a quantifier whose
body binds any other element field, reaches a nested one, holds no structure predicate
or more than one, or is not a plain conjunction -- and every ``forall``, since the index
shows the elements that match, never that every element does.

An occurrence names its reaction by ID, so the index is built only over a corpus that
states each reaction once; two files carrying one reaction would have each copy's
structures answering for the other.

The grammar bounds what a query can cost (one pass and a sort), and ``search``
takes a wall-clock timeout that interrupts the final query. Name resolution, the
one-time library and index builds, screening, and verification all run before the timer
starts: each is bounded by the corpus rather than by the query, so a slow one is slow
for every caller and shows up in the logs rather than in a timeout. The two builds have
separate triggers -- the library on the first substructure predicate, the index on the
first query it takes a clause of -- and a search wanting both pays both, upwards of a
minute over the whole corpus, on top of whatever timeout it was given.
"""

import array
import collections
import contextlib
import dataclasses
import glob
import math
import re
import threading
import time
from collections.abc import Callable, Iterable, Iterator, Sequence
from typing import Any, Self

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, DataStructs
from rdkit.Chem import rdSubstructLibrary

from ord_schema import artifacts, projection, resolvers, structures
from ord_schema.agent import pivot, query
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Structures read per batch while building the library. Bounds peak memory during the
# build without making the per-batch Python overhead matter.
_BUILD_BATCH = 50_000

# Match sets kept for reuse. Each is a bitmap over the corpus ID space -- 2 MB at ORD's
# scale -- so this trades a few tens of megabytes for the repeat of a query that costs
# seconds.
_CACHED_MATCHES = 16

# How much memory the materialized column sets may hold between them, when the caller
# states no budget. The whole projection is 18.46 GB in memory at ORD's scale and no
# default holds it, so every budget is a partial cache and the question is only which
# columns fit. This holds any single one of the large ones -- workups is 5.08 GB,
# inputs 3.93 GB, outcomes 3.04 GB -- which is what a query pairing a structure clause
# with a projection one needs. A corpus asked for more says so on every query it costs,
# naming this argument; see _warn_refused.
#
# Held to after a set is built rather than before, since what one costs is known only
# once it exists: the peak is this plus the largest single set, and building a set costs
# its size whether or not it is kept.
_NARROW_BUDGET_BYTES = 4 * 1024**3

# Top-level projection columns, which is the granularity a compiled query names them at.
# A narrower table holding only the ones a query mentions answers it identically, so the
# query is compiled against it as written.
_TOP_LEVEL = tuple(projection.SCHEMA.names)


# What a structure predicate's answer depends on: the operation, whether the string
# is read as SMARTS or as SMILES, the string, and the similarity threshold.
_MatchKey = tuple[str, bool, str, float | None]


def _reads_as_smarts(parameter: query.StructureParameter) -> bool:
    """Returns whether the predicate's string is read as SMARTS rather than SMILES.

    Which parser reads a string is part of what the string means -- ``C1=CC=CC=C1`` is
    benzene as SMILES and six aliphatic carbons as SMARTS -- so the query molecule and
    the key a match set is cached under have to make this call the same way.

    Args:
        parameter: The structure predicate.

    Returns:
        True for a substructure predicate stating its own pattern, which is the only
        form a query writes in SMARTS; a resolved compound is always SMILES.
    """
    return parameter.pattern is not None and parameter.op == "substructure"


class PairingError(ValueError):
    """The projections and structures artifacts do not form a consistent corpus."""


def _format_bytes(value: float) -> str:
    """Returns a byte count in the largest unit it reaches, for a person to read.

    Args:
        value: Bytes.

    Returns:
        The figure and its unit, e.g. ``3.3 GB``. A budget stated in gigabytes says
        nothing about a table that costs kilobytes, and a warning naming both has to
        make the two comparable.
    """
    for unit, scale in (("GB", 1024**3), ("MB", 1024**2), ("kB", 1024)):
        if abs(value) >= scale:
            return f"{value / scale:.1f} {unit}"
    return f"{value:.0f} B"


def _resolve_with_resolvers(name: str) -> str:
    """Resolves a compound name to SMILES through ord_schema.resolvers."""
    smiles, _ = resolvers.resolve_name("name", name)
    return smiles


def _sql_string(value: str) -> str:
    """Returns ``value`` as a SQL string literal, single quotes escaped."""
    escaped = value.replace("'", "''")
    return f"'{escaped}'"


def _sql_strings(values: Iterable[str]) -> str:
    """Returns ``values`` as a SQL list literal, single quotes escaped."""
    return "[" + ", ".join(_sql_string(value) for value in values) + "]"


def _max_structure_id(path: str) -> int | None:
    """Returns the largest ``structure_id`` a projection carries.

    Read from Parquet footer statistics, so it decodes no column data however large the
    projection is.

    Args:
        path: Path to a projection.

    Returns:
        The maximum across every ``structure_id`` column, or None when the projection
        holds no compound with a structure.
    """
    with pq.ParquetFile(path) as projected:
        metadata = projected.metadata
        columns = [
            index
            for index in range(metadata.num_columns)
            if metadata.schema.column(index).name == "structure_id"
        ]
        largest = None
        for group in range(metadata.num_row_groups):
            for index in columns:
                statistics = metadata.row_group(group).column(index).statistics
                if statistics is not None and statistics.has_min_max:
                    value = statistics.max
                    largest = value if largest is None else max(largest, value)
    return largest


def _index(pattern: str, artifact: str, require_current: bool) -> dict[str, str]:
    """Returns the artifacts matching ``pattern``, keyed by source dataset.

    Keyed by the source hash rather than the filename because that is what an
    artifact *is*: two files pair when they restate the same dataset, whatever
    they are called or wherever they sit. Keying on the basename would let two
    files in different directories collapse onto one another silently.

    Args:
        pattern: Glob matching the artifact files.
        artifact: Artifact name every match must hold.
        require_current: Refuse artifacts not written by the current versions.

    Returns:
        A mapping from source dataset hash to path.

    Raises:
        PairingError: If a match holds another artifact, is stale under
            ``require_current``, or restates a dataset another match already did.
    """
    index: dict[str, str] = {}
    for path in sorted(glob.glob(pattern, recursive=True)):
        stamps = artifacts.load_stamps(path)
        if stamps.artifact != artifact:
            raise PairingError(f"{path} is a {stamps.artifact}, not a {artifact}")
        if require_current and not artifacts.stamps_are_current(stamps, artifact):
            raise PairingError(f"{path} is stale; derive it again first")
        if stamps.source_md5 in index:
            raise PairingError(
                f"{path} and {index[stamps.source_md5]} are both {artifact} "
                "artifacts of the same source dataset; which one answers a query "
                "would be arbitrary"
            )
        index[stamps.source_md5] = path
    return index


def _pair(
    projection_pattern: str, structures_pattern: str, require_current: bool
) -> list[tuple[str, str]]:
    """Returns (projection, structures) path pairs, verified by their stamps.

    Args:
        projection_pattern: Glob matching the projection artifacts.
        structures_pattern: Glob matching the structures artifacts.
        require_current: Refuse artifacts not written by the current versions.

    Returns:
        One pair per source dataset, ordered by the projection's source hash.

    Raises:
        PairingError: If no projection matches, or if either side states a dataset the
            other does not, since a projection's IDs index its own partner's molecules
            and nobody else's.
    """
    projections = _index(projection_pattern, projection.ARTIFACT, require_current)
    structure_files = _index(structures_pattern, structures.ARTIFACT, require_current)
    if not projections:
        raise PairingError(f"no projections matched: {projection_pattern}")
    unpaired = projections.keys() ^ structure_files.keys()
    if unpaired:
        orphans = sorted((projections | structure_files)[key] for key in unpaired)
        raise PairingError(
            "projections and structures artifacts do not pair up; these have no "
            f"counterpart derived from the same source dataset: {orphans}"
        )
    return [(projections[key], structure_files[key]) for key in sorted(projections)]


# A structure that did not parse still gets a library entry, so every structure ID
# resolves to one and the mapping needs no special case. An empty molecule fingerprints
# to no bits, and a query with no atoms is refused, so it can never match.
_UNPARSEABLE = Chem.Mol().ToBinary()
_NO_BITS = DataStructs.ExplicitBitVect(structures.PATTERN_FP_SIZE)

# The one element field the occurrence index carries besides the structure. A query
# binding anything else to the element runs against the projection instead.
_INDEXED_FIELD = "reaction_role"


def _indexed_paths(schema: pa.Schema = projection.SCHEMA) -> dict[str, str]:
    """Returns the paths an occurrence index covers, mapped to their traversals.

    Taken from the same schema walk the structures artifact is built from and resolved
    through the compiler, so the index reaches the elements the projection does by
    construction rather than by a list here agreeing with one there.

    Every path the schema can carry a structure at has to come out covered, which is
    what lets the build check what it indexed against what the corpus says it holds. A
    path this walk dropped would make a whole corpus look like one that lost structures.

    Args:
        schema: The projection schema.

    Returns:
        Dotted query paths to the DuckDB expression yielding their elements, for every
        level whose elements carry a structure.

    Raises:
        ValueError: If a structure-bearing level carries no ``_INDEXED_FIELD`` or does
            not resolve to a repeated expression, so the index cannot cover it; or if
            the schema carries no structure at all. Raised at import, naming the path,
            because the answer is to change this module rather than to retry.
    """
    paths: dict[str, str] = {}
    for _, path, dtype in structures.structure_levels(schema):
        resolved = query.resolve(path, schema=schema, allow_internal=True)
        if _INDEXED_FIELD not in [field.name for field in dtype]:
            raise ValueError(
                f"{path} carries a structure but no {_INDEXED_FIELD}, so the "
                "occurrence index cannot cover it and a structure sitting only "
                "there would look like one the index lost"
            )
        if not resolved.repeated:
            raise ValueError(
                f"{path} carries a structure but resolves to one element rather than "
                "a repeated expression, which the build cannot unnest"
            )
        paths[path] = resolved.expression
    if not paths:
        raise ValueError(
            "the schema carries no structure at any path, so an occurrence index has "
            "nothing to index and its build would assemble no SQL at all"
        )
    return paths


INDEXED_PATHS = _indexed_paths()


def _index_condition(
    path: str,
    fields: dict[str, str],
    allocate: Callable[[], str],
) -> str | None:
    """Returns a row condition standing for an element quantifier, or None.

    The compiler's ``ElementIndex``: it has already read the quantifier's body down
    to a structure predicate on the element's own ``smiles`` and the field
    equalities beside it, and this decides whether an occurrence row can carry that
    question. A row holds the path, the corpus-wide structure ID, and one element
    field, so a query binding any other field, or asking about a path the index does
    not cover, is left to the elements.

    What comes back is a semi-join rather than a whole query, which is what lets the
    rest of the query stay exactly what it was: the aggregate, the ordering, the
    limit, and every clause beside this one compile unchanged, and this condition
    narrows the same rows they are about. A reaction is one row of ``reactions``, so
    the join cannot multiply rows and needs no DISTINCT.

    Args:
        path: The path the quantifier ranges over.
        fields: Element fields the body requires, mapped to their literals.
        allocate: Names the structure parameter the match set binds under.

    Returns:
        The condition, or None to leave the quantifier compiled over the elements.
    """
    if path not in INDEXED_PATHS or set(fields) - {_INDEXED_FIELD}:
        return None
    conditions = [
        f"occurrence.path = {_sql_string(path)}",
        f"get_bit(CAST(${allocate()} AS BITSTRING), occurrence.global_id::INTEGER) = 1",
    ]
    role = fields.get(_INDEXED_FIELD)
    if role is not None:
        conditions.append(f"occurrence.{_INDEXED_FIELD} = {_sql_string(role)}")
    # S608: every fragment is a schema-derived path, a compiler-issued parameter
    # name, or an escaped literal. The inner relation is aliased so that the outer
    # reaction_id and the inner one cannot be confused for each other.
    return (
        "reaction_id IN (SELECT occurrence.reaction_id "  # noqa: S608
        f"FROM occurrences AS occurrence WHERE {' AND '.join(conditions)})"
    )


def _group(entry_of: array.array, count: int) -> tuple[array.array, array.array]:
    """Inverts structure-to-entry into entry-to-structures.

    Args:
        entry_of: The library entry each structure ID resolves to, indexed by ID.
        count: How many entries the library holds.

    Returns:
        ``(members, starts)``, where the structure IDs sharing entry ``i`` are
        ``members[starts[i]:starts[i + 1]]``. Counting sort, so one pass each way and
        no per-entry list object.
    """
    starts = array.array("I", bytes(4 * (count + 1)))
    for entry in entry_of:
        starts[entry + 1] += 1
    for index in range(1, count + 1):
        starts[index] += starts[index - 1]
    members = array.array("I", bytes(4 * len(entry_of)))
    cursor = starts[:-1]
    for global_id, entry in enumerate(entry_of):
        members[cursor[entry]] = global_id
        cursor[entry] += 1
    return members, starts


def _mentioned(sql: str) -> frozenset[str]:
    """Returns the top-level projection columns a compiled query names.

    Read back off the SQL rather than walked from the query, because the SQL is
    what has to resolve: a column named there and missing from the table is a
    catalog error, and one named only in the query is nothing at all.

    Matched as a word outside a string literal and not behind a dot. A top-level
    column is always where a path starts, while a field reached through one is
    always qualified -- an element's ``e0.smiles`` and a reaction's own ``smiles``
    are different columns spelled alike, and matching the qualified form would
    materialize the reaction-level column for a query that reads only the element's.
    Literals are stripped for the same reason: the semi-joins carry path literals
    like ``'inputs.components'``, and a name inside one is never a column reference.

    Args:
        sql: The compiled query.

    Returns:
        The columns the query mentions, always including ``reaction_id``.
    """
    # An escaped quote inside a literal is two quotes, so this eats those too.
    outside = re.sub(r"'(?:[^']|'')*'", "''", sql)
    mentioned = {
        column
        for column in _TOP_LEVEL
        if re.search(rf"(?<![.\w]){re.escape(column)}\b", outside)
    }
    mentioned.add("reaction_id")
    return frozenset(mentioned)


def _memory_bytes(cursor: duckdb.DuckDBPyConnection) -> int:
    """Returns the bytes DuckDB holds in its in-memory tables.

    The database-wide figure, so the cost of one table is the difference across
    creating it. ``duckdb_tables`` reports a row *count* rather than a size, which
    is why this asks the memory accounting instead.

    Args:
        cursor: Any cursor on the corpus connection.

    Returns:
        Bytes held, across every table in the database.
    """
    row = cursor.execute(
        "SELECT coalesce(sum(memory_usage_bytes), 0) FROM duckdb_memory() "
        "WHERE tag = 'IN_MEMORY_TABLE'"
    ).fetchone()
    return int(row[0]) if row is not None and row[0] is not None else 0


def _run_with_timeout(
    cursor: duckdb.DuckDBPyConnection,
    sql: str,
    parameters: dict[str, Any],
    timeout_seconds: float,
) -> pa.Table:
    """Runs ``sql``, interrupting it if it outlasts ``timeout_seconds``.

    ``Timer.cancel`` only sets a flag, so a timer that has already passed its own
    check fires anyway. The lock makes the interrupt and the teardown exclusive, so
    it either lands while this query is running or not at all -- and in particular
    never after the caller closes the cursor.

    Args:
        cursor: The cursor to run on, and the one an expired timer interrupts.
        sql: The compiled query.
        parameters: Values to bind.
        timeout_seconds: Wall-clock bound.

    Returns:
        The result as an Arrow table.

    Raises:
        TimeoutError: If the query is interrupted by the timer.
    """
    lock = threading.Lock()
    running = True

    def interrupt() -> None:
        with lock:
            if running:
                cursor.interrupt()

    timer = threading.Timer(timeout_seconds, interrupt)
    timer.start()
    try:
        return cursor.execute(sql, parameters).to_arrow_table()
    except duckdb.InterruptException as error:
        raise TimeoutError(f"query exceeded {timeout_seconds} seconds") from error
    finally:
        timer.cancel()
        with lock:
            running = False


@dataclasses.dataclass(frozen=True)
class _Held:
    """Something worth keeping in memory, and the statement that builds it.

    One cache and one budget cover both kinds, so eviction can weigh a column set
    against a pivot rather than each starving the other in its own pool.

    Attributes:
        key: Identity in the cache. Carries the kind, so a pivot over a path and a
            column set naming it cannot collide.
        select: The query ``CREATE TABLE`` wraps.
        description: How the entry is named in a log line or a warning.
    """

    key: tuple[str, ...]
    select: str
    description: str


@dataclasses.dataclass
class _Narrow:
    """A materialized table holding some of the projection's top-level columns.

    Attributes:
        name: The table's name in the catalog.
        held: What it costs to keep, in bytes.
        readers: Searches currently reading it. An entry with any is passed over by
            eviction, since a search holds only the name until it runs.
    """

    name: str
    held: int
    readers: int


class Corpus:
    """A searchable pairing of projections with their structures artifacts.

    Opens one DuckDB connection over both artifact sets and publishes the relation a
    compiled query runs against. Use as a context manager, or call ``close``.
    """

    def __init__(
        self,
        projection_pattern: str,
        structures_pattern: str,
        *,
        resolver: Callable[[str], str] | None = None,
        threads: int = -1,
        require_current: bool = True,
        narrow_budget_bytes: int | None = None,
        pivots_dir: str | None = None,
    ) -> None:
        """Pairs the artifacts, checks their stamps, and prepares the relations.

        Args:
            projection_pattern: Glob matching the projection files.
            structures_pattern: Glob matching the structures artifacts.
            resolver: Maps a compound name to SMILES. Defaults to
                ``ord_schema.resolvers``, which calls external services; inject
                something local for tests or offline use.
            threads: Threads RDKit uses to match a substructure query; -1 means one
                per core. These are real threads with the GIL released, so concurrent
                searches overlap instead of serializing. The count is per search, not
                per corpus, so a server expecting several at once should divide the
                core count by that concurrency rather than leaving this at -1.
            require_current: Refuse artifacts not written by the current library,
                artifact, and RDKit versions. The screen's completeness guarantee
                assumes the query and the artifact fingerprint identically, so turning
                this off is only for reading a corpus known to match anyway.
            narrow_budget_bytes: How much memory the materialized column sets may
                hold between them; 4 GB by default. Worth stating on a machine that can
                spend more, since a column set the budget refuses is read from the
                Parquet files on every query that names it -- over ORD, seconds against
                tenths of a second. The whole projection is 18.46 GB in memory, so a
                corpus asked about most of it wants most of that. Zero materializes
                nothing at all, which is what bounds a small machine: a set is measured
                by building it, so any other figure has a peak of itself plus the
                largest set a query names, whether or not that set is then kept.
            pivots_dir: Directory holding derived pivot artifacts, one subdirectory per
                repeated level. Given one, a quantifier over a level reads the artifact
                instead of unnesting the projection to build it, which is minutes per
                level over ORD and is paid by whichever query asks first. A level with
                no subdirectory is still built in process, so a partial set of
                artifacts is a partial speedup rather than a missing answer.

        Raises:
            PairingError: If either pattern matches nothing, a file on one side has no
                partner on the other, a pair's stamps disagree about the source
                dataset, or -- with ``require_current`` -- any artifact is stale.
            ValueError: If ``narrow_budget_bytes`` is negative, which no amount of
                memory is; zero is allowed, and materializes nothing.
        """
        self._resolver = resolver if resolver is not None else _resolve_with_resolvers
        if narrow_budget_bytes is not None and narrow_budget_bytes < 0:
            raise ValueError(
                f"narrow_budget_bytes is {narrow_budget_bytes}, which is not an amount "
                "of memory; pass zero to materialize nothing"
            )
        self._narrow_budget = (
            narrow_budget_bytes
            if narrow_budget_bytes is not None
            else _NARROW_BUDGET_BYTES
        )
        # Said at open because the alternative is saying it nowhere: a process killed
        # for holding too much leaves no log of its own, and what it thought it was
        # allowed to hold is the first thing worth knowing afterwards.
        logger.info(
            "materialized column sets may hold %s", _format_bytes(self._narrow_budget)
        )
        self._threads = threads
        # Built on first substructure query; see _library. The library holds one entry
        # per distinct molecule, and _members with _starts map an entry back to the
        # structure IDs sharing it.
        self._substructure_library: rdSubstructLibrary.SubstructLibrary | None = None
        self._members = array.array("I")
        self._starts = array.array("I")
        self._library_lock = threading.Lock()
        # Built on the first query that can use it; see _occurrences. A build that finds
        # the corpus unindexable keeps its reason here, since the corpus cannot change
        # while this object is open and rebuilding costs a pass per path to reach it.
        self._occurrences_built = False
        self._refused: str | None = None
        self._occurrences_lock = threading.Lock()
        # Recent match sets, most recently used last; see _matches.
        self._matched: collections.OrderedDict[_MatchKey, str] = (
            collections.OrderedDict()
        )
        # Match sets a thread is computing, so callers asking the same question while it
        # runs wait for that answer instead of repeating it; see _matches.
        self._matching: dict[_MatchKey, threading.Event] = {}
        self._matches_lock = threading.Lock()
        # Materialized column sets, most recently used last; see _narrowed_table. The
        # serial names them, and never repeats, so a build that fails partway leaves no
        # name for the next one to collide with.
        self._narrowed: collections.OrderedDict[tuple[str, ...], _Narrow] = (
            collections.OrderedDict()
        )
        # Column sets that came out larger than the whole budget, and what each cost.
        # Kept because the projection they are built from does not change while this
        # object is open, so building one again costs a pass over the corpus to reach
        # the same refusal -- and because every query naming one is slow for a reason
        # worth restating; see _warn_refused.
        self._narrow_refused: dict[tuple[str, ...], int] = {}
        self._narrow_serial = 0
        self._narrow_lock = threading.Lock()
        # Held across a materialization, so the two memory readings bracket one table
        # and no search waits on a build to be told what the cache already holds.
        self._narrow_build_lock = threading.Lock()
        pairs = _pair(projection_pattern, structures_pattern, require_current)
        self._pivots_dir = pivots_dir
        self._require_current = require_current
        # The projections a pivot artifact has to have been derived from, and the views
        # already published over the artifacts that were; see _pivot_view.
        self._projections = [pair[0] for pair in pairs]
        self._pivot_views: dict[str, str | None] = {}
        self._pivot_lock = threading.Lock()
        self._connection = duckdb.connect()
        try:
            self._total, self._searchable = self._prepare(pairs)
        except Exception:
            # Nothing else holds the connection yet, and the caller has no object to
            # close one from if __init__ does not return.
            self._connection.close()
            raise

    def _prepare(self, pairs: list[tuple[str, str]]) -> tuple[int, int]:
        """Publishes the relations, and returns the total and searchable row counts."""
        offsets = []
        total = 0
        for projected, structured in pairs:
            with pq.ParquetFile(structured) as artifact:
                count = artifact.metadata.num_rows
            # The IDs the projection carries have to land inside the partner's rows.
            # Stamps cannot see this: an artifact rederived from a rewritten
            # projection keeps the same source hash, so a short one pairs cleanly and
            # then aliases its neighbor's molecules -- in range for get_bit, wrong
            # about the chemistry, and silent. Read from footer statistics, so it
            # costs no column data.
            largest = _max_structure_id(projected)
            if largest is not None and largest >= count:
                raise PairingError(
                    f"{projected} carries structure_id {largest} but {structured} "
                    f"holds only {count} structures, so its IDs would join to another "
                    "dataset's molecules; derive the structures artifact from this "
                    "projection again"
                )
            offsets.append((projected, structured, total))
            total += count
        registered = pa.table(
            {
                "projection_filename": [offset[0] for offset in offsets],
                "structures_filename": [offset[1] for offset in offsets],
                query.STRUCTURE_OFFSET: pa.array(
                    [offset[2] for offset in offsets], type=pa.int64()
                ),
            }
        )
        # A registered object belongs to the connection that registered it, so a
        # per-search cursor raises a catalog error on any view built over one. Copying
        # it into the catalog puts the views below within reach of every cursor.
        self._connection.register("registered_offsets", registered)
        self._connection.execute(
            "CREATE TABLE structure_offsets AS SELECT * FROM registered_offsets"
        )
        self._connection.unregister("registered_offsets")
        # Inlined rather than bound because DuckDB refuses parameters in DDL. The
        # paths came from this process's own glob, and the quoting still escapes them.
        projection_files = _sql_strings(pair[0] for pair in pairs)
        structure_files = _sql_strings(pair[1] for pair in pairs)
        # S608: the fragments are this module's constants and this process's own
        # glob results, quoted with their single quotes escaped.
        self._connection.execute(
            f"""
            CREATE VIEW {query.TABLE} AS
            SELECT p.* EXCLUDE (filename), o.{query.STRUCTURE_OFFSET}
            FROM read_parquet({projection_files}, filename=true) p
            JOIN structure_offsets o ON p.filename = o.projection_filename
            """  # noqa: S608
        )
        self._connection.execute(
            f"""
            CREATE VIEW corpus_structures AS
            SELECT s.* EXCLUDE (filename),
                   (s.structure_id + o.{query.STRUCTURE_OFFSET})::BIGINT AS global_id
            FROM read_parquet({structure_files}, filename=true) s
            JOIN structure_offsets o ON s.filename = o.structures_filename
            """  # noqa: S608
        )
        # Both views join on a path that Python read and DuckDB re-resolved. Anything
        # that makes the two disagree -- a glob metacharacter in a directory name is
        # the reachable one, since read_parquet globs each element it is handed --
        # drops rows from the join rather than failing, leaving a corpus that answers
        # every query with silence. The structures side is counted against a total
        # taken from the footers; the reactions side has no such total to check
        # against, so a projection whose path failed the same way goes unnoticed here
        # and shows up only as reactions nobody can find.
        counts = self._connection.execute(
            "SELECT count(*), count(pattern_fp) FROM corpus_structures"
        ).fetchone()
        assert counts is not None  # An aggregate over any relation returns one row.
        joined, searchable = counts
        if joined != total:
            raise PairingError(
                f"the structures artifacts hold {total} rows but only {joined} joined "
                "to their offsets; a path did not survive read_parquet"
            )
        return total, searchable

    def __enter__(self) -> Self:
        """Returns the corpus itself; closing on exit is the whole protocol."""
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Closes the connection."""
        self.close()

    def close(self) -> None:
        """Closes the connection."""
        self._connection.close()

    def _query_molecule(
        self, parameter: query.StructureParameter, resolve: Callable[[str], str]
    ) -> Chem.Mol:
        """Returns the query molecule for a structure predicate.

        Args:
            parameter: The predicate to build a molecule for.
            resolve: Maps a compound name to SMILES.

        Returns:
            A pattern from SMARTS for a substructure predicate stated as one; a
            molecule from SMILES otherwise (including every resolved compound).

        Raises:
            ValueError: If a resolved compound's SMILES does not parse -- the
                resolver's contract is RDKit-readable SMILES, so this is a resolver
                bug, not a query error.
        """
        if _reads_as_smarts(parameter):
            assert parameter.pattern is not None  # What reading it as SMARTS means.
            return Chem.MolFromSmarts(parameter.pattern)  # Validated at model time.
        if parameter.pattern is not None:
            smiles = parameter.pattern
        else:
            assert parameter.compound is not None  # One of the two is always set.
            smiles = resolve(parameter.compound)
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(
                f"resolver returned unreadable SMILES {smiles!r} "
                f"for {parameter.compound!r}"
            )
        return molecule

    def _library(self) -> rdSubstructLibrary.SubstructLibrary:
        """Returns the substructure library, building it on first use.

        Built lazily because it costs seconds and gigabytes, and a corpus asked only
        for similarity or scalar queries never needs it.

        Structures are deduplicated per dataset, so a molecule used by several datasets
        occupies a row in each and would otherwise be matched once per copy. The library
        holds one entry per distinct SMILES -- 1,435,426 of the corpus's 2,016,224 rows,
        so 29% of the matching disappears -- and ``_members`` maps an entry back to
        every structure ID sharing its molecule. Two structures with one SMILES have one
        molecule and one fingerprint, both derived from that SMILES by the RDKit the
        stamps pin, so which row is read does not matter.

        Concurrent first searches serialize on the build rather than each building a
        copy, which at corpus scale is about 1.5 GB apiece. Whoever waits needs the
        finished library anyway.

        Returns:
            A library over every distinct molecule in the corpus, screened by the
            pattern fingerprints the artifact already stores.

        Raises:
            PairingError: If the corpus IDs are not exactly ``0`` to ``total - 1``, so
                that the position an ID is assumed to occupy holds another structure; if
                the library and the entry table disagree about how many distinct
                molecules there are; or if a structure holds one of its derived columns
                without the other.
        """
        with self._library_lock:
            if self._substructure_library is not None:
                return self._substructure_library
            start = time.perf_counter()
            molecules = rdSubstructLibrary.CachedMolHolder()
            patterns = rdSubstructLibrary.PatternHolder(structures.PATTERN_FP_SIZE)
            cursor = self._connection.cursor()
            # Which library entry each structure ID resolves to, in ID order. An array
            # of unsigned ints rather than a list, which at corpus scale is 8 MB against
            # something closer to 60.
            entry_of = array.array("I")
            entries: dict[str, int] = {}
            position = 0
            try:
                reader = cursor.execute(
                    "SELECT global_id, smiles, mol_binary, pattern_fp "
                    "FROM corpus_structures ORDER BY global_id"
                ).to_arrow_reader(_BUILD_BATCH)
                for batch in reader:
                    identifiers = batch.column("global_id").to_pylist()
                    smiles_values = batch.column("smiles").to_pylist()
                    blobs = batch.column("mol_binary").to_pylist()
                    fingerprints = batch.column("pattern_fp").to_pylist()
                    for global_id, smiles, blob, fingerprint in zip(
                        identifiers, smiles_values, blobs, fingerprints, strict=True
                    ):
                        # An ID is read as the position it occupies in this pass, which
                        # it is only while the IDs are every integer from zero, each
                        # once. A duplicate or a gap slides every later structure one
                        # position along, which a row count cannot see.
                        if global_id != position:
                            raise PairingError(
                                f"the corpus states structure ID {global_id} where "
                                f"{position} was expected, so its IDs are not one "
                                "unbroken run and every later structure sits at a "
                                "position other than its own ID; derive the structures "
                                "artifacts again"
                            )
                        # A structures artifact writes the derived columns together or
                        # not at all. One without the other means the row came from
                        # somewhere else, and taking either branch would be a guess:
                        # dropping a live fingerprint makes the structure unmatchable
                        # while it still counts as searchable, and a missing one
                        # reaches RDKit as a type error naming no row.
                        if (blob is None) != (fingerprint is None):
                            raise PairingError(
                                f"structure {global_id} has "
                                f"{'no ' if blob is None else 'a '}molecule but "
                                f"{'no ' if fingerprint is None else 'a '}fingerprint; "
                                "a structures artifact writes both or neither, so this "
                                "one was not written by this library"
                            )
                        entry = entries.get(smiles)
                        if entry is None:
                            entry = entries[smiles] = len(entries)
                            if blob is None:
                                molecules.AddBinary(_UNPARSEABLE)
                                patterns.AddFingerprint(_NO_BITS)
                            else:
                                molecules.AddBinary(blob)
                                patterns.AddFingerprint(
                                    DataStructs.CreateFromBinaryText(fingerprint)
                                )
                        entry_of.append(entry)
                        position += 1
            finally:
                cursor.close()
            library = rdSubstructLibrary.SubstructLibrary(molecules, patterns)
            if position != self._total:
                raise PairingError(
                    f"the library covers {position} structures but the corpus has "
                    f"{self._total}, so a structure ID would fall outside the run this "
                    "read and resolve to no entry"
                )
            # The holders and the entry table are filled in the same branch, so a
            # divergence means one of them missed a row -- after which every entry above
            # the divergence maps to a neighbor's structures: in range, wrong molecule,
            # and invisible to a count of either one alone.
            if len(library) != len(entries):
                raise PairingError(
                    f"the library holds {len(library)} molecules against "
                    f"{len(entries)} distinct SMILES, so an entry index does not name "
                    "the molecule it was built from"
                )
            self._members, self._starts = _group(entry_of, len(entries))
            logger.info(
                "built a substructure library over %d distinct molecules covering "
                "%d structures in %.1fs",
                len(library),
                position,
                time.perf_counter() - start,
            )
            self._substructure_library = library
            return self._substructure_library

    def _occurrences(self) -> None:
        """Builds the occurrence index, once.

        One row per structure occurrence, carrying the corpus-wide ID, the path it sat
        at, and the element's own ``reaction_role``. Keeping the role beside the
        structure is what preserves element binding: "pyridine as the solvent" stays a
        condition on one row rather than an intersection of two reaction sets, which
        over-returns.

        Costs one pass over the projections per indexed path, so it is built when a
        query first wants it rather than at open. Concurrent first queries serialize
        here and share the result.

        Raises:
            PairingError: If the corpus states a reaction more than once, or the index
                does not reach every structure the corpus holds. An unreached structure
                is one whose reactions no routed query can find, and the answer comes
                back empty rather than wrong -- which reads like an answer rather than
                like a corpus that cannot be searched. The reason is kept and re-raised,
                since a corpus does not change under an open ``Corpus`` and rebuilding
                is several passes to reach the same refusal.
        """
        with self._occurrences_lock:
            if self._refused is not None:
                raise PairingError(self._refused)
            if self._occurrences_built:
                return
            start = time.perf_counter()
            # A path's elements are already a list, so unnest yields one row each. The
            # offset comes from the row's own file, which the reactions view carries.
            selects = "\nUNION ALL\n".join(
                f"""
                SELECT reaction_id, {_sql_string(path)} AS path,
                       (element.structure_id + {query.STRUCTURE_OFFSET})::UINTEGER
                           AS global_id,
                       element.{_INDEXED_FIELD} AS {_INDEXED_FIELD}
                FROM {query.TABLE}, unnest({expression}) AS unnested(element)
                WHERE element.structure_id IS NOT NULL
                """  # noqa: S608
                for path, expression in INDEXED_PATHS.items()
            )
            # Its own cursor: this runs while searches are in flight, and the shared
            # connection holds their results.
            cursor = self._connection.cursor()
            try:
                # An occurrence names the reaction it sat in by ID, so the semi-join
                # answers for every row carrying that ID. Two files stating the same
                # reaction pair cleanly -- what pairing checks is that each projection
                # has its structures, not that the corpus states a reaction once -- and
                # each one's structures would then answer the other's queries. Costs a
                # scan of one column against the several this build already runs.
                duplicated = cursor.execute(
                    f"SELECT reaction_id, count(*) FROM {query.TABLE} "  # noqa: S608
                    "GROUP BY reaction_id HAVING count(*) > 1 "
                    "ORDER BY reaction_id LIMIT 1"
                ).fetchone()
                if duplicated is not None:
                    self._refused = (
                        f"reaction {duplicated[0]} is stated by {duplicated[1]} rows "
                        "of the corpus, and the occurrence index finds reactions by "
                        "ID, so one copy's structures would answer for the other; "
                        "glob one artifact per reaction"
                    )
                    raise PairingError(self._refused)
                # Under the build lock, so the only other thing that puts a table
                # in this database does not do it while a materialization is measuring
                # what one costs -- a difference of two readings taken across everyone.
                # OR REPLACE, so a build interrupted after the table exists is a build
                # the next query repeats rather than one that collides with itself
                # forever. S608: the fragments are this module's own schema walk and
                # the compiler's traversals, not anything a query supplies.
                with self._narrow_build_lock:
                    cursor.execute(f"CREATE OR REPLACE TABLE occurrences AS {selects}")
                counts = dict(
                    cursor.execute(
                        "SELECT path, count(*) FROM occurrences GROUP BY path"
                    ).fetchall()
                )
                reached = cursor.execute(
                    "SELECT count(DISTINCT global_id) FROM occurrences"
                ).fetchone()
                assert reached is not None  # An aggregate returns one row.
                # A structures artifact holds one row per distinct structure its
                # projection uses, so every ID it states is one some element carries and
                # the index has to find all of them. Counting them is what catches a
                # single traversal reaching nothing, which a total over every path
                # cannot: the other paths carry the total and the dead one answers its
                # queries with silence.
                if reached[0] != self._total:
                    self._refused = (
                        f"the occurrence index reached {reached[0]} of the corpus's "
                        f"{self._total} structures over {sorted(INDEXED_PATHS)}; "
                        "either the projections are not the schema this walk was built "
                        "from, or their rows did not survive the filename join, so the "
                        "reactions holding the rest cannot be found"
                    )
                    raise PairingError(self._refused)
                self._occurrences_built = True
            finally:
                cursor.close()
            # Per path, so how the occurrences are distributed is visible rather than
            # only their total. A path holding none is normal -- most corpora have no
            # authentic standards -- which is why the refusal above counts structures
            # rather than paths.
            logger.info(
                "indexed %d structure occurrences in %.1fs: %s",
                sum(counts.values()),
                time.perf_counter() - start,
                ", ".join(f"{path} {counts.get(path, 0)}" for path in INDEXED_PATHS),
            )

    def _substructure_ids(
        self, parameter: query.StructureParameter, resolve: Callable[[str], str]
    ) -> list[int]:
        """Screens and verifies a substructure predicate; returns global IDs.

        RDKit screens and matches in one call, across its own threads with the GIL
        released, so concurrent searches overlap here instead of serializing.

        Args:
            parameter: The predicate to match.
            resolve: Maps a compound name to SMILES.

        Returns:
            The corpus-wide structure IDs whose molecules contain the query.
        """
        molecule = self._query_molecule(parameter, resolve)
        library = self._library()
        # maxResults defaults to 1000, which would silently truncate: a broad pattern
        # matches hundreds of thousands of ORD's distinct molecules.
        matched = library.GetMatches(
            molecule, numThreads=self._threads, maxResults=len(library) or 1
        )
        # A library entry is a molecule, not a structure: every ID sharing that molecule
        # matched too, and the answer is stated in IDs.
        identifiers: list[int] = []
        for entry in matched:
            identifiers.extend(
                self._members[self._starts[entry] : self._starts[entry + 1]]
            )
        return identifiers

    def _similarity_ids(
        self,
        cursor: duckdb.DuckDBPyConnection,
        parameter: query.StructureParameter,
        resolve: Callable[[str], str],
    ) -> list[int]:
        """Screens a similarity predicate; the fingerprint is the whole answer."""
        molecule = self._query_molecule(parameter, resolve)
        fingerprint = structures.morgan_fingerprint(molecule)
        blob = DataStructs.BitVectToBinaryText(fingerprint)
        popcount = fingerprint.GetNumOnBits()
        threshold = parameter.threshold
        assert threshold is not None  # A similarity predicate always carries one.
        # Tanimoto >= t bounds the candidate popcount to [t * q, q / t], which the
        # artifact's morgan_popcount column answers without touching the fingerprints.
        rows = cursor.execute(
            """
            SELECT global_id FROM corpus_structures
            WHERE morgan_fp IS NOT NULL
              AND morgan_popcount BETWEEN $lo AND $hi
              AND bit_count(CAST(morgan_fp AS BITSTRING)
                            & CAST($q AS BITSTRING))::DOUBLE
                  / ($n + morgan_popcount
                     - bit_count(CAST(morgan_fp AS BITSTRING) & CAST($q AS BITSTRING)))
                  >= $t
            """,
            {
                "q": blob,
                "n": popcount,
                "lo": math.ceil(threshold * popcount),
                "hi": math.floor(popcount / threshold),
                "t": threshold,
            },
        ).fetchall()
        return [row[0] for row in rows]

    def _bitmap(self, matched: Sequence[int]) -> str:
        """Returns the match set as a bitmap over the corpus-wide ID space."""
        bits = bytearray(b"0" * max(self._total, 1))
        for global_id in matched:
            bits[global_id] = ord("1")
        return bits.decode()

    def _materialize(self, held: _Held) -> _Narrow | None:
        """Returns the table ``held`` describes, building it if the cache lacks one.

        The caller owns a read on whatever comes back and has to release it, which
        ``_narrowed_table`` does; an entry with a reader is never evicted.

        What the cache holds is read under the short lock the bookkeeping takes, and
        only a build waits on another build. A search whose columns are already
        materialized is answered while one is running, which is the common case on a
        corpus serving several questions at once.

        Args:
            held: What to keep, and the statement that builds it.

        Returns:
            The entry, or None when the budget is zero or this costs more than the
            cache may hold in total. Nothing is kept either way and the projection
            answers directly; something refused as too large is remembered, so the next
            query wanting it does not build it again to throw it away again.
        """
        if not self._narrow_budget:
            # The one budget that can be settled without building, since what a set
            # costs is known only once it exists. A small container spends nothing and
            # answers from Parquet.
            return None
        with self._narrow_lock:
            settled, entry, refused = self._cached(held.key)
        if settled:
            # Warned outside the lock: it is the one every search takes to read the
            # cache, and a log handler is not something to hold it through.
            if refused is not None:
                self._warn_refused(held, refused)
            return entry
        with self._narrow_build_lock:
            # Asked again under the build lock: whoever held it may have been building
            # exactly this, and a second table of it would be memory held twice for one
            # answer.
            with self._narrow_lock:
                settled, entry, refused = self._cached(held.key)
            if settled:
                if refused is not None:
                    self._warn_refused(held, refused)
                return entry
            return self._build(held)

    def _cached(self, key: tuple[str, ...]) -> tuple[bool, _Narrow | None, int | None]:
        """Returns whether the cache settles ``key``, and how.

        Called with ``_narrow_lock`` held. A hit takes the read on the caller's behalf,
        since an entry released back to the cache between the lookup and the read could
        be evicted in between.

        Args:
            key: Identity of the entry, as ``_Held.key`` gives it.

        Returns:
            Whether the cache settles the question, the entry for a hit, and what it
            cost where it was refused as too large -- which the caller says out loud
            once it is no longer holding the lock.
        """
        cached = self._narrowed.get(key)
        if cached is not None:
            self._narrowed.move_to_end(key)
            cached.readers += 1
            return True, cached, None
        refused = self._narrow_refused.get(key)
        if refused is not None:
            return True, None, refused
        return False, None, None

    def _warn_refused(self, held: _Held, cost: int) -> None:
        """Says that a query is about to be answered the slow way, and what to change.

        Warned on every query that wants the entry rather than once when it was first
        refused, because the query is slow every time and whoever is asking why is
        reading the log now, not the one line printed when the corpus opened.

        Args:
            held: What was refused, and how to name it.
            cost: What building it took when it was measured.
        """
        logger.warning(
            "materializing %s takes %s, over this corpus's budget of %s, so every "
            "query wanting it reads the Parquet files instead -- seconds rather than "
            "tenths of a second at ORD's scale. Open the Corpus with a larger "
            "narrow_budget_bytes if the machine has the memory to spare.",
            held.description,
            _format_bytes(cost),
            _format_bytes(self._narrow_budget),
        )

    def _build(self, held: _Held) -> _Narrow | None:
        """Materializes ``held``, recording its cost or that it costs too much.

        Called with ``_narrow_build_lock`` held and ``_narrow_lock`` free, so the two
        memory readings bracket this table and nothing else: every other statement that
        puts a table in this database -- another materialization, or the index build --
        takes the same lock.

        Args:
            held: What to keep, and the statement that builds it.

        Returns:
            The entry, held once for the caller, or None if it exceeds the budget.
        """
        with self._narrow_lock:
            # Never reused, so a build that fails between the CREATE and the entry
            # cannot leave a name the next attempt collides with.
            self._narrow_serial += 1
            name = f"narrow_{self._narrow_serial}"
        start = time.perf_counter()
        # Its own cursor, for the same reason the index build takes one: searches
        # are in flight, and the shared connection holds their results.
        cursor = self._connection.cursor()
        try:
            # Said before the pass rather than after it, because the pass is what
            # takes the time and whoever is watching a query hang wants to know what it
            # is waiting for while it waits.
            logger.info("building %s", held.description)
            before = _memory_bytes(cursor)
            # The select was assembled by whoever built the _Held, which is where
            # its fragments are accounted for.
            cursor.execute(f"CREATE TABLE {name} AS {held.select}")
            cost = max(_memory_bytes(cursor) - before, 0)
            if cost > self._narrow_budget:
                # Too big to keep, and keeping it is the only reason to build it.
                cursor.execute(f"DROP TABLE {name}")
                with self._narrow_lock:
                    self._narrow_refused[held.key] = cost
                self._warn_refused(held, cost)
                return None
        except Exception:
            # A table nobody tracks is memory nobody frees, and the failure the
            # caller sees has to be the one that happened -- so the cleanup swallows
            # everything, and says what it swallowed. A DROP that fails here leaves
            # memory the cache will never account for, which is worth a line even
            # though it is not what the caller asked about.
            try:
                cursor.execute(f"DROP TABLE IF EXISTS {name}")
            except Exception:
                logger.exception("could not drop %s after a failed build", name)
            raise
        finally:
            cursor.close()
        logger.info(
            "materialized %s as %s, %s in %.1fs",
            held.description,
            name,
            _format_bytes(cost),
            time.perf_counter() - start,
        )
        entry = _Narrow(name=name, held=cost, readers=1)
        with self._narrow_lock:
            self._narrowed[held.key] = entry
            self._evict()
        return entry

    def _evict(self) -> None:
        """Drops materialized tables, least recently used first, down to the budget.

        Called with ``_narrow_lock`` held. An entry a search is reading is passed over
        rather than dropped: the search bound the table's name into its SQL and reads it
        only after resolving names and matching structures, so dropping it there would
        fail a query that had already been answered correctly. Passing over every
        candidate leaves the cache above its budget until a reader finishes, which is
        the direction that keeps answers right.
        """
        for key in list(self._narrowed):
            if (
                sum(entry.held for entry in self._narrowed.values())
                <= self._narrow_budget
            ):
                return
            entry = self._narrowed[key]
            if entry.readers:
                continue
            cursor = self._connection.cursor()
            try:
                cursor.execute(f"DROP TABLE {entry.name}")
            except duckdb.Error:
                # Bookkeeping, not the caller's query: the entry goes either way, since
                # keeping one whose table may be gone would answer from a name that
                # does not resolve.
                logger.exception("could not drop the materialized %s", entry.name)
            finally:
                cursor.close()
            del self._narrowed[key]
            logger.info("evicted the materialized %s", entry.name)

    @contextlib.contextmanager
    def _narrowed_table(self, columns: frozenset[str]) -> Iterator[str | None]:
        """Yields a table holding only ``columns``, held against eviction while read.

        Reading a handful of columns out of a 442-leaf projection spread over 53 files
        costs mostly per-file overhead, which is the same whatever the query asks for.
        Paying it once and answering from memory afterwards is worth a materialization
        for the columns a corpus is actually asked about: a temperature filter falls
        from 1.24s to 0.21s, a group-by on stirring type from 0.75s to 0.003s.

        The rows are the same rows; their order is not. Neither relation orders anything
        a query did not ask to be ordered, and they do not agree on the accident, so a
        caller wanting a particular order has to say so -- and a caller wanting a
        particular *subset*, which is what a ``limit`` with no ordering asks for, is
        asking a question neither relation answers the same way twice.

        Args:
            columns: Top-level projection columns the query names.

        Yields:
            The table's name, or None when nothing is materialized -- a budget of zero,
            or a set costing more memory than the whole cache is allowed -- in which
            case the projection answers directly.
        """
        selected = ", ".join([*sorted(columns), query.STRUCTURE_OFFSET])
        entry = self._materialize(
            _Held(
                key=("columns", *sorted(columns)),
                # S608: the names come from the projection schema and this module.
                select=f"SELECT {selected} FROM {query.TABLE}",  # noqa: S608
                description=str(sorted(columns)),
            )
        )
        if entry is None:
            yield None
            return
        try:
            yield entry.name
        finally:
            with self._narrow_lock:
                entry.readers -= 1

    def _pivot_view(self, path: str) -> str | None:
        """Returns a view over the pivot artifacts for ``path``, or None if none exist.

        Published once per level and remembered, including the answer that there is
        nothing to publish, so a corpus without artifacts globs the directory once
        rather than on every query.

        Args:
            path: The repeated level, as the query grammar names it.

        Returns:
            The view's name, or None to leave the level to be built in process.

        Raises:
            PairingError: If the artifacts for this level were not derived from the
                projections this corpus reads -- a pivot of another corpus answers
                confidently and wrongly, since its reaction IDs resolve against the
                wrong rows -- or, with ``require_current``, if any of them is stale.
        """
        if self._pivots_dir is None:
            return None
        with self._pivot_lock:
            if path in self._pivot_views:
                return self._pivot_views[path]
            # Recursive, because a pivot tree mirrors the projections it was derived
            # from: a sharded corpus puts them under <level>/<shard>/ rather than
            # directly under the level.
            files = sorted(
                glob.glob(f"{self._pivots_dir}/{path}/**/*.parquet", recursive=True)
            )
            if not files:
                self._pivot_views[path] = None
                return None
            self._check_pivots(path, files)
            name = pivot.table_name(path)
            # S608: the name comes from this module's walk of the schema, and the paths
            # from this process's own glob, quoted with their single quotes escaped.
            self._connection.execute(
                f"CREATE OR REPLACE VIEW {name} AS "  # noqa: S608
                f"SELECT * FROM read_parquet({_sql_strings(files)})"
            )
            logger.info("read %d pivot artifacts for %s", len(files), path)
            self._pivot_views[path] = name
            return name

    def _check_pivots(self, path: str, files: list[str]) -> None:
        """Refuses pivot artifacts that do not belong to this corpus.

        A pivot names its reactions by ID, and an ID resolves against whatever rows the
        projections hold. Artifacts derived from another corpus therefore answer a
        quantifier confidently and wrongly rather than failing, which is why the check
        is on the source datasets rather than on the file count.

        Args:
            path: The repeated level the files claim to hold.
            files: The artifacts found for it.

        Raises:
            PairingError: If a file is not a pivot over ``path``, if the set of source
                datasets differs from the projections', or -- with ``require_current``
                -- if any artifact is stale.
        """
        wanted = {artifacts.load_stamps(name).source_md5 for name in self._projections}
        found = set()
        for name in files:
            stamps = artifacts.load_stamps(name)
            if stamps.artifact != pivot.ARTIFACT:
                raise PairingError(
                    f"{name} is a {stamps.artifact}, not a {pivot.ARTIFACT}"
                )
            held = pivot.pivot_path(name)
            if held != path:
                raise PairingError(
                    f"{name} holds the pivot over {held}, but sits where {path} is "
                    "read from, so a quantifier would be answered by the wrong level"
                )
            if self._require_current and not artifacts.stamps_are_current(
                stamps, pivot.ARTIFACT
            ):
                raise PairingError(f"{name} is a stale {pivot.ARTIFACT}")
            found.add(stamps.source_md5)
        if found != wanted:
            raise PairingError(
                f"the {pivot.ARTIFACT} artifacts for {path} were derived from "
                f"{len(found)} source datasets and the projections from {len(wanted)}, "
                f"and {len(wanted - found)} of the projections' are missing; a pivot "
                "of another corpus names reactions this one does not hold"
            )

    @contextlib.contextmanager
    def _pivoted_table(self, path: str) -> Iterator[str | None]:
        """Yields a pivot over ``path``, held against eviction while read.

        One row per element of the level, carrying the ordinals that say which element
        it was and the element's own fields with the repeated ones removed. A quantifier
        answered from it is a semi-join rather than a decode of the whole nested column,
        which over ORD is milliseconds against seconds.

        Building it costs a pass that unnests every level down to this one, so the first
        query wanting a path waits for it. That cost belongs in a derived artifact
        rather than in a query, and until there is one it is announced rather than
        hidden.

        Args:
            path: The repeated level, as the query grammar names it.

        Yields:
            The table's name, or None when the schema has no such level or the budget
            refuses it -- in which case the quantifier compiles over the elements.
        """
        level = pivot.LEVELS.get(path)
        if level is None:
            yield None
            return
        published = self._pivot_view(path)
        if published is not None:
            # A view over artifacts costs no memory to hold and nothing to release, so
            # there is no reader to take and no budget to spend.
            yield published
            return
        entry = self._materialize(
            _Held(
                key=("pivot", path),
                select=pivot.select(level, query.TABLE),
                description=f"the pivot over {path}",
            )
        )
        if entry is None:
            yield None
            return
        try:
            yield entry.name
        finally:
            with self._narrow_lock:
                entry.readers -= 1

    def _matches(
        self,
        cursor: duckdb.DuckDBPyConnection,
        parameter: query.StructureParameter,
        resolve: Callable[[str], str],
    ) -> str:
        """Returns a structure predicate's match set as a bitmap, reusing a recent one.

        The chemistry depends on the query molecule, the operation, and the threshold,
        so the same predicate asked twice has the same answer -- and asking it costs a
        pass over every molecule in the corpus. A compound is keyed by what it resolved
        to rather than by its name, so two names for one molecule share an entry and a
        resolver that starts answering differently is a different key.

        One question is answered once even while it is still being answered: a caller
        arriving during another thread's pass waits for that pass rather than starting
        a second, so a burst of identical requests costs one. The wait is per predicate,
        so unrelated searches still overlap. A pass that raises publishes nothing and a
        waiter behind it takes the matching over, which leaves the error with the caller
        that hit it rather than with everyone queued behind.

        Args:
            cursor: The cursor a similarity screen runs on.
            parameter: The predicate to evaluate.
            resolve: Maps a compound name to SMILES.

        Returns:
            The bitmap over corpus-wide structure IDs.

        Raises:
            ValueError: If the predicate names neither a pattern nor a compound, or if
                a resolved compound's SMILES does not parse.
            PairingError: If the library does not come out one entry per distinct
                molecule over an unbroken run of IDs.
        """
        if parameter.pattern is not None:
            pattern = parameter.pattern
        elif parameter.compound is not None:
            pattern = resolve(parameter.compound)
        else:
            raise ValueError(
                f"structure parameter {parameter.name!r} names neither a pattern nor a "
                "compound, which the grammar does not allow"
            )
        # Which parser reads the string is part of what the string means, so it is part
        # of the key -- and it is asked of the same function the query molecule comes
        # from, so the two cannot come apart. Resolvers answer in the Kekule form a
        # SMARTS pattern is written in, so a name and a pattern do collide.
        key = (parameter.op, _reads_as_smarts(parameter), pattern, parameter.threshold)
        while True:
            with self._matches_lock:
                cached = self._matched.get(key)
                if cached is not None:
                    # Reinserted so the eviction below drops what has gone longest
                    # unused rather than what was computed longest ago.
                    self._matched.move_to_end(key)
                    logger.info(
                        "%s %r reused a cached match set", parameter.op, pattern
                    )
                    return cached
                computing = self._matching.get(key)
                if computing is None:
                    self._matching[key] = threading.Event()
                    break
            logger.info("%s %r is already running; waiting", parameter.op, pattern)
            computing.wait()
        try:
            if parameter.op == "substructure":
                matched = self._substructure_ids(parameter, resolve)
            else:
                matched = self._similarity_ids(cursor, parameter, resolve)
            logger.info(
                "%s %r matched %d of %d structures (%d unsearchable)",
                parameter.op,
                pattern,
                len(matched),
                self._searchable,
                self._total - self._searchable,
            )
            bitmap = self._bitmap(matched)
            with self._matches_lock:
                self._matched[key] = bitmap
                while len(self._matched) > _CACHED_MATCHES:
                    self._matched.popitem(last=False)
        finally:
            # Woken after the result is published, so a waiter finds it on the next turn
            # of the loop; a failed pass leaves nothing there and the waiter recomputes.
            with self._matches_lock:
                waiters = self._matching.pop(key, None)
            if waiters is not None:
                waiters.set()
        return bitmap

    def search(
        self, request: query.Query, *, timeout_seconds: float | None = None
    ) -> pa.Table:
        """Compiles and runs a query, returning the result as an Arrow table.

        Runs on its own cursor, so concurrent searches sharing this corpus do not read
        each other's results. Two searches wait on each other at the one-time library
        and index builds, at a structure predicate one of them is already matching, and
        at a materialization, which is one lock for all column sets: the cost of a table
        is read as the difference two memory readings make, and two builds interleaved
        would each be charged the other's.

        A quantifier the occurrence index can answer becomes a condition on the row
        rather than a filter over the elements, and the rest of the query -- the clauses
        beside it, the aggregate, the ordering, the limit -- compiles exactly as it
        would have. The screening and verification are the compiler's either way, so an
        indexed query and a projection one differ only in how they reach the reactions
        holding a match. Which happened is logged, because it is the one thing about a
        search that the result does not show.

        A query whose named columns fit the cache's budget is compiled a second time
        against a table holding only those columns, which for an indexed structure
        query is little more than the reaction ID. The narrow table is held for the
        life of the search, since the query names it seconds before it reads it and
        eviction would otherwise drop it in between.

        Args:
            request: The query to run.
            timeout_seconds: Wall-clock bound on the final query. Name resolution, the
                library and index builds, screening, and verification run before the
                timer starts and are not counted, so a first substructure search can
                take seconds longer than this allows. None runs unbounded.

        Returns:
            The selected columns: ``reaction_id`` for a plain query, the group and
            measure columns for an aggregated one.

        Raises:
            query.QueryError: If the query does not compile.
            ValueError: If a compound name cannot be resolved.
            PairingError: If the corpus IDs do not come out one unbroken run, or the
                occurrence index refuses the corpus -- a reaction stated by more than
                one row, or a structure no indexed path reaches.
            TimeoutError: If the query exceeds ``timeout_seconds``.
        """
        # Records into this search rather than onto the corpus: two searches share one
        # Corpus, and whether the index answered is a fact about one of them.
        indexing = False

        def index(
            path: str, fields: dict[str, str], allocate: Callable[[], str]
        ) -> str | None:
            nonlocal indexing
            condition = _index_condition(path, fields, allocate)
            indexing = indexing or condition is not None
            return condition

        # Cached across the whole search: a compound named in both a value and a
        # structure predicate is one external lookup, not two.
        resolved: dict[str, str] = {}

        def resolve(name: str) -> str:
            if name not in resolved:
                resolved[name] = self._resolver(name)
            return resolved[name]

        with contextlib.ExitStack() as reading:
            pivoted: dict[str, str | None] = {}

            def pivot_table(path: str) -> str | None:
                # Built while compiling rather than after, because the budget decides
                # whether the table exists at all and the SQL names it either way. The
                # answer is remembered so the second compilation, against the narrow
                # relation, cannot route a quantifier differently from the first.
                if path not in pivoted:
                    pivoted[path] = reading.enter_context(self._pivoted_table(path))
                return pivoted[path]

            compiled = query.compile_query(request, index=index, pivot=pivot_table)
            if indexing:
                # Built after compiling and before running: the condition names the
                # table, so the table has to exist by the time the query does. Logged
                # after the build, so a build that refuses the corpus does not leave a
                # line claiming the index answered a query nothing ran.
                self._occurrences()
                logger.info("the occurrence index answers part of this query")
            else:
                logger.info("the projection answers this query")
            # Read off the compiled SQL, so the second compilation is the same query
            # against a different relation -- no rewriting of SQL text, which a query
            # whose own string literal named the relation would otherwise corrupt.
            narrow = reading.enter_context(
                self._narrowed_table(_mentioned(compiled.sql))
            )
            if narrow is not None:
                compiled = query.compile_query(
                    request, table=narrow, index=index, pivot=pivot_table
                )
            cursor = self._connection.cursor()
            try:
                parameters: dict[str, Any] = {
                    name: resolve(name) for name in compiled.compounds
                }
                for parameter in compiled.structures:
                    parameters[parameter.name] = self._matches(
                        cursor, parameter, resolve
                    )
                if timeout_seconds is None:
                    return cursor.execute(compiled.sql, parameters).to_arrow_table()
                return _run_with_timeout(
                    cursor, compiled.sql, parameters, timeout_seconds
                )
            finally:
                cursor.close()
