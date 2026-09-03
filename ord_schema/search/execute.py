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

"""Running a compiled query against a corpus of projections and structures base.

The compiler (:mod:`ord_schema.search.query`) leaves two kinds of work undone on
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

Reading the projection at all costs a re-parse of every file's footer, and over a
442-leaf schema in 53 files that parse is most of what a query spends however few leaves
it goes on to read. The corpus keeps the parsed footers instead; see ``_cache_footers``.

The grammar bounds what a query can cost (one pass and a sort), and ``search``
takes a wall-clock timeout over the whole call, bounding when the caller is answered
rather than when work stops; see ``_Deadline`` for which phases can be interrupted and
which are only checked. Both builds happen
at open unless the corpus was given ``warm=False``, which leaves the index to the first
query taking a clause of it and the library to the first substructure predicate. A
search that pays for either pays upwards of a minute over the whole corpus, on top of
whatever timeout it was given -- unless an ``occurrences_dir`` covers every indexed
path, where the index is a view over those files and there is nothing to build.
"""

import array
import collections
import contextlib
import dataclasses
import functools
import glob
import hashlib
import math
import pathlib
import threading
import time
from collections.abc import Callable, Iterable, Iterator, Sequence
from typing import Any, Self

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, DataStructs
from rdkit.Chem import rdSubstructLibrary

from ord_schema import resolvers
from ord_schema.artifacts import (
    base,
    occurrences,
    pivot,
    projection,
    structures,
)
from ord_schema.logging import get_logger
from ord_schema.search import query

logger = get_logger(__name__)

# Structures read per batch while building the library. Bounds peak memory during the
# build without making the per-batch Python overhead matter.
_BUILD_BATCH = 50_000

# Match sets kept for reuse. Each is a bitmap over the corpus ID space -- 2 MB at ORD's
# scale -- so this trades a few tens of megabytes for the repeat of a query that costs
# seconds.
_CACHED_MATCHES = 16

# How much memory the pivots built in process may hold between them, when the caller
# states no budget. The ones worth the most are 0.39 to 4.18 GB at ORD's scale, so this
# holds any single one of them and every budget is a partial cache; the question is only
# which fit. A corpus asked for more says so on every query it costs, naming this
# argument; see _warn_refused. A pivot read from an artifact costs nothing here, and
# answers within tens of milliseconds of one built, so a deployment holding much of this
# is one that has not derived them.
#
# Held to after a pivot is built rather than before, since what one costs is known only
# once it exists: the peak is this plus the largest single pivot, and building one costs
# its size whether or not it is kept.
_PIVOT_BUDGET_BYTES = 4 * 1024**3


# What the process holds beside DuckDB's own accounting while the occurrence index
# builds: about 1.4 GB resident before it starts -- the interpreter, RDKit, the paired
# footers -- and up to 2.8 GB more during it. A cgroup cap leaving less than this above
# the memory_limit is one the kernel ends at, and a kill arrives as neither an exception
# nor a log line. Measured in a container over the full corpus: an 8 GiB cap took
# DuckDB's default 6.3 GiB limit and was killed 85s in, while the same build under a
# 6.5 GiB limit and a 12 GiB cap finished, peaking at 9.3 GB resident.
_PROCESS_HEADROOM_BYTES = 4 * 1024**3

# Where a container states its memory cap, cgroup v2 first. Read at open rather than
# taken as given, since DuckDB sizes its default limit from the same cap and the two
# together are what has to fit. A cap can be stated at any level from the process's own
# cgroup up to the mount root, and a container in a systemd slice sits under several, so
# the path this process occupies is read from /proc rather than assumed to be the root,
# and every level of it is consulted.
_CGROUP_V2_ROOT = pathlib.Path("/sys/fs/cgroup")
_CGROUP_V1_MEMORY_ROOT = pathlib.Path("/sys/fs/cgroup/memory")
_PROC_SELF_CGROUP = pathlib.Path("/proc/self/cgroup")


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


# The predicates the artifact answers with one equality: the column each compares, and
# the function that puts a query molecule into that column's terms.
_HASH_MATCHES: dict[str, tuple[str, Callable[[Chem.Mol], str | None]]] = {
    "same_compound": ("mol_hash", structures.mol_hash),
    "same_parent": ("parent_hash", structures.parent_hash),
}


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


def _canonical(smiles: str) -> str:
    """Returns the SMILES as the corpus spells it.

    A resolver answers in whatever form its service uses -- PubChem returns pyridine
    as ``C1=CC=NC=C1`` -- while the projection stores what RDKit canonicalizes,
    ``c1ccncc1``. Comparing the two as strings matches nothing and says so by
    returning no rows, which is the worst way for a query to be wrong. Structure
    predicates are unaffected: they go through RDKit rather than string equality.

    Args:
        smiles: What the resolver returned.

    Returns:
        The canonical form, or the input unchanged if RDKit cannot parse it, since a
        name that resolves to something unparseable is the resolver's problem to report
        and not this function's to hide.
    """
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(molecule) if molecule is not None else smiles


def sql_string(value: str) -> str:
    """Returns ``value`` as a SQL string literal, single quotes escaped."""
    escaped = value.replace("'", "''")
    return f"'{escaped}'"


def _sql_paths(paths: Iterable[str]) -> str:
    """Returns ``paths`` as a SQL list literal naming each file and no other.

    read_parquet reads its arguments as globs, so a path is not simply itself: a
    directory named with a bracket becomes a character class that matches a sibling
    instead, and one with a star or a question mark matches its neighbors as well as
    itself. Either way the query answers off files nobody asked for, which is the kind
    of wrong that arrives as data rather than as an error.

    Args:
        paths: Paths to files, as this process found them.

    Returns:
        A SQL list literal, single quotes escaped and glob characters spelled as the
        one-character classes that match them literally.
    """
    return "[" + ", ".join(sql_string(glob.escape(path)) for path in paths) + "]"


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


def _index(
    pattern: str, artifact: str, schema: pa.Schema, require_current: bool
) -> dict[str, tuple[str, base.Stamps]]:
    """Returns the artifacts matching ``pattern``, keyed by source dataset.

    Keyed by the source hash rather than the filename because that is what an
    artifact *is*: two files pair when they restate the same dataset, whatever
    they are called or wherever they sit. Keying on the basename would let two
    files in different directories collapse onto one another silently.

    Args:
        pattern: Glob matching the artifact files.
        artifact: Artifact name every match must hold.
        schema: The columns this library's version of the artifact carries.
        require_current: Refuse artifacts not written by the current versions.

    Returns:
        A mapping from source dataset hash to the path and the stamps read off it. The
        stamps come back rather than being dropped because reading them opens a footer,
        and the corpus fingerprint wants exactly the same ones.

    Raises:
        PairingError: If a match holds another artifact, lacks a column this library
            reads, is stale under ``require_current``, or restates a dataset another
            match already did.
    """
    index: dict[str, tuple[str, base.Stamps]] = {}
    for path in sorted(glob.glob(pattern, recursive=True)):
        stamps = base.load_stamps(path)
        if stamps.artifact != artifact:
            raise PairingError(f"{path} is a {stamps.artifact}, not a {artifact}")
        missing = base.missing_columns(path, schema)
        if missing:
            # The stamps say nothing about a file's columns, so this is the only check
            # that sees a schema change. Without it a file predating a column is read
            # as a corpus member and fails deep in DuckDB -- as a binder error on the
            # query that wants the column, or, where some files have it and some do
            # not, as a schema mismatch that takes down every structure query.
            raise PairingError(
                f"{path} is a {artifact} artifact without {missing}, which this "
                "library reads; derive it again first"
            )
        if require_current and not base.stamps_are_current(stamps, artifact):
            raise PairingError(f"{path} is stale; derive it again first")
        if stamps.source_md5 in index:
            raise PairingError(
                f"{path} and {index[stamps.source_md5][0]} are both {artifact} "
                "artifacts of the same source dataset; which one answers a query "
                "would be arbitrary"
            )
        index[stamps.source_md5] = (path, stamps)
    return index


def _fingerprint(stamps: Iterable[base.Stamps]) -> str:
    """Returns a short digest naming the artifacts a corpus opened.

    Args:
        stamps: The stamps of every artifact in the corpus, in any order. Both halves
            of each pair belong here: the structures artifact answers every structure
            predicate and is derived separately, so a corpus whose structures were
            rebuilt can answer differently with its projections untouched.

    Returns:
        Sixteen hex characters, sorted so the same corpus opened through a different
        glob digests the same. The whole stamp goes in rather than the source hash
        alone, because two artifacts built from identical sources by different library
        versions can disagree too.
    """
    digest = hashlib.sha256()
    # Sorted as text rather than as tuples: a stamp field is optional, and comparing
    # None against a string to order two artifacts raises.
    for value in sorted(repr(dataclasses.astuple(one)) for one in stamps):
        digest.update(value.encode())
    return digest.hexdigest()[:16]


def _pair(
    projection_pattern: str, structures_pattern: str, require_current: bool
) -> tuple[list[tuple[str, str, str]], str]:
    """Returns (projection, structures, source) triples, verified by their stamps.

    Args:
        projection_pattern: Glob matching the projection base.
        structures_pattern: Glob matching the structures base.
        require_current: Refuse artifacts not written by the current versions.

    Returns:
        One triple per source dataset, ordered by the source hash, carrying that hash
        beside the paths -- every artifact derived from the dataset names it, so it is
        what pairs a pivot with the offset this corpus gave its projection. Also the
        fingerprint over the stamps this already read to verify them.

    Raises:
        PairingError: If no projection matches, if either side lacks a column this
            library reads, or if either side states a dataset the other does not, since
            a projection's IDs index its own partner's molecules and nobody else's.
    """
    projections = _index(
        projection_pattern, projection.ARTIFACT, projection.SCHEMA, require_current
    )
    structure_files = _index(
        structures_pattern, structures.ARTIFACT, structures.SCHEMA, require_current
    )
    if not projections:
        raise PairingError(f"no projections matched: {projection_pattern}")
    unpaired = projections.keys() ^ structure_files.keys()
    if unpaired:
        orphans = sorted((projections | structure_files)[key][0] for key in unpaired)
        raise PairingError(
            "projections and structures artifacts do not pair up; these have no "
            f"counterpart derived from the same source dataset: {orphans}"
        )
    pairs = [
        (projections[key][0], structure_files[key][0], key)
        for key in sorted(projections)
    ]
    fingerprint = _fingerprint(
        stamps
        for index in (projections, structure_files)
        for _, stamps in index.values()
    )
    return pairs, fingerprint


# A structure that did not parse still gets a library entry, so every structure ID
# resolves to one and the mapping needs no special case. An empty molecule fingerprints
# to no bits, and a query with no atoms is refused, so it can never match.
_UNPARSEABLE = Chem.Mol().ToBinary()
_NO_BITS = DataStructs.ExplicitBitVect(structures.PATTERN_FP_SIZE)

# The one element field the occurrence index carries besides the structure. A query
# binding anything else to the element runs against the projection instead. Taken from
# the artifact, so the column the index filters on and the column an occurrences
# artifact carries cannot come apart.
_INDEXED_FIELD = occurrences.INDEXED_FIELD


def _indexed_paths(schema: pa.Schema = projection.SCHEMA) -> dict[str, str]:
    """Returns the paths an occurrence index covers, mapped to their traversals.

    Which paths those are is ``occurrences.indexed_paths``, so the compiler routes a
    quantifier to the index exactly where an artifact carries the rows to answer it. Two
    walks of the schema would be two answers waiting to disagree, and the disagreement
    is silent: a path routed here that no artifact holds is a quantifier answered from
    nothing.

    What this adds is the traversal, which only the compiler knows: the expression
    unnesting the level out of the projection, for the paths no pivot artifact covers.

    Args:
        schema: The projection schema.

    Returns:
        Dotted query paths to the DuckDB expression yielding their elements, for every
        level whose elements carry a structure.

    Raises:
        ValueError: If a structure-bearing path does not resolve to a repeated
            expression, so the build cannot unnest it, or if
            ``occurrences.indexed_paths`` refuses the schema. Raised at import, naming
            the path, because the answer is to change one of these modules, not to
            retry.
    """
    paths: dict[str, str] = {}
    for path in occurrences.indexed_paths(schema):
        resolved = query.resolve(path, schema=schema, allow_internal=True)
        if not resolved.repeated:
            raise ValueError(
                f"{path} carries a structure but resolves to one element rather than "
                "a repeated expression, which the build cannot unnest"
            )
        paths[path] = resolved.expression
    return paths


INDEXED_PATHS = _indexed_paths()


def _occurrences_from_projection(path: str, expression: str) -> str:
    """Returns the SELECT unnesting ``path``'s occurrences out of the projection.

    Used where no pivot artifact covers the path's level. A path's elements are already
    a list, so unnest yields one row each, and the offset comes from the row's own file,
    which the reactions view carries.

    Args:
        path: The indexed path.
        expression: The traversal yielding its elements, from ``INDEXED_PATHS``.

    Returns:
        The SELECT.
    """
    # S608: the traversal is the compiler's own, resolved against the projection schema,
    # and the path is a key of INDEXED_PATHS, quoted.
    return f"""
        SELECT reaction_id, {sql_string(path)} AS path,
               (element.structure_id + {query.STRUCTURE_OFFSET})::UINTEGER AS global_id,
               element.{_INDEXED_FIELD} AS {_INDEXED_FIELD}
        FROM {query.TABLE}, unnest({expression}) AS unnested(element)
        WHERE element.structure_id IS NOT NULL
        """  # noqa: S608


def _warn_when_the_bound_was_reached(
    request: query.Query, compiled: query.Compiled, returned: int
) -> None:
    """Says when this corpus's bound, not the query's, may have cut the answer short.

    Said after the query rather than while compiling it, because a bound applied is not
    an answer cut short: a query the bound never reached is the ordinary case, and
    warning about it on every generous limit is how a reader learns to skip the line
    that matters. A result that comes back exactly full is the most a row count can say
    -- it may hold every match and end there -- and nothing in the table says which.

    An aggregated query is the louder case. A truncated list of reactions is a sample of
    the matching ones; a truncated list of groups is part of a distribution read as the
    whole of it, and where the query ordered by nothing it is an arbitrary part.

    Args:
        request: The query as asked, whose own ``limit`` is not this corpus's doing.
        compiled: What it compiled to, carrying the limit that reached the SQL.
        returned: Rows the query returned.
    """
    if compiled.limit is None or compiled.limit == request.limit:
        return
    if returned < compiled.limit:
        return
    if request.aggregate is not None:
        logger.warning(
            "this query grouped into at least the %d rows the corpus bounds it at, so "
            "the groups it returned are %s of the distribution rather than all of it; "
            "ask for fewer groups, or raise max_rows",
            compiled.limit,
            "an arbitrary part" if not request.order_by else "the leading part",
        )
    else:
        logger.warning(
            "this query returned the %d rows the corpus bounds it at, so there may be "
            "matches it did not return; raise max_rows to see them",
            compiled.limit,
        )


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
        f"occurrence.path = {sql_string(path)}",
        f"get_bit(CAST(${allocate()} AS BITSTRING), occurrence.global_id::INTEGER) = 1",
    ]
    role = fields.get(_INDEXED_FIELD)
    if role is not None:
        conditions.append(f"occurrence.{_INDEXED_FIELD} = {sql_string(role)}")
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


def _cache_footers(connection: duckdb.DuckDBPyConnection) -> None:
    """Keeps the Parquet footers parsed between queries.

    A scan re-reads and re-parses the footer of every file it touches, and the
    projection's footers are large: 442 leaves over 53 files, with statistics per row
    group. Over ORD that parse is most of what a query costs, and it is charged again on
    every query however little column data the query goes on to read -- a temperature
    filter falls from 0.76s to 0.10s once the footers are held, a group-by on stirring
    type from 0.73s to 0.055s.

    Roughly 200 MB across the whole corpus, and bounded by the files rather than by what
    is asked of them, which is what makes it worth spending unconditionally where the
    gigabytes a materialized pivot costs are weighed against a budget. DuckDB re-reads a
    file that has been rewritten, so an artifact replaced under an open corpus is read
    as it stands.

    Set globally rather than on this connection: a search runs on its own cursor, which
    is its own session, and a session-scoped setting would leave every search paying the
    parse.

    Args:
        connection: The corpus connection, on the database every cursor shares.
    """
    connection.execute("SET GLOBAL parquet_metadata_cache=true")


def _cgroup_relative_path(controller: str) -> str:
    """Returns the cgroup this process occupies, relative to that hierarchy's mount.

    Args:
        controller: A cgroup v1 controller name, or the empty string for the v2
            hierarchy, whose ``/proc/self/cgroup`` line names no controller.

    Returns:
        The path with no leading slash, empty where /proc says nothing about this
        hierarchy -- which includes every platform without cgroups.
    """
    try:
        text = _PROC_SELF_CGROUP.read_text()
    except OSError:
        return ""
    for line in text.splitlines():
        # Lines read "hierarchy:controllers:path", with controllers empty under v2.
        _, _, rest = line.partition(":")
        controllers, _, path = rest.partition(":")
        named = controllers.split(",") if controllers else []
        matched = controller in named if named else not controller
        if matched:
            return path.lstrip("/")
    return ""


def _stated_cap_bytes(path: pathlib.Path) -> int | None:
    """Returns the cap one cgroup file states, or None where it states none.

    Args:
        path: A ``memory.max`` or ``memory.limit_in_bytes`` file, which need not exist.

    Returns:
        The figure in bytes. None covers an unreadable file, ``max``, and the
        top-of-range number that is how cgroup v1 spells the same thing.
    """
    try:
        text = path.read_text().strip()
    except OSError:
        return None
    if text == "max":
        return None
    try:
        cap = int(text)
    except ValueError:
        return None
    return None if cap >= 2**60 else cap


def _container_cap_bytes() -> int | None:
    """Returns the memory this process may hold before the kernel intervenes, in bytes.

    A cap can be stated at any level between the process's own cgroup and the mount
    root, and the one that ends the process is the smallest of them, so every level is
    read rather than the root alone.

    Returns:
        The smallest cap that applies, or None where none does -- a machine outside a
        container, a cgroup hierarchy that states no limit anywhere above this process,
        or any platform without cgroups.
    """
    caps = []
    for root, filename, controller in (
        (_CGROUP_V2_ROOT, "memory.max", ""),
        (_CGROUP_V1_MEMORY_ROOT, "memory.limit_in_bytes", "memory"),
    ):
        directory = root
        directories = [root]
        for part in pathlib.PurePosixPath(_cgroup_relative_path(controller)).parts:
            directory = directory / part
            directories.append(directory)
        caps.extend(
            cap
            for cap in (_stated_cap_bytes(d / filename) for d in directories)
            if cap is not None
        )
    return min(caps) if caps else None


def _setting_bytes(value: str) -> int | None:
    """Returns a DuckDB size setting in bytes.

    Args:
        value: A setting as DuckDB prints it, e.g. ``6.3 GiB``.

    Returns:
        The figure in bytes, or None if it does not read as a size, which is what
        ``memory_limit`` says when nothing bounds it.
    """
    scales = {"B": 1, "KiB": 1024, "MiB": 1024**2, "GiB": 1024**3, "TiB": 1024**4}
    number, _, unit = value.strip().partition(" ")
    scale = scales.get(unit)
    if scale is None:
        return None
    try:
        return int(float(number) * scale)
    except ValueError:
        return None


def _warn_when_the_cap_leaves_no_headroom(
    connection: duckdb.DuckDBPyConnection,
) -> None:
    """Warns when DuckDB's limit leaves the rest of the process too little to run in.

    DuckDB reads the cgroup, so a container gets a limit sized to its own cap rather
    than to the host -- but that limit bounds DuckDB's buffers, not the process, and
    the default takes about 80% of the cap for them. What the corpus holds beside them
    is what the kernel kills for.

    Args:
        connection: The corpus connection.
    """
    cap = _container_cap_bytes()
    if cap is None:
        return
    setting = connection.execute("SELECT current_setting('memory_limit')").fetchone()
    limit = _setting_bytes(setting[0]) if setting is not None else None
    if limit is None or cap - limit >= _PROCESS_HEADROOM_BYTES:
        return
    logger.warning(
        "this container may hold %s and DuckDB's memory_limit is %s, leaving %s for "
        "everything else the process holds, which building the occurrence index has "
        "been measured to want %s of. Reaching the cap is a kill rather than an "
        "exception, so open the Corpus with a memory_limit that leaves that much, or "
        "give the container more.",
        _format_bytes(cap),
        _format_bytes(limit),
        _format_bytes(cap - limit),
        _format_bytes(_PROCESS_HEADROOM_BYTES),
    )


class _Deadline:
    """When a search has to be finished, and what it may still spend.

    A search is several phases -- compiling, which may build a pivot; the occurrence
    index; name resolution; screening and verifying a structure predicate; and the query
    itself -- and only the last is a DuckDB statement a timer can interrupt. The rest is
    RDKit with the GIL released, an external resolver, or a build another search waits
    on. So this bounds when the caller gets an answer or an error rather than when work
    stops: each phase is checked as it finishes, and the query gets whatever is left.

    A shared build is deliberately not interrupted. The occurrence index, the
    substructure library, and a materialized pivot are built once for the corpus under a
    lock, so a timer that killed one would fail every search queued behind it -- callers
    that asked for no bound included.
    """

    def __init__(self, seconds: float | None) -> None:
        """Starts the clock.

        Args:
            seconds: Wall-clock bound, or None for no deadline.
        """
        self._at = None if seconds is None else time.monotonic() + seconds
        self._seconds = seconds

    def remaining(self) -> float | None:
        """Returns the seconds left, or None where there is no deadline."""
        if self._at is None:
            return None
        return self._at - time.monotonic()

    def check(self, phase: str) -> None:
        """Raises if the deadline has passed.

        Args:
            phase: What has just finished, named in the error so the caller learns
                which part of the search was the slow one.

        Raises:
            TimeoutError: If the deadline has passed.
        """
        left = self.remaining()
        if left is not None and left <= 0:
            raise TimeoutError(
                f"the search exceeded {self._seconds} seconds while {phase}; it ran to "
                "completion because that phase cannot be interrupted"
            )


def _run_with_timeout(
    cursor: duckdb.DuckDBPyConnection,
    sql: str,
    parameters: dict[str, Any],
    timeout_seconds: float,
) -> pa.Table:
    """Runs ``sql``, interrupting it if it outlasts ``timeout_seconds``.

    The one phase of a search that is interruptible, and the only one whose work is this
    caller's alone: the cursor is the search's own, so interrupting it fails nothing
    else.

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
        TimeoutError: If the query is interrupted by the timer, if ``timeout_seconds``
            has already run out before it starts, or if the query answers after the
            bound because its interrupt was lost.
    """
    # Refused rather than armed. A timer set for a bound already spent fires while the
    # cursor is still idle, and DuckDB clears an interrupt nothing was running when it
    # arrived -- so the query would then run unbounded and answer after the deadline it
    # was given. Decided from the same value that would arm the timer, so no time passes
    # between the test and the arming for the bound to expire in.
    if timeout_seconds <= 0:
        raise TimeoutError(
            "the search spent its whole bound before the query could start"
        )
    lock = threading.Lock()
    running = True

    def interrupt() -> None:
        with lock:
            if running:
                cursor.interrupt()

    timer = threading.Timer(timeout_seconds, interrupt)
    started = time.monotonic()
    timer.start()
    try:
        answered = cursor.execute(sql, parameters).to_arrow_table()
    except duckdb.InterruptException as error:
        raise TimeoutError(f"query exceeded {timeout_seconds} seconds") from error
    finally:
        timer.cancel()
        with lock:
            running = False
    # The interrupt is not guaranteed to land. A timer armed for a very short bound can
    # fire while the cursor is still idle, and DuckDB clears an interrupt that arrived
    # with nothing running -- measured, a bound of a nanosecond loses it about one run
    # in five -- after which the query goes on unbounded. So the answer is checked
    # against the clock rather than trusted because no interrupt arrived: the bound is a
    # promise about when the caller is answered, and one kept only when a race goes the
    # right way is not a promise.
    elapsed = time.monotonic() - started
    if elapsed > timeout_seconds:
        raise TimeoutError(
            f"the query answered after {elapsed:.3f} seconds, past the "
            f"{timeout_seconds} it was given; its interrupt did not land"
        )
    return answered


@dataclasses.dataclass
class _Pivot:
    """A pivot over one repeated level, materialized in memory.

    Attributes:
        name: The table's name in the catalog.
        held: What it costs to keep, in bytes.
        readers: Searches currently reading it. A pivot with any is passed over by
            eviction, since a search holds only the name until it runs.
    """

    name: str
    held: int
    readers: int


class Corpus:
    """A searchable pairing of projections with their structures base.

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
        pivot_budget_bytes: int | None = None,
        pivots_dir: str | None = None,
        derive_pivots: bool = False,
        require_pivots: bool = False,
        occurrences_dir: str | None = None,
        require_occurrences: bool = False,
        memory_limit: str | None = None,
        warm: bool = True,
        max_rows: int | None = None,
    ) -> None:
        """Pairs the artifacts, checks their stamps, and prepares the relations.

        Args:
            projection_pattern: Glob matching the projection files.
            structures_pattern: Glob matching the structures base.
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
            pivot_budget_bytes: How much memory the pivots built in process may hold
                between them; 4 GB by default. Worth stating on a machine that can spend
                more, since a quantifier over a level the budget refuses is answered by
                unnesting the projection, which is seconds rather than milliseconds. The
                four levels that answer the most cost 9.21 GB built, so a corpus asked
                about all of them and given no artifacts wants most of that. Zero
                materializes nothing at all, which is what bounds a small machine: a
                pivot is measured by building it, so any other figure has a peak of
                itself plus the largest pivot a query wants, kept or not.
            pivots_dir: Directory holding derived pivot artifacts, one subdirectory per
                repeated level. Given one, a quantifier over a level reads the artifact
                instead of unnesting the projection to build it, which is minutes per
                level over ORD and is paid by whichever query asks first. Worth pointing
                at wherever a deployment can: the four levels that answer the most are
                514 MB as artifacts against 9.21 GB built, and answer within tens of
                milliseconds of the built ones. A level with no subdirectory is still
                built in process, so a partial set of artifacts is a partial speedup
                rather than a missing answer.
            require_pivots: Refuse a corpus whose ``pivots_dir`` does not hold every
                level, rather than building the missing ones in memory on the query
                that first wants one. Needs a ``pivots_dir``. What it costs is a footer
                read per artifact at startup; what it buys is that no query pays for an
                unnest nobody meant to run. It leaves ``pivot_budget_bytes`` with
                nothing to govern -- every level is read rather than built -- so the two
                together are a contradiction rather than a belt and braces.
            derive_pivots: Write the artifact for a level that has none, rather than
                building the level in memory. Needs a writable ``pivots_dir``. The pass
                costs what building in memory costs, and once, since the next process
                over the same directory reads what this one wrote; it also holds a
                projection at a time rather than the whole level, so a level too large
                for the budget is still derivable. Off by default: a deployment reading
                artifacts someone else derives should not start writing them because a
                directory was mistyped, and a corpus-scale derivation belongs in
                ``ord_schema.artifacts.scripts.derive_pivots``. Needs ``warm`` off,
                since warming reads the levels the occurrence index covers and a
                derivation has to precede the read of the level it completes.
            occurrences_dir: Directory holding derived occurrence artifacts, one
                subdirectory per indexed path. Given one covering every path, the
                occurrence index is a view over those files rather than a table built
                at open: the corpus stops holding the 1.19 GB the table costs over ORD,
                stops wanting the 5-6.5 GB of DuckDB memory and 16-25 GB of temporary
                files to build it, and pays instead about 1.28x per query -- less on
                the narrow paths, where reading one path's files beats filtering across
                every occurrence in the corpus. A path with no artifacts is read from a
                pivot or unnested from the projection as before, and one uncovered path
                materializes the whole index, since a view would repeat that traversal
                on every query; ``require_occurrences`` refuses that rather than
                falling back to it. Derived by
                ``ord_schema.artifacts.scripts.derive_occurrences``.
            require_occurrences: Refuse a corpus whose ``occurrences_dir`` does not
                cover every indexed path, rather than materializing the index because
                one path was missing. Needs an ``occurrences_dir``. What it costs is a
                footer read per artifact at startup; what it buys is that a deployment
                sized for the view is not handed the table, which wants several
                gigabytes it was not given and takes the container down rather than
                answering slowly.
            memory_limit: What DuckDB may hold, spelled as it spells sizes: ``6500MiB``,
                ``8GB``. Left unset, DuckDB takes about 80% of the machine, or of the
                container's cap, which it does read. That default suits a process that
                is mostly DuckDB, and this one is not: RDKit, the interpreter, and the
                index build's own allocations sit outside that accounting, and together
                have been measured at up to 2.8 GB during a build. Exceeding a cap is a
                kill rather than an exception, so a container wants this set low enough
                to leave that much. Over ORD, ``6500MiB`` under a 12 GiB cap builds the
                occurrence index, where DuckDB's own default under an 8 GiB cap is
                killed partway.
            warm: Build at open what a structure query would otherwise build on the
                query itself, on by default: the occurrence index, then the substructure
                library. An ``occurrences_dir`` covering every indexed path leaves the
                index nothing to build, so what this then pays for it is the scan that
                checks it reaches every structure the corpus holds. Built, the index
                costs a read of the pivot artifacts covering the indexed paths and a
                pass over the projections for every path they do not
                cover -- over ORD the difference between a second and a minute -- and it
                wants 5-6.5 GB of DuckDB memory to build and holds 1.19 GB afterwards.
                The library is about 8s and 2.2 GB on top. Together they are most of
                what a container is sized against, which is the reason to spend them
                where falling short is a failed deployment rather than a failed request;
                see ``check_index``. A corpus asked only scalar queries needs neither,
                and one asked scalar and similarity queries needs the index without the
                library -- ``warm=False`` beside ``check_index`` builds that pair.
            max_rows: Most rows a search may return, applied whether or not the query
                asked for a limit and clamping one that asks for more. The grammar
                leaves ``limit`` optional, so a query without one returns every matching
                reaction -- at ORD's scale millions of rows, built into an Arrow table
                in this process, from a query that reads as ordinary. A caller serving
                queries it did not write wants this set. None leaves every query's own
                limit, or none, exactly as stated.

        Raises:
            PairingError: If either pattern matches nothing, a file on one side has no
                partner on the other, a pair's stamps disagree about the source
                dataset, or -- with ``require_current`` -- any artifact is stale. With
                ``warm``, also whatever building the occurrence index and the
                substructure library refuse: a corpus stating a reaction twice, one
                whose index does not reach every structure it holds, one whose IDs are
                not every integer from zero, or pivot artifacts over an indexed level
                that belong to another corpus, all of which otherwise surface on the
                first query that wants what refused them rather than at open.
            duckdb.OutOfMemoryException: With ``warm``, if DuckDB is held below what the
                occurrence index wants to build; see ``check_index`` for that floor, and
                for the container cap that kills the process rather than raising.
            ValueError: If ``require_pivots`` is set beside a ``pivot_budget_bytes``,
                which budgets for a build it rules out. If ``derive_pivots`` is set
                beside ``warm``, which reads the levels the index covers before the
                derivation can complete them. If ``max_rows`` is less than one, which is
                a bound no query can satisfy rather than a small one. If
                ``pivot_budget_bytes`` is negative, which no amount of memory is; zero
                is allowed, and materializes nothing; if ``derive_pivots`` or
                ``require_pivots`` is set without a ``pivots_dir``; or if
                ``require_occurrences`` is set without an ``occurrences_dir``.
        """
        self._resolver = resolver if resolver is not None else _resolve_with_resolvers
        if max_rows is not None and max_rows <= 0:
            raise ValueError(
                f"max_rows is {max_rows}, which no query can return; leave it unset to "
                "bound nothing"
            )
        self._max_rows = max_rows
        if pivot_budget_bytes is not None and pivot_budget_bytes < 0:
            raise ValueError(
                f"pivot_budget_bytes is {pivot_budget_bytes}, which is not an amount "
                "of memory; pass zero to materialize nothing"
            )
        if require_pivots and pivot_budget_bytes is not None:
            raise ValueError(
                "require_pivots leaves no pivot to build in memory, so a "
                "pivot_budget_bytes beside it budgets for something that cannot "
                "happen; pass one or the other"
            )
        # Zero where every level is read: nothing reaches the builder to be budgeted.
        self._pivot_budget = (
            0
            if require_pivots
            else (
                pivot_budget_bytes
                if pivot_budget_bytes is not None
                else _PIVOT_BUDGET_BYTES
            )
        )
        # Said at open because the alternative is saying it nowhere: a process killed
        # for holding too much leaves no log of its own, and what it thought it was
        # allowed to hold is the first thing worth knowing afterwards.
        if require_pivots:
            logger.info("every pivot is read from %s, and none is built", pivots_dir)
        else:
            logger.info(
                "pivots built in process may hold %s", _format_bytes(self._pivot_budget)
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
        # Pivots built in process, most recently used last; see _pivoted_table. The
        # serial names them, and never repeats, so a build that fails partway leaves no
        # name for the next one to collide with.
        self._pivoted: collections.OrderedDict[str, _Pivot] = collections.OrderedDict()
        # Levels whose pivot came out larger than the whole budget, and what each cost.
        # Kept because the projection they are built from does not change while this
        # object is open, so building one again costs a pass over the corpus to reach
        # the same refusal -- and because every query over one is slow for a reason
        # worth restating; see _warn_refused.
        self._too_large: dict[str, int] = {}
        self._pivot_serial = 0
        self._pivots_lock = threading.Lock()
        # Held across a materialization, so the two memory readings bracket one table
        # and no search waits on a build to be told what the cache already holds.
        self._build_lock = threading.Lock()
        if derive_pivots and pivots_dir is None:
            raise ValueError(
                "derive_pivots was set without a pivots_dir, so there is nowhere to "
                "write what it would derive"
            )
        if require_pivots and pivots_dir is None:
            raise ValueError(
                "require_pivots was set without a pivots_dir, so there is nowhere to "
                "look for what it would require"
            )
        if require_occurrences and occurrences_dir is None:
            raise ValueError(
                "require_occurrences was set without an occurrences_dir, so there is "
                "nowhere to look for what it would require"
            )
        if derive_pivots and warm:
            # The index reads the artifacts over the levels it indexes as files, so
            # warming checks a level derive_pivots means to complete before the
            # derivation that completes it has run -- and a set covering some of the
            # projections is exactly what the check refuses; see _derive_pivot. The
            # ordering that would work derives every indexed level at open, minutes
            # charged to a corpus that asked for artifacts as they were wanted.
            raise ValueError(
                "derive_pivots was set beside warm, which reads the levels the "
                "occurrence index covers before the derivation can complete them; open "
                "with warm=False and let a query or check_index build the index once "
                "the artifacts are there"
            )
        pairs, self._fingerprint = _pair(
            projection_pattern, structures_pattern, require_current
        )
        self._pivots_dir = pivots_dir
        self._occurrences_dir = occurrences_dir
        self._derive_pivots = derive_pivots
        self._require_current = require_current
        # The pattern the derivation reads, which is what mirrors the projections'
        # layout under a level's subdirectory; see _derive_pivot.
        self._projection_pattern = projection_pattern
        self._derive_lock = threading.Lock()
        # The projections a pivot artifact has to have been derived from, and the views
        # already published over the artifacts that were; see _pivot_view.
        self._projections = [pair[0] for pair in pairs]
        self._pivot_views: dict[str, str] = {}
        self._views_lock = threading.Lock()
        # The artifacts found for each level, with the source each names, under their
        # own lock: the occurrence build reads them as files while a query reads the
        # same level as a view, and the lock holds a footer read per artifact, which is
        # not something to hold the lock a published view is looked up under.
        self._pivot_artifacts_found: dict[str, dict[str, str] | None] = {}
        # The occurrence artifacts found for each indexed path, under the same lock and
        # for the same reason.
        self._occurrence_artifacts_found: dict[str, dict[str, str] | None] = {}
        self._artifacts_lock = threading.Lock()
        # What a dataset's structure IDs add to reach a corpus-wide one, keyed by the
        # source every artifact derived from it names; filled by _prepare.
        self._offset_of_source: dict[str, int] = {}
        self._connection = duckdb.connect(
            config={"memory_limit": memory_limit} if memory_limit is not None else {}
        )
        try:
            _cache_footers(self._connection)
            _warn_when_the_cap_leaves_no_headroom(self._connection)
            self._total, self._searchable = self._prepare(pairs)
            if require_pivots:
                self._require_every_pivot()
            if require_occurrences:
                self._require_every_occurrence()
            if warm:
                self._warm()
        except Exception:
            # Nothing else holds the connection yet, and the caller has no object to
            # close one from if __init__ does not return.
            self._connection.close()
            raise

    def _warm(self) -> None:
        """Builds at open everything a structure query would otherwise build on demand.

        The index before the library, so the library's gigabytes are not resident
        through the build that wants the most memory of anything here.

        Raises:
            PairingError: Whatever the index or the library refuses; see
                ``_occurrences`` and ``_library``.
        """
        self._occurrences()
        self._library()

    def _prepare(self, pairs: list[tuple[str, str, str]]) -> tuple[int, int]:
        """Publishes the relations, and returns the total and searchable row counts."""
        offsets = []
        total = 0
        stated = 0
        for projected, structured, source in pairs:
            with pq.ParquetFile(structured) as artifact:
                count = artifact.metadata.num_rows
            with pq.ParquetFile(projected) as artifact:
                stated += artifact.metadata.num_rows
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
            # Keyed by the source rather than by the file, so an artifact derived from
            # this projection finds the offset wherever it is filed; see
            # ``_pivot_offsets``.
            self._offset_of_source[source] = total
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
        projection_files = _sql_paths(pair[0] for pair in pairs)
        structure_files = _sql_paths(pair[1] for pair in pairs)
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
        # every query with silence. Each side is counted against the total its own
        # footers state, which is where the loss is visible: a reaction the view lost is
        # one no query can find, and the occurrence index cannot stand in for the count,
        # since a path read from a pivot artifact reaches the structures whatever the
        # view holds.
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
        reactions = self._connection.execute(
            f"SELECT count(*) FROM {query.TABLE}"  # noqa: S608
        ).fetchone()
        assert reactions is not None  # An aggregate over any relation returns one row.
        if reactions[0] != stated:
            raise PairingError(
                f"the projections hold {stated} reactions but only {reactions[0]} "
                "joined to their offsets; a path did not survive read_parquet"
            )
        return total, searchable

    @property
    def fingerprint(self) -> str:
        """Returns a short digest naming the artifacts this corpus opened.

        A question log records this beside the query rather than the rows the query
        returned, so an old record stays reproducible: the same query against the same
        fingerprint answers the same way, and a different fingerprint says why it does
        not.
        """
        return self._fingerprint

    def _require_every_pivot(self) -> None:
        """Refuses a corpus whose pivots_dir does not hold every level, current.

        ``check_pivots`` reports what is held; this insists. A level with no artifacts
        is otherwise built in process on the query that first wants it, which is
        minutes of unnesting charged to whoever asked and no error anywhere -- the
        failure a caller naming a pivots_dir is least likely to be expecting, since
        naming one says the artifacts are meant to answer.

        Raises:
            PairingError: If any level has no artifacts. Wrong provenance, a level
                filed under another, and staleness are raised by the publishing itself.
        """
        missing = [
            path for path in sorted(pivot.LEVELS) if self._pivot_view(path) is None
        ]
        if missing:
            raise PairingError(
                f"{self._pivots_dir} holds no pivot artifacts for {missing}, which "
                "this corpus was opened requiring; derive them first with "
                "ord_schema.artifacts.scripts.derive_pivots"
            )

    def check_pivots(self) -> dict[str, int]:
        """Publishes every pivot artifact this corpus can read, and returns the counts.

        Meant for startup. Artifacts are otherwise found by the first query that wants
        a level, so a tree derived from another corpus, filed under the wrong level, or
        left stale is a ``PairingError`` raised at whatever hour that query arrives.
        Calling this makes the same refusal happen while a server is still starting,
        where it is a failed deployment rather than a failed request.

        Publishing is the check: a level whose artifacts are wrong cannot be published,
        and one that is fine is left ready, so the first real query does not pay the
        glob. Cheap either way -- a view over Parquet reads footers, not rows.

        Reports what is held rather than arranging for it: this reaches every level the
        schema has, and ``derive_pivots`` here would be 39 unnests of the projection at
        startup for a deployment that asks about two.

        Returns:
            The number of artifacts found per level, for the levels that have any.
            Empty when the corpus was opened without ``pivots_dir``, which is not an
            error: every level is then built in process.

        Raises:
            PairingError: If any level's artifacts were not derived from this corpus's
                projections, hold a different level than the one they are filed under,
                or -- with ``require_current`` -- are stale.
        """
        if self._pivots_dir is None:
            return {}
        found = {}
        for path in sorted(pivot.LEVELS):
            if self._pivot_view(path) is not None:
                found[path] = len(pivot.artifact_paths(self._pivots_dir, path))
        logger.info(
            "%d of %d levels are held as artifacts: %s",
            len(found),
            len(pivot.LEVELS),
            ", ".join(f"{path} {count}" for path, count in found.items()) or "none",
        )
        return found

    def check_index(self) -> dict[str, int]:
        """Builds the occurrence index and returns what it reached, per path.

        The sibling of ``check_pivots``, and for the same reason: a corpus opened
        without ``warm`` builds the index on the first query that can spend it, so a
        corpus it refuses -- a reaction stated twice, a structure no indexed path
        reaches -- fails that query rather than the deployment. It also has a memory
        floor, and meeting that one the other way is worse still. Over ORD the build
        wants 5-6.5 GB of DuckDB memory and writes 16-25 GB of temporary files getting
        there; below that it raises ``duckdb.OutOfMemoryException`` rather than running
        slowly, since a block it cannot pin is not one it can spill.

        That exception is the good outcome, and under a cgroup cap it is only reached
        with room to spare above ``memory_limit``: the build holds up to 2.8 GB outside
        DuckDB's accounting, and reaching the cap is a kill, which arrives as neither an
        exception nor a log line. DuckDB's default limit is about 80% of the cap, which
        does not leave it -- an 8 GiB container is killed 85s in, where a 12 GiB one
        holding DuckDB to 6500MiB finishes at 9.3 GB resident. ``Corpus`` warns at open
        when a cap it can read leaves too little.

        None of that applies to an index published as a view, which an
        ``occurrences_dir`` covering every indexed path gets: there is nothing to
        materialize, so the floor is the scan this call makes rather than a build. Over
        ORD that is 0.13s and 0.32 GiB resident, against 2.9s and 1.92 GiB for the table
        read from pivots.

        A corpus opened with ``warm`` has paid whichever of those it owes already, and
        the call is then its counts. Leaving it cold suits one asked only for scalar or
        similarity queries, which never spends what a built index costs: about a minute
        over ORD where no artifact covers a path, and 1.19 GB resident once built.

        Returns:
            The number of occurrences found per indexed path, including the paths that
            hold none -- most corpora record no authentic standards, and a zero there is
            ordinary rather than a sign the walk missed something.

        Raises:
            PairingError: If the corpus states a reaction more than once, or the index
                does not reach every structure the corpus holds.
        """
        self._occurrences()
        cursor = self._connection.cursor()
        try:
            counts = dict(
                cursor.execute(
                    "SELECT path, count(*) FROM occurrences GROUP BY path"
                ).fetchall()
            )
        finally:
            cursor.close()
        return {path: counts.get(path, 0) for path in INDEXED_PATHS}

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

        Built at open under ``warm``, and otherwise on the first substructure query: it
        costs seconds and gigabytes, and a corpus asked only for similarity or scalar
        queries never needs it.

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

    def _occurrences_from_pivot(self, path: str) -> str | None:
        """Returns the SELECT reading ``path``'s occurrences from a pivot, or None.

        A pivot holds one row per element of the nearest repeated level, so an indexed
        path that descends from one through singular struct fields is read from that
        level's artifact rather than unnested out of the projection.

        The corpus-wide structure ID is the element's own plus its dataset's offset,
        which is a running total over the corpus's pairs and so belongs to no artifact:
        the same pivot read beside a different set of datasets needs a different
        offset. It is joined here from the file each row came from, keyed by the source
        the file and its projection both name.

        Args:
            path: The indexed path.

        Returns:
            The SELECT, or None where no repeated level covers the path, or where the
            directory holds no artifact for the level that does.

        Raises:
            PairingError: If the artifacts do not pair with this corpus; see
                ``_pivot_artifacts``.
        """
        reached = pivot.reach(path)
        if reached is None:
            return None
        level, remainder, _ = reached
        sources = self._pivot_artifacts(level.path)
        if sources is None:
            return None
        offsets = self._artifact_offsets(
            f"{pivot.table_name(level.path)}_offsets", sources
        )
        element = ".".join(["x", pivot.ELEMENT, *remainder])
        # S608: the level and its remainder come from this module's walk of the schema,
        # and the paths from this process's own glob, quoted with their single quotes
        # escaped.
        return f"""
            SELECT x.{pivot.REACTION_ID}, {sql_string(path)} AS path,
                   ({element}.structure_id + o.{query.STRUCTURE_OFFSET})::UINTEGER
                       AS global_id,
                   {element}.{_INDEXED_FIELD} AS {_INDEXED_FIELD}
            FROM read_parquet({_sql_paths(sources)}, filename=true) x
            JOIN {offsets} o ON x.filename = o.artifact_filename
            WHERE {element}.structure_id IS NOT NULL
            """  # noqa: S608

    def _artifact_offsets(self, published: str, sources: dict[str, str]) -> str:
        """Publishes what each artifact's rows add to reach a corpus-wide ID.

        An artifact and the projection it descends from restate one source dataset, so
        they carry the same ``source_md5`` however many derivations apart they sit, and
        that hash is what pairs an artifact with the offset its projection was given.
        Pairing on the filename instead would rest on the two trees being laid out
        alike, which nothing enforces.

        Args:
            published: The relation to create, named for the level or path it serves so
                that a retry after a refusal replaces these rows rather than leaving a
                relation behind per attempt.
            sources: The source each artifact names, from the check that accepted them.

        Returns:
            ``published``.
        """
        rows = []
        for name, source in sources.items():
            # The checks require these artifacts' sources to be exactly the
            # projections', and the offsets are keyed by those same hashes.
            offset = self._offset_of_source.get(source)
            assert offset is not None, f"{name} names a source the corpus does not hold"
            rows.append((name, offset))
        registered = pa.table(
            {
                "artifact_filename": [row[0] for row in rows],
                query.STRUCTURE_OFFSET: pa.array(
                    [row[1] for row in rows], type=pa.int64()
                ),
            }
        )
        # Under the build lock, for the reason _build states: a materialization brackets
        # its table with two memory readings taken across everyone. S608: the name comes
        # from this module's walk of the schema.
        with self._build_lock:
            self._connection.register("registered_artifact_offsets", registered)
            self._connection.execute(
                f"CREATE OR REPLACE TABLE {published} AS "  # noqa: S608
                "SELECT * FROM registered_artifact_offsets"
            )
            self._connection.unregister("registered_artifact_offsets")
        return published

    def _occurrence_artifacts(self, path: str) -> dict[str, str] | None:
        """Returns the artifacts holding ``path``'s occurrences, checked, or None.

        Globbed and checked once per path and remembered, including the answer that
        there are none, so a corpus without artifacts reads the directory once rather
        than on every query.

        Args:
            path: The indexed path, as the query grammar names it.

        Returns:
            The source dataset each artifact was derived from, keyed by its path, or
            None where the corpus has no directory or the directory holds none for this
            path.

        Raises:
            PairingError: If the artifacts do not belong to this corpus; see
                ``_check_occurrences``.
        """
        if self._occurrences_dir is None:
            return None
        with self._artifacts_lock:
            if path in self._occurrence_artifacts_found:
                return self._occurrence_artifacts_found[path]
            files = [
                str(found)
                for found in occurrences.artifact_paths(self._occurrences_dir, path)
            ]
            if not files:
                self._occurrence_artifacts_found[path] = None
                return None
            sources = self._check_occurrences(path, files)
            self._occurrence_artifacts_found[path] = sources
            return sources

    def _check_occurrences(self, path: str, files: list[str]) -> dict[str, str]:
        """Refuses occurrence artifacts that do not belong to this corpus.

        The check ``_check_pivots`` makes, for the same reason: an occurrence names its
        reaction by ID and its structure by a dataset-local ``structure_id``, and both
        resolve against whatever this corpus holds. Artifacts derived from another
        corpus therefore answer a quantifier confidently and wrongly rather than
        failing.

        Every indexed path writes the same three columns, so where a file sits says
        nothing about what it holds and the stamped path is the only thing that does. A
        file under the wrong one would answer that path's quantifiers with another
        level's elements.

        Args:
            path: The indexed path the files claim to hold.
            files: The artifacts found for it.

        Returns:
            The source dataset each artifact was derived from, keyed by its path, which
            is what pairs an artifact with the offset its projection was given; see
            ``_artifact_offsets``.

        Raises:
            PairingError: If a file is not an occurrences artifact, if it is stamped
                with another path, if two artifacts restate one source dataset, if the
                set of source datasets differs from the projections', or -- with
                ``require_current`` -- if any artifact is stale.
        """
        wanted = {base.load_stamps(name).source_md5 for name in self._projections}
        sources: dict[str, str] = {}
        derived_from: dict[str, str] = {}
        for name in files:
            stamps = base.load_stamps(name)
            if stamps.artifact != occurrences.ARTIFACT:
                raise PairingError(
                    f"{name} is a {stamps.artifact}, not an {occurrences.ARTIFACT} "
                    "artifact"
                )
            held = occurrences.occurrence_path(name)
            if held != path:
                raise PairingError(
                    f"{name} holds the occurrences at {held}, but sits where {path} is "
                    "read from, so a quantifier over it would be answered by another "
                    "path's elements"
                )
            if self._require_current and not base.stamps_are_current(
                stamps, occurrences.ARTIFACT
            ):
                raise PairingError(f"{name} is a stale {occurrences.ARTIFACT} artifact")
            # Two artifacts of one dataset pass the set comparison below and are both
            # read, so every occurrence at the path is stated twice. A quantifier cannot
            # see it -- a semi-join returns a reaction once however many rows name it --
            # but the structure count that decides whether the index reaches the corpus
            # is taken over distinct IDs, so it cannot see it either, and check_index
            # reports the doubled count as the corpus's own.
            if stamps.source_md5 in derived_from:
                raise PairingError(
                    f"{name} and {derived_from[stamps.source_md5]} are both "
                    f"{occurrences.ARTIFACT} artifacts at {path} of the same source "
                    "dataset, so its occurrences would each be read twice"
                )
            derived_from[stamps.source_md5] = name
            sources[name] = stamps.source_md5
        if set(sources.values()) != wanted:
            found = set(sources.values())
            raise PairingError(
                f"the {occurrences.ARTIFACT} artifacts at {path} were derived from "
                f"{len(found)} source datasets and the projections from {len(wanted)}, "
                f"and {len(wanted - found)} of the projections' are missing; an "
                "artifact of another corpus names reactions this one does not hold"
            )
        return sources

    def _occurrences_from_artifact(self, path: str) -> str | None:
        """Returns the SELECT reading ``path``'s occurrences from artifacts, or None.

        The artifact is already one row per occurrence, so the read adds only the
        corpus-wide ID: the element's own plus its dataset's offset, which is a running
        total over the corpus's pairs and so belongs to no artifact. It is joined here
        from the file each row came from, keyed by the source that file and its
        projection both name.

        Args:
            path: The indexed path.

        Returns:
            The SELECT, or None where the corpus has no artifacts for the path.

        Raises:
            PairingError: If the artifacts do not belong to this corpus; see
                ``_check_occurrences``.
        """
        sources = self._occurrence_artifacts(path)
        if sources is None:
            return None
        offsets = self._artifact_offsets(
            f"occurrences_{path.replace('.', '_')}_offsets", sources
        )
        # The artifact holds no row for an element carrying no structure, so nothing
        # here filters: what a pivot needs a WHERE for, the derivation already did.
        # S608: the path is a key of INDEXED_PATHS and the filenames come from this
        # process's own glob, quoted with their single quotes escaped.
        return f"""
            SELECT x.reaction_id, {sql_string(path)} AS path,
                   (x.structure_id + o.{query.STRUCTURE_OFFSET})::UINTEGER AS global_id,
                   x.{_INDEXED_FIELD}
            FROM read_parquet({_sql_paths(sources)}, filename=true) x
            JOIN {offsets} o ON x.filename = o.artifact_filename
            """  # noqa: S608

    def _require_every_occurrence(self) -> None:
        """Refuses a corpus whose occurrences_dir does not hold every indexed path.

        A path with no artifact is otherwise read from a pivot or unnested out of the
        projection, and either leaves the index a materialized table -- 1.19 GB held
        over ORD, and 5-6.5 GB of DuckDB memory wanted to build it. A deployment sized
        for the view is sized for neither, so the fallback that keeps a corpus
        answering is the one that gets the container killed, which is the failure a
        caller naming an occurrences_dir is least likely to be expecting.

        Raises:
            PairingError: If any indexed path has no artifacts. Wrong provenance, a
                path filed under another, and staleness are raised by the check itself.
        """
        missing = [
            path
            for path in sorted(INDEXED_PATHS)
            if self._occurrence_artifacts(path) is None
        ]
        if missing:
            raise PairingError(
                f"{self._occurrences_dir} holds no {occurrences.ARTIFACT} artifacts "
                f"for {missing}, which this corpus was opened requiring; derive them "
                "first with ord_schema.artifacts.scripts.derive_occurrences"
            )

    def _occurrences(self) -> None:
        """Builds the occurrence index, once.

        One row per structure occurrence, carrying the corpus-wide ID, the path it sat
        at, and the element's own ``reaction_role``. Keeping the role beside the
        structure is what preserves element binding: "pyridine as the solvent" stays a
        condition on one row rather than an intersection of two reaction sets, which
        over-returns.

        Costs a read of the pivot artifacts covering the indexed paths and one pass over
        the projections per path they do not cover, charged at open under ``warm`` and
        otherwise to the query that first wants it. Concurrent first queries serialize
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
            # Three ways to reach a path's occurrences, cheapest first. An occurrences
            # artifact is the rows themselves and needs only the offset joined on. A
            # level's pivot holds one row per element, which is a scan of a few hundred
            # megabytes against a traversal of every reaction. Failing both, the
            # projection is unnested. Decided per path, so a directory holding some and
            # not others is a partial speedup rather than a missing index.
            reads, from_artifacts, from_pivots = [], [], []
            for path, expression in INDEXED_PATHS.items():
                reading = self._occurrences_from_artifact(path)
                if reading is not None:
                    from_artifacts.append(path)
                else:
                    reading = self._occurrences_from_pivot(path)
                    if reading is not None:
                        from_pivots.append(path)
                    else:
                        reading = _occurrences_from_projection(path, expression)
                reads.append(reading)
            selects = "\nUNION ALL\n".join(reads)
            # A view where every path is an artifact, and a table otherwise. Read per
            # query the artifacts cost 1.28x the built table summed over every path and
            # pattern -- and less on the narrow paths, where reading one path's files
            # beats filtering path = ? across every occurrence in the corpus -- against
            # the 1.19 GB the table holds and the 5-6.5 GB it wants to build. A path
            # read from a pivot or unnested out of the projection is the traversal this
            # index exists to avoid, and a view would repeat it on every query rather
            # than pay it once, so one uncovered path materializes all of them.
            held = "VIEW" if len(from_artifacts) == len(INDEXED_PATHS) else "TABLE"
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
                # Under the build lock, which every statement that puts a relation in
                # this database takes, so none of them lands while a materialization is
                # measuring what one costs -- a difference of two readings taken across
                # everyone.
                # OR REPLACE, so a build interrupted after the relation exists is a
                # build the next query repeats rather than one that collides with
                # itself forever. The kind cannot change under an open corpus, since
                # what it turns on is which artifacts the directory held at open.
                # S608: the fragments are this module's own schema walk and the
                # compiler's traversals, not anything a query supplies.
                with self._build_lock:
                    cursor.execute(f"CREATE OR REPLACE {held} occurrences AS {selects}")
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
                        f"{self._total} structures over {sorted(INDEXED_PATHS)}, so "
                        "the reactions holding the rest cannot be found; the "
                        "projections are not the schema this walk was built from"
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
                "indexed %d structure occurrences in %.1fs as a %s, %d of %d paths "
                "read from occurrence artifacts and %d from pivot artifacts: %s",
                sum(counts.values()),
                time.perf_counter() - start,
                held.lower(),
                len(from_artifacts),
                len(INDEXED_PATHS),
                len(from_pivots),
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

    def _hashed_ids(
        self,
        cursor: duckdb.DuckDBPyConnection,
        parameter: query.StructureParameter,
        resolve: Callable[[str], str],
    ) -> list[int]:
        """Matches on one of the artifact's hash columns; returns global IDs.

        The artifact already holds every structure's hashes, so the whole answer is one
        equality against a column -- no screen, no verification, and no chemistry
        beyond hashing the query molecule the way the artifact hashed its own.

        Args:
            cursor: The cursor the match runs on.
            parameter: The predicate to match, whose ``op`` chooses the column.
            resolve: Maps a compound name to SMILES.

        Returns:
            The corpus-wide structure IDs sharing the query molecule's hash.

        Raises:
            ValueError: If a resolved compound's SMILES does not parse, or if RDKit
                refuses to hash the query molecule -- either leaves the question
                unaskable rather than quietly unanswerable.
        """
        column, hasher = _HASH_MATCHES[parameter.op]
        molecule = self._query_molecule(parameter, resolve)
        hashed = hasher(molecule)
        if hashed is None:
            raise ValueError(
                f"could not hash the query molecule for {parameter.name!r}; the corpus "
                "compares compounds by hash, so there is nothing to compare it to"
            )
        # S608: the column comes from this module's own table, keyed by an op the
        # grammar holds to a Literal; the molecule reaches the query as a parameter.
        rows = cursor.execute(
            f"SELECT global_id FROM corpus_structures WHERE {column} = $h",  # noqa: S608
            {"h": hashed},
        ).fetchall()
        return [row[0] for row in rows]

    def _bitmap(self, matched: Sequence[int]) -> str:
        """Returns the match set as a bitmap over the corpus-wide ID space."""
        bits = bytearray(b"0" * max(self._total, 1))
        for global_id in matched:
            bits[global_id] = ord("1")
        return bits.decode()

    def _materialize(self, path: str, select: str) -> _Pivot | None:
        """Returns the pivot over ``path``, building it if the cache lacks one.

        The caller owns a read on whatever comes back and has to release it, which
        ``_pivoted_table`` does; a pivot with a reader is never evicted.

        What the cache holds is read under the short lock the bookkeeping takes, and
        only a build waits on another build. A search over a level already materialized
        is answered while one is running, which is the common case on a corpus serving
        several questions at once.

        Args:
            path: The repeated level, as the query grammar names it.
            select: The query ``CREATE TABLE`` wraps.

        Returns:
            The pivot, or None when the budget is zero or this costs more than the cache
            may hold in total. Nothing is kept either way and the quantifier compiles
            over the elements; a level refused as too large is remembered, so the next
            query wanting it does not build it again to throw it away again.
        """
        if not self._pivot_budget:
            # The one budget that can be settled without building, since what a pivot
            # costs is known only once it exists. A small container spends nothing and
            # answers from Parquet.
            return None
        with self._pivots_lock:
            settled, entry, refused = self._cached(path)
        if settled:
            # Warned outside the lock: it is the one every search takes to read the
            # cache, and a log handler is not something to hold it through.
            if refused is not None:
                self._warn_refused(path, refused)
            return entry
        with self._build_lock:
            # Asked again under the build lock: whoever held it may have been building
            # exactly this, and a second table of it would be memory held twice for one
            # answer.
            with self._pivots_lock:
                settled, entry, refused = self._cached(path)
            if settled:
                if refused is not None:
                    self._warn_refused(path, refused)
                return entry
            return self._build(path, select)

    def _cached(self, path: str) -> tuple[bool, _Pivot | None, int | None]:
        """Returns whether the cache settles ``path``, and how.

        Called with ``_pivots_lock`` held. A hit takes the read on the caller's behalf,
        since a pivot released back to the cache between the lookup and the read could
        be evicted in between.

        Args:
            path: The repeated level, as the query grammar names it.

        Returns:
            Whether the cache settles the question, the pivot for a hit, and what it
            cost where it was refused as too large -- which the caller says out loud
            once it is no longer holding the lock.
        """
        cached = self._pivoted.get(path)
        if cached is not None:
            self._pivoted.move_to_end(path)
            cached.readers += 1
            return True, cached, None
        refused = self._too_large.get(path)
        if refused is not None:
            return True, None, refused
        return False, None, None

    def _warn_refused(self, path: str, cost: int) -> None:
        """Says that a query is about to be answered the slow way, and what to change.

        Warned on every query over the level rather than once when it was first refused,
        because the query is slow every time and whoever is asking why is reading the
        log now, not the one line printed when the corpus opened.

        Args:
            path: The level whose pivot was refused.
            cost: What building it took when it was measured.
        """
        logger.warning(
            "the pivot over %s takes %s, over this corpus's budget of %s, so every "
            "query over that level unnests the projection instead. Derive the pivot as "
            "an artifact, or open the Corpus with a larger pivot_budget_bytes if the "
            "machine has the memory to spare.",
            path,
            _format_bytes(cost),
            _format_bytes(self._pivot_budget),
        )

    def _build(self, path: str, select: str) -> _Pivot | None:
        """Materializes the pivot over ``path``, or records that it costs too much.

        Called with ``_build_lock`` held and ``_pivots_lock`` free, so the two memory
        readings bracket this table and nothing else: every other statement that puts a
        table in this database -- another materialization, or the index build -- takes
        the same lock.

        Args:
            path: The repeated level, as the query grammar names it.
            select: The query ``CREATE TABLE`` wraps.

        Returns:
            The pivot, held once for the caller, or None if it exceeds the budget.
        """
        with self._pivots_lock:
            # Never reused, so a build that fails between the CREATE and the entry
            # cannot leave a name the next attempt collides with.
            self._pivot_serial += 1
            name = f"pivoted_{self._pivot_serial}"
        start = time.perf_counter()
        # Its own cursor, for the same reason the index build takes one: searches
        # are in flight, and the shared connection holds their results.
        cursor = self._connection.cursor()
        try:
            # Said before the pass rather than after it, because the pass is what
            # takes the time and whoever is watching a query hang wants to know what it
            # is waiting for while it waits.
            logger.info("building the pivot over %s", path)
            before = _memory_bytes(cursor)
            # The select comes from ord_schema.artifacts.pivot, which is where its
            # fragments are accounted for.
            cursor.execute(f"CREATE TABLE {name} AS {select}")
            cost = max(_memory_bytes(cursor) - before, 0)
            if cost > self._pivot_budget:
                # Too big to keep, and keeping it is the only reason to build it.
                cursor.execute(f"DROP TABLE {name}")
                with self._pivots_lock:
                    self._too_large[path] = cost
                self._warn_refused(path, cost)
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
            "materialized the pivot over %s as %s, %s in %.1fs",
            path,
            name,
            _format_bytes(cost),
            time.perf_counter() - start,
        )
        entry = _Pivot(name=name, held=cost, readers=1)
        with self._pivots_lock:
            self._pivoted[path] = entry
            self._evict()
        return entry

    def _evict(self) -> None:
        """Drops materialized pivots, least recently used first, down to the budget.

        Called with ``_pivots_lock`` held. A pivot a search is reading is passed over
        rather than dropped: the search bound the table's name into its SQL and reads it
        only after resolving names and matching structures, so dropping it there would
        fail a query that had already been answered correctly. Passing over every
        candidate leaves the cache above its budget until a reader finishes, which is
        the direction that keeps answers right.
        """
        for path in list(self._pivoted):
            if (
                sum(entry.held for entry in self._pivoted.values())
                <= self._pivot_budget
            ):
                return
            entry = self._pivoted[path]
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
            del self._pivoted[path]
            logger.info("evicted the pivot over %s, held as %s", path, entry.name)

    def _pivot_view(self, path: str) -> str | None:
        """Returns a view over the pivot artifacts for ``path``, or None if none exist.

        Published once per level, over whichever artifacts ``_pivot_artifacts`` found
        and checked, and remembered under its own lock so a query over a published level
        waits on nothing.

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
        sources = self._pivot_artifacts(path)
        if sources is None:
            return None
        files = list(sources)
        with self._views_lock:
            published = self._pivot_views.get(path)
            if published is not None:
                return published
            name = pivot.table_name(path)
            # S608: the name comes from this module's walk of the schema, and the paths
            # from this process's own glob, quoted with their single quotes escaped.
            self._connection.execute(
                f"CREATE OR REPLACE VIEW {name} AS "  # noqa: S608
                f"SELECT * FROM read_parquet({_sql_paths(files)})"
            )
            logger.info("read %d pivot artifacts for %s", len(files), path)
            self._pivot_views[path] = name
            return name

    def _pivot_artifacts(self, path: str) -> dict[str, str] | None:
        """Returns the artifacts holding the pivot over ``path``, checked, or None.

        Globbed and checked once per level and remembered, including the answer that
        there are none, so a corpus without artifacts reads the directory once rather
        than on every query. The stamps the check reads are kept with the paths, since
        the source an artifact names is what pairs it with an offset and reading the
        footers again would be one open per artifact per level.

        Args:
            path: The repeated level, as the query grammar names it.

        Returns:
            The source dataset each artifact was derived from, keyed by its path, or
            None where the directory holds none for this level.

        Raises:
            PairingError: If the artifacts were not derived from the projections this
                corpus reads, or -- with ``require_current`` -- if any is stale.
        """
        if self._pivots_dir is None:
            return None
        with self._artifacts_lock:
            if path in self._pivot_artifacts_found:
                return self._pivot_artifacts_found[path]
            files = [
                str(found) for found in pivot.artifact_paths(self._pivots_dir, path)
            ]
            if not files:
                self._pivot_artifacts_found[path] = None
                return None
            sources = self._check_pivots(path, files)
            self._pivot_artifacts_found[path] = sources
            return sources

    def _check_pivots(self, path: str, files: list[str]) -> dict[str, str]:
        """Refuses pivot artifacts that do not belong to this corpus.

        A pivot names its reactions by ID, and an ID resolves against whatever rows the
        projections hold. Artifacts derived from another corpus therefore answer a
        quantifier confidently and wrongly rather than failing, which is why the check
        is on the source datasets rather than on the file count.

        Args:
            path: The repeated level the files claim to hold.
            files: The artifacts found for it.

        Returns:
            The source dataset each artifact was derived from, keyed by its path, which
            is what pairs an artifact with the offset its projection was given; see
            ``_pivot_offsets``.

        Raises:
            PairingError: If a file is not a pivot over ``path``, if two artifacts
                restate one source dataset, if the set of source datasets differs from
                the projections', or -- with ``require_current`` -- if any artifact is
                stale.
        """
        wanted = {base.load_stamps(name).source_md5 for name in self._projections}
        sources: dict[str, str] = {}
        derived_from: dict[str, str] = {}
        for name in files:
            stamps = base.load_stamps(name)
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
            if self._require_current and not base.stamps_are_current(
                stamps, pivot.ARTIFACT
            ):
                raise PairingError(f"{name} is a stale {pivot.ARTIFACT}")
            # Two artifacts of one dataset pass the set comparison below and are both
            # read, so every element of the level is stated twice. A quantifier cannot
            # see it -- a semi-join returns a reaction once however many rows name it --
            # but the occurrence index counts those rows, and check_index reports the
            # count.
            if stamps.source_md5 in derived_from:
                raise PairingError(
                    f"{name} and {derived_from[stamps.source_md5]} are both "
                    f"{pivot.ARTIFACT} artifacts over {path} of the same source "
                    "dataset, so the level's elements would each be read twice"
                )
            derived_from[stamps.source_md5] = name
            sources[name] = stamps.source_md5
        if set(sources.values()) != wanted:
            found = set(sources.values())
            raise PairingError(
                f"the {pivot.ARTIFACT} artifacts for {path} were derived from "
                f"{len(found)} source datasets and the projections from {len(wanted)}, "
                f"and {len(wanted - found)} of the projections' are missing; a pivot "
                "of another corpus names reactions this one does not hold"
            )
        return sources

    def _derive_pivot(self, path: str) -> str | None:
        """Writes the pivot artifacts for ``path`` and publishes them.

        The same derivation ``derive_pivots`` runs offline, aimed at this corpus's own
        projections: one artifact per projection, stamped, and skipped where one is
        already current -- so a run interrupted partway is finished by the next rather
        than started again. It reads a projection at a time and streams the artifact's
        rows out a batch at a time, but holds the level's ancestors whole on the way --
        under DuckDB's own memory limit rather than the pivot budget, which bounds what
        this corpus holds in process and not what a derivation spends.

        Derived before the level is read rather than after refusing to read it. What an
        interrupted run leaves behind is a set covering some of the projections, which
        is exactly what ``_pivot_view`` refuses; checking first would turn the
        interruption into a level that can never be completed from here, since the pass
        that would complete it is the one behind the refusal. The cost of that ordering
        is one ``derive_tree`` pass per level per process where every artifact is
        already current, which reads footers and writes nothing.

        One derivation at a time, and the published view is looked up again under the
        lock: whoever held it may have been deriving exactly this, and a second pass
        would rewrite the artifacts the first one just wrote.

        Args:
            path: The repeated level, as the query grammar names it.

        Returns:
            The view's name, or None if the derivation wrote nothing -- which means the
            pattern reaches no projection this corpus reads, and the level is left to be
            built in process.

        Raises:
            PairingError: If what is there once the derivation has run still does not
                pair with this corpus -- a stranger's artifacts, or a level filed under
                another one's name, neither of which a pass over the projections
                displaces. ``_pivot_view`` decides.
            ValueError: If a projection is stale, since a pivot derived from one would
                claim a provenance it does not have. Reachable only with
                ``require_current`` off, which is what let the corpus open at all.
        """
        assert self._pivots_dir is not None  # Checked when derive_pivots was accepted.
        with self._derive_lock:
            with self._views_lock:
                published = self._pivot_views.get(path)
            if published is not None:
                return published
            start = time.perf_counter()
            logger.info("deriving the pivot artifacts for %s", path)
            written, skipped, _ = base.derive_tree(
                self._projection_pattern,
                str(pathlib.Path(self._pivots_dir) / path),
                artifact=pivot.ARTIFACT,
                write=functools.partial(pivot.write_pivot, level_path=path),
                # The same columns derive_pivots declares, so a tree built by one and
                # extended by the other agrees on which artifacts are current.
                schema=pivot.schema(pivot.LEVELS[path]),
                parent_artifact=projection.ARTIFACT,
            )
            logger.info(
                "derived %d pivot artifacts for %s in %.1fs (%d already current)",
                written,
                path,
                time.perf_counter() - start,
                skipped,
            )
            if not written and not skipped:
                return None
            # A glob that found nothing is remembered, so it has to be forgotten before
            # the same call can find what this just wrote.
            with self._artifacts_lock:
                self._pivot_artifacts_found.pop(path, None)
            return self._pivot_view(path)

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
        if self._derive_pivots:
            # The same pass, spent where it survives the process rather than dying with
            # it -- and against no budget, since what comes back is a view over Parquet.
            # It reads the level as well as writing it, because a set an interrupted run
            # left partial is one that has to be completed before it can be read.
            published = self._derive_pivot(path)
        else:
            published = self._pivot_view(path)
        if published is not None:
            # A view over artifacts costs no memory to hold and nothing to release, so
            # there is no reader to take and no budget to spend.
            yield published
            return
        entry = self._materialize(path, pivot.select(level, query.TABLE))
        if entry is None:
            yield None
            return
        try:
            yield entry.name
        finally:
            with self._pivots_lock:
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
            cursor: The cursor a similarity screen or a hash match runs on.
            parameter: The predicate to evaluate.
            resolve: Maps a compound name to SMILES.

        Returns:
            The bitmap over corpus-wide structure IDs.

        Raises:
            ValueError: If the predicate names neither a pattern nor a compound, if a
                resolved compound's SMILES does not parse, or if RDKit refuses to hash
                the query molecule of a predicate that matches on a hash.
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
            elif parameter.op in _HASH_MATCHES:
                matched = self._hashed_ids(cursor, parameter, resolve)
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
        at a materialization, which is one lock for all pivots: the cost of a table is
        read as the difference two memory readings make, and two builds interleaved
        would each be charged the other's.

        A quantifier the occurrence index can answer becomes a condition on the row
        rather than a filter over the elements, and the rest of the query -- the clauses
        beside it, the aggregate, the ordering, the limit -- compiles exactly as it
        would have. The screening and verification are the compiler's either way, so an
        indexed query and a projection one differ only in how they reach the reactions
        holding a match. Which happened is logged, because it is the one thing about a
        search that the result does not show.

        A pivot the query routes a quantifier to is held for the life of the search,
        since the SQL names it seconds before it reads it and eviction would otherwise
        drop it in between.

        Args:
            request: The query to run.
            timeout_seconds: Wall-clock bound on the whole call, not on the final
                query alone. The query is interrupted; a phase already running is not,
                and raises when it finishes instead -- so a search can take longer than
                this and still report the overrun. ``_Deadline`` says which phases those
                are and why. None runs unbounded.

        Returns:
            The selected columns: ``reaction_id`` for a plain query, the group and
            measure columns for an aggregated one. At most this corpus's ``max_rows``,
            where it has one: an answer that reached the bound is logged as possibly
            cut, and the table itself does not say. ``Compiled.limit`` is what actually
            reached the SQL, for a caller that wants to check rather than read a log.

        Raises:
            query.QueryError: If the query does not compile.
            ValueError: If a compound name cannot be resolved.
            PairingError: If the corpus IDs do not come out one unbroken run, or the
                occurrence index refuses the corpus -- a reaction stated by more than
                one row, or a structure no indexed path reaches.
            TimeoutError: If the search exceeds ``timeout_seconds``, naming the phase
                that was running when the bound passed.
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
                resolved[name] = _canonical(self._resolver(name))
            return resolved[name]

        deadline = _Deadline(timeout_seconds)
        with contextlib.ExitStack() as reading:
            pivoted: dict[str, str | None] = {}

            def pivot_table(path: str) -> str | None:
                # Built while compiling rather than after, because the budget decides
                # whether the table exists at all and the SQL names it either way. The
                # answer is remembered so a level two quantifiers range over is one
                # build and one read, held once for the whole search.
                if path not in pivoted:
                    pivoted[path] = reading.enter_context(self._pivoted_table(path))
                return pivoted[path]

            compiled = query.compile_query(
                request, index=index, pivot=pivot_table, max_rows=self._max_rows
            )
            # Compiling is cheap, but pivot_table above is not: a level with no artifact
            # is unnested out of the projection here, which is minutes over a corpus.
            deadline.check("compiling the query")
            if compiled.limit != request.limit:
                logger.info("the corpus bounds this query at %d rows", compiled.limit)
            if indexing:
                # Built after compiling and before running: the condition names the
                # table, so the table has to exist by the time the query does. Logged
                # after the build, so a build that refuses the corpus does not leave a
                # line claiming the index answered a query nothing ran.
                self._occurrences()
                deadline.check("building the occurrence index")
                logger.info("the occurrence index answers part of this query")
            else:
                logger.info("the projection answers this query")
            cursor = self._connection.cursor()
            try:
                parameters: dict[str, Any] = {
                    name: resolve(name) for name in compiled.compounds
                }
                # An external service, so how long it takes is not this corpus's to
                # predict; a resolver that hangs is the likeliest way a search overruns
                # without a single slow query in it.
                deadline.check("resolving compound names")
                for parameter in compiled.structures:
                    parameters[parameter.name] = self._matches(
                        cursor, parameter, resolve
                    )
                    # Per predicate rather than after all of them, so a query with two
                    # does not pay for the second once the first has spent the budget.
                    # Named by what the caller wrote rather than by parameter.name,
                    # which is the compiler's own placeholder and identifies nothing
                    # outside the SQL it was allocated for.
                    asked = parameter.pattern or parameter.compound
                    deadline.check(f"matching {parameter.op} {asked!r}")
                left = deadline.remaining()
                if left is None:
                    answered = cursor.execute(compiled.sql, parameters).to_arrow_table()
                else:
                    answered = _run_with_timeout(cursor, compiled.sql, parameters, left)
                _warn_when_the_bound_was_reached(request, compiled, answered.num_rows)
                return answered
            finally:
                cursor.close()
