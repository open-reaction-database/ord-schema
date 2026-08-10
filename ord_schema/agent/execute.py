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
element's own ``reaction_role`` -- turns that scan into a filter over a narrow table.
Keeping the role beside the structure is what preserves element binding: "pyridine as
the solvent" stays a condition on a single row rather than an intersection of two
reaction sets, which over-returns. A query the index cannot answer, because it reaches
an element field the index does not carry or needs a projection column to group or sort
by, runs against the projection; both paths screen and verify through the same compiler,
so they differ only in how they reach the reactions.

The grammar bounds what a query can cost (one pass and a sort), and ``search``
takes a wall-clock timeout that interrupts the final query. Name resolution, the
one-time library and index builds, screening, and verification all run before the timer
starts: each is bounded by the corpus rather than by the query, so a slow one is slow
for every caller and shows up in the logs rather than in a timeout. The two builds have
separate triggers -- the library on the first substructure predicate, the index on the
first query the planner routes -- and a search wanting both pays both, upwards of a
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
from ord_schema.agent import query
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Structures read per batch while building the library. Bounds peak memory during the
# build without making the per-batch Python overhead matter.
_BUILD_BATCH = 50_000

# Match sets kept for reuse. Each is a bitmap over the corpus ID space -- 2 MB at ORD's
# scale -- so this trades a few tens of megabytes for the repeat of a query that costs
# seconds.
_CACHED_MATCHES = 16

# How much memory the materialized column sets may hold between them. Held to by
# evicting the least recently used, so a corpus asked many different questions settles
# on the columns it is actually asked about. Enforced after a table is built rather than
# before, since what one costs is known only once it exists, so the peak is this plus
# the largest single set; a set larger than the whole budget is dropped again unkept.
_NARROW_BUDGET_BYTES = 2 * 1024**3

# Top-level projection columns, which is the granularity a compiled query names them at.
# A narrower table holding only the ones a query mentions answers it identically, so the
# query is compiled against it as written.
_TOP_LEVEL = tuple(projection.SCHEMA.names)


# What a structure predicate's answer depends on: the operation, whether the string
# is read as SMARTS or as SMILES, the string, and the similarity threshold.
_MatchKey = tuple[str, bool, str, float | None]


class PairingError(ValueError):
    """The projections and structures artifacts do not form a consistent corpus."""


def _resolve_with_resolvers(name: str) -> str:
    """Resolves a compound name to SMILES through ord_schema.resolvers."""
    smiles, _ = resolvers.resolve_name("name", name)
    return smiles


def _sql_strings(values: Iterable[str]) -> str:
    """Returns ``values`` as a SQL list literal, single quotes escaped."""
    quoted = ", ".join("'" + value.replace("'", "''") + "'" for value in values)
    return f"[{quoted}]"


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

    Walked from the schema and resolved through the compiler, so the index reaches the
    same elements the projection does by construction rather than by a list here
    agreeing with one there.

    Args:
        schema: The projection schema.

    Returns:
        Dotted query paths to the DuckDB expression yielding their elements, for every
        repeated level whose elements carry both a structure and ``_INDEXED_FIELD``.
    """
    paths: dict[str, str] = {}

    def walk(dtype: pa.DataType, prefix: str) -> None:
        if pa.types.is_struct(dtype):
            names = [field.name for field in dtype]
            if "structure_id" in names and _INDEXED_FIELD in names:
                # The schema walk speaks Parquet; the grammar has no wrapper levels.
                path = prefix.replace(".list.element", "").replace(
                    ".key_value.value", ""
                )
                resolved = query.resolve(path, allow_internal=True)
                if resolved.repeated:
                    paths[path] = resolved.expression
            for child in dtype:
                walk(child.type, f"{prefix}.{child.name}")
        elif pa.types.is_list(dtype):
            walk(dtype.value_type, f"{prefix}.list.element")
        elif pa.types.is_map(dtype):
            walk(dtype.item_type, f"{prefix}.key_value.value")

    for field in schema:
        walk(field.type, field.name)
    return paths


INDEXED_PATHS = _indexed_paths()


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


@dataclasses.dataclass(frozen=True)
class _IndexTerms:
    """What a predicate the occurrence index can answer asks of one row.

    Attributes:
        role: The role the element must hold, or None if the query does not bind one.
            Absence is not a wildcard the index applies -- it is a query that named no
            role, so the condition is left off.
    """

    role: str | None


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

        Raises:
            PairingError: If either pattern matches nothing, a file on one side has no
                partner on the other, a pair's stamps disagree about the source
                dataset, or -- with ``require_current`` -- any artifact is stale.
        """
        self._resolver = resolver if resolver is not None else _resolve_with_resolvers
        self._threads = threads
        # Built on first substructure query; see _library. The library holds one entry
        # per distinct molecule, and _members with _starts map an entry back to the
        # structure IDs sharing it.
        self._substructure_library: rdSubstructLibrary.SubstructLibrary | None = None
        self._members = array.array("I")
        self._starts = array.array("I")
        self._library_lock = threading.Lock()
        # Built on the first query that can use it; see _occurrences.
        self._occurrences_built = False
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
        self._narrowed: collections.OrderedDict[frozenset[str], _Narrow] = (
            collections.OrderedDict()
        )
        self._narrow_serial = 0
        self._narrow_lock = threading.Lock()
        pairs = self._pair(projection_pattern, structures_pattern, require_current)
        self._connection = duckdb.connect()
        try:
            self._total, self._searchable = self._prepare(pairs)
        except Exception:
            # Nothing else holds the connection yet, and the caller has no object to
            # close one from if __init__ does not return.
            self._connection.close()
            raise

    @staticmethod
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

    @staticmethod
    def _pair(
        projection_pattern: str, structures_pattern: str, require_current: bool
    ) -> list[tuple[str, str]]:
        """Returns (projection, structures) path pairs, verified by their stamps."""
        projections = Corpus._index(
            projection_pattern, projection.ARTIFACT, require_current
        )
        structure_files = Corpus._index(
            structures_pattern, structures.ARTIFACT, require_current
        )
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
        if parameter.pattern is not None and parameter.op == "substructure":
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

    @staticmethod
    def _index_terms(where: query.Predicate) -> _IndexTerms | None:
        """Returns what an element predicate asks of one occurrence row, or None.

        Args:
            where: The predicate inside a quantifier, relative to the element.

        Returns:
            The terms, if every clause is a structure predicate on the element's own
            structure or an equality on the indexed field, or None if any clause reaches
            something the index does not carry.
        """
        clauses = where.clauses if isinstance(where, query.And) else [where]
        structures_seen = 0
        role: str | None = None
        for clause in clauses:
            if isinstance(clause, (query.Substructure, query.Similarity)):
                if clause.path != "smiles":
                    return None
                structures_seen += 1
            elif (
                isinstance(clause, query.Comparison)
                and clause.op == "eq"
                and clause.path == _INDEXED_FIELD
                and isinstance(clause.value.literal, str)
                and role is None
            ):
                role = clause.value.literal
            else:
                return None
        # Two structure predicates on one element are two bitmaps, and one pass over one
        # column cannot apply both to the same row.
        if structures_seen != 1:
            return None
        return _IndexTerms(role=role)

    @staticmethod
    def _plan(request: query.Query) -> query.Compiled | None:
        """Returns a compiled occurrence-index query, or None to use the projection.

        An occurrence row carries a structure and one element field, so the index
        answers exactly one shape: a bare ``exists`` at an indexed path whose body is a
        conjunction of one structure predicate on the element's own ``smiles`` and at
        most one ``reaction_role`` equality, with nothing aggregated and nothing
        ordered. Every other query -- a ``forall``, a negation or a disjunction,
        another element field, a scalar outside the quantifier, a group-by column, an
        ordering -- reaches something an occurrence row does not hold, and the
        projection answers it.

        A ``limit`` rides along. Neither relation orders what a query did not ask to be
        ordered, so a limited query with no ``order_by`` selects some rows of the match
        set on both paths, and which ones differs between them.

        Args:
            request: The query to plan.

        Returns:
            A compiled query over ``occurrences``, or None if the projection has to
            answer it.

        Raises:
            query.QueryError: If the query does not compile. Planning compiles it to
                reach the structure parameter, so a malformed query fails here rather
                than on the projection path.
        """
        where = request.where
        if (
            request.aggregate is not None
            or request.order_by
            or not isinstance(where, query.Quantifier)
            or where.op != "exists"
            or where.path not in INDEXED_PATHS
        ):
            return None
        terms = Corpus._index_terms(where.where)
        if terms is None:
            return None
        compiled = query.compile_query(request)
        # The bitmap and the molecule behind it are the compiler's, so the two paths
        # screen and verify identically and differ only in how they reach the reactions
        # holding a match.
        if len(compiled.structures) != 1:
            return None
        conditions = [
            f"path = '{where.path}'",
            f"get_bit(CAST(${compiled.structures[0].name} AS BITSTRING), "
            "global_id::INTEGER) = 1",
        ]
        if terms.role is not None:
            escaped = terms.role.replace("'", "''")
            conditions.append(f"{_INDEXED_FIELD} = '{escaped}'")
        limit = f" LIMIT {request.limit}" if request.limit is not None else ""
        # S608: every fragment is a schema-derived path, a compiler-issued parameter
        # name, or an escaped literal.
        sql = (
            "SELECT DISTINCT reaction_id FROM occurrences "  # noqa: S608
            f"WHERE {' AND '.join(conditions)}{limit}"
        )
        # Which relation answered is the one thing a result does not show, and the two
        # are supposed to agree, so a disagreement is only ever debugged from here.
        logger.info("the occurrence index answers %s", where.path)
        return dataclasses.replace(compiled, sql=sql)

    def _occurrences(self) -> None:
        """Builds the occurrence index, once, if it is not already built.

        One row per structure occurrence, carrying the corpus-wide ID, the path it sat
        at, and the element's own ``reaction_role``. Keeping the role beside the
        structure is what preserves element binding: "pyridine as the solvent" stays a
        condition on one row rather than an intersection of two reaction sets, which
        over-returns.

        Costs one pass over the projections per indexed path, so it is built when a
        query first wants it rather than at open. Concurrent first queries serialize
        here and share the result.

        Raises:
            PairingError: If a corpus holding structures indexes none of them. Every
                path the projection can carry a structure at is indexed, so an empty
                index means a traversal reached nothing -- and an empty index answers
                every structure query with "no matches", which reads like an answer
                rather than like a corpus that cannot be searched.
        """
        with self._occurrences_lock:
            if self._occurrences_built:
                return
            start = time.perf_counter()
            # A path's elements are already a list, so unnest yields one row each. The
            # offset comes from the row's own file, which the reactions view carries.
            selects = "\nUNION ALL\n".join(
                f"""
                SELECT reaction_id, '{path}' AS path,
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
                # OR REPLACE, so a build interrupted after the table exists is a build
                # the next query repeats rather than one that collides with itself
                # forever. S608: the fragments are this module's own schema walk and
                # the compiler's traversals, not anything a query supplies.
                cursor.execute(f"CREATE OR REPLACE TABLE occurrences AS {selects}")
                indexed = cursor.execute(
                    "SELECT path, count(*) FROM occurrences GROUP BY path"
                ).fetchall()
                counts = dict(indexed)
                total = sum(counts.values())
                if self._total and not total:
                    raise PairingError(
                        f"the occurrence index came out empty over {self._total} "
                        f"structures; none of {sorted(INDEXED_PATHS)} reached an "
                        "element, so the projections are not the schema this walk was "
                        "built from"
                    )
                self._occurrences_built = True
            finally:
                cursor.close()
            # Per path, so one that reaches nothing is visible here rather than only in
            # the answers it fails to contribute to.
            logger.info(
                "indexed %d structure occurrences in %.1fs: %s",
                total,
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

    @staticmethod
    def _mentioned(sql: str) -> frozenset[str]:
        """Returns the top-level projection columns a compiled query names.

        Read back off the SQL rather than walked from the query, because the SQL is
        what has to resolve: a column named there and missing from the table is a
        catalog error, and one named only in the query is nothing at all.

        Matched as a word anywhere in the SQL, so this errs wide: a nested field sharing
        a leaf name with a top-level column names both -- ``smiles`` is a component's
        field and a reaction's own column -- and a name inside a string literal counts
        too. Each costs a column that gets materialized and never read.

        Args:
            sql: The compiled query.

        Returns:
            The columns the query mentions, always including ``reaction_id``.
        """
        mentioned = {
            column
            for column in _TOP_LEVEL
            if re.search(rf"\b{re.escape(column)}\b", sql)
        }
        mentioned.add("reaction_id")
        return frozenset(mentioned)

    @staticmethod
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

    def _materialize(self, columns: frozenset[str]) -> _Narrow | None:
        """Returns a held table of ``columns``, building it if the cache lacks one.

        The caller owns a read on whatever comes back and has to release it, which
        ``_narrowed_table`` does; an entry with a reader is never evicted.

        Args:
            columns: Top-level projection columns the query names.

        Returns:
            The entry, or None if this column set costs more than the cache may hold in
            total, in which case nothing is kept and the projection answers directly.
        """
        with self._narrow_lock:
            cached = self._narrowed.get(columns)
            if cached is not None:
                self._narrowed.move_to_end(columns)
                cached.readers += 1
                return cached
            # Never reused, so a build that fails between the CREATE and the entry
            # cannot leave a name the next attempt collides with.
            self._narrow_serial += 1
            name = f"narrow_{self._narrow_serial}"
            selected = ", ".join([*sorted(columns), query.STRUCTURE_OFFSET])
            start = time.perf_counter()
            # Its own cursor, for the same reason the index build takes one: searches
            # are in flight, and the shared connection holds their results.
            cursor = self._connection.cursor()
            try:
                before = self._memory_bytes(cursor)
                # S608: the names come from the projection schema and this module.
                cursor.execute(
                    f"CREATE TABLE {name} AS "  # noqa: S608
                    f"SELECT {selected} FROM {query.TABLE}"
                )
                # Measured under the lock, so no other materialization moves the figure
                # between the two reads. An index build can, which costs an
                # overstatement and an eviction, never a wrong answer.
                held = max(self._memory_bytes(cursor) - before, 0)
                if held > _NARROW_BUDGET_BYTES:
                    # Too big to keep, and keeping it is the only reason to build it.
                    cursor.execute(f"DROP TABLE {name}")
                    logger.info(
                        "not caching %s: %.1f GB exceeds the budget",
                        sorted(columns),
                        held / 1024**3,
                    )
                    return None
            except Exception:
                # A table nobody tracks is memory nobody frees, and the failure the
                # caller sees has to be the one that happened.
                with contextlib.suppress(duckdb.Error):
                    cursor.execute(f"DROP TABLE IF EXISTS {name}")
                raise
            finally:
                cursor.close()
            logger.info(
                "materialized %s as %s, %.2f GB in %.1fs",
                sorted(columns),
                name,
                held / 1024**3,
                time.perf_counter() - start,
            )
            entry = _Narrow(name=name, held=held, readers=1)
            self._narrowed[columns] = entry
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
        for held in list(self._narrowed):
            if sum(entry.held for entry in self._narrowed.values()) <= (
                _NARROW_BUDGET_BYTES
            ):
                return
            entry = self._narrowed[held]
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
            del self._narrowed[held]
            logger.info("evicted the materialized %s", sorted(held))

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
            The table's name, or None if materializing it would cost more memory than
            the whole cache is allowed, in which case the projection answers directly.
        """
        entry = self._materialize(columns)
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
        # of the key: _query_molecule reads a stated substructure pattern as SMARTS and
        # everything else as SMILES, and one text is two query molecules across them --
        # C1=CC=CC=C1 is benzene as SMILES and six aliphatic carbons as SMARTS.
        # Resolvers answer in that Kekule form, so a name and a pattern do collide.
        from_smarts = parameter.pattern is not None and parameter.op == "substructure"
        key = (parameter.op, from_smarts, pattern, parameter.threshold)
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

        A query the occurrence index can answer runs against it instead of the
        projection, which is the same answer reached without scanning every reaction.
        The screening and verification are the compiler's either way, so the two paths
        differ only in how they find the reactions holding a match. Which one answered
        is logged, because it is the one thing about a search that the result does not
        show.

        A query the index declines is compiled a second time, against a table holding
        only the columns the first compilation named. The narrow table is held for the
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
            PairingError: If the corpus IDs do not come out one unbroken run, or a
                corpus holding structures indexes none of them.
            TimeoutError: If the query exceeds ``timeout_seconds``.
        """
        indexed = self._plan(request)
        compiled = indexed if indexed is not None else query.compile_query(request)
        # Cached across the whole search: a compound named in both a value and a
        # structure predicate is one external lookup, not two.
        resolved: dict[str, str] = {}

        def resolve(name: str) -> str:
            if name not in resolved:
                resolved[name] = self._resolver(name)
            return resolved[name]

        with contextlib.ExitStack() as reading:
            if indexed is not None:
                self._occurrences()
            else:
                logger.info("the projection answers this query")
                # The index answers by itself; everything else reads the projection, and
                # reads it faster from a table holding only the columns it names. Which
                # ones those are is read off the compiled SQL, so the second compilation
                # is over the same query against a different relation -- no rewriting of
                # SQL text, which a query whose own string literal named the relation
                # would otherwise corrupt.
                narrow = reading.enter_context(
                    self._narrowed_table(self._mentioned(compiled.sql))
                )
                if narrow is not None:
                    compiled = query.compile_query(request, table=narrow)
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
                return self._run_with_timeout(
                    cursor, compiled.sql, parameters, timeout_seconds
                )
            finally:
                cursor.close()

    @staticmethod
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
