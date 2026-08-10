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
by, runs against the projection unchanged; both paths screen and verify through the same
compiler, so they differ only in how they reach the reactions.

The grammar bounds what a query can cost (one pass and a sort), and ``search``
takes a wall-clock timeout that interrupts the final query. Name resolution, the
one-time library and index builds, screening, and verification all run before the timer
starts: each is bounded by the corpus rather than by the query, so a slow one is slow
for every caller and shows up in the logs rather than in a timeout. The first
substructure search therefore takes both builds -- upwards of a minute over the whole
corpus, the index being four passes over the projections -- on top of whatever timeout
it was given.
"""

import dataclasses
import glob
import math
import threading
import time
from collections.abc import Callable, Iterable, Sequence
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


# An ID whose structure did not parse still occupies a slot, so a library index is a
# corpus-wide structure ID and needs no translation. An empty molecule fingerprints to
# no bits, and a query with no atoms is refused, so it can never match.
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
        # Built on first substructure query; see _library.
        self._substructure_library: rdSubstructLibrary.SubstructLibrary | None = None
        self._library_lock = threading.Lock()
        # Built on the first query that can use it; see _occurrences.
        self._occurrences_built = False
        self._occurrences_lock = threading.Lock()
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
        for similarity or scalar queries never needs it. Streamed in ID order so a
        library index *is* a corpus-wide structure ID: a structure that did not parse
        still occupies its slot, holding an empty molecule that nothing matches.

        Concurrent first searches serialize on the build rather than each building a
        copy, which at corpus scale is about 2 GB apiece. Whoever waits needs the
        finished library anyway.

        Returns:
            A library over every structure in the corpus, screened by the pattern
            fingerprints the artifact already stores.

        Raises:
            PairingError: If the corpus IDs are not exactly ``0`` to ``total - 1``, so
                that a library index would name a different structure, or if a
                structure holds one of its derived columns without the other.
        """
        with self._library_lock:
            if self._substructure_library is not None:
                return self._substructure_library
            start = time.perf_counter()
            molecules = rdSubstructLibrary.CachedMolHolder()
            patterns = rdSubstructLibrary.PatternHolder(structures.PATTERN_FP_SIZE)
            cursor = self._connection.cursor()
            position = 0
            try:
                reader = cursor.execute(
                    "SELECT global_id, mol_binary, pattern_fp FROM corpus_structures "
                    "ORDER BY global_id"
                ).to_arrow_reader(_BUILD_BATCH)
                for batch in reader:
                    identifiers = batch.column("global_id").to_pylist()
                    blobs = batch.column("mol_binary").to_pylist()
                    fingerprints = batch.column("pattern_fp").to_pylist()
                    for global_id, blob, fingerprint in zip(
                        identifiers, blobs, fingerprints, strict=True
                    ):
                        # Sorting alone does not make an index an ID; it does that only
                        # while the IDs are every integer from zero, each once. A
                        # duplicate or a gap slides every later structure onto a
                        # neighbor's index, which a row count cannot see.
                        if global_id != position:
                            raise PairingError(
                                f"the corpus states structure ID {global_id} where "
                                f"{position} was expected, so its IDs are not one "
                                "unbroken run and a library index would name a "
                                "different structure; derive the structures artifacts "
                                "again"
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
                        if blob is None:
                            molecules.AddBinary(_UNPARSEABLE)
                            patterns.AddFingerprint(_NO_BITS)
                        else:
                            molecules.AddBinary(blob)
                            patterns.AddFingerprint(
                                DataStructs.CreateFromBinaryText(fingerprint)
                            )
                        position += 1
            finally:
                cursor.close()
            library = rdSubstructLibrary.SubstructLibrary(molecules, patterns)
            if len(library) != self._total:
                raise PairingError(
                    f"the library holds {len(library)} structures but the corpus has "
                    f"{self._total}; its indices would not be structure IDs"
                )
            logger.info(
                "built a substructure library over %d structures in %.1fs",
                len(library),
                time.perf_counter() - start,
            )
            self._substructure_library = library
            return self._substructure_library

    @staticmethod
    def _index_terms(where: query.Predicate) -> tuple[bool, str | None] | None:
        """Returns whether an element predicate is one the index holds.

        Args:
            where: The predicate inside a quantifier, relative to the element.

        Returns:
            ``(has_structure, role)`` if every clause is a structure predicate on the
            element's own structure or an equality on the indexed field, or None if any
            clause reaches something the index does not carry.
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
        return True, role

    @staticmethod
    def _plan(request: query.Query) -> query.Compiled | None:
        """Returns a compiled occurrence-index query, or None to use the projection.

        The index carries a structure and one field per element, so it answers exactly
        the shape it holds: some element at an indexed path contains a structure, and
        optionally equals a role. Everything else -- another element field, a scalar
        outside the quantifier, a group-by column, an ordering -- lives only in the
        projection, and asking the index for it would answer a different question.

        Args:
            request: The query to plan.

        Returns:
            A compiled query over ``occurrences``, or None if the projection has to
            answer it.
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
        _, role = terms
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
        if role is not None:
            escaped = role.replace("'", "''")
            conditions.append(f"{_INDEXED_FIELD} = '{escaped}'")
        limit = f" LIMIT {request.limit}" if request.limit is not None else ""
        # S608: every fragment is a schema-derived path, a compiler-issued parameter
        # name, or an escaped literal.
        sql = (
            "SELECT DISTINCT reaction_id FROM occurrences "  # noqa: S608
            f"WHERE {' AND '.join(conditions)}{limit}"
        )
        return dataclasses.replace(compiled, sql=sql)

    def _occurrences(self) -> None:
        """Builds the occurrence index, once, if it is not already built.

        One row per structure occurrence, carrying the corpus-wide ID, the path it sat
        at, and the element's own ``reaction_role``. Keeping the role beside the
        structure is what preserves element binding: "pyridine as the solvent" stays a
        condition on one row rather than an intersection of two reaction sets, which
        over-returns.

        Costs a full pass over the projections, so it is built when a query first wants
        it rather than at open. Concurrent first queries serialize here and share the
        result.
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
            # S608: the fragments are this module's own schema walk and the compiler's
            # traversals, not anything a query supplies.
            self._connection.execute(f"CREATE TABLE occurrences AS {selects}")
            count = self._connection.execute(
                "SELECT count(*) FROM occurrences"
            ).fetchone()
            assert count is not None  # An aggregate over any relation returns one row.
            logger.info(
                "indexed %d structure occurrences in %.1fs",
                count[0],
                time.perf_counter() - start,
            )
            self._occurrences_built = True

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
        # matches hundreds of thousands of structures in ORD.
        return list(
            library.GetMatches(
                molecule, numThreads=self._threads, maxResults=len(library) or 1
            )
        )

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

    def search(
        self, request: query.Query, *, timeout_seconds: float | None = None
    ) -> pa.Table:
        """Compiles and runs a query, returning the result as an Arrow table.

        Runs on its own cursor, so concurrent searches sharing this corpus do not read
        each other's results. They still serialize on the library build, and only on
        that: the first substructure search builds it while the others wait.

        A query the occurrence index can answer runs against it instead of the
        projection, which is the same answer reached without scanning every reaction.
        The screening and verification are the compiler's either way, so the two paths
        differ only in how they find the reactions holding a match.

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
            PairingError: If the library does not come out one entry per structure.
            TimeoutError: If the query exceeds ``timeout_seconds``.
        """
        indexed = self._plan(request)
        compiled = indexed if indexed is not None else query.compile_query(request)
        if indexed is not None:
            self._occurrences()
        # Cached across the whole search: a compound named in both a value and a
        # structure predicate is one external lookup, not two.
        resolved: dict[str, str] = {}

        def resolve(name: str) -> str:
            if name not in resolved:
                resolved[name] = self._resolver(name)
            return resolved[name]

        cursor = self._connection.cursor()
        try:
            parameters: dict[str, Any] = {
                name: resolve(name) for name in compiled.compounds
            }
            for parameter in compiled.structures:
                if parameter.op == "substructure":
                    matched = self._substructure_ids(parameter, resolve)
                else:
                    matched = self._similarity_ids(cursor, parameter, resolve)
                unsearchable = self._total - self._searchable
                logger.info(
                    "%s %r matched %d of %d structures (%d unsearchable)",
                    parameter.op,
                    parameter.pattern
                    if parameter.pattern is not None
                    else parameter.compound,
                    len(matched),
                    self._searchable,
                    unsearchable,
                )
                parameters[parameter.name] = self._bitmap(matched)
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
