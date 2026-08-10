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

The grammar bounds what a query can cost (one pass and a sort), and ``search``
takes a wall-clock timeout that interrupts the final query. Name resolution, screening,
and verification run before the timer starts: each is bounded by the corpus rather than
by the query, so a slow one is slow for every caller and shows up in the logs rather
than in a timeout.
"""

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
                per core. These are real threads with the GIL released, so a server
                may share one corpus across concurrent requests.
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
        self._substruct_library: rdSubstructLibrary.SubstructLibrary | None = None
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
        self._connection.register(
            "structure_offsets",
            pa.table(
                {
                    "projection_filename": [offset[0] for offset in offsets],
                    "structures_filename": [offset[1] for offset in offsets],
                    query.STRUCTURE_OFFSET: pa.array(
                        [offset[2] for offset in offsets], type=pa.int64()
                    ),
                }
            ),
        )
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
        # every query with silence. Counting once here catches it whatever the cause.
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

        Returns:
            A library over every structure in the corpus, screened by the pattern
            fingerprints the artifact already stores.
        """
        if self._substruct_library is None:
            start = time.perf_counter()
            molecules = rdSubstructLibrary.CachedMolHolder()
            patterns = rdSubstructLibrary.PatternHolder(structures.PATTERN_FP_SIZE)
            reader = self._connection.execute(
                "SELECT mol_binary, pattern_fp FROM corpus_structures "
                "ORDER BY global_id"
            ).to_arrow_reader(_BUILD_BATCH)
            for batch in reader:
                blobs = batch.column("mol_binary").to_pylist()
                fingerprints = batch.column("pattern_fp").to_pylist()
                for blob, fingerprint in zip(blobs, fingerprints, strict=True):
                    if blob is None:
                        molecules.AddBinary(_UNPARSEABLE)
                        patterns.AddFingerprint(_NO_BITS)
                    else:
                        molecules.AddBinary(blob)
                        patterns.AddFingerprint(
                            DataStructs.CreateFromBinaryText(fingerprint)
                        )
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
            self._substruct_library = library
        return self._substruct_library

    def _substructure_ids(
        self, parameter: query.StructureParameter, resolve: Callable[[str], str]
    ) -> list[int]:
        """Screens and verifies a substructure predicate; returns global IDs.

        RDKit screens and matches in one call, across its own threads with the GIL
        released, so a server handling concurrent requests neither forks nor
        serializes on this.
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
        self, parameter: query.StructureParameter, resolve: Callable[[str], str]
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
        rows = self._connection.execute(
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

        Args:
            request: The query to run.
            timeout_seconds: Wall-clock bound on the final query. Name resolution,
                screening, and verification run before the timer starts and are not
                counted. None runs unbounded.

        Returns:
            The selected columns: ``reaction_id`` for a plain query, the group and
            measure columns for an aggregated one.

        Raises:
            query.QueryError: If the query does not compile.
            ValueError: If a compound name cannot be resolved.
            TimeoutError: If the query exceeds ``timeout_seconds``.
        """
        compiled = query.compile_query(request)
        # Cached across the whole search: a compound named in both a value and a
        # structure predicate is one external lookup, not two.
        resolved: dict[str, str] = {}

        def resolve(name: str) -> str:
            if name not in resolved:
                resolved[name] = self._resolver(name)
            return resolved[name]

        parameters: dict[str, Any] = {
            name: resolve(name) for name in compiled.compounds
        }
        for parameter in compiled.structures:
            if parameter.op == "substructure":
                matched = self._substructure_ids(parameter, resolve)
            else:
                matched = self._similarity_ids(parameter, resolve)
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
            return self._connection.execute(compiled.sql, parameters).to_arrow_table()
        return self._run_with_timeout(compiled.sql, parameters, timeout_seconds)

    def _run_with_timeout(
        self, sql: str, parameters: dict[str, Any], timeout_seconds: float
    ) -> pa.Table:
        """Runs ``sql``, interrupting it if it outlasts ``timeout_seconds``.

        ``Timer.cancel`` only sets a flag, so a timer that has already passed its own
        check fires anyway. The lock makes the interrupt and the teardown exclusive: it
        either lands while this query owns the connection, or not at all. Without it a
        late interrupt reaches whatever query runs next and reports that one as having
        timed out.

        Args:
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
                    self._connection.interrupt()

        timer = threading.Timer(timeout_seconds, interrupt)
        timer.start()
        try:
            return self._connection.execute(sql, parameters).to_arrow_table()
        except duckdb.InterruptException as error:
            raise TimeoutError(f"query exceeded {timeout_seconds} seconds") from error
        finally:
            timer.cancel()
            with lock:
                running = False
