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
  element's ``structure_id``; the chemistry -- a fingerprint screen over the structures
  artifact, then exact subgraph verification from the serialized molecules -- runs
  here, and the verified match set binds as a ``BITSTRING`` parameter. Similarity has
  no verification step: Tanimoto is defined on the Morgan fingerprint, so the screen is
  the whole answer.

Structure ids are dataset-local, so a corpus-wide bitmap needs an offset per file.
``Corpus`` pairs every projection with its structures artifact, refuses a pair whose
stamps disagree about the source dataset -- the two files are two statements of one
derivation, and a mismatched pair would join ids to the wrong molecules silently --
and publishes a ``reactions`` relation carrying each row's offset in the column the
compiled SQL expects.

The grammar bounds what a query can cost (one pass and a sort), so the remaining
runaway risk is corpus size times a slow filter; ``search`` takes a wall-clock timeout
that interrupts the query rather than trusting it.
"""

import glob
import math
import pathlib
import threading
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Self

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from joblib import Parallel, delayed
from rdkit import Chem, DataStructs

from ord_schema import artifacts, projection, resolvers, structures
from ord_schema.agent import query
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Structures verified per joblib task. Large enough that pickling blobs is amortized,
# small enough that a dozen workers stay busy on a typical survivor set.
_VERIFY_CHUNK = 4096


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


def _verify_chunk(
    pattern: str, from_smarts: bool, blobs: Sequence[bytes]
) -> list[bool]:
    """Returns whether each serialized molecule contains the query as a subgraph.

    A module-level function so joblib can ship it to worker processes; the query is
    rebuilt from its string per chunk because a mol does not cross the boundary.

    Args:
        pattern: The query, as SMARTS or SMILES.
        from_smarts: Whether ``pattern`` is SMARTS.
        blobs: Serialized molecules (``mol_binary`` column values).

    Returns:
        One flag per blob.
    """
    parsed = Chem.MolFromSmarts(pattern) if from_smarts else Chem.MolFromSmiles(pattern)
    return [
        Chem.Mol(bytes(blob)).HasSubstructMatch(parsed)  # ty: ignore[no-matching-overload]
        for blob in blobs
    ]


def _pattern_blob(molecule: Chem.Mol) -> tuple[bytes, int]:
    """Returns the packed pattern fingerprint and its popcount."""
    fingerprint = Chem.PatternFingerprint(molecule, fpSize=structures.PATTERN_FP_SIZE)
    return DataStructs.BitVectToBinaryText(fingerprint), fingerprint.GetNumOnBits()


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
        n_jobs: int = 1,
        require_current: bool = True,
    ) -> None:
        """Pairs the artifacts, checks their stamps, and prepares the relations.

        Args:
            projection_pattern: Glob matching the projection files.
            structures_pattern: Glob matching the structures artifacts.
            resolver: Maps a compound name to SMILES. Defaults to
                ``ord_schema.resolvers``, which calls external services; inject
                something local for tests or offline use.
            n_jobs: Worker processes for substructure verification. 1 verifies inline.
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
        self._n_jobs = n_jobs
        pairs = self._pair(projection_pattern, structures_pattern, require_current)
        self._connection = duckdb.connect()
        self._total = self._prepare(pairs)

    @staticmethod
    def _pair(
        projection_pattern: str, structures_pattern: str, require_current: bool
    ) -> list[tuple[str, str]]:
        """Returns (projection, structures) path pairs, verified by their stamps."""
        projections = {
            pathlib.Path(path).name: path
            for path in glob.glob(projection_pattern, recursive=True)
        }
        structure_files = {
            pathlib.Path(path).name: path
            for path in glob.glob(structures_pattern, recursive=True)
        }
        if not projections:
            raise PairingError(f"no projections matched: {projection_pattern}")
        unpaired = sorted(projections.keys() ^ structure_files.keys())
        if unpaired:
            raise PairingError(
                f"projections and structures artifacts do not pair up; only one side "
                f"has {unpaired}"
            )
        pairs = []
        for name in sorted(projections):
            projected, structured = projections[name], structure_files[name]
            left = artifacts.load_stamps(projected)
            right = artifacts.load_stamps(structured)
            for stamps, path, artifact in (
                (left, projected, projection.ARTIFACT),
                (right, structured, structures.ARTIFACT),
            ):
                if stamps.artifact != artifact:
                    raise PairingError(
                        f"{path} is a {stamps.artifact}, not a {artifact}"
                    )
                if require_current and not artifacts.stamps_are_current(
                    stamps, artifact
                ):
                    raise PairingError(f"{path} is stale; derive it again first")
            if left.source_md5 != right.source_md5:
                raise PairingError(
                    f"{projected} and {structured} were derived from different "
                    "sources; their structure ids do not agree"
                )
            pairs.append((projected, structured))
        return pairs

    def _prepare(self, pairs: list[tuple[str, str]]) -> int:
        """Publishes the relations and returns the corpus-wide structure count."""
        offsets = []
        total = 0
        for projected, structured in pairs:
            with pq.ParquetFile(structured) as artifact:
                count = artifact.metadata.num_rows
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
        return total

    def __enter__(self) -> Self:
        """Returns the corpus itself; closing on exit is the whole protocol."""
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Closes the connection."""
        self.close()

    def close(self) -> None:
        """Closes the connection."""
        self._connection.close()

    def _query_molecule(self, parameter: query.StructureParameter) -> Chem.Mol:
        """Returns the query molecule for a structure predicate.

        Args:
            parameter: The predicate to build a molecule for.

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
            smiles = self._resolver(parameter.compound)
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(
                f"resolver returned unreadable SMILES {smiles!r} "
                f"for {parameter.compound!r}"
            )
        return molecule

    def _substructure_ids(self, parameter: query.StructureParameter) -> list[int]:
        """Screens and verifies a substructure predicate; returns global ids."""
        molecule = self._query_molecule(parameter)
        blob, popcount = _pattern_blob(molecule)
        survivors = self._connection.execute(
            """
            SELECT global_id, mol_binary FROM corpus_structures
            WHERE pattern_fp IS NOT NULL
              AND bit_count(CAST(pattern_fp AS BITSTRING) & CAST($q AS BITSTRING)) = $n
            """,
            {"q": blob, "n": popcount},
        ).fetchall()
        if not survivors:
            return []
        pattern = parameter.pattern if parameter.op == "substructure" else None
        from_smarts = pattern is not None
        if not from_smarts:
            pattern = Chem.MolToSmiles(molecule)
        chunks = [
            survivors[start : start + _VERIFY_CHUNK]
            for start in range(0, len(survivors), _VERIFY_CHUNK)
        ]
        verified = Parallel(n_jobs=self._n_jobs)(
            delayed(_verify_chunk)(pattern, from_smarts, [row[1] for row in chunk])
            for chunk in chunks
        )
        return [
            row[0]
            for chunk, flags in zip(chunks, verified, strict=True)
            for row, match in zip(chunk, flags, strict=True)
            if match
        ]

    def _similarity_ids(self, parameter: query.StructureParameter) -> list[int]:
        """Screens a similarity predicate; the fingerprint is the whole answer."""
        molecule = self._query_molecule(parameter)
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
        """Returns the match set as a bitmap over the corpus-wide id space."""
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
            timeout_seconds: Wall-clock bound on the final query (structure screening
                and verification are bounded by the corpus, not by the query, and are
                not counted). None runs unbounded.

        Returns:
            The selected columns: ``reaction_id`` for a plain query, the group and
            measure columns for an aggregated one.

        Raises:
            query.QueryError: If the query does not compile.
            ValueError: If a compound name cannot be resolved.
            TimeoutError: If the query exceeds ``timeout_seconds``.
        """
        compiled = query.compile_query(request)
        parameters: dict[str, Any] = {}
        for name in compiled.compounds:
            parameters[name] = self._resolver(name)
        for parameter in compiled.structures:
            if parameter.op == "substructure":
                matched = self._substructure_ids(parameter)
            else:
                matched = self._similarity_ids(parameter)
            logger.info(
                "%s %r matched %d of %d structures",
                parameter.op,
                parameter.pattern or parameter.compound,
                len(matched),
                self._total,
            )
            parameters[parameter.name] = self._bitmap(matched)
        if timeout_seconds is None:
            return self._connection.execute(compiled.sql, parameters).to_arrow_table()
        timer = threading.Timer(timeout_seconds, self._connection.interrupt)
        timer.start()
        try:
            return self._connection.execute(compiled.sql, parameters).to_arrow_table()
        except duckdb.InterruptException as error:
            raise TimeoutError(f"query exceeded {timeout_seconds} seconds") from error
        finally:
            timer.cancel()
