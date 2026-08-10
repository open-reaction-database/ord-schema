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

"""A structure-search artifact: one row per distinct structure in a projection.

Structure search runs in two steps that this artifact serves in turn. A **screen**
tests the query's fingerprint bits against every row -- a bitwise scan a query engine
runs in milliseconds, with no chemistry library in sight -- and is complete but not
exact: it never misses a true match, and passes false ones. **Verification** then runs
real subgraph isomorphism over the survivors, from the serialized molecule carried
here, which skips re-parsing SMILES and runs about five times faster for it. Neither
step needs an index; the reason a structure search elsewhere needs one is the cost of
scanning rows one at a time, and a columnar scan of every distinct structure in the
corpus is cheap.

The artifact is deduplicated, which is what makes both steps affordable: the corpus
holds roughly fifteen component occurrences per distinct structure, so screening and
verification run against the distinct set and their results fan back out through the
projection. The join is ``structure_id`` -- assigned by ``write_projection`` in
first-seen order, carried on every compound there that records a structure and numbered
``0..n-1`` here, so a match set travels back into a query as a bitmap indexed by ID.
The two files are two statements of one derivation and are meaningful only together:
IDs are not stable across builds, nothing outside the artifacts should record them, and
a projection rewritten in place needs its structures artifact rederived with it -- the
currency stamps name the source dataset rather than the projection file, so the skip
check cannot see that the pairing changed.

The screen's completeness is an invariant, not a hope: if query ``Q`` is a subgraph of
target ``T``, every bit of ``Q``'s pattern fingerprint is set in ``T``'s, so
``bit_count(fp & q) = popcount(q)`` never rejects a true match -- measured across 221
query/target pairs, including explicit-hydrogen queries, with no exception. What a
query layer still has to handle is that these molecules are built from SMILES and so
hold their hydrogens implicitly: a SMARTS naming one as its own atom matches nothing,
which :mod:`ord_schema.agent.query` rewrites rather than runs. Similarity needs no
verification at all: Tanimoto
is defined on the Morgan fingerprint, so the screen *is* the answer, and
``morgan_popcount`` bounds it (``popcount(B)`` must lie within ``[t * popcount(A),
popcount(A) / t]`` for Tanimoto ``>= t``).

A structure whose SMILES RDKit cannot parse keeps its row -- the ID space must stay
dense -- with every derived column null. It can never match a structure query, and the
count of such rows is logged per dataset.
"""

import os
from collections.abc import Iterator
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from ord_schema import artifacts, atomic_io, projection
from ord_schema.logging import get_logger

logger = get_logger(__name__)

ARTIFACT = "structures"

# Fingerprint parameters. Screen precision plateaus by 2048 bits -- the residual false
# positives are the gap between feature containment and subgraph isomorphism, which no
# width closes -- so wider fingerprints buy storage, not selectivity.
PATTERN_FP_SIZE = 2048
MORGAN_FP_SIZE = 2048
MORGAN_RADIUS = 2

_MORGAN = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS, fpSize=MORGAN_FP_SIZE
)


def morgan_fingerprint(molecule: Chem.Mol) -> DataStructs.ExplicitBitVect:
    """Returns the Morgan fingerprint with the parameters this artifact stores.

    The query side of a similarity search must fingerprint identically to the artifact
    or the comparison is meaningless, so this is the one place the parameters live.

    Args:
        molecule: The molecule to fingerprint.

    Returns:
        A ``MORGAN_FP_SIZE``-bit fingerprint at radius ``MORGAN_RADIUS``.
    """
    return _MORGAN.GetFingerprint(molecule)


# Named so a published file says how its fingerprints were made; the reader recovers
# the parameters from the footer rather than from this library's history.
_META_FINGERPRINT = "ord.fingerprint"

SCHEMA = pa.schema(
    [
        # The join key: every compound in the sibling projection with a ``smiles``
        # carries one, and rows here are numbered 0..n-1 in the same first-seen order.
        pa.field("structure_id", pa.uint32(), nullable=False),
        pa.field("smiles", pa.string(), nullable=False),
        pa.field(
            "pattern_fp",
            pa.binary(),
            metadata={_META_FINGERPRINT: f"rdkit-pattern;fpSize={PATTERN_FP_SIZE}"},
        ),
        pa.field(
            "morgan_fp",
            pa.binary(),
            metadata={
                _META_FINGERPRINT: (
                    f"morgan;radius={MORGAN_RADIUS};fpSize={MORGAN_FP_SIZE}"
                )
            },
        ),
        pa.field("morgan_popcount", pa.uint32()),
        # Chem.Mol(mol_binary) reconstructs the molecule without re-parsing SMILES,
        # which is what verification spends most of its time on otherwise.
        pa.field("mol_binary", pa.binary()),
    ]
)


def _structure_columns(schema: pa.Schema = projection.SCHEMA) -> list[str]:
    """Returns the projection's (smiles, structure_id) leaves as Parquet column paths.

    Walked from the schema rather than listed by hand, so a compound field added
    upstream is read from wherever it lands -- components, products, workup inputs,
    authentic standards -- without anyone updating a list here.

    Args:
        schema: The projection schema.

    Returns:
        Parquet physical column paths, in schema order.
    """
    paths: list[str] = []

    def walk(dtype: pa.DataType, prefix: str) -> None:
        if pa.types.is_struct(dtype):
            if "structure_id" in [field.name for field in dtype]:
                paths.append(f"{prefix}.smiles")
                paths.append(f"{prefix}.structure_id")
            for child in dtype:
                walk(child.type, f"{prefix}.{child.name}")
        elif pa.types.is_list(dtype):
            walk(dtype.value_type, f"{prefix}.list.element")
        elif pa.types.is_map(dtype):
            walk(dtype.item_type, f"{prefix}.key_value.value")

    for field in schema:
        walk(field.type, field.name)
    return paths


def _pairs(value: Any) -> Iterator[tuple[str | None, int | None]]:
    """Yields every (smiles, structure_id) pair within a pruned projection row.

    Args:
        value: A projection row read through ``_structure_columns``, or any value
            within one: structs arrive as dicts, lists as lists, and maps as lists of
            ``(key, value)`` tuples.

    Yields:
        One pair per compound struct, however deeply nested.
    """
    if isinstance(value, dict):
        if "structure_id" in value:
            yield value.get("smiles"), value["structure_id"]
        for child in value.values():
            yield from _pairs(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _pairs(child)


def _collect(source: str | os.PathLike[str]) -> list[str]:
    """Reads the (smiles, structure_id) mapping back out of a projection.

    Args:
        source: Path to the projection.

    Returns:
        The distinct SMILES, indexed by ``structure_id``.

    Raises:
        ValueError: If ``source``'s schema does not hold the structure columns this
            library's projection defines, or if the pairs do not state a single dense
            one-to-one mapping -- an ID bound to two SMILES, one SMILES under two IDs,
            a compound with one half of the pair, or a gap in the ID space. Each means
            ``source`` was not written by this library's ``write_projection``, and an
            artifact derived from it would join wrongly rather than fail.
    """
    smiles_by_id: dict[int, str] = {}
    id_by_smiles: dict[str, int] = {}
    columns = _structure_columns()
    with pq.ParquetFile(source) as projected:
        # Checked against the file's own schema because read_row_group silently drops
        # requested columns the file does not have: an old-schema projection would
        # otherwise state no pairs at all and derive an empty artifact that stamps
        # itself current.
        available = _structure_columns(projected.schema_arrow)
        if available != columns:
            raise ValueError(
                f"{source}: structure columns disagree with this library's projection "
                f"schema (missing {sorted(set(columns) - set(available))}, unexpected "
                f"{sorted(set(available) - set(columns))}); derive the projection "
                "again first"
            )
        for row_group in range(projected.num_row_groups):
            table = projected.read_row_group(row_group, columns=columns)
            for row in table.to_pylist():
                for smiles, structure_id in _pairs(row):
                    if smiles is None and structure_id is None:
                        continue  # A compound whose structure the source never stated.
                    if smiles is None or structure_id is None:
                        raise ValueError(
                            f"{source}: compound with smiles={smiles!r} but "
                            f"structure_id={structure_id!r}; the projection assigns "
                            "both or neither"
                        )
                    known = smiles_by_id.setdefault(structure_id, smiles)
                    if known != smiles:
                        raise ValueError(
                            f"{source}: structure_id {structure_id} is both "
                            f"{known!r} and {smiles!r}"
                        )
                    known_id = id_by_smiles.setdefault(smiles, structure_id)
                    if known_id != structure_id:
                        raise ValueError(
                            f"{source}: {smiles!r} is both structure_id {known_id} "
                            f"and {structure_id}; the artifact would hold one "
                            "structure twice"
                        )
    if set(smiles_by_id) != set(range(len(smiles_by_id))):
        missing = sorted(set(range(len(smiles_by_id))) - set(smiles_by_id))[:5]
        raise ValueError(
            f"{source}: structure IDs are not dense; first missing {missing}"
        )
    return [smiles_by_id[structure_id] for structure_id in range(len(smiles_by_id))]


def _row(structure_id: int, smiles: str) -> dict[str, Any]:
    """Featurizes one structure into a row matching ``SCHEMA``.

    Args:
        structure_id: The row's ID, equal to its position in the artifact.
        smiles: The structure, as the projection derived it.

    Returns:
        A dict keyed by column name; every derived column is None when RDKit cannot
        parse ``smiles``.
    """
    row: dict[str, Any] = {
        "structure_id": structure_id,
        "smiles": smiles,
        "pattern_fp": None,
        "morgan_fp": None,
        "morgan_popcount": None,
        "mol_binary": None,
    }
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return row
    morgan = morgan_fingerprint(mol)
    row["pattern_fp"] = DataStructs.BitVectToBinaryText(
        Chem.PatternFingerprint(mol, fpSize=PATTERN_FP_SIZE)
    )
    row["morgan_fp"] = DataStructs.BitVectToBinaryText(morgan)
    row["morgan_popcount"] = morgan.GetNumOnBits()
    row["mol_binary"] = mol.ToBinary()
    return row


def is_current(path: str | os.PathLike[str], source_md5: str) -> bool:
    """Returns whether ``path`` is a structures artifact of ``source_md5`` by us.

    Delegates to ``artifacts.is_current``, which requires the artifact name, the source
    content, the library version, and the artifact version to all match. A missing or
    unreadable file is not current.

    Args:
        path: Path to check. Need not exist.
        source_md5: Hash of the source dataset the artifact should restate, which is
            the dataset's own even though the artifact derives from a projection.

    Returns:
        Whether the file is a structures artifact this library would write today.
    """
    return artifacts.is_current(path, ARTIFACT, source_md5)


def write_structures(
    source: str | os.PathLike[str],
    output: str | os.PathLike[str],
    *,
    source_md5: str | None = None,
    source_dataset_id: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 10_000,
) -> int:
    """Derives a structures artifact from a projection and writes it.

    Only the smiles and structure_id leaves are decoded, one projection row group at a
    time. The output is published atomically, so a failure partway leaves any existing
    artifact untouched.

    Args:
        source: Path to the projection whose structures to featurize.
        output: Path to write the artifact to.
        source_md5: Hash of the *source dataset* to stamp, if the caller already read
            one. Taken from the projection's own stamps when omitted, so the artifact
            names the dataset it reflects rather than the file it read.
        source_dataset_id: Source dataset ID to stamp, if the caller already read one.
            Taken from the projection's stamps when omitted.
        compression: Parquet compression codec.
        row_group_size: Rows per output row group.

    Returns:
        Number of rows written: the number of distinct structures in the projection.

    Raises:
        ValueError: If ``source`` is not a projection, or does not state a single dense
            (smiles, structure_id) mapping.
    """
    parent = artifacts.load_stamps(source)
    if parent.artifact != projection.ARTIFACT:
        raise ValueError(
            f"{source} is a {parent.artifact}, not a {projection.ARTIFACT}; the "
            "structures artifact featurizes a projection"
        )
    # derive_tree refuses stale parents, but this writer is public and its output
    # inherits the dataset hash: an artifact derived from a stale projection would
    # claim a provenance it does not have and nothing would ever mark it stale again.
    if not artifacts.stamps_are_current(parent, projection.ARTIFACT):
        raise ValueError(
            f"{source} is a stale {projection.ARTIFACT}; derive it again first"
        )
    if source_md5 is None:
        source_md5 = parent.source_md5
    if source_dataset_id is None:
        source_dataset_id = parent.source_dataset_id
    stamps = artifacts.current_stamps(ARTIFACT, source_dataset_id, source_md5)
    schema = SCHEMA.with_metadata(artifacts.to_metadata(stamps))
    all_smiles = _collect(source)
    unparseable = 0
    with (
        atomic_io.atomic_path(output) as temp_path,
        pq.ParquetWriter(temp_path, schema, compression=compression) as writer,
    ):
        # A projection with no structures leaves the loop empty; the writer's close
        # still produces a valid zero-row file carrying the schema and stamps.
        for start in range(0, len(all_smiles), row_group_size):
            batch = [
                _row(structure_id, smiles)
                for structure_id, smiles in enumerate(
                    all_smiles[start : start + row_group_size], start=start
                )
            ]
            unparseable += sum(row["mol_binary"] is None for row in batch)
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))
    if unparseable:
        # A warning because any nonzero count is anomalous: the projection's SMILES
        # are RDKit-canonical, so a re-parse failure means version skew or an RDKit
        # bug. Per dataset rather than per row so the count is diffable between runs.
        logger.warning("%s: %d structures do not parse", source, unparseable)
    return len(all_smiles)
