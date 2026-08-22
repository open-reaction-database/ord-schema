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
runs in milliseconds without a chemistry library -- and is complete but not exact: it
never misses a true match, and passes false ones. **Verification** then runs
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
which :mod:`ord_schema.search.query` rewrites rather than runs. Similarity needs no
verification at all: Tanimoto is defined on the Morgan fingerprint, so the screen *is*
the answer, and ``morgan_popcount`` bounds it (``popcount(B)`` must lie within ``[t *
popcount(A), popcount(A) / t]`` for Tanimoto ``>= t``).

Beside the SMILES the projection derived, a row carries a ``mol_hash``: an identity
for the molecule that ignores how it was drawn. Equality on the stored spelling is
exact and answers the wrong question surprisingly often -- acetic acid and acetate, an
amine and its ammonium, a 2-pyridone and its 2-hydroxypyridine tautomer are each one
reagent written two ways, and each compares unequal.

It is RDKit's registration hash, taken of the uncharged molecule. Uncharging is what
makes the scheme insensitive to protonation as well as to tautomer: one of its layers
is the molecular formula, and a formula states the charge. Registration is the right
API rather than a bare tautomer hash because it normalizes what a bare hash does not --
atom-map labels, unnecessary explicit hydrogens, enhanced stereo, S-group data -- so a
compound written with atom maps is the same compound. It is a hash rather than a
canonical tautomer because canonicalizing one costs what the corpus does not have:
measured over ORD's own molecules, RDKit's tautomer enumerator runs at ~500 structures
a second against ~2,600 for this hash, which is an hour against thirteen minutes over
the two million distinct structures here, before the per-file parallelism the build
already has. The stored SMILES beside it is what a reader looks at; the hash is only
ever compared.

Fragments are left alone, so sodium acetate stays distinct from acetic acid. The
question that counts a reagent and the salt it was sold as as one is a different one,
and it gets its own column: a ``parent_hash``, the same hash taken of the molecule with
its recognized counterions removed. Two columns rather than one choice between them,
because either answer is wrong for the other question -- potassium and caesium carbonate
are one reagent to a query about carbonate and two to a query about the caesium -- and a
corpus-wide column is the wrong place to guess which was meant. See ``salt_parent`` for
where stripping stops.

The hash is RDKit's and moves with it, which is why ``base`` stamps the RDKit version
and ``is_current`` refuses an artifact built by another one. Without that stamp an
RDKit upgrade would silently change what a stored column means while every query kept
running.

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
from rdkit.Chem import RegistrationHash, SaltRemover, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

from ord_schema import atomic_io
from ord_schema.artifacts import base, projection
from ord_schema.logging import get_logger

logger = get_logger(__name__)

ARTIFACT = "structures"

# Screen precision plateaus by 2048 bits -- the residual false positives are the gap
# between feature containment and subgraph isomorphism, which no width closes -- so
# wider fingerprints buy storage, not selectivity.
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
        # What an equality compares when the question is "the same compound" rather
        # than "the same drawing of one"; see the module docstring.
        pa.field("mol_hash", pa.string()),
        # The same, of the molecule without its counterions, for the question that
        # counts a reagent and the salt it was sold as as one.
        pa.field("parent_hash", pa.string()),
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


def structure_levels(
    schema: pa.Schema = projection.SCHEMA,
) -> list[tuple[str, str, pa.DataType]]:
    """Returns the levels of ``schema`` whose elements carry a structure.

    Walked rather than listed by hand, so a compound field added upstream is found
    wherever it lands -- components, products, workup inputs, authentic standards --
    without anyone updating a list. Every caller that has to agree about where a
    structure can sit reads this one walk, since two walks are two answers waiting to
    disagree.

    Args:
        schema: Schema to walk; the projection schema this library defines by default,
            or a written file's ``schema_arrow`` to ask what that file holds.

    Returns:
        One ``(column, path, dtype)`` per structure-bearing struct, in schema order:
        ``column`` is the Parquet physical prefix, carrying the ``list.element`` and
        ``key_value.value`` segments a file has; ``path`` is the same level as the query
        grammar names it, which has no wrapper segments; and ``dtype`` is the struct, so
        a caller can ask what else the element carries.
    """
    levels: list[tuple[str, str, pa.DataType]] = []

    def walk(dtype: pa.DataType, column: str, path: str) -> None:
        if pa.types.is_struct(dtype):
            if "structure_id" in [field.name for field in dtype]:
                levels.append((column, path, dtype))
            for child in dtype:
                walk(child.type, f"{column}.{child.name}", f"{path}.{child.name}")
        elif pa.types.is_list(dtype):
            walk(dtype.value_type, f"{column}.list.element", path)
        elif pa.types.is_map(dtype):
            walk(dtype.item_type, f"{column}.key_value.value", path)

    for field in schema:
        walk(field.type, field.name, field.name)
    return levels


def _structure_columns(schema: pa.Schema = projection.SCHEMA) -> list[str]:
    """Returns the projection's (smiles, structure_id) leaves as Parquet column paths.

    Args:
        schema: Schema to read the leaves from; the projection schema this library
            defines by default, or a written file's ``schema_arrow`` to ask what that
            file holds.

    Returns:
        Parquet physical column paths, in schema order.
    """
    return [
        f"{column}.{leaf}"
        for column, _, _ in structure_levels(schema)
        for leaf in ("smiles", "structure_id")
    ]


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


# Built once: constructing an uncharger per molecule dominates the work it does, and
# it carries no per-molecule state.
_UNCHARGER = rdMolStandardize.Uncharger()
# The registration hash's tautomer-insensitive scheme drops the two layers that spell a
# structure exactly and keeps the formula, the tautomer hash, and the S-group data.
_SCHEME = RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS


# The standard counterion list: chloride, bromide, sodium, potassium, TFA and the rest
# of what a reagent is sold as. Not the largest fragment, which is the other way to
# spell "parent" and answers palladium acetate with acetic acid.
_SALTS = SaltRemover.SaltRemover()


def _molecule(structure: Chem.Mol | str) -> Chem.Mol | None:
    """Returns a molecule from either a molecule or SMILES, or None if it does not read.

    Args:
        structure: A parsed molecule, or SMILES to read one from.

    Returns:
        The molecule, or None where RDKit cannot parse the string.
    """
    return Chem.MolFromSmiles(structure) if isinstance(structure, str) else structure


def salt_parent(molecule: Chem.Mol) -> Chem.Mol:
    """Returns the molecule without the counterions it was sold as.

    Sodium acetate becomes acetate and triethylamine hydrochloride becomes
    triethylamine, so a reagent recorded with a spectator counterion and the same
    reagent recorded without one come out the same.

    Stripping stops where the molecule *is* the salt: nothing is removed unless what
    survives still holds carbon, which leaves sodium hydride whole rather than turning
    it into hydrogen, and palladium acetate whole rather than into palladium. That rule
    is why this is not RDKit's fragment parent, which keeps the largest fragment and
    would answer palladium acetate with acetic acid.

    Args:
        molecule: A parsed molecule.

    Returns:
        The stripped molecule, or the original where stripping would leave no carbon.
    """
    stripped = _SALTS.StripMol(molecule, dontRemoveEverything=True)
    if stripped.GetNumAtoms() and any(
        atom.GetAtomicNum() == 6 for atom in stripped.GetAtoms()
    ):
        return stripped
    return molecule


def parent_hash(structure: Chem.Mol | str) -> str | None:
    """Returns the hash of a molecule's salt parent.

    What ``mol_hash`` ignores, this ignores too, and counterions besides. Two reagents
    whose parents differ still differ: potassium and sodium carbonate share a parent and
    hash alike, while sodium carbonate and sodium acetate do not.

    Args:
        structure: A parsed molecule, or SMILES to read one from.

    Returns:
        The hash, or None where RDKit cannot read or hash the structure.
    """
    molecule = _molecule(structure)
    return None if molecule is None else mol_hash(salt_parent(molecule))


def mol_hash(structure: Chem.Mol | str) -> str | None:
    """Returns an identity for a molecule that ignores tautomer and protonation state.

    Two drawings of one reagent hash the same: acetic acid and acetate, an amine and
    its ammonium, a 2-pyridone and its 2-hydroxypyridine tautomer, a compound written
    with atom-map labels and the same compound without them. Fragments are left alone,
    so sodium acetate hashes differently from acetic acid, and so is stereochemistry,
    so enantiomers stay apart.

    The molecule is uncharged before it is hashed, which is what makes the scheme
    insensitive to protonation as well: one of its layers is the molecular formula, and
    a formula states the charge.

    Args:
        structure: A parsed molecule, or SMILES to read one from.

    Returns:
        The hash, or None where RDKit cannot read or hash the structure, which leaves
        the column null rather than storing a value the query side would not reproduce.
    """
    molecule = _molecule(structure)
    if molecule is None:
        return None
    try:
        layers = RegistrationHash.GetMolLayers(
            _UNCHARGER.uncharge(molecule), enable_tautomer_hash_v2=True
        )
        return RegistrationHash.GetMolHash(layers, _SCHEME)
    except (
        Chem.AtomValenceException,
        Chem.KekulizeException,
        RuntimeError,
        ValueError,
    ):
        return None


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
        "mol_hash": None,
        "parent_hash": None,
        "pattern_fp": None,
        "morgan_fp": None,
        "morgan_popcount": None,
        "mol_binary": None,
    }
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return row
    row["mol_hash"] = mol_hash(mol)
    row["parent_hash"] = parent_hash(mol)
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

    Delegates to ``base.is_current``, which requires the artifact name, the source
    content, the library version, the artifact version, and the RDKit version to all
    match. A missing or unreadable file is not current.

    Args:
        path: Path to check. Need not exist.
        source_md5: Hash of the source dataset the artifact should restate, which is
            the dataset's own even though the artifact derives from a projection.

    Returns:
        Whether the file is a structures artifact this library would write today.
    """
    return base.is_current(path, ARTIFACT, source_md5)


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
        compression: Parquet codec, any name ``pq.ParquetWriter`` accepts.
        row_group_size: Rows per output row group, which is also how many structures
            are featurized and held in memory at a time.

    Returns:
        Number of rows written: the number of distinct structures in the projection.

    Raises:
        ValueError: If ``source`` is not a projection, or does not state a single dense
            (smiles, structure_id) mapping.
    """
    parent = base.load_stamps(source)
    if parent.artifact != projection.ARTIFACT:
        raise ValueError(
            f"{source} is a {parent.artifact}, not a {projection.ARTIFACT}; the "
            "structures artifact featurizes a projection"
        )
    # derive_tree refuses stale parents, but this writer is public and its output
    # inherits the dataset hash: an artifact derived from a stale projection would
    # claim a provenance it does not have and nothing would ever mark it stale again.
    if not base.stamps_are_current(parent, projection.ARTIFACT):
        raise ValueError(
            f"{source} is a stale {projection.ARTIFACT}; derive it again first"
        )
    if source_md5 is None:
        source_md5 = parent.source_md5
    if source_dataset_id is None:
        source_dataset_id = parent.source_dataset_id
    stamps = base.current_stamps(ARTIFACT, source_dataset_id, source_md5)
    schema = SCHEMA.with_metadata(base.to_metadata(stamps))
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
