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

"""Tests for ord_schema.artifacts.structures."""

import pathlib
from importlib import metadata

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem, DataStructs

from ord_schema import parquet
from ord_schema.artifacts import base, projection, structures
from ord_schema.proto import reaction_pb2


def _reaction(
    reaction_id: str = "ord-0001", *, input_smiles: str = "c1ccccc1"
) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    component = reaction.inputs["a"].components.add()
    component.identifiers.add(type="SMILES", value=input_smiles)
    outcome = reaction.outcomes.add()
    outcome.products.add().identifiers.add(type="SMILES", value="Cc1ccccc1")
    return reaction


def _project(tmp_path, reactions=None) -> pathlib.Path:
    """Writes a source dataset and its projection, which is what structures reads."""
    source = tmp_path / "source.parquet"
    with parquet.DatasetWriter(
        source, name="test", description="test", dataset_id="ord_dataset-1"
    ) as writer:
        writer.write_all(reactions or [_reaction()])
    path = tmp_path / "projection.parquet"
    projection.write_projection(source, path)
    return path


def test_structure_columns_cover_every_compound_location():
    columns = structures._structure_columns()
    # Compounds sit in four places; the two beyond inputs and products are easy to
    # overlook, which is why the columns are walked from the schema rather than listed.
    for prefix in (
        "inputs.key_value.value.components",
        "workups.list.element.input.components",
        "outcomes.list.element.products",
        "outcomes.list.element.products.list.element.measurements.list.element"
        ".authentic_standard",
    ):
        assert f"{prefix}.list.element.smiles" in columns or (
            f"{prefix}.smiles" in columns
        ), prefix


def test_write_structures_round_trips(tmp_path):
    source = _project(tmp_path)
    output = tmp_path / "structures.parquet"
    assert structures.write_structures(source, output) == 2
    table = pq.read_table(output)
    assert table.schema.names == structures.SCHEMA.names
    assert table.column("structure_id").to_pylist() == [0, 1]
    assert table.column("smiles").to_pylist() == ["c1ccccc1", "Cc1ccccc1"]


def test_shared_structures_are_deduplicated(tmp_path):
    # Benzene appears in every reaction; the artifact holds it once, under the ID the
    # projection stamped everywhere it occurs.
    reactions = [_reaction(f"ord-{i:04d}") for i in range(5)]
    source = _project(tmp_path, reactions)
    output = tmp_path / "structures.parquet"
    assert structures.write_structures(source, output) == 2
    table = pq.read_table(output)
    assert table.column("structure_id").to_pylist() == [0, 1]
    assert table.column("smiles").to_pylist() == ["c1ccccc1", "Cc1ccccc1"]


def test_row_groups_keep_ids_aligned(tmp_path):
    # Position in the file is the ID, so each output row group must continue the
    # sequence and each row's derived columns must be its own molecule's --
    # misalignment here is silent join corruption, not an error.
    reaction = _reaction()
    component = reaction.inputs["b"].components.add()
    component.identifiers.add(type="SMILES", value="CCO")
    source = _project(tmp_path, [reaction])
    output = tmp_path / "structures.parquet"
    assert structures.write_structures(source, output, row_group_size=1) == 3
    with pq.ParquetFile(output) as artifact:
        assert artifact.num_row_groups == 3
    table = pq.read_table(output)
    assert table.column("structure_id").to_pylist() == [0, 1, 2]
    for row in table.to_pylist():
        mol = Chem.Mol(row["mol_binary"])
        assert Chem.MolToSmiles(mol) == Chem.CanonSmiles(row["smiles"]), row["smiles"]
        fingerprint = DataStructs.CreateFromBinaryText(row["morgan_fp"])
        assert row["morgan_popcount"] == fingerprint.GetNumOnBits()


def test_the_screen_is_complete_and_verification_is_exact(tmp_path):
    # Toluene contains benzene: every query bit must be set in the target fingerprint
    # (the screen never rejects a true match), and the serialized molecule must verify
    # identically to one re-parsed from SMILES.
    source = _project(tmp_path)
    output = tmp_path / "structures.parquet"
    structures.write_structures(source, output)
    rows = pq.read_table(output).to_pylist()
    target = next(row for row in rows if row["smiles"] == "Cc1ccccc1")
    query = Chem.MolFromSmarts("c1ccccc1")
    query_fp = DataStructs.CreateFromBinaryText(
        DataStructs.BitVectToBinaryText(
            Chem.PatternFingerprint(query, fpSize=structures.PATTERN_FP_SIZE)
        )
    )
    target_fp = DataStructs.CreateFromBinaryText(target["pattern_fp"])
    assert (query_fp & target_fp).GetNumOnBits() == query_fp.GetNumOnBits()
    mol = Chem.Mol(target["mol_binary"])
    assert mol.HasSubstructMatch(query)
    assert Chem.MolToSmiles(mol) == Chem.CanonSmiles("Cc1ccccc1")
    # Containment is directional: benzene's row must not pass a toluene query, or the
    # screen is comparing something other than the fingerprints.
    benzene_fp = DataStructs.CreateFromBinaryText(
        next(row for row in rows if row["smiles"] == "c1ccccc1")["pattern_fp"]
    )
    assert (target_fp & benzene_fp).GetNumOnBits() < target_fp.GetNumOnBits()
    # A query with a generic atom exercises RDKit's query-fingerprint generalization,
    # the part of PatternFingerprint most likely to shift under an upgrade.
    generic = Chem.MolFromSmarts("*c1ccccc1")
    generic_fp = DataStructs.CreateFromBinaryText(
        DataStructs.BitVectToBinaryText(
            Chem.PatternFingerprint(generic, fpSize=structures.PATTERN_FP_SIZE)
        )
    )
    assert (generic_fp & target_fp).GetNumOnBits() == generic_fp.GetNumOnBits()


def test_morgan_popcount_counts_the_morgan_fingerprint(tmp_path):
    source = _project(tmp_path)
    output = tmp_path / "structures.parquet"
    structures.write_structures(source, output)
    for row in pq.read_table(output).to_pylist():
        fingerprint = DataStructs.CreateFromBinaryText(row["morgan_fp"])
        assert row["morgan_popcount"] == fingerprint.GetNumOnBits()
        assert row["morgan_popcount"] > 0


def test_an_unparseable_structure_keeps_its_row_with_null_columns():
    row = structures._row(7, "not-a-smiles")
    assert row["structure_id"] == 7
    assert row["smiles"] == "not-a-smiles"
    for column in (
        "mol_hash",
        "pattern_fp",
        "morgan_fp",
        "morgan_popcount",
        "mol_binary",
    ):
        assert row[column] is None, column


def test_structures_from_workups_and_standards_are_included(tmp_path):
    reaction = _reaction()
    workup = reaction.workups.add()
    workup.input.components.add().identifiers.add(type="SMILES", value="CCO")
    measurement = reaction.outcomes[0].products[0].measurements.add()
    measurement.authentic_standard.identifiers.add(type="SMILES", value="CCN")
    source = _project(tmp_path, [reaction])
    output = tmp_path / "structures.parquet"
    assert structures.write_structures(source, output) == 4
    smiles = pq.read_table(output).column("smiles").to_pylist()
    assert "CCO" in smiles
    assert "CCN" in smiles


def test_a_dataset_with_no_structures_writes_an_empty_artifact(tmp_path):
    # The empty file is the point: it carries the schema and stamps, so the derive
    # driver can tell "done, nothing to featurize" from "never derived".
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    component = reaction.inputs["a"].components.add()
    component.identifiers.add(type="NAME", value="mystery")
    reaction.outcomes.add()
    source = _project(tmp_path, [reaction])
    output = tmp_path / "structures.parquet"
    assert structures.write_structures(source, output) == 0
    table = pq.read_table(output)
    assert table.num_rows == 0
    assert table.schema.names == structures.SCHEMA.names
    stamps = base.load_stamps(output)
    assert stamps.artifact == structures.ARTIFACT
    assert structures.is_current(output, stamps.source_md5)


def test_stamps_name_the_source_dataset_not_the_projection(tmp_path):
    source_path = tmp_path / "source.parquet"
    with parquet.DatasetWriter(
        source_path, name="test", description="test", dataset_id="ord_dataset-1"
    ) as writer:
        writer.write_all([_reaction()])
    projected = tmp_path / "projection.parquet"
    projection.write_projection(source_path, projected)
    output = tmp_path / "structures.parquet"
    structures.write_structures(projected, output)
    stamps = base.load_stamps(output)
    assert stamps.artifact == structures.ARTIFACT
    assert stamps.source_md5 == parquet.DatasetView(source_path).md5()
    assert stamps.source_dataset_id == "ord_dataset-1"
    assert stamps.ord_schema_version == metadata.version("ord-schema")
    assert structures.is_current(output, stamps.source_md5)


def test_a_source_dataset_is_refused(tmp_path):
    source = tmp_path / "source.parquet"
    with parquet.DatasetWriter(
        source, name="test", description="test", dataset_id="ord_dataset-1"
    ) as writer:
        writer.write_all([_reaction()])
    with pytest.raises(ValueError, match="not a derived artifact"):
        structures.write_structures(source, tmp_path / "structures.parquet")


def _corrupt_projection(tmp_path, rows) -> pathlib.Path:
    """Writes ``rows`` under a projection's schema and stamps, bypassing the writer."""
    stamps = base.current_stamps(projection.ARTIFACT, "ord_dataset-1", "0" * 32)
    schema = projection.SCHEMA.with_metadata(base.to_metadata(stamps))
    path = tmp_path / "projection.parquet"
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path)
    return path


def test_a_smiles_without_an_id_is_refused(tmp_path):
    # message_row without a mapping leaves IDs null; deriving from such a file would
    # join wrongly rather than fail, so the mismatch is an error here.
    row = projection.message_row(_reaction())
    row["reaction_id"] = "ord-0001"
    path = _corrupt_projection(tmp_path, [row])
    with pytest.raises(ValueError, match="both or neither"):
        structures.write_structures(path, tmp_path / "structures.parquet")


def test_conflicting_ids_are_refused(tmp_path):
    ids_a = {"c1ccccc1": 0, "Cc1ccccc1": 1}
    ids_b = {"c1ccccc1": 1, "Cc1ccccc1": 0}
    row_a = projection.message_row(_reaction("ord-0001"), ids_a)
    row_b = projection.message_row(_reaction("ord-0002"), ids_b)
    row_a["reaction_id"], row_b["reaction_id"] = "ord-0001", "ord-0002"
    path = _corrupt_projection(tmp_path, [row_a, row_b])
    with pytest.raises(ValueError, match="is both"):
        structures.write_structures(path, tmp_path / "structures.parquet")


def test_one_smiles_under_two_ids_is_refused(tmp_path):
    # The reverse defect of conflicting IDs: the mapping must be one-to-one in both
    # directions, or the artifact holds duplicate rows and an inflated distinct count.
    ids_a = {"c1ccccc1": 0, "Cc1ccccc1": 1}
    ids_b = {"c1ccccc1": 2, "Cc1ccccc1": 1}
    row_a = projection.message_row(_reaction("ord-0001"), ids_a)
    row_b = projection.message_row(_reaction("ord-0002"), ids_b)
    row_a["reaction_id"], row_b["reaction_id"] = "ord-0001", "ord-0002"
    path = _corrupt_projection(tmp_path, [row_a, row_b])
    with pytest.raises(ValueError, match="one structure twice"):
        structures.write_structures(path, tmp_path / "structures.parquet")


def test_a_gap_in_the_id_space_is_refused(tmp_path):
    row = projection.message_row(_reaction(), {"c1ccccc1": 3, "Cc1ccccc1": 4})
    row["reaction_id"] = "ord-0001"
    path = _corrupt_projection(tmp_path, [row])
    with pytest.raises(ValueError, match="not dense"):
        structures.write_structures(path, tmp_path / "structures.parquet")


def test_an_unparseable_structure_survives_the_write_path(tmp_path):
    # Cannot arise from a same-version projection -- its SMILES are RDKit-canonical --
    # so it is staged through _corrupt_projection: the row must keep its position (the
    # ID space stays dense) rather than being cleaned out of the batch.
    ids: dict[str, int] = {}
    row = projection.message_row(_reaction(), ids)
    (component,) = dict(row["inputs"])["a"]["components"]
    component["smiles"] = "not-a-smiles"  # Keeps the ID benzene was assigned.
    row["reaction_id"] = "ord-0001"
    path = _corrupt_projection(tmp_path, [row])
    output = tmp_path / "structures.parquet"
    assert structures.write_structures(path, output) == 2
    rows = pq.read_table(output).to_pylist()
    bad = next(r for r in rows if r["smiles"] == "not-a-smiles")
    assert bad["structure_id"] == ids["c1ccccc1"]
    for column in (
        "mol_hash",
        "pattern_fp",
        "morgan_fp",
        "morgan_popcount",
        "mol_binary",
    ):
        assert bad[column] is None, column
    good = next(r for r in rows if r["smiles"] != "not-a-smiles")
    assert good["mol_binary"] is not None


def test_a_projection_missing_the_id_columns_is_refused(tmp_path):
    # read_row_group silently drops requested columns a file does not have, so an
    # old-schema projection would otherwise state no pairs at all and derive an empty
    # artifact that stamps itself current.
    stamps = base.current_stamps(projection.ARTIFACT, "ord_dataset-1", "0" * 32)
    component = pa.struct([pa.field("smiles", pa.string())])
    reaction_input = pa.struct([pa.field("components", pa.list_(component))])
    schema = pa.schema(
        [
            pa.field("reaction_id", pa.string()),
            pa.field("inputs", pa.map_(pa.string(), reaction_input)),
        ]
    ).with_metadata(base.to_metadata(stamps))
    path = tmp_path / "projection.parquet"
    pq.write_table(
        pa.Table.from_pylist(
            [{"reaction_id": "ord-0001", "inputs": []}], schema=schema
        ),
        path,
    )
    with pytest.raises(ValueError, match="derive the projection again"):
        structures.write_structures(path, tmp_path / "structures.parquet")


def test_a_stale_projection_is_refused(tmp_path):
    # The output inherits the dataset hash, so an artifact derived from a stale
    # projection would claim a provenance it does not have and never read stale again.
    stamps = base.Stamps(
        artifact=projection.ARTIFACT,
        source_dataset_id="ord_dataset-1",
        source_md5="0" * 32,
        ord_schema_version="0.0.0",
        artifact_version=base.ARTIFACT_VERSION,
        rdkit_version="0000.00.0",
    )
    schema = projection.SCHEMA.with_metadata(base.to_metadata(stamps))
    path = tmp_path / "projection.parquet"
    pq.write_table(pa.Table.from_pylist([], schema=schema), path)
    with pytest.raises(ValueError, match="stale"):
        structures.write_structures(path, tmp_path / "structures.parquet")


@pytest.mark.parametrize(
    ("label", "left", "right"),
    [
        ("protonation", "CC(=O)O", "CC(=O)[O-]"),
        ("ammonium", "CCN", "CC[NH3+]"),
        ("keto-enol", "CC(=O)C", "CC(O)=C"),
        ("amide-imidic acid", "CC(=O)N", "CC(O)=N"),
        ("lactam-lactim", "O=c1cccc[nH]1", "Oc1ccccn1"),
        ("kekule", "C1=CC=NC=C1", "c1ccncc1"),
    ],
)
def test_one_compound_drawn_two_ways_hashes_the_same(label, left, right):
    del label
    assert structures.mol_hash(left) == structures.mol_hash(right)


@pytest.mark.parametrize(
    ("label", "left", "right"),
    [
        ("different molecules", "CCO", "CC(=O)O"),
        ("different heteroatom", "c1ccncc1", "c1ccccc1"),
        # A salt is a different reagent rather than a different drawing of one, and
        # every rule that collapses the two mangles something else: keeping the largest
        # fragment turns Pd(OAc)2 into acetic acid, and stripping recognized
        # counterions turns NaH into hydrogen.
        ("salt", "CC(=O)O", "CC(=O)[O-].[Na+]"),
        ("cocrystal", "c1ccccc1", "c1ccccc1.Cc1ccccc1"),
    ],
)
def test_different_compounds_hash_differently(label, left, right):
    del label
    assert structures.mol_hash(left) != structures.mol_hash(right)


def test_an_unhashable_smiles_is_none_rather_than_a_guess():
    # A column holding a spelling the query side would not reproduce is worse than a
    # null one, which simply never matches.
    assert structures.mol_hash("not a molecule") is None


def test_the_artifact_carries_a_hash_beside_every_readable_structure(tmp_path):
    source = _project(tmp_path, [_reaction(input_smiles="CC(=O)[O-]")])
    output = tmp_path / "structures.parquet"
    structures.write_structures(source, output)
    rows = pq.read_table(output).to_pylist()
    # The acetate is stored as drawn and hashes to what acetic acid hashes to, which is
    # the whole point: the corpus keeps its own spelling and still answers for it.
    acetate = next(row for row in rows if "[O-]" in row["smiles"])
    assert acetate["mol_hash"] == structures.mol_hash("CC(=O)O")
    assert all(row["mol_hash"] for row in rows)
