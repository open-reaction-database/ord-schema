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

"""Tests for ord_schema.structures."""

import pathlib
from importlib import metadata

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem, DataStructs

from ord_schema import artifacts, parquet, projection, structures
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
    # Benzene appears in every reaction; the artifact holds it once, under the id the
    # projection stamped everywhere it occurs.
    reactions = [_reaction(f"ord-{i:04d}") for i in range(5)]
    source = _project(tmp_path, reactions)
    output = tmp_path / "structures.parquet"
    assert structures.write_structures(source, output) == 2


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
    for column in ("pattern_fp", "morgan_fp", "morgan_popcount", "mol_binary"):
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
    stamps = artifacts.load_stamps(output)
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
    stamps = artifacts.current_stamps(projection.ARTIFACT, "ord_dataset-1", "0" * 32)
    schema = projection.SCHEMA.with_metadata(artifacts.to_metadata(stamps))
    path = tmp_path / "projection.parquet"
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path)
    return path


def test_a_smiles_without_an_id_is_refused(tmp_path):
    # message_row without a mapping leaves ids null; deriving from such a file would
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
    # The reverse defect of conflicting ids: the mapping must be one-to-one in both
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
