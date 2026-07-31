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

"""Tests for ord_schema.projection."""

from typing import cast

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from google.protobuf.descriptor import Descriptor

from ord_schema import artifacts, parquet, projection
from ord_schema.proto import dataset_pb2, reaction_pb2


def _reaction(reaction_id: str = "ord-0001") -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction()
    reaction.reaction_id = reaction_id
    component = reaction.inputs["toluene"].components.add()
    component.identifiers.add(type="SMILES", value="Cc1ccccc1")
    component.identifiers.add(type="NAME", value="toluene")
    component.amount.mass.value = 1.0
    component.amount.mass.units = reaction_pb2.Mass.GRAM
    outcome = reaction.outcomes.add()
    product = outcome.products.add()
    product.identifiers.add(type="SMILES", value="Cc1ccccc1O")
    return reaction


def _dataset_path(tmp_path, reactions, dataset_id="ord_dataset-1"):
    dataset = dataset_pb2.Dataset(
        dataset_id=dataset_id,
        name="test",
        description="test dataset",
        reactions=reactions,
    )
    path = tmp_path / "source.parquet"
    parquet.save_dataset(dataset, path)
    return path


def _leaf_count(data_type: pa.DataType) -> int:
    if pa.types.is_struct(data_type):
        return sum(_leaf_count(field.type) for field in data_type)
    if pa.types.is_list(data_type):
        return _leaf_count(data_type.value_type)
    if pa.types.is_map(data_type):
        return _leaf_count(data_type.item_type) + 1
    return 1


# Schema generation


def test_schema_covers_every_reaction_field():
    descriptor = cast(Descriptor, reaction_pb2.Reaction.DESCRIPTOR)
    expected = {field.name for field in descriptor.fields}
    assert expected <= set(projection.SCHEMA.names)


def test_schema_adds_a_collapsed_smiles_at_reaction_and_component_level():
    assert "smiles" in projection.SCHEMA.names
    components = (
        projection.SCHEMA.field("inputs").type.item_type.field("components").type
    )
    assert "smiles" in [field.name for field in components.value_type]


def test_united_fields_collapse_to_a_float_named_for_its_unit():
    temperature = projection.SCHEMA.field("conditions").type.field("temperature").type
    assert temperature.field("setpoint_kelvin").type == pa.float64()
    assert "setpoint" not in [field.name for field in temperature]


def test_enums_project_as_strings():
    stirring = projection.SCHEMA.field("conditions").type.field("stirring").type
    assert stirring.field("type").type == pa.string()


def test_repeated_messages_project_as_lists_and_maps_as_maps():
    assert pa.types.is_list(projection.SCHEMA.field("outcomes").type)
    assert pa.types.is_map(projection.SCHEMA.field("inputs").type)


def test_schema_is_deep_but_finite():
    # The projection only terminates because the message graph is acyclic; if that ever
    # stops holding, build_schema raises rather than recursing forever.
    assert sum(_leaf_count(field.type) for field in projection.SCHEMA) > 300


def test_build_schema_rejects_a_recursive_message():
    descriptor = cast(Descriptor, reaction_pb2.Reaction.DESCRIPTOR)
    with pytest.raises(ValueError, match="reaches itself"):
        # Seed the stack with Reaction so the first visit looks like a revisit.
        projection._struct_fields(descriptor, frozenset({descriptor.full_name}))


# Row projection


def test_unset_scalars_are_null_rather_than_the_proto_default():
    row = projection.reaction_row(reaction_pb2.Reaction())
    assert row["reaction_id"] is None
    assert row["provenance"] is None


def test_zero_is_preserved_where_the_source_says_zero():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = 0.0
    setpoint.units = reaction_pb2.Temperature.KELVIN
    row = projection.reaction_row(reaction)
    assert row["conditions"]["temperature"]["setpoint_kelvin"] == 0.0


def test_enum_values_project_as_their_names():
    reaction = _reaction()
    reaction.conditions.stirring.type = reaction_pb2.StirringConditions.STIR_BAR
    row = projection.reaction_row(reaction)
    assert row["conditions"]["stirring"]["type"] == "STIR_BAR"


def test_temperature_converts_to_kelvin():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = 25.0
    setpoint.units = reaction_pb2.Temperature.CELSIUS
    row = projection.reaction_row(reaction)
    assert row["conditions"]["temperature"]["setpoint_kelvin"] == pytest.approx(298.15)


def test_amount_mass_converts_to_grams():
    reaction = _reaction()
    mass = reaction.inputs["toluene"].components[0].amount.mass
    mass.value, mass.units = 1000.0, reaction_pb2.Mass.MILLIGRAM
    row = projection.reaction_row(reaction)
    component = row["inputs"][0][1]["components"][0]
    assert component["amount"]["mass_grams"] == pytest.approx(1.0)


def test_a_unit_the_resolver_cannot_read_is_null_not_a_wrong_number():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = 25.0  # units left UNSPECIFIED
    row = projection.reaction_row(reaction)
    assert row["conditions"]["temperature"]["setpoint_kelvin"] is None


def test_structural_identifiers_collapse_into_smiles():
    row = projection.reaction_row(_reaction())
    component = row["inputs"][0][1]["components"][0]
    assert component["smiles"] == "Cc1ccccc1"
    assert [identifier["type"] for identifier in component["identifiers"]] == ["NAME"]


def test_non_structural_identifiers_survive_the_collapse():
    reaction = _reaction()
    component = reaction.inputs["toluene"].components[0]
    component.identifiers.add(type="CAS_NUMBER", value="108-88-3")
    row = projection.reaction_row(reaction)
    kept = {
        identifier["type"]
        for identifier in row["inputs"][0][1]["components"][0]["identifiers"]
    }
    assert kept == {"NAME", "CAS_NUMBER"}


def test_a_name_only_compound_keeps_its_name_and_has_no_smiles():
    reaction = reaction_pb2.Reaction(reaction_id="ord-0002")
    component = reaction.inputs["quench"].components.add()
    component.identifiers.add(type="NAME", value="ice water")
    row = projection.reaction_row(reaction)
    projected = row["inputs"][0][1]["components"][0]
    assert projected["smiles"] is None
    assert projected["identifiers"][0]["value"] == "ice water"


def test_reaction_smiles_is_generated_when_the_source_records_none():
    row = projection.reaction_row(_reaction())
    assert row["smiles"] is not None
    assert ">>" in row["smiles"]


def test_empty_repeated_fields_are_null_rather_than_empty_lists():
    row = projection.reaction_row(reaction_pb2.Reaction(reaction_id="ord-0003"))
    assert row["outcomes"] is None
    assert row["observations"] is None


# Writing


def test_write_projection_round_trips(tmp_path):
    source = _dataset_path(tmp_path, [_reaction("ord-1"), _reaction("ord-2")])
    output = tmp_path / "projection.parquet"
    assert projection.write_projection(source, output) == 2
    table = pq.read_table(output)
    assert table.num_rows == 2
    assert table.column("reaction_id").to_pylist() == ["ord-1", "ord-2"]


def test_write_projection_stamps_the_footer(tmp_path):
    source = _dataset_path(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    stamps = artifacts.load_stamps(output)
    assert stamps.artifact == projection.ARTIFACT
    assert stamps.artifact_version == artifacts.ARTIFACT_VERSION
    assert stamps.source_dataset_id == "ord_dataset-1"


def test_is_current_tracks_source_content(tmp_path):
    source = _dataset_path(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    source_md5, _ = parquet.streaming_md5(source)
    assert projection.is_current(output, source_md5)
    assert not projection.is_current(output, "0" * 32)


def test_is_current_rejects_a_different_artifact(tmp_path):
    source = _dataset_path(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    source_md5, _ = parquet.streaming_md5(source)
    assert not artifacts.is_current(output, "view", source_md5)


def test_is_current_is_false_for_a_missing_file(tmp_path):
    assert not projection.is_current(tmp_path / "absent.parquet", "0" * 32)


def test_is_current_is_false_for_a_source_dataset(tmp_path):
    source = _dataset_path(tmp_path, [_reaction()])
    source_md5, _ = parquet.streaming_md5(source)
    assert not projection.is_current(source, source_md5)


def test_a_failed_write_leaves_an_existing_projection_intact(tmp_path, monkeypatch):
    source = _dataset_path(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    original = output.read_bytes()

    def _boom(*args, **kwargs):
        raise RuntimeError("projection failed")

    monkeypatch.setattr(projection, "reaction_row", _boom)
    with pytest.raises(RuntimeError, match="projection failed"):
        projection.write_projection(source, output)
    assert output.read_bytes() == original
    assert list(tmp_path.glob("*.tmp")) == []
