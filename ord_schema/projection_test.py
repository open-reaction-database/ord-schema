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

import pathlib
from types import SimpleNamespace
from typing import cast

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from google.protobuf.descriptor import Descriptor, FieldDescriptor
from rdkit import Chem

from ord_schema import artifacts, parquet, projection
from ord_schema.proto import reaction_pb2


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


def _source(
    tmp_path, reactions, dataset_id="ord_dataset-1", row_group_size=1000
) -> pathlib.Path:
    path = tmp_path / "source.parquet"
    with parquet.DatasetWriter(
        path,
        name="test",
        description="test dataset",
        dataset_id=dataset_id,
        row_group_size=row_group_size,
    ) as writer:
        writer.write_all(reactions)
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


def test_united_fields_collapse_to_floats_named_for_their_unit():
    temperature = projection.SCHEMA.field("conditions").type.field("temperature").type
    assert temperature.field("setpoint_kelvin").type == pa.float64()
    assert temperature.field("setpoint_precision_kelvin").type == pa.float64()
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
    row = projection.message_row(reaction_pb2.Reaction())
    assert row["reaction_id"] is None
    assert row["provenance"] is None


def test_zero_is_preserved_where_the_source_says_zero():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = 0.0
    setpoint.units = reaction_pb2.Temperature.KELVIN
    row = projection.message_row(reaction)
    assert row["conditions"]["temperature"]["setpoint_kelvin"] == 0.0


def test_enum_values_project_as_their_names():
    reaction = _reaction()
    reaction.conditions.stirring.type = reaction_pb2.StirringConditions.STIR_BAR
    row = projection.message_row(reaction)
    assert row["conditions"]["stirring"]["type"] == "STIR_BAR"


def test_temperature_converts_to_kelvin():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = 25.0
    setpoint.units = reaction_pb2.Temperature.CELSIUS
    row = projection.message_row(reaction)
    assert row["conditions"]["temperature"]["setpoint_kelvin"] == pytest.approx(298.15)


def test_precision_converts_with_its_value():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value, setpoint.precision = 77.0, 1.8
    setpoint.units = reaction_pb2.Temperature.FAHRENHEIT
    temperature = projection.message_row(reaction)["conditions"]["temperature"]
    assert temperature["setpoint_kelvin"] == pytest.approx(298.15)
    # Fahrenheit is an offset scale, so the precision converts by the ratio alone.
    assert temperature["setpoint_precision_kelvin"] == pytest.approx(1.0)


def test_unrecorded_precision_is_null_not_zero():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = 25.0
    setpoint.units = reaction_pb2.Temperature.CELSIUS
    temperature = projection.message_row(reaction)["conditions"]["temperature"]
    assert temperature["setpoint_kelvin"] == pytest.approx(298.15)
    assert temperature["setpoint_precision_kelvin"] is None


def test_precision_is_null_where_its_value_cannot_be_converted():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value, setpoint.precision = 25.0, 0.5  # units left UNSPECIFIED
    temperature = projection.message_row(reaction)["conditions"]["temperature"]
    assert temperature["setpoint_kelvin"] is None
    assert temperature["setpoint_precision_kelvin"] is None


def test_amount_mass_converts_to_grams():
    reaction = _reaction()
    mass = reaction.inputs["toluene"].components[0].amount.mass
    mass.value, mass.units = 1000.0, reaction_pb2.Mass.MILLIGRAM
    row = projection.message_row(reaction)
    component = row["inputs"][0][1]["components"][0]
    assert component["amount"]["mass_grams"] == pytest.approx(1.0)


def test_a_unit_the_resolver_cannot_read_is_null_not_a_wrong_number():
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = 25.0  # units left UNSPECIFIED
    row = projection.message_row(reaction)
    assert row["conditions"]["temperature"]["setpoint_kelvin"] is None


def test_structural_identifiers_collapse_into_smiles():
    row = projection.message_row(_reaction())
    component = row["inputs"][0][1]["components"][0]
    assert component["smiles"] == "Cc1ccccc1"
    assert [identifier["type"] for identifier in component["identifiers"]] == ["NAME"]


def test_non_structural_identifiers_survive_the_collapse():
    reaction = _reaction()
    component = reaction.inputs["toluene"].components[0]
    component.identifiers.add(type="CAS_NUMBER", value="108-88-3")
    row = projection.message_row(reaction)
    kept = {
        identifier["type"]
        for identifier in row["inputs"][0][1]["components"][0]["identifiers"]
    }
    assert kept == {"NAME", "CAS_NUMBER"}


_TOLUENE_MOLBLOCK = Chem.MolToMolBlock(Chem.MolFromSmiles("Cc1ccccc1"))


@pytest.mark.parametrize(
    ("identifier_type", "value"),
    [
        ("INCHI", "InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3"),
        ("MOLBLOCK", _TOLUENE_MOLBLOCK),
        ("CXSMILES", "Cc1ccccc1"),
    ],
)
def test_every_structural_compound_type_collapses(identifier_type, value):
    reaction = reaction_pb2.Reaction(reaction_id="ord-0005")
    component = reaction.inputs["x"].components.add()
    component.identifiers.add(type=identifier_type, value=value)
    component.identifiers.add(type="NAME", value="toluene")
    projected = projection.message_row(reaction)["inputs"][0][1]["components"][0]
    assert projected["smiles"] == "Cc1ccccc1"
    assert [identifier["type"] for identifier in projected["identifiers"]] == ["NAME"]


def test_reaction_structural_identifiers_collapse_too():
    # smiles is generated rather than read from the recorded identifier, so the two can
    # differ; the source column stays authoritative, as it does for a dropped MOLBLOCK.
    reaction = _reaction()
    reaction.identifiers.add(type="REACTION_CXSMILES", value="C>>CO |f:0.1|")
    reaction.identifiers.add(type="REACTION_TYPE", value="oxidation")
    row = projection.message_row(reaction)
    assert row["smiles"] is not None
    assert [identifier["type"] for identifier in row["identifiers"]] == [
        "REACTION_TYPE"
    ]


def test_a_product_compound_collapses_like_an_input_component():
    row = projection.message_row(_reaction())
    assert row["outcomes"][0]["products"][0]["smiles"] == "Cc1ccccc1O"


def test_a_malformed_identifier_does_not_hide_a_later_valid_one():
    # smiles_from_compound stops at the first structural identifier, so without the
    # fallback this compound projects no structure though it records one.
    reaction = reaction_pb2.Reaction(reaction_id="ord-0007")
    component = reaction.inputs["x"].components.add()
    component.identifiers.add(type="SMILES", value="not a smiles")
    component.identifiers.add(
        type="INCHI", value="InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3"
    )
    projected = projection.message_row(reaction)["inputs"][0][1]["components"][0]
    assert projected["smiles"] == "Cc1ccccc1"


def test_repeated_scalars_project_as_a_list_and_empty_as_null():
    # A different branch from repeated messages, and the only one carrying bare scalars.
    details = reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails()
    assert projection.message_row(details)["eic_masses"] is None
    details.eic_masses.extend([1.5, 2.5])
    assert projection.message_row(details)["eic_masses"] == [1.5, 2.5]


def test_a_structural_identifier_that_fails_to_collapse_is_kept():
    # Dropping it as well would report the compound as having no recorded structure.
    reaction = reaction_pb2.Reaction(reaction_id="ord-0004")
    component = reaction.inputs["mystery"].components.add()
    component.identifiers.add(type="SMILES", value="not a smiles")
    projected = projection.message_row(reaction)["inputs"][0][1]["components"][0]
    assert projected["smiles"] is None
    assert projected["identifiers"] == [
        {"type": "SMILES", "value": "not a smiles", "details": None}
    ]


def test_a_name_only_compound_keeps_its_name_and_has_no_smiles():
    reaction = reaction_pb2.Reaction(reaction_id="ord-0002")
    component = reaction.inputs["quench"].components.add()
    component.identifiers.add(type="NAME", value="ice water")
    row = projection.message_row(reaction)
    projected = row["inputs"][0][1]["components"][0]
    assert projected["smiles"] is None
    assert projected["identifiers"][0]["value"] == "ice water"


def test_reaction_smiles_is_generated_when_the_source_records_none():
    row = projection.message_row(_reaction())
    assert row["smiles"] is not None
    assert ">>" in row["smiles"]


def test_a_generated_reaction_smiles_is_complete_or_null():
    # Generating around an unreadable reactant would put a SMILES missing one in the
    # same column as a complete one, with nothing to tell them apart.
    reaction = _reaction()
    mystery = reaction.inputs["toluene"].components.add()
    mystery.identifiers.add(type="NAME", value="mystery reactant")
    assert projection.message_row(reaction)["smiles"] is None


def test_an_agent_recorded_only_by_name_does_not_null_the_reaction_smiles():
    # Agents are excluded from the generated SMILES, so an unreadable one is silent.
    reaction = _reaction()
    catalyst = reaction.inputs["cat"].components.add(reaction_role="CATALYST")
    catalyst.identifiers.add(type="NAME", value="ItBu")
    assert projection.message_row(reaction)["smiles"] is not None


def test_empty_repeated_fields_are_null_rather_than_empty_lists():
    row = projection.message_row(reaction_pb2.Reaction(reaction_id="ord-0003"))
    assert row["outcomes"] is None
    assert row["observations"] is None


def test_a_reaction_id_that_disagrees_with_its_column_is_an_error():
    with pytest.raises(ValueError, match="disagrees"):
        projection.reaction_row(_reaction("ord-in-message"), "ord-in-column")


@pytest.mark.parametrize("column", ["", None])
def test_a_reaction_id_column_with_no_value_is_an_error(column):
    # None is required to travel as itself rather than as a defaulted "no column was
    # read", which would make a null row skip the check entirely.
    with pytest.raises(ValueError, match="reaction_id column is"):
        projection.reaction_row(_reaction("ord-in-message"), column)
    with pytest.raises(ValueError, match="reaction_id column is"):
        projection.reaction_row(reaction_pb2.Reaction(), column)


def test_a_matching_column_projects_the_message_id():
    # Checked, never read: DatasetView.md5 does not hash the column, so a value taken
    # from it would sit outside the staleness hash.
    row = projection.reaction_row(_reaction("ord-0001"), "ord-0001")
    assert row["reaction_id"] == "ord-0001"


# Writing


def test_write_projection_round_trips(tmp_path):
    source = _source(tmp_path, [_reaction("ord-1"), _reaction("ord-2")])
    output = tmp_path / "projection.parquet"
    assert projection.write_projection(source, output) == 2
    table = pq.read_table(output)
    assert table.num_rows == 2
    assert table.column("reaction_id").to_pylist() == ["ord-1", "ord-2"]


def test_write_projection_reads_the_reaction_id_column(tmp_path):
    # Only a file built outside DatasetWriter can disagree, since the writer fills the
    # column from the message; the projection has to notice rather than pick one.
    path = tmp_path / "inconsistent.parquet"
    blob = _reaction("ord-in-message").SerializeToString(deterministic=True)
    schema = pq.read_schema(_source(tmp_path, [_reaction()]))
    pq.write_table(
        pa.table({"reaction_id": ["ord-in-column"], "reaction": [blob]}, schema=schema),
        path,
    )
    with pytest.raises(ValueError, match=f"{path}.*disagrees"):
        projection.write_projection(path, tmp_path / "projection.parquet")


def test_write_projection_rejects_a_null_reaction_id_column(tmp_path):
    # A producer declaring the column nullable can put a null in it, which
    # iter_reactions hands through as None; the blob's ID must not stand in for it.
    path = tmp_path / "null_id.parquet"
    reaction = _reaction("ord-in-message")
    schema = pa.schema(
        [
            pa.field("reaction_id", pa.string(), nullable=True),
            pa.field("reaction", pa.binary(), nullable=False),
        ]
    ).with_metadata(pq.read_schema(_source(tmp_path, [_reaction()])).metadata)
    pq.write_table(
        pa.table(
            {
                "reaction_id": [None],
                "reaction": [reaction.SerializeToString(deterministic=True)],
            },
            schema=schema,
        ),
        path,
    )
    with pytest.raises(ValueError, match="reaction_id column is None"):
        projection.write_projection(path, tmp_path / "projection.parquet")


def test_write_projection_covers_every_row_group(tmp_path):
    # Every other test fits in one row group, so a regression that stopped after the
    # first would truncate real datasets to their first row_group_size reactions.
    source = _source(
        tmp_path, [_reaction(f"ord-{index}") for index in range(5)], row_group_size=2
    )
    assert parquet.DatasetView(source).num_row_groups == 3
    output = tmp_path / "projection.parquet"
    assert projection.write_projection(source, output) == 5
    assert pq.read_table(output).column("reaction_id").to_pylist() == [
        f"ord-{index}" for index in range(5)
    ]


def test_write_projection_stamps_the_footer(tmp_path):
    source = _source(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    stamps = artifacts.load_stamps(output)
    assert stamps.artifact == projection.ARTIFACT
    assert stamps.artifact_version == artifacts.ARTIFACT_VERSION
    assert stamps.source_dataset_id == "ord_dataset-1"


def test_is_current_tracks_source_content(tmp_path):
    source = _source(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    assert projection.is_current(output, parquet.DatasetView(source).md5())
    assert not projection.is_current(output, "0" * 32)


def test_is_current_rejects_a_different_artifact(tmp_path):
    source = _source(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    assert not artifacts.is_current(output, "view", parquet.DatasetView(source).md5())


def test_is_current_is_false_for_a_missing_file(tmp_path):
    assert not projection.is_current(tmp_path / "absent.parquet", "0" * 32)


def test_is_current_is_false_for_a_source_dataset(tmp_path):
    source = _source(tmp_path, [_reaction()])
    assert not projection.is_current(source, parquet.DatasetView(source).md5())


def test_a_failed_write_leaves_an_existing_projection_intact(tmp_path, monkeypatch):
    source = _source(tmp_path, [_reaction()])
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


def test_precision_without_a_value_is_null():
    # The pair has to agree: a null setpoint beside a stated uncertainty reads as a
    # measurement nobody took but everybody bounded.
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.precision = 0.5
    setpoint.units = reaction_pb2.Temperature.CELSIUS
    temperature = projection.message_row(reaction)["conditions"]["temperature"]
    assert temperature["setpoint_kelvin"] is None
    assert temperature["setpoint_precision_kelvin"] is None


def test_enum_fields_carry_their_members_in_metadata():
    # Arrow has no type that records an enum's choices: a dictionary column carries only
    # the values a batch contained, and DuckDB reads it back as VARCHAR either way.
    # Metadata does survive, nested to any depth, so the published file says what the
    # spellings are without anyone holding this library.
    components = (
        projection.SCHEMA.field("inputs").type.item_type.field("components").type
    )
    members = projection.enum_members(components.value_type.field("reaction_role"))
    assert members is not None
    assert members[:2] == ("UNSPECIFIED", "REACTANT")
    assert projection.enum_members(projection.SCHEMA.field("reaction_id")) is None


def test_enum_metadata_survives_a_parquet_round_trip(tmp_path):
    source = _source(tmp_path, [_reaction()])
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    schema = pq.read_schema(output)
    components = schema.field("inputs").type.item_type.field("components").type
    members = projection.enum_members(components.value_type.field("reaction_role"))
    assert members is not None
    assert members[:2] == ("UNSPECIFIED", "REACTANT")


def test_a_repeated_united_field_is_refused():
    # A tripwire for a schema change upstream, so it needs its own proof: a repeated
    # united field takes the same branch as a singular one, and expanding it would
    # publish two scalar columns holding one measurement per reaction where the source
    # holds several. No such field exists today, hence the synthetic descriptor.
    field = SimpleNamespace(
        name="durations",
        full_name="ord.Fake.durations",
        message_type=SimpleNamespace(name="Time"),
        label=FieldDescriptor.LABEL_REPEATED,
    )
    descriptor = SimpleNamespace(name="Fake", full_name="ord.Fake", fields=[field])
    with pytest.raises(ValueError, match="repeated united message"):
        projection._struct_fields(cast(Descriptor, descriptor), frozenset())
