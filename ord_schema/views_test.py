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

"""Tests for ord_schema.views."""

import pathlib
from importlib import metadata
from typing import NamedTuple

import pyarrow.parquet as pq
import pytest

from ord_schema import artifacts, parquet, projection, views
from ord_schema.proto import dataset_pb2, reaction_pb2


def _compound(smiles: str | None = None, *, name: str | None = None):
    compound = reaction_pb2.Compound()
    if smiles is not None:
        compound.identifiers.add(type="SMILES", value=smiles)
    if name is not None:
        compound.identifiers.add(type="NAME", value=name)
    return compound


def _reaction(reaction_id: str = "ord-0001") -> reaction_pb2.Reaction:
    """A minimal reaction: one benzene input, one toluene product."""
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.identifiers.add(type="REACTION_SMILES", value="c1ccccc1>>Cc1ccccc1")
    reaction.inputs["a"].components.append(_compound("c1ccccc1"))
    outcome = reaction.outcomes.add()
    outcome.products.add().identifiers.add(type="SMILES", value="Cc1ccccc1")
    return reaction


def _dataset(*reactions, dataset_id="ord_dataset-00000000000000000000000000000000"):
    return dataset_pb2.Dataset(
        dataset_id=dataset_id,
        name="test",
        description="desc",
        reactions=list(reactions) or [_reaction()],
    )


def test_reaction_row_projects_the_documented_columns():
    row = views.reaction_row(projection.message_row(_reaction()))
    assert set(row) == set(views.SCHEMA.names)
    assert row["reaction_id"] == "ord-0001"
    assert row["reaction_smiles"] == "c1ccccc1>>Cc1ccccc1"
    assert row["input_smiles"] == ["c1ccccc1"]
    assert row["output_smiles"] == ["Cc1ccccc1"]


def test_reaction_row_is_null_where_the_source_is_silent():
    row = views.reaction_row(projection.message_row(_reaction()))
    for column in (
        "yield_percent",
        "temperature_kelvin",
        "reaction_time_seconds",
        "doi",
        "patent",
    ):
        assert row[column] is None, column


def test_reaction_smiles_is_generated_when_not_stored():
    reaction = _reaction()
    del reaction.identifiers[:]
    row = views.reaction_row(projection.message_row(reaction))
    assert row["reaction_smiles"] == "c1ccccc1>>Cc1ccccc1"


def test_a_recorded_reaction_smiles_keeps_its_atom_mapping():
    # The reason to keep the deposited string rather than regenerate: generation builds
    # from the components and cannot reconstruct a mapping.
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    reaction.identifiers.add(
        type="REACTION_SMILES",
        value="[CH3:1][OH:2].[C:3](=O)O>CCO.[Na+]>[CH3:1][O:2][C:3]=O",
    )
    value = views.reaction_row(projection.message_row(reaction))["reaction_smiles"]
    for atom_map in (":1]", ":2]", ":3]"):
        assert atom_map in value
    # ...with the agents normalized away.
    assert value.count(">") == 2
    assert "[Na+]" not in value


def test_agents_are_left_out_of_the_generated_reaction_smiles():
    # An empty agent block is idiomatic, and excluding agents means a ligand recorded
    # only by name cannot decide whether the reaction gets a SMILES at all.
    reaction = _reaction()
    del reaction.identifiers[:]
    catalyst = reaction.inputs["cat"].components.add(reaction_role="CATALYST")
    catalyst.identifiers.add(type="SMILES", value="[Pd]")
    assert (
        views.reaction_row(projection.message_row(reaction))["reaction_smiles"]
        == "c1ccccc1>>Cc1ccccc1"
    )


def test_a_ligand_recorded_only_by_name_does_not_null_the_reaction_smiles():
    reaction = _reaction()
    del reaction.identifiers[:]
    ligand = reaction.inputs["cat"].components.add(reaction_role="CATALYST")
    ligand.identifiers.add(type="NAME", value="ItBu")
    assert (
        views.reaction_row(projection.message_row(reaction))["reaction_smiles"]
        == "c1ccccc1>>Cc1ccccc1"
    )


def test_an_unreadable_reactant_nulls_a_generated_reaction_smiles():
    # Strict where it counts: a SMILES missing a reactant describes another reaction.
    reaction = _reaction()
    del reaction.identifiers[:]  # Nothing recorded, so it has to be generated.
    reaction.inputs["b"].components.add().identifiers.add(
        type="NAME", value="mystery reactant"
    )
    assert (
        views.reaction_row(projection.message_row(reaction))["reaction_smiles"] is None
    )


def test_enhanced_stereochemistry_survives_generation():
    # Component CXSMILES blocks cannot be string-joined into a reaction; the reaction is
    # assembled from parsed components so RDKit re-emits one block with fixed indices.
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    component = reaction.inputs["a"].components.add()
    component.identifiers.add(type="CXSMILES", value="C[C@H](N)O |o1:1|")
    reaction.outcomes.add().products.add().identifiers.add(
        type="SMILES", value="C[C@H](N)OC"
    )
    assert (
        "|o1:"
        in views.reaction_row(projection.message_row(reaction))["reaction_smiles"]
    )


def test_plain_reaction_smiles_is_used_when_that_is_all_there_is():
    assert (
        views.reaction_row(projection.message_row(_reaction()))["reaction_smiles"]
        == "c1ccccc1>>Cc1ccccc1"
    )


def test_component_cxsmiles_is_preferred_over_plain_smiles():
    reaction = _reaction()
    del reaction.inputs["a"]
    compound = reaction_pb2.Compound()
    compound.identifiers.add(type="SMILES", value="C[C@H](N)O")
    compound.identifiers.add(type="CXSMILES", value="C[C@H](N)O |o1:1|")
    reaction.inputs["a"].components.append(compound)
    assert views.reaction_row(projection.message_row(reaction))["input_smiles"] == [
        "C[C@H](N)O |o1:1|"
    ]


def test_component_smiles_are_canonicalized():
    reaction = _reaction()
    del reaction.inputs["a"]
    reaction.inputs["a"].components.append(_compound("C1=CC=CC=C1"))
    assert views.reaction_row(projection.message_row(reaction))["input_smiles"] == [
        "c1ccccc1"
    ]


def test_component_smiles_derive_from_other_structural_identifiers():
    reaction = _reaction()
    del reaction.inputs["a"]
    compound = reaction_pb2.Compound()
    compound.identifiers.add(type="INCHI", value="InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")
    reaction.inputs["a"].components.append(compound)
    assert views.reaction_row(projection.message_row(reaction))["input_smiles"] == [
        "c1ccccc1"
    ]


def test_unparseable_component_smiles_are_skipped():
    reaction = _reaction()
    reaction.inputs["a"].components.append(_compound("not-a-smiles"))
    assert views.reaction_row(projection.message_row(reaction))["input_smiles"] == [
        "c1ccccc1"
    ]


def test_reaction_smiles_is_null_when_it_cannot_be_generated():
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    reaction.inputs["a"].components.append(_compound(name="benzene"))
    assert (
        views.reaction_row(projection.message_row(reaction))["reaction_smiles"] is None
    )


def test_input_smiles_follow_sorted_keys_not_insertion_order():
    # Realistic keys on purpose: protobuf map iteration for single-character keys
    # happens to come out sorted, so "a"/"b"/"c" would pass with or without the sort.
    keys = (("m0_base", "C"), ("m1_solvent", "CC"), ("m2_reagent", "CCC"))
    forward, backward = _reaction(), _reaction()
    del forward.inputs["a"], backward.inputs["a"]
    for key, smiles in keys:
        forward.inputs[key].components.append(_compound(smiles))
    for key, smiles in reversed(keys):
        backward.inputs[key].components.append(_compound(smiles))
    expected = ["C", "CC", "CCC"]
    assert (
        views.reaction_row(projection.message_row(forward))["input_smiles"] == expected
    )
    assert (
        views.reaction_row(projection.message_row(backward))["input_smiles"] == expected
    )


def test_components_without_smiles_are_skipped():
    reaction = _reaction()
    reaction.inputs["a"].components.append(_compound(name="mystery solvent"))
    assert views.reaction_row(projection.message_row(reaction))["input_smiles"] == [
        "c1ccccc1"
    ]


def test_yield_takes_the_largest_measurement_on_any_product():
    reaction = _reaction()
    outcome = reaction.outcomes[0]
    for value in (12.0, 87.5, 40.0):
        measurement = outcome.products[0].measurements.add(type="YIELD")
        measurement.percentage.value = value
    assert views.reaction_row(projection.message_row(reaction))[
        "yield_percent"
    ] == pytest.approx(87.5)


def test_non_yield_measurements_are_ignored():
    reaction = _reaction()
    measurement = reaction.outcomes[0].products[0].measurements.add(type="PURITY")
    measurement.percentage.value = 99.0
    assert views.reaction_row(projection.message_row(reaction))["yield_percent"] is None


def test_yield_spans_every_outcome():
    # output_smiles already covers every outcome; reading the yield off the first only
    # would report null for a yield the source states, in a column where null is
    # documented to mean the source is silent.
    reaction = _reaction()
    reaction.outcomes[0].products[0].measurements.add(
        type="YIELD"
    ).percentage.value = 10.0
    second = reaction.outcomes.add()
    second.products.add().measurements.add(type="YIELD").percentage.value = 90.0
    assert views.reaction_row(projection.message_row(reaction))[
        "yield_percent"
    ] == pytest.approx(90.0)


def test_a_yield_recorded_only_on_a_later_outcome_is_not_null():
    reaction = _reaction()
    second = reaction.outcomes.add()
    second.products.add().measurements.add(type="YIELD").percentage.value = 85.0
    assert views.reaction_row(projection.message_row(reaction))[
        "yield_percent"
    ] == pytest.approx(85.0)


def test_a_non_finite_yield_is_skipped_rather_than_winning():
    # NaN compares False against everything, so leaving it in the running maximum makes
    # the published value depend on the order measurements happen to be recorded.
    for order in ([float("nan"), 91.0], [91.0, float("nan")]):
        reaction = _reaction()
        for value in order:
            reaction.outcomes[0].products[0].measurements.add(
                type="YIELD"
            ).percentage.value = value
        assert views.reaction_row(projection.message_row(reaction))[
            "yield_percent"
        ] == pytest.approx(91.0)


def test_a_yield_percentage_with_no_value_is_null_not_zero():
    # An unset value reads as 0.0 through protobuf, and a fabricated zero yield reads
    # as a failed reaction rather than an unrecorded one.
    reaction = _reaction()
    measurement = reaction.outcomes[0].products[0].measurements.add(type="YIELD")
    measurement.percentage.precision = 1.0
    assert views.reaction_row(projection.message_row(reaction))["yield_percent"] is None


def test_a_yield_the_source_records_as_zero_is_kept():
    # The complement of the test above, and the case that breaks first if presence
    # handling regresses anywhere on the path: Percentage.value carries explicit
    # presence, so a reaction reported as 0% has to reach the column as 0.0 rather
    # than collapse into the null that means "unrecorded".
    reaction = _reaction()
    measurement = reaction.outcomes[0].products[0].measurements.add(type="YIELD")
    measurement.percentage.value = 0.0
    assert views.reaction_row(projection.message_row(reaction))["yield_percent"] == 0.0


def test_a_yield_not_recorded_as_a_percentage_is_null():
    # Reading float_value would assume a scale the source never stated.
    reaction = _reaction()
    measurement = reaction.outcomes[0].products[0].measurements.add(type="YIELD")
    measurement.float_value.value = 72.0
    assert views.reaction_row(projection.message_row(reaction))["yield_percent"] is None


@pytest.mark.parametrize(
    ("units_enum", "value", "expected"),
    [
        (reaction_pb2.Temperature.KELVIN, 300.0, 300.0),
        (reaction_pb2.Temperature.CELSIUS, 25.0, 298.15),
        (reaction_pb2.Temperature.FAHRENHEIT, 77.0, 298.15),
    ],
)
def test_temperature_converts_to_kelvin(units_enum, value, expected):
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value, setpoint.units = value, units_enum
    row = views.reaction_row(projection.message_row(reaction))
    assert row["temperature_kelvin"] == pytest.approx(expected)


@pytest.mark.parametrize(
    ("units_enum", "value", "expected"),
    [
        (reaction_pb2.Time.SECOND, 90.0, 90.0),
        (reaction_pb2.Time.MINUTE, 2.0, 120.0),
        (reaction_pb2.Time.HOUR, 1.5, 5400.0),
        (reaction_pb2.Time.DAY, 1.0, 86400.0),
    ],
)
def test_reaction_time_converts_to_seconds(units_enum, value, expected):
    reaction = _reaction()
    reaction_time = reaction.outcomes[0].reaction_time
    reaction_time.value, reaction_time.units = value, units_enum
    assert views.reaction_row(projection.message_row(reaction))[
        "reaction_time_seconds"
    ] == pytest.approx(expected)


def test_measurements_with_unspecified_units_read_as_null():
    reaction = _reaction()
    reaction.conditions.temperature.setpoint.value = 300.0
    assert (
        views.reaction_row(projection.message_row(reaction))["temperature_kelvin"]
        is None
    )


def test_provenance_columns():
    reaction = _reaction()
    reaction.provenance.doi = "10.1000/example"
    reaction.provenance.patent = "US1234567"
    row = views.reaction_row(projection.message_row(reaction))
    assert row["doi"] == "10.1000/example"
    assert row["patent"] == "US1234567"


class _Projected(NamedTuple):
    """A projection to derive a view from, and the source hash it should stamp."""

    path: pathlib.Path
    source_md5: str


def _project(tmp_path, dataset=None, *, name="ds.parquet", **kwargs) -> _Projected:
    """Writes a source dataset and its projection, which is what a view reads."""
    source = tmp_path / name
    parquet.save_dataset(dataset or _dataset(), str(source), **kwargs)
    path = tmp_path / f"projection-{name}"
    projection.write_projection(source, path)
    return _Projected(path, parquet.DatasetView(source).md5())


def test_write_view_round_trip(tmp_path):
    reactions = [_reaction(f"ord-{i:04d}") for i in range(5)]
    source = _project(tmp_path, _dataset(*reactions))
    output = str(tmp_path / "view.parquet")
    assert views.write_view(source.path, output) == 5
    table = pq.read_table(output)
    assert table.num_rows == 5
    assert table.schema.names == views.SCHEMA.names
    assert table.column("reaction_id").to_pylist() == [f"ord-{i:04d}" for i in range(5)]
    assert table.column("input_smiles").to_pylist() == [["c1ccccc1"]] * 5


def test_write_view_spans_multiple_row_groups(tmp_path):
    reactions = [_reaction(f"ord-{i:04d}") for i in range(10)]
    source = _project(tmp_path, _dataset(*reactions), row_group_size=3)
    assert pq.ParquetFile(source.path).num_row_groups > 1
    output = str(tmp_path / "view.parquet")
    assert views.write_view(source.path, output) == 10
    assert pq.read_table(output).column("reaction_id").to_pylist() == [
        f"ord-{i:04d}" for i in range(10)
    ]


def test_write_view_stamps_the_footer(tmp_path):
    source = _project(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source.path, output)
    stamps = artifacts.load_stamps(output)
    assert stamps.artifact == views.ARTIFACT
    assert stamps.source_md5 == source.source_md5
    assert stamps.source_dataset_id == "ord_dataset-00000000000000000000000000000000"
    assert stamps.ord_schema_version == metadata.version("ord-schema")
    assert stamps.artifact_version == artifacts.ARTIFACT_VERSION


def test_a_view_is_not_current_as_a_projection(tmp_path):
    # One shared artifact version means the name is what separates the two.
    source = _project(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source.path, output)
    assert not artifacts.is_current(output, "projection", source.source_md5)


def test_is_current_tracks_source_content(tmp_path):
    source = _project(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source.path, output)
    assert views.is_current(output, source.source_md5)
    assert not views.is_current(output, "0" * 32)


def test_is_current_tracks_the_artifact_version(tmp_path, monkeypatch):
    source = _project(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source.path, output)
    monkeypatch.setattr(artifacts, "ARTIFACT_VERSION", "next")
    assert not views.is_current(output, source.source_md5)


def test_is_current_is_false_for_a_missing_view(tmp_path):
    assert not views.is_current(str(tmp_path / "absent.parquet"), "0" * 32)


def test_the_written_schema_is_the_documented_one(tmp_path):
    # Asserting against views.SCHEMA would only prove the writer used it. Consumers bind
    # to these names and types, so they are spelled out -- as Parquet reports them,
    # which renames the list element from "item" to "element" on the way through.
    source = _project(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source.path, output)
    schema = pq.read_schema(output)
    assert [(field.name, str(field.type)) for field in schema] == [
        ("reaction_id", "string"),
        ("reaction_smiles", "string"),
        ("input_smiles", "list<element: string>"),
        ("output_smiles", "list<element: string>"),
        ("yield_percent", "float"),
        ("temperature_kelvin", "float"),
        ("reaction_time_seconds", "float"),
        ("doi", "string"),
        ("patent", "string"),
    ]
    assert not schema.field("reaction_id").nullable


def test_write_view_round_trips_a_measurement_through_parquet(tmp_path):
    # float32 is the source's own precision -- proto Percentage.value is a float -- but
    # it is consumer-visible, so the value a query sees is pinned here.
    reaction = _reaction()
    reaction.outcomes[0].products[0].measurements.add(
        type="YIELD"
    ).percentage.value = 87.65
    source = _project(tmp_path, _dataset(reaction))
    output = str(tmp_path / "view.parquet")
    views.write_view(source.path, output)
    assert pq.read_table(output).column("yield_percent").to_pylist() == [
        pytest.approx(87.65, rel=1e-6)
    ]


def test_temperature_reads_the_setpoint_not_an_achieved_measurement(tmp_path):
    # One of several plausible readings, and the one the column is named for.
    del tmp_path
    reaction = _reaction()
    reaction.conditions.temperature.setpoint.value = 300.0
    reaction.conditions.temperature.setpoint.units = reaction_pb2.Temperature.KELVIN
    achieved = reaction.conditions.temperature.measurements.add()
    achieved.temperature.value = 350.0
    achieved.temperature.units = reaction_pb2.Temperature.KELVIN
    assert views.reaction_row(projection.message_row(reaction))[
        "temperature_kelvin"
    ] == pytest.approx(300.0)


def test_a_product_without_a_readable_structure_is_skipped(tmp_path):
    # Symmetric with the input side: a None inside the list would be a worse shape for a
    # consumer than a shorter list.
    del tmp_path
    reaction = _reaction()
    reaction.outcomes[0].products.add().identifiers.add(type="NAME", value="mystery")
    assert views.reaction_row(projection.message_row(reaction))["output_smiles"] == [
        "Cc1ccccc1"
    ]


def test_write_view_leaves_an_existing_view_intact_on_failure(tmp_path, monkeypatch):
    source = _project(tmp_path)
    output = tmp_path / "view.parquet"
    views.write_view(source.path, output)
    before = output.read_bytes()

    def _boom(*args, **kwargs):
        raise RuntimeError("derivation failed")

    monkeypatch.setattr(views, "reaction_row", _boom)
    # The source and the reaction have to be in the message: a corpus run has thousands
    # of datasets and millions of reactions, and neither is in the traceback otherwise.
    with pytest.raises(ValueError, match=r"ds\.parquet: ord-0001: derivation failed"):
        views.write_view(source.path, output)
    assert output.read_bytes() == before
    # The temp sibling must not survive either.
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "ds.parquet",
        "projection-ds.parquet",
        "view.parquet",
    ]


def test_components_identified_only_by_cxsmiles_are_kept():
    reaction = _reaction()
    del reaction.inputs["a"]
    compound = reaction_pb2.Compound()
    compound.identifiers.add(type="CXSMILES", value="c1ccccc1 |f:0.1|")
    reaction.inputs["a"].components.append(compound)
    assert views.reaction_row(projection.message_row(reaction))["input_smiles"] == [
        "c1ccccc1"
    ]


def test_every_projection_column_read_is_a_leaf_of_the_projection_schema(tmp_path):
    # pq.read_row_group drops a column path it cannot find and raises nothing, so a
    # leaf renamed in the projection would surface as a KeyError on the first reaction
    # that happens to populate it -- or not at all, for a path _nested short-circuits
    # past. Checked against the schema instead, where a rename fails immediately.
    path = tmp_path / "empty.parquet"
    pq.write_table(projection.SCHEMA.empty_table(), path)
    with pq.ParquetFile(path) as parquet_file:
        paths = {
            parquet_file.schema.column(i).path for i in range(len(parquet_file.schema))
        }
    assert set(views._PROJECTION_COLUMNS) <= paths


def test_write_view_populates_every_column_through_parquet(tmp_path):
    # The unit tests above narrow an in-memory row holding every projection column;
    # write_view narrows a column-pruned Parquet read. Only a round trip that populates
    # all nine columns exercises the paths in _PROJECTION_COLUMNS.
    reaction = _reaction()
    reaction.conditions.temperature.setpoint.value = 25.0
    reaction.conditions.temperature.setpoint.units = reaction_pb2.Temperature.CELSIUS
    outcome = reaction.outcomes[0]
    outcome.reaction_time.value = 2.0
    outcome.reaction_time.units = reaction_pb2.Time.HOUR
    outcome.products[0].measurements.add(type="YIELD").percentage.value = 87.5
    reaction.provenance.doi = "10.1000/example"
    reaction.provenance.patent = "US1234567"
    source = _project(tmp_path, _dataset(reaction))
    output = tmp_path / "view.parquet"
    assert views.write_view(source.path, output) == 1
    row = pq.read_table(output).to_pylist()[0]
    assert row == {
        "reaction_id": "ord-0001",
        "reaction_smiles": "c1ccccc1>>Cc1ccccc1",
        "input_smiles": ["c1ccccc1"],
        "output_smiles": ["Cc1ccccc1"],
        "yield_percent": pytest.approx(87.5),
        "temperature_kelvin": pytest.approx(298.15),
        "reaction_time_seconds": pytest.approx(7200.0),
        "doi": "10.1000/example",
        "patent": "US1234567",
    }


def test_reaction_time_takes_the_first_outcome_that_records_one(tmp_path):
    # Documented in the module docstring, and untested until now: mutating the guard so
    # a later outcome wins left the whole suite green.
    reaction = _reaction()
    reaction.outcomes[0].reaction_time.value = 1.0
    reaction.outcomes[0].reaction_time.units = reaction_pb2.Time.HOUR
    later = reaction.outcomes.add()
    later.reaction_time.value = 5.0
    later.reaction_time.units = reaction_pb2.Time.HOUR
    row = views.reaction_row(projection.message_row(reaction))
    assert row["reaction_time_seconds"] == pytest.approx(3600.0)
