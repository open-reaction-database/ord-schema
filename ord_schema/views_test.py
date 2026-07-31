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

from importlib import metadata

import pyarrow.parquet as pq
import pytest

from ord_schema import parquet, views
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
    row = views.reaction_row("ord-0001", _reaction())
    assert set(row) == set(views.SCHEMA.names)
    assert row["reaction_id"] == "ord-0001"
    assert row["reaction_smiles"] == "c1ccccc1>>Cc1ccccc1"
    assert row["input_smiles"] == ["c1ccccc1"]
    assert row["output_smiles"] == ["Cc1ccccc1"]


def test_reaction_row_is_null_where_the_source_is_silent():
    row = views.reaction_row("ord-0001", _reaction())
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
    row = views.reaction_row("ord-0001", reaction)
    assert row["reaction_smiles"] == "c1ccccc1>>Cc1ccccc1"


_CXSMILES = "CC(=O)O.CCO>>CC(=O)OCC |f:0.1|"


def test_stored_reaction_cxsmiles_is_used_whole():
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    reaction.identifiers.add(type="REACTION_CXSMILES", value=_CXSMILES)
    assert views.reaction_row("ord-0001", reaction)["reaction_smiles"] == _CXSMILES


def test_cxsmiles_is_preferred_over_plain_smiles():
    # CXSMILES is a superset, so the richer identifier wins regardless of the order
    # the two appear in.
    reaction = _reaction()
    reaction.identifiers.add(type="REACTION_CXSMILES", value=_CXSMILES)
    assert views.reaction_row("ord-0001", reaction)["reaction_smiles"] == _CXSMILES


def test_plain_reaction_smiles_is_used_when_that_is_all_there_is():
    assert (
        views.reaction_row("ord-0001", _reaction())["reaction_smiles"]
        == "c1ccccc1>>Cc1ccccc1"
    )


def test_trailing_whitespace_is_preserved_rather_than_truncated():
    padded = "CC(=O)O.CCO>>CC(=O)OCC\xa0\xa0"
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    reaction.identifiers.add(type="REACTION_SMILES", value=padded)
    assert views.reaction_row("ord-0001", reaction)["reaction_smiles"] == padded


def test_component_cxsmiles_is_preferred_over_plain_smiles():
    reaction = _reaction()
    del reaction.inputs["a"]
    compound = reaction_pb2.Compound()
    compound.identifiers.add(type="SMILES", value="C[C@H](N)O")
    compound.identifiers.add(type="CXSMILES", value="C[C@H](N)O |o1:1|")
    reaction.inputs["a"].components.append(compound)
    assert views.reaction_row("x", reaction)["input_smiles"] == ["C[C@H](N)O |o1:1|"]


def test_component_smiles_are_canonicalized():
    reaction = _reaction()
    del reaction.inputs["a"]
    reaction.inputs["a"].components.append(_compound("C1=CC=CC=C1"))
    assert views.reaction_row("x", reaction)["input_smiles"] == ["c1ccccc1"]


def test_component_smiles_derive_from_other_structural_identifiers():
    reaction = _reaction()
    del reaction.inputs["a"]
    compound = reaction_pb2.Compound()
    compound.identifiers.add(type="INCHI", value="InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")
    reaction.inputs["a"].components.append(compound)
    assert views.reaction_row("x", reaction)["input_smiles"] == ["c1ccccc1"]


def test_unparseable_component_smiles_are_skipped():
    reaction = _reaction()
    reaction.inputs["a"].components.append(_compound("not-a-smiles"))
    assert views.reaction_row("x", reaction)["input_smiles"] == ["c1ccccc1"]


def test_reaction_smiles_is_null_when_it_cannot_be_generated():
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    reaction.inputs["a"].components.append(_compound(name="benzene"))
    assert views.reaction_row("ord-0001", reaction)["reaction_smiles"] is None


def test_input_smiles_follow_sorted_keys_not_insertion_order():
    forward, backward = _reaction(), _reaction()
    del forward.inputs["a"], backward.inputs["a"]
    for key, smiles in (("a", "C"), ("b", "CC"), ("c", "CCC")):
        forward.inputs[key].components.append(_compound(smiles))
    for key, smiles in (("c", "CCC"), ("b", "CC"), ("a", "C")):
        backward.inputs[key].components.append(_compound(smiles))
    expected = ["C", "CC", "CCC"]
    assert views.reaction_row("x", forward)["input_smiles"] == expected
    assert views.reaction_row("x", backward)["input_smiles"] == expected


def test_components_without_smiles_are_skipped():
    reaction = _reaction()
    reaction.inputs["a"].components.append(_compound(name="mystery solvent"))
    assert views.reaction_row("x", reaction)["input_smiles"] == ["c1ccccc1"]


def test_yield_takes_the_largest_measurement_on_any_product():
    reaction = _reaction()
    outcome = reaction.outcomes[0]
    for value in (12.0, 87.5, 40.0):
        measurement = outcome.products[0].measurements.add(type="YIELD")
        measurement.percentage.value = value
    assert views.reaction_row("x", reaction)["yield_percent"] == pytest.approx(87.5)


def test_non_yield_measurements_are_ignored():
    reaction = _reaction()
    measurement = reaction.outcomes[0].products[0].measurements.add(type="PURITY")
    measurement.percentage.value = 99.0
    assert views.reaction_row("x", reaction)["yield_percent"] is None


def test_yield_comes_from_the_first_outcome():
    reaction = _reaction()
    reaction.outcomes[0].products[0].measurements.add(
        type="YIELD"
    ).percentage.value = 10.0
    second = reaction.outcomes.add()
    second.products.add().measurements.add(type="YIELD").percentage.value = 90.0
    assert views.reaction_row("x", reaction)["yield_percent"] == pytest.approx(10.0)


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
    row = views.reaction_row("x", reaction)
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
    assert views.reaction_row("x", reaction)["reaction_time_seconds"] == pytest.approx(
        expected
    )


def test_measurements_with_unspecified_units_read_as_null():
    reaction = _reaction()
    reaction.conditions.temperature.setpoint.value = 300.0
    assert views.reaction_row("x", reaction)["temperature_kelvin"] is None


def test_provenance_columns():
    reaction = _reaction()
    reaction.provenance.doi = "10.1000/example"
    reaction.provenance.patent = "US1234567"
    row = views.reaction_row("x", reaction)
    assert row["doi"] == "10.1000/example"
    assert row["patent"] == "US1234567"


def _write_source(tmp_path, dataset=None, *, name="ds.parquet", **kwargs):
    path = str(tmp_path / name)
    parquet.save_dataset(dataset or _dataset(), path, **kwargs)
    return path


def test_write_view_round_trip(tmp_path):
    reactions = [_reaction(f"ord-{i:04d}") for i in range(5)]
    source = _write_source(tmp_path, _dataset(*reactions))
    output = str(tmp_path / "view.parquet")
    assert views.write_view(source, output) == 5
    table = pq.read_table(output)
    assert table.num_rows == 5
    assert table.schema.names == views.SCHEMA.names
    assert table.column("reaction_id").to_pylist() == [f"ord-{i:04d}" for i in range(5)]
    assert table.column("input_smiles").to_pylist() == [["c1ccccc1"]] * 5


def test_write_view_spans_multiple_row_groups(tmp_path):
    reactions = [_reaction(f"ord-{i:04d}") for i in range(10)]
    source = _write_source(tmp_path, _dataset(*reactions), row_group_size=3)
    assert parquet.num_row_groups(source) > 1
    output = str(tmp_path / "view.parquet")
    assert views.write_view(source, output) == 10
    assert pq.read_table(output).column("reaction_id").to_pylist() == [
        f"ord-{i:04d}" for i in range(10)
    ]


def test_write_view_stamps_the_footer(tmp_path):
    source = _write_source(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source, output)
    stamps = views.load_stamps(output)
    source_md5, _ = parquet.streaming_md5(source)
    assert stamps.source_md5 == source_md5
    assert stamps.source_dataset_id == "ord_dataset-00000000000000000000000000000000"
    assert stamps.ord_schema_version == metadata.version("ord-schema")
    assert stamps.view_version == views.VIEW_VERSION


def test_load_stamps_rejects_a_source_dataset(tmp_path):
    source = _write_source(tmp_path)
    with pytest.raises(ValueError, match="not a view"):
        views.load_stamps(source)


def test_is_current_tracks_source_content(tmp_path):
    source = _write_source(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source, output)
    source_md5, _ = parquet.streaming_md5(source)
    assert views.is_current(output, source_md5)
    assert not views.is_current(output, "0" * 32)


def test_is_current_tracks_the_view_version(tmp_path, monkeypatch):
    source = _write_source(tmp_path)
    output = str(tmp_path / "view.parquet")
    views.write_view(source, output)
    source_md5, _ = parquet.streaming_md5(source)
    monkeypatch.setattr(views, "VIEW_VERSION", f"{views.VIEW_VERSION}-next")
    assert not views.is_current(output, source_md5)


def test_is_current_is_false_for_a_missing_view(tmp_path):
    assert not views.is_current(str(tmp_path / "absent.parquet"), "0" * 32)


def test_write_view_leaves_an_existing_view_intact_on_failure(tmp_path, monkeypatch):
    source = _write_source(tmp_path)
    output = tmp_path / "view.parquet"
    views.write_view(source, output)
    before = output.read_bytes()

    def _boom(*args, **kwargs):
        raise RuntimeError("derivation failed")

    monkeypatch.setattr(views, "reaction_row", _boom)
    with pytest.raises(RuntimeError, match="derivation failed"):
        views.write_view(source, output)
    assert output.read_bytes() == before
    # The temp sibling must not survive either.
    assert sorted(p.name for p in tmp_path.iterdir()) == ["ds.parquet", "view.parquet"]


def test_components_identified_only_by_cxsmiles_are_kept():
    reaction = _reaction()
    del reaction.inputs["a"]
    compound = reaction_pb2.Compound()
    compound.identifiers.add(type="CXSMILES", value="c1ccccc1 |f:0.1|")
    reaction.inputs["a"].components.append(compound)
    assert views.reaction_row("x", reaction)["input_smiles"] == ["c1ccccc1"]
