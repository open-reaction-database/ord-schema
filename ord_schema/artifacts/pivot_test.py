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

"""Tests for ord_schema.artifacts.pivot."""

import logging
import pathlib

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import base, pivot, projection
from ord_schema.proto import dataset_pb2, reaction_pb2


def test_covers_the_levels_a_query_quantifies_over():
    levels = pivot.repeated_levels()
    for path in (
        "workups",
        "inputs",
        "inputs.components",
        "outcomes",
        "outcomes.products",
        "outcomes.products.measurements",
        "conditions.temperature.measurements",
    ):
        assert path in levels


def test_ordinals_accumulate_down_the_path():
    levels = pivot.repeated_levels()
    assert levels["workups"].ordinals == ("workup_index",)
    assert levels["outcomes"].ordinals == ("outcome_index",)
    assert levels["outcomes.products"].ordinals == ("outcome_index", "product_index")
    assert levels["outcomes.products.measurements"].ordinals == (
        "outcome_index",
        "product_index",
        "measurement_index",
    )
    # A map is a repeated level too, so its components sit one ordinal deeper.
    assert levels["inputs.components"].ordinals == ("input_index", "component_index")


def test_element_type_drops_repeated_fields_and_keeps_structs():
    element = pivot.repeated_levels()["outcomes.products"].element_type
    names = [field.name for field in element]
    assert "isolated_color" in names
    assert "texture" in names
    assert "measurements" not in names
    assert "identifiers" not in names
    assert "analyses" not in names
    assert [field.name for field in element.field("texture").type] == [
        "type",
        "details",
    ]


def test_a_repeated_field_inside_a_kept_struct_is_dropped_too():
    element = pivot.repeated_levels()["outcomes.products.measurements"].element_type
    standard = element.field("authentic_standard").type
    names = [field.name for field in standard]
    assert "smiles" in names
    assert "identifiers" not in names
    assert "preparations" not in names


def test_the_structure_id_is_kept():
    # A pivot can answer a structure predicate, which needs the corpus-wide ID the
    # schema marks internal for a model's benefit rather than this walk's.
    element = pivot.repeated_levels()["inputs.components"].element_type
    assert "structure_id" in [field.name for field in element]


def test_a_level_whose_elements_are_all_repeated_is_not_covered():
    schema = pa.schema(
        [
            pa.field(
                "things",
                pa.list_(pa.struct([pa.field("bits", pa.list_(pa.int32()))])),
            )
        ]
    )
    assert pivot.repeated_levels(schema) == {}


def test_a_list_of_scalars_is_not_a_level():
    schema = pa.schema([pa.field("names", pa.list_(pa.string()))])
    assert pivot.repeated_levels(schema) == {}


def test_enum_metadata_survives_pruning():
    # The compiler validates an enum comparison against the members on the field, so a
    # pruned type that dropped them would refuse a query the projection accepts.
    element = pivot.repeated_levels()["inputs.components"].element_type
    field = element.field("reaction_role")
    assert field.metadata is not None


def test_table_name_is_a_bare_identifier():
    assert pivot.table_name("outcomes.products") == "pivot_outcomes_products"
    assert pivot.table_name("workups") == "pivot_workups"


def test_ordinals_that_would_collide_are_refused():
    # Two segments singularizing to one name would give a row two columns of that name,
    # and the second would silently win.
    schema = pa.schema(
        [
            pa.field(
                "analyses",
                pa.list_(
                    pa.struct(
                        [
                            pa.field(
                                "analysis",
                                pa.list_(pa.struct([pa.field("value", pa.int32())])),
                            )
                        ]
                    )
                ),
            )
        ]
    )
    with pytest.raises(ValueError, match="analysis_index"):
        pivot.repeated_levels(schema)


@pytest.fixture(scope="module")
def projected(tmp_path_factory) -> pathlib.Path:
    """A written projection with two outcomes on one reaction and none on another."""
    root = tmp_path_factory.mktemp("pivot")
    split = reaction_pb2.Reaction(reaction_id="ord-pv01")
    first = split.outcomes.add()
    first.products.add(isolated_color="white", is_desired_product=False)
    second = split.outcomes.add()
    second.products.add(isolated_color="yellow", is_desired_product=True)
    second.products.add(isolated_color="red", is_desired_product=False)
    source = root / "ord_dataset-pv.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-pv",
            name="test",
            description="test",
            reactions=[split, reaction_pb2.Reaction(reaction_id="ord-pv02")],
        ),
        str(source),
    )
    output = root / "projection.parquet"
    projection.write_projection(source, output)
    return output


def test_a_written_pivot_carries_its_rows_and_its_ordinals(projected, tmp_path):
    output = tmp_path / "products.parquet"
    written = pivot.write_pivot(projected, output, level_path="outcomes.products")
    assert written == 3
    table = pq.read_table(output)
    assert table.column_names == [
        "reaction_id",
        "outcome_index",
        "product_index",
        "element",
    ]
    rows = [
        (
            row["reaction_id"],
            row["outcome_index"],
            row["product_index"],
            row["element"]["isolated_color"],
        )
        for row in table.to_pylist()
    ]
    assert sorted(rows) == [
        ("ord-pv01", 1, 1, "white"),
        ("ord-pv01", 2, 1, "yellow"),
        ("ord-pv01", 2, 2, "red"),
    ]


def test_a_written_pivot_is_stamped_with_its_level(projected, tmp_path):
    output = tmp_path / "products.parquet"
    pivot.write_pivot(projected, output, level_path="outcomes.products")
    stamps = base.load_stamps(output)
    assert stamps.artifact == pivot.ARTIFACT
    assert stamps.source_dataset_id == "ord_dataset-pv"
    assert stamps.source_md5 == base.load_stamps(projected).source_md5
    assert base.stamps_are_current(stamps, pivot.ARTIFACT)
    assert pivot.pivot_path(output) == "outcomes.products"


def test_a_level_with_no_elements_still_writes_a_readable_file(projected, tmp_path):
    # A reader globs the level it wants; a file that is missing and a file that is
    # empty are different things, and only the second says the build ran.
    output = tmp_path / "workups.parquet"
    assert pivot.write_pivot(projected, output, level_path="workups") == 0
    assert pq.read_table(output).num_rows == 0
    assert pivot.pivot_path(output) == "workups"


def test_a_level_with_no_elements_is_never_unnested(projected, tmp_path, caplog):
    # The saving: discovering that a level holds nothing by unnesting it costs a full
    # pass over every ancestor, and most levels of this schema hold nothing.
    with caplog.at_level(logging.INFO, logger="ord_schema.artifacts.pivot"):
        assert (
            pivot.write_pivot(projected, tmp_path / "w.parquet", level_path="workups")
            == 0
        )
    assert any("no elements at workups" in record.message for record in caplog.records)


@pytest.mark.parametrize("level_path", sorted(pivot.LEVELS))
def test_counting_a_level_agrees_with_unnesting_it(projected, tmp_path, level_path):
    # The count decides whether the unnest runs at all, so a count that disagreed with
    # it would write an empty artifact over a level that holds rows.
    connection = duckdb.connect()
    try:
        connection.read_parquet(str(projected)).create_view("reactions")
        # S608: the expression comes from the module's own walk of the schema.
        row = connection.execute(
            f"SELECT {pivot.count_expression(pivot.LEVELS[level_path])} "  # noqa: S608
            "FROM reactions"
        ).fetchone()
        assert row is not None
        counted = row[0]
    finally:
        connection.close()
    written = pivot.write_pivot(
        projected, tmp_path / "level.parquet", level_path=level_path
    )
    assert counted == written


def test_a_path_that_is_not_a_level_is_refused(projected, tmp_path):
    with pytest.raises(ValueError, match="not a repeated level"):
        pivot.write_pivot(
            projected, tmp_path / "x.parquet", level_path="conditions.temperature"
        )


def test_pivoting_something_that_is_not_a_projection_is_refused(projected, tmp_path):
    output = tmp_path / "products.parquet"
    pivot.write_pivot(projected, output, level_path="outcomes.products")
    with pytest.raises(ValueError, match="not a projection"):
        pivot.write_pivot(
            output, tmp_path / "again.parquet", level_path="outcomes.products"
        )


def test_reach_finds_a_level_that_is_the_path_itself():
    reached = pivot.reach("outcomes.products")
    assert reached is not None
    level, remainder, dtype = reached
    assert level.path == "outcomes.products"
    assert remainder == ()
    assert dtype is level.element_type


def test_reach_descends_to_a_singular_struct_under_a_level():
    # One authentic standard per measurement rather than a list of its own, so the
    # level it ranges over is the measurements, and the pivot over those carries it.
    reached = pivot.reach("outcomes.products.measurements.authentic_standard")
    assert reached is not None
    level, remainder, dtype = reached
    assert level.path == "outcomes.products.measurements"
    assert remainder == ("authentic_standard",)
    assert "smiles" in [field.name for field in dtype]


def test_reach_declines_a_path_leaving_the_pruned_element():
    # measurements is a repeated field of a product, so it is not on the pruned type;
    # it is a level of its own, and reach finds that level rather than descending.
    reached = pivot.reach("outcomes.products.measurements")
    assert reached is not None
    level, remainder, _ = reached
    assert level.path == "outcomes.products.measurements"
    assert remainder == ()
    # A repeated field whose elements are scalars is no level of its own -- there is no
    # element struct to pivot -- and pruning removed it from the measurement, so
    # descending to it finds nothing either way.
    assert (
        pivot.reach("outcomes.products.measurements.mass_spec_details.eic_masses")
        is None
    )
    # Whereas a repeated field of structs is a level, reached as itself rather than by
    # descending into the measurement that holds it.
    nested = pivot.reach(
        "outcomes.products.measurements.authentic_standard.identifiers"
    )
    assert nested is not None
    standard, remainder, _ = nested
    assert standard.path == (
        "outcomes.products.measurements.authentic_standard.identifiers"
    )
    assert remainder == ()


def test_reach_declines_a_path_with_no_repeated_ancestor():
    assert pivot.reach("conditions.temperature") is None
    assert pivot.reach("reaction_id") is None
