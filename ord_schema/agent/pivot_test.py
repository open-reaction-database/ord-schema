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

"""Tests for ord_schema.agent.pivot."""

import pyarrow as pa
import pytest

from ord_schema.agent import pivot


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
