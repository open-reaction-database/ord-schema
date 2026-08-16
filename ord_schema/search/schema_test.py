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

"""Tests for ord_schema.search.schema."""

from importlib import resources

import pyarrow as pa
import pytest

from ord_schema.artifacts import projection
from ord_schema.search import schema


def test_scalars():
    rendered = schema.describe(
        pa.schema([pa.field("reaction_id", pa.string()), pa.field("n", pa.int64())])
    )
    assert rendered == "reaction_id: VARCHAR\nn: BIGINT"


def test_struct_children_are_indented():
    rendered = schema.describe(
        pa.schema(
            [
                pa.field(
                    "provenance",
                    pa.struct([pa.field("doi", pa.string())]),
                )
            ]
        )
    )
    assert rendered == "provenance: STRUCT\n  doi: VARCHAR"


def test_list_of_scalars_names_its_element_type():
    rendered = schema.describe(pa.schema([pa.field("names", pa.list_(pa.string()))]))
    assert rendered == "names: LIST<VARCHAR>"


def test_list_of_structs_expands():
    rendered = schema.describe(
        pa.schema(
            [
                pa.field(
                    "identifiers",
                    pa.list_(pa.struct([pa.field("value", pa.string())])),
                )
            ]
        )
    )
    assert rendered == "identifiers: LIST<STRUCT>\n  value: VARCHAR"


def test_map_names_its_key_type_and_expands_its_value():
    rendered = schema.describe(
        pa.schema(
            [
                pa.field(
                    "inputs",
                    pa.map_(pa.string(), pa.struct([pa.field("n", pa.int32())])),
                )
            ]
        )
    )
    assert rendered == "inputs: MAP<VARCHAR, STRUCT>\n  n: INTEGER"


def test_nesting_compounds():
    rendered = schema.describe(
        pa.schema(
            [
                pa.field(
                    "inputs",
                    pa.map_(
                        pa.string(),
                        pa.struct(
                            [
                                pa.field(
                                    "components",
                                    pa.list_(
                                        pa.struct([pa.field("smiles", pa.string())])
                                    ),
                                )
                            ]
                        ),
                    ),
                )
            ]
        )
    )
    assert rendered == (
        "inputs: MAP<VARCHAR, STRUCT>\n  components: LIST<STRUCT>\n    smiles: VARCHAR"
    )


def test_unnamed_leaf_type_raises():
    # A leaf the projection cannot currently hold, so the description would silently
    # hand a model an Arrow spelling no DuckDB query can use.
    with pytest.raises(KeyError):
        schema.describe(pa.schema([pa.field("when", pa.timestamp("s"))]))


def test_projection_schema_renders():
    rendered = schema.describe()
    lines = rendered.splitlines()
    assert len(lines) == sum(1 for _ in _fields(projection.SCHEMA))
    # Unit-carrying names are what let the prompt say nothing about units, and the
    # indentation is the path: conditions -> temperature/pressure -> setpoint.
    assert "conditions: STRUCT" in lines
    assert "    setpoint_kelvin: DOUBLE" in lines
    # The collapsed structural identifier is a top-level column, not nested.
    assert "smiles: VARCHAR" in lines
    assert "reaction_id: VARCHAR" in lines


def _fields(schema_or_type):
    """Yields every renderable field reachable from a schema or nested type.

    Internal fields (structure_id) are skipped to match describe(), which hides them:
    they join artifacts together rather than stating a fact a model may query.
    """
    for field in schema_or_type:
        if projection.is_internal(field):
            continue
        yield field
        data_type = field.type
        if pa.types.is_list(data_type):
            data_type = data_type.value_type
        elif pa.types.is_map(data_type):
            data_type = data_type.item_type
        if pa.types.is_struct(data_type):
            yield from _fields(data_type)


def test_the_snapshot_matches_what_describe_renders():
    # The rendered schema is what a model is told it may query, so a proto field added
    # upstream silently rewrites the prompt. Checking a snapshot in puts that change in
    # the diff, where it can be read, rather than leaving it to be discovered from a
    # model's behavior. The snapshot is a review artifact, never a source: describe()
    # generates what actually ships, and this test is what keeps the two honest.
    path = resources.files("ord_schema.search") / "projection_schema.txt"
    assert path.read_text(encoding="utf-8") == schema.describe() + "\n", (
        "the projection schema changed; regenerate the snapshot with\n"
        "  python -c 'from ord_schema.search import schema; "
        'open("ord_schema/search/projection_schema.txt", "w").write(schema.describe() '
        "+ chr(10))'"
    )


def test_enum_columns_carry_their_members():
    lines = schema.describe().splitlines()
    role = next(line for line in lines if line.strip().startswith("reaction_role:"))
    assert "REACTANT" in role
    assert "SOLVENT" in role
    # A plain string column is not decorated, so the parentheses mean something.
    assert next(line for line in lines if line.startswith("reaction_id:")).endswith(
        "VARCHAR"
    )


def test_only_a_field_carrying_members_is_annotated():
    rendered = schema.describe(
        pa.schema(
            [
                pa.field("kind", pa.string(), metadata={b"ord.enum": b"A,B"}),
                pa.field("free", pa.string()),
            ]
        )
    )
    assert rendered == "kind: VARCHAR  (A | B)\nfree: VARCHAR"


def test_internal_columns_are_hidden():
    # structure_id joins the projection to its structures artifact; its values are not
    # stable across builds, so a model must not learn the name to compare against.
    assert "structure_id" not in schema.describe()
