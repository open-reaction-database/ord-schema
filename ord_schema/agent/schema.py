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

"""The projection schema, rendered for a model that has to write SQL against it.

A model cannot query columns it has not been told about, so the whole schema goes in the
prompt. It is small enough to: 442 leaves render as 537 lines and roughly 5,400 tokens,
which is why translation can stay a single forced tool call instead of a retrieval loop
over column metadata.

An enum projects as its value *name*, so each enum column carries its members beside
it. Without them a model has the column and no way to learn the spelling it must compare
against, which is a guess it will sometimes get wrong silently. They are read from the
field's own metadata, so this stays a rendering of the schema and nothing else.

The rendering is generated from ``projection.SCHEMA`` rather than written by hand,
because that schema is itself generated from the proto descriptors -- a field added
upstream becomes a column with nobody deciding it should, and a hand-written description
would silently fall behind it. ``ord_schema.agent.sql`` validates against the same
schema object, so the columns a model is told about and the columns its query is checked
against cannot disagree.

Types are named in DuckDB's vocabulary rather than Arrow's, since DuckDB is the dialect
the model writes. Units are already in the column names -- ``setpoint_kelvin``,
``mass_grams`` -- so the tree carries them without a word of prose, and nothing has to
explain that temperatures are kelvin.
"""

import pyarrow as pa

from ord_schema import projection

# Arrow leaf types the projection can hold, in DuckDB's names. Deliberately not a
# fallback to ``str(dtype)``: a type the projection gains later should surface here as a
# KeyError rather than reach a prompt as an Arrow spelling no DuckDB query can use.
_SCALARS: dict[pa.DataType, str] = {
    pa.string(): "VARCHAR",
    pa.binary(): "BLOB",
    pa.bool_(): "BOOLEAN",
    pa.float32(): "FLOAT",
    pa.float64(): "DOUBLE",
    pa.int32(): "INTEGER",
    pa.int64(): "BIGINT",
    pa.uint32(): "UINTEGER",
    pa.uint64(): "UBIGINT",
}

_INDENT = "  "


def _scalar_name(dtype: pa.DataType) -> str:
    """Returns the DuckDB name for a leaf type.

    Args:
        dtype: Arrow type of a leaf column.

    Returns:
        The DuckDB type name.

    Raises:
        KeyError: If the projection holds a leaf type this module does not name.
    """
    return _SCALARS[dtype]


def _render(field: pa.Field, depth: int, lines: list[str]) -> None:
    """Appends the lines describing one field, and recurses into its children.

    A container is named by its header (``LIST<STRUCT>``) with its children indented
    beneath, so the query that reaches a leaf can be read off the indentation. An enum
    leaf is followed by the values it may hold, read from the field's own metadata.

    Args:
        field: Field to describe.
        depth: Indentation depth.
        lines: Accumulator to append to.
    """
    indent = _INDENT * depth
    name, dtype = field.name, field.type
    if pa.types.is_struct(dtype):
        lines.append(f"{indent}{name}: STRUCT")
        for child in dtype:
            _render(child, depth + 1, lines)
    elif pa.types.is_list(dtype):
        value = dtype.value_type
        if pa.types.is_struct(value):
            lines.append(f"{indent}{name}: LIST<STRUCT>")
            for child in value:
                _render(child, depth + 1, lines)
        else:
            lines.append(f"{indent}{name}: LIST<{_scalar_name(value)}>")
    elif pa.types.is_map(dtype):
        key = _scalar_name(dtype.key_type)
        value = dtype.item_type
        if pa.types.is_struct(value):
            lines.append(f"{indent}{name}: MAP<{key}, STRUCT>")
            for child in value:
                _render(child, depth + 1, lines)
        else:
            lines.append(f"{indent}{name}: MAP<{key}, {_scalar_name(value)}>")
    else:
        members = projection.enum_members(field)
        values = f"  ({' | '.join(members)})" if members else ""
        lines.append(f"{indent}{name}: {_scalar_name(dtype)}{values}")


def describe(schema: pa.Schema = projection.SCHEMA) -> str:
    """Returns ``schema`` as an indented type tree.

    Args:
        schema: Schema to render; the projection schema by default.

    Returns:
        One line per field, two spaces of indentation per level of nesting, with
        container fields naming their shape and their children indented beneath, and
        each enum column followed by the values it may hold. The indentation gives the
        path to a leaf, not the syntax for reaching it: a child of a ``LIST`` or ``MAP``
        needs a list function rather than a dot.

    Raises:
        KeyError: If the schema holds a leaf type this module does not name.
    """
    lines: list[str] = []
    for field in schema:
        _render(field, 0, lines)
    return "\n".join(lines)
