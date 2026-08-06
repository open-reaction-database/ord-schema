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

"""A queryable projection of the whole Reaction proto.

Source datasets store reactions as opaque wire-format bytes, so any question about them
costs a full deserialization. A projection restates the *same* reactions as nested
Parquet -- message to ``STRUCT``, repeated field to ``LIST``, map to ``MAP`` -- so a
query engine can read one leaf without touching the rest.

**Nothing is dropped.** Every field of every message reachable from ``Reaction``
appears, because a field left out is a question nobody can ask, and the point of the
artifact is to support questions nobody enumerated in advance. The schema is generated
from the proto descriptors rather than written by hand, so a new field appears here
without anyone deciding it is worth carrying.

That works only because the message graph is acyclic: no message reaches itself, so the
recursion terminates. ``build_schema`` enforces this rather than assuming it -- a
recursive message added upstream would otherwise recurse forever at import time.

Two normalizations are applied, and only two. Both cost no query and remove a real trap:

* **United messages become one canonical float.** A ``{value, precision, units}``
  triple projected verbatim means ``WHERE temperature > 350`` silently misses every row
  recorded in Celsius. Each becomes a single column named for its unit --
  ``setpoint_kelvin`` -- so a question in other units stays expressible while the
  mixed-unit comparison that would quietly return the wrong rows does not.
* **Structural identifiers collapse to one canonical ``smiles``.** ``SMILES``,
  ``CXSMILES``, ``INCHI``, and ``MOLBLOCK`` all answer "what is this molecule," so the
  projection answers it once rather than making every consumer re-implement the
  preference order. They collapse only when that succeeds: an identifier RDKit cannot
  read stays in the list, since dropping it as well would report a compound as having no
  recorded structure when the source recorded one.

Every other identifier is kept, as a list. ``NAME`` alone covers compounds that no
structural identifier reaches, and a compound may carry several of them, so pivoting
identifiers into named scalar fields would silently drop data.

Anything else -- precision, the original units, the identifier types that were collapsed
-- remains in the source ``reaction`` column, which stays authoritative for byte-exact
round-tripping. The projection promises that every field is *readable*, not that it is
reproducible from the projection alone.

.. warning::
   Query repeated levels with list lambdas, not ``UNNEST`` in the ``FROM`` clause.
   ``UNNEST`` materializes the exploded intermediate and is 27-200x slower on identical
   answers; over the full corpus the difference is 0.9 seconds against not finishing.

   .. code-block:: sql

      -- fast
      WHERE len(list_filter(
              flatten(list_transform(map_values(inputs), i -> i.components)),
              c -> len(list_filter(c.identifiers,
                                   x -> x.type = 'NAME' AND x.value = 'THF')) > 0)) > 0

      -- same answer, does not finish
      FROM p, UNNEST(map_values(inputs)) t(i), UNNEST(i.components) u(c),
              UNNEST(c.identifiers) v(x)
      WHERE x.type = 'NAME' AND x.value = 'THF'
"""

import os
from typing import Any, cast

import pyarrow as pa
import pyarrow.parquet as pq
from google.protobuf.descriptor import Descriptor, FieldDescriptor
from google.protobuf.message import Message

import ord_schema
from ord_schema import artifacts, atomic_io, message_helpers, parquet, units
from ord_schema.proto import reaction_pb2

ARTIFACT = "projection"

# Canonical unit per united message type, and the suffix its column carries. Keyed by
# message name because the same target applies wherever the message appears. Only types
# reachable from Reaction are listed; Concentration is united but no field uses it.
_CANONICAL_UNITS: dict[str, tuple[str, str]] = {
    "Current": ("A", "amperes"),
    "FlowRate": ("mL/min", "milliliters_per_minute"),
    "Length": ("cm", "centimeters"),
    "Mass": ("g", "grams"),
    "Moles": ("mol", "moles"),
    "Pressure": ("kPa", "kilopascals"),
    "Temperature": ("K", "kelvin"),
    "Time": ("s", "seconds"),
    "Voltage": ("V", "volts"),
    "Volume": ("L", "liters"),
    "Wavelength": ("nm", "nanometers"),
}

# Messages that get a collapsed `smiles` field, and the identifier types it replaces.
_STRUCTURAL_COMPOUND_TYPES = frozenset(
    {
        reaction_pb2.CompoundIdentifier.SMILES,
        reaction_pb2.CompoundIdentifier.CXSMILES,
        reaction_pb2.CompoundIdentifier.INCHI,
        reaction_pb2.CompoundIdentifier.MOLBLOCK,
    }
)
_STRUCTURAL_REACTION_TYPES = frozenset(
    {
        reaction_pb2.ReactionIdentifier.REACTION_SMILES,
        reaction_pb2.ReactionIdentifier.REACTION_CXSMILES,
    }
)
_COMPOUND_MESSAGES = frozenset({"Compound", "ProductCompound"})

_ARROW_SCALARS: dict[int, pa.DataType] = {
    FieldDescriptor.TYPE_DOUBLE: pa.float64(),
    FieldDescriptor.TYPE_FLOAT: pa.float32(),
    FieldDescriptor.TYPE_INT64: pa.int64(),
    FieldDescriptor.TYPE_UINT64: pa.uint64(),
    FieldDescriptor.TYPE_INT32: pa.int32(),
    FieldDescriptor.TYPE_UINT32: pa.uint32(),
    FieldDescriptor.TYPE_BOOL: pa.bool_(),
    FieldDescriptor.TYPE_STRING: pa.string(),
    FieldDescriptor.TYPE_BYTES: pa.binary(),
    # Enums project as their value name rather than their number: readable in a query,
    # and stable across renumbering in a way the integer is not.
    FieldDescriptor.TYPE_ENUM: pa.string(),
}

_RESOLVER = units.UnitResolver()


def _canonical_unit(field: FieldDescriptor) -> tuple[str, str] | None:
    """Returns ``field``'s (target unit, column suffix), or None if it is not united."""
    if field.message_type is None:
        return None
    return _CANONICAL_UNITS.get(field.message_type.name)


def column_name(field: FieldDescriptor) -> str:
    """Returns the column name for ``field``, carrying the unit where one applies."""
    canonical = _canonical_unit(field)
    if canonical is not None:
        return f"{field.name}_{canonical[1]}"
    return field.name


def _field_type(field: FieldDescriptor, stack: frozenset[str]) -> pa.DataType:
    if _canonical_unit(field) is not None:
        return pa.float64()
    message_type = field.message_type
    if message_type is not None and message_type.GetOptions().map_entry:
        key_field = message_type.fields_by_name["key"]
        value_field = message_type.fields_by_name["value"]
        return pa.map_(_field_type(key_field, stack), _field_type(value_field, stack))
    if message_type is not None:
        inner = pa.struct(_struct_fields(message_type, stack))
    else:
        inner = _ARROW_SCALARS[field.type]
    if field.label == FieldDescriptor.LABEL_REPEATED:
        return pa.list_(inner)
    return inner


def _struct_fields(descriptor: Descriptor, stack: frozenset[str]) -> list[pa.Field]:
    if descriptor.full_name in stack:
        raise ValueError(
            f"{descriptor.full_name} reaches itself; a total projection of a recursive "
            "message does not terminate"
        )
    stack = stack | {descriptor.full_name}
    fields = []
    if descriptor.name in _COMPOUND_MESSAGES or descriptor.name == "Reaction":
        fields.append(pa.field("smiles", pa.string()))
    fields.extend(
        pa.field(column_name(field), _field_type(field, stack))
        for field in descriptor.fields
    )
    return fields


def build_schema() -> pa.Schema:
    """Returns the Arrow schema projected from the Reaction descriptor.

    Returns:
        One top-level field per Reaction field, plus the collapsed ``smiles``.

    Raises:
        ValueError: If the message graph contains a cycle, which a total projection
            cannot represent.
    """
    # The protobuf stubs type DESCRIPTOR as optional; on a generated class it never is.
    descriptor = cast(Descriptor, reaction_pb2.Reaction.DESCRIPTOR)
    return pa.schema(_struct_fields(descriptor, frozenset()))


SCHEMA = build_schema()


def _canonical_value(message: ord_schema.UnitMessage, target: str) -> float | None:
    """Returns a united message's value in ``target`` units, or None if unreadable."""
    if not message.HasField("value") or not message.units:
        return None
    try:
        return _RESOLVER.convert(message, target).value
    except (KeyError, ValueError):
        # A unit the resolver cannot convert reads as null rather than as a wrong
        # number in the wrong unit.
        return None


def _smiles(message: Message) -> str | None:
    """Returns the canonical SMILES for a compound or reaction, or None."""
    if cast(Descriptor, message.DESCRIPTOR).name == "Reaction":
        try:
            return (
                message_helpers.get_reaction_smiles(
                    cast(reaction_pb2.Reaction, message), generate_if_missing=True
                )
                or None
            )
        except ValueError:
            return None
    try:
        return (
            message_helpers.smiles_from_compound(cast(reaction_pb2.Compound, message))
            or None
        )
    except ValueError:
        return None


def _kept_identifiers(
    descriptor: Descriptor, values: Any, *, collapsed: bool
) -> list[Any]:
    """Returns identifiers the collapsed ``smiles`` column does not already carry.

    Args:
        descriptor: Descriptor of the message holding ``values``.
        values: The message's ``identifiers`` field.
        collapsed: Whether ``smiles`` holds the structure these identifiers record. A
            failed collapse keeps every identifier, so the structure the source did
            record stays readable.

    Returns:
        The identifiers to project.
    """
    if not collapsed:
        return list(values)
    structural = (
        _STRUCTURAL_COMPOUND_TYPES
        if descriptor.name in _COMPOUND_MESSAGES
        else _STRUCTURAL_REACTION_TYPES
    )
    return [value for value in values if value.type not in structural]


def _enum_name(field: FieldDescriptor, number: int) -> str:
    return field.enum_type.values_by_number[number].name


def message_row(message: Message) -> dict[str, Any]:
    """Projects ``message`` to a dict matching its struct type in ``SCHEMA``.

    Unset fields are None rather than the proto default, so a consumer can tell "the
    source is silent" from "the source says zero".

    Args:
        message: Any message reachable from Reaction.

    Returns:
        A dict keyed by projected column name.
    """
    descriptor = cast(Descriptor, message.DESCRIPTOR)
    row: dict[str, Any] = {}
    smiles = None
    if descriptor.name in _COMPOUND_MESSAGES or descriptor.name == "Reaction":
        smiles = _smiles(message)
        row["smiles"] = smiles
    for field in descriptor.fields:
        name = column_name(field)
        canonical = _canonical_unit(field)
        if canonical is not None:
            row[name] = _canonical_value(getattr(message, field.name), canonical[0])
            continue
        message_type = field.message_type
        if message_type is not None and message_type.GetOptions().map_entry:
            value_field = message_type.fields_by_name["value"]
            items = getattr(message, field.name).items()
            if value_field.message_type is not None:
                row[name] = [(key, message_row(value)) for key, value in items] or None
            else:
                row[name] = list(items) or None
            continue
        if field.label == FieldDescriptor.LABEL_REPEATED:
            values = getattr(message, field.name)
            if field.name == "identifiers":
                values = _kept_identifiers(
                    descriptor, values, collapsed=smiles is not None
                )
            if message_type is not None:
                row[name] = [message_row(value) for value in values] or None
            elif field.type == FieldDescriptor.TYPE_ENUM:
                row[name] = [_enum_name(field, value) for value in values] or None
            else:
                row[name] = list(values) or None
            continue
        if message_type is not None:
            row[name] = (
                message_row(getattr(message, field.name))
                if message.HasField(field.name)
                else None
            )
            continue
        try:
            present = message.HasField(field.name)
        except ValueError:
            # Proto3 scalars without explicit presence: absent and default are the same
            # wire state, so treat the default as absent.
            present = bool(getattr(message, field.name))
        if not present:
            row[name] = None
        elif field.type == FieldDescriptor.TYPE_ENUM:
            row[name] = _enum_name(field, getattr(message, field.name))
        else:
            row[name] = getattr(message, field.name)
    return row


def reaction_row(
    reaction: reaction_pb2.Reaction, reaction_id: str | None = None
) -> dict[str, Any]:
    """Projects a Reaction to a dict matching ``SCHEMA``.

    Args:
        reaction: Reaction to project.
        reaction_id: The source's ``reaction_id`` column value, if the caller read one.
            That column is what identifies the row in the source -- it is the key
            ``DatasetView.get_reaction`` looks up and the one a join has to match -- so
            it wins over a message that records no ID of its own.

    Returns:
        A dict keyed by projected column name.

    Raises:
        ValueError: If ``reaction_id`` and the message record different IDs, which
            leaves no safe answer to project.
    """
    row = message_row(reaction)
    if reaction_id is None or reaction_id == reaction.reaction_id:
        return row
    if reaction.reaction_id:
        raise ValueError(
            f"reaction_id column {reaction_id!r} disagrees with the reaction's own "
            f"{reaction.reaction_id!r}"
        )
    row["reaction_id"] = reaction_id
    return row


def is_current(path: str | os.PathLike[str], source_md5: str) -> bool:
    """Returns whether ``path`` is a projection of ``source_md5`` by this library."""
    return artifacts.is_current(path, ARTIFACT, source_md5)


def write_projection(
    source: parquet.DatasetView,
    output: str | os.PathLike[str],
    *,
    source_md5: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 1000,
) -> int:
    """Projects a source Parquet dataset and writes it.

    Reactions are read and written one source row group at a time, so peak memory is
    bounded by the largest row group rather than the dataset. The output is published
    atomically, so a failure partway leaves any existing projection untouched.

    Args:
        source: View of the source Parquet dataset.
        output: Path to write the projection to.
        source_md5: Source hash to stamp, if the caller already computed one. Hashed
            here when omitted, which costs a full pass over the source.
        compression: Parquet compression codec.
        row_group_size: Rows per output row group.

    Returns:
        Number of rows written.
    """
    if source_md5 is None:
        source_md5 = source.md5()
    stamps = artifacts.stamps(ARTIFACT, source.dataset_id or None, source_md5)
    schema = SCHEMA.with_metadata(artifacts.to_metadata(stamps))
    rows = 0
    with (
        atomic_io.atomic_path(output) as temp_path,
        pq.ParquetWriter(temp_path, schema, compression=compression) as writer,
    ):
        for row_group in range(source.num_row_groups):
            batch = [
                reaction_row(reaction, reaction_id)
                for reaction_id, reaction in source.iter_reactions(row_group=row_group)
            ]
            rows += len(batch)
            writer.write_table(
                pa.Table.from_pylist(batch, schema=schema),
                row_group_size=row_group_size,
            )
    return rows
