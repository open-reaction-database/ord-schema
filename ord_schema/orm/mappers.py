# Copyright 2022 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Table mappings for Reaction protos.

Notes:
    * Foreign keys to the `reaction` table are done using the `id` column, not the ORD reaction ID (`reaction_id`).
      However, the `reaction_id` column is used when the reaction ID is specifically called for, as with crude inputs.
    * When a message type is used in multiple places, use Single Table Inheritance; see
      https://docs.sqlalchemy.org/en/14/orm/inheritance.html#single-table-inheritance.
    * The naming might be a bit confusing: classes for single-use and parent protos have the same name
      as the corresponding message, while classes for multi-use protos (child classes) are named as
      _<ContainingClass><Attribute>. For example, Reaction.identifiers uses the ReactionIdentifier class,
      while Reaction.inputs uses the _ReactionInputs class (since _ReactionInput is used in multiple places).
      The effect is that all tables storing proto information are named after their corresponding proto
      messages, with extra child tables for storing relationships.
    * Parent and child class names are prepended with underscores; e.g. _ReactionInput. Every parent class has an
      associated with_polymorphic wrapper for query construction, named without the leading underscore;
      e.g. ReactionInput.
    * Only message types are allowed for repeated/mapped values in the ORM (not scalar types). Specifically:
        * MassSpecMeasurementDetails.eic_masses is converted from repeated float to repeated FloatValue.
        * Dataset.reaction_ids is converted from repeated string to repeated ReactionId.
    * Source.id was renamed in the ORM to Source.vendor_id to avoid conflicting with the `id` column.
"""
from __future__ import annotations

from collections import defaultdict
from operator import attrgetter
from typing import Any, Mapping, Optional, Type

from google.protobuf.descriptor import Descriptor, FieldDescriptor
from google.protobuf.message import Message
from inflection import underscore
from sqlalchemy import Boolean, Column, Enum, Float, Integer, ForeignKey, LargeBinary, Text
from sqlalchemy.orm import relationship

import ord_schema.orm.structure  # pylint: disable=unused-import
from ord_schema import message_helpers
from ord_schema.orm import Base
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2




def get_message_type(full_name: str) -> Any:
    if not full_name.startswith("ord."):
        raise ValueError(f"Unexpected message type: {full_name}")
    op = attrgetter(full_name[4:])
    try:
        return op(reaction_pb2)
    except AttributeError:
        return op(dataset_pb2)


def get_parents(message_type: Type[Message]) -> dict[Type[Message], list[tuple[Type[Message], str, bool]]]:
    """Returns the parent message types for each message type."""
    parents = defaultdict(list)
    for child, parent, field_name, unique in _get_message_contexts(message_type.DESCRIPTOR, None, None, None):
        if parent is not None:
            parents[get_message_type(child)].append((get_message_type(parent), field_name, unique))
    parents[message_type] = []  # Add the base type.
    return parents


def _get_message_contexts(
    descriptor: Descriptor, parent: str | None, field_name: str | None, unique: bool | None
) -> set[tuple[str, str | None, str | None, bool | None]]:
    """Returns the set of contexts for each message type."""
    counts = {(descriptor.full_name, parent, field_name, unique)}
    for field in descriptor.fields:
        if field.type == FieldDescriptor.TYPE_MESSAGE:
            if set(field.message_type.fields_by_name.keys()) == {"key", "value"}:
                # Check for maps.
                field_message_type = field.message_type.fields_by_name["value"].message_type
            else:
                field_message_type = field.message_type
            counts |= _get_message_contexts(
                field_message_type,
                parent=descriptor.full_name,
                field_name=field.name,
                unique=field.label != FieldDescriptor.LABEL_REPEATED,
            )
    return counts


def build_mappers(message_type: Type[Message]) -> dict[Type[Message], Type]:
    mappers = {}
    parents = get_parents(message_type)
    for message_type in parents.keys():
        mappers[message_type] = _build_mapper(message_type, parents=parents)
    return mappers


_FIELD_TYPES = {
    FieldDescriptor.TYPE_BOOL: Boolean,
    FieldDescriptor.TYPE_BYTES: LargeBinary,
    FieldDescriptor.TYPE_DOUBLE: Float,
    FieldDescriptor.TYPE_FLOAT: Float,
    FieldDescriptor.TYPE_INT32: Integer,
    FieldDescriptor.TYPE_INT64: Integer,
    FieldDescriptor.TYPE_STRING: Text,
}


def _build_mapper(
    message_type: Type[Message], parents: dict[Type[Message], list[tuple[Type[Message], str, bool]]]
) -> Type:
    table_name = underscore(message_type.DESCRIPTOR.name)
    attrs = {"__tablename__": table_name, "id": Column(Integer, primary_key=True)}
    add_key = False
    for parent_type, field_name, unique in parents[message_type]:
        parent_table_name = underscore(parent_type.DESCRIPTOR.name)
        attrs[f"{parent_table_name}_id"] = Column(
            Integer,
            ForeignKey(f"{parent_table_name}.id", ondelete="CASCADE"),
            nullable=len(parents[message_type]) > 1,  # Use multiple foreign keys instead of inheritance.
            unique=unique,
        )
        if isinstance(getattr(parent_type(), field_name), Mapping):
            add_key = True
    if add_key:
        attrs["key"] = Column(Text)  # Map key.
    for field in message_type.DESCRIPTOR.fields:
        if field.name == "eic_masses":
            attrs[field.name] = relationship("_EicMass")
        elif field.name == "reaction_ids":
            attrs[field.name] = relationship("_ReactionId")
        elif field.type == FieldDescriptor.TYPE_MESSAGE:
            if set(field.message_type.fields_by_name.keys()) == {"key", "value"}:
                field_message_type = field.message_type.fields_by_name["value"].message_type
            else:
                field_message_type = field.message_type
            attrs[field.name] = relationship(
                field_message_type.name, uselist=field.label == FieldDescriptor.LABEL_REPEATED
            )
        elif field.type == FieldDescriptor.TYPE_ENUM:
            attrs[field.name] = Column(Enum(*field.enum_type.values_by_name.keys()), name=field.enum_type.name)
        else:
            attrs[field.name] = Column(_FIELD_TYPES[field.type])
    if message_type == dataset_pb2.Dataset:
        attrs["dataset_id"] = Column(Text, nullable=False, unique=True)
    if message_type == reaction_pb2.Reaction:
        attrs["proto"] = Column(LargeBinary, nullable=False)
        attrs["reaction_id"] = Column(Text, nullable=False, unique=True)
    if message_type in {reaction_pb2.Compound, reaction_pb2.ProductCompound}:
        attrs["structure"] = relationship("Structure", uselist=False)
    return type(message_type.DESCRIPTOR.name, (Base,), attrs)


class _ReactionId(Base):
    """Mapper for dataset_pb2.Dataset.reaction_ids."""

    __tablename__ = "reaction_id"
    id = Column(Integer, primary_key=True)
    dataset_id = Column(Integer, ForeignKey("dataset.id", ondelete="CASCADE"), nullable=False)
    reaction_id = Column(Text, ForeignKey("reaction.reaction_id", ondelete="CASCADE"), nullable=False)


class _EicMass(Base):
    """Mapper for reaction_pb2.ProductMeasurement.MassSpecMeasurement.eic_masses."""

    __tablename__ = "eic_mass"
    id = Column(Integer, primary_key=True)
    mass_spec_measurement_details_id = Column(
        Integer, ForeignKey("mass_spec_measurement_details.id", ondelete="CASCADE"), nullable=False
    )
    value = Column(Float)


MAPPERS: dict[Type[Message], Type] = build_mappers(dataset_pb2.Dataset)
MESSAGE_TYPES: dict[Type, Type[Message]] = {value: key for key, value in MAPPERS.items()}


def from_proto(  # pylint: disable=too-many-branches
    message: Message, mapper: Optional[Type[Base]] = None, key: Optional[str] = None
) -> Base:
    """Converts a protobuf message into an ORM object.

    Args:
        message: Protobuf message.
        mapper: ORM mapper class. For top-level protos like Dataset and Reaction this can be left as None; it must
            be provided for Child subclasses to properly handle polymorphism.
        key: Map key (we store maps as rows of (key, value) tuples).

    Returns:
        ORM object.
    """
    if mapper is None:
        mapper = MAPPERS[type(message)]
    kwargs = {}
    if key is not None:
        kwargs["key"] = key
    if type(message) == reaction_pb2.Reaction:
        kwargs["proto"] = message.SerializeToString()
    for field, value in message.ListFields():
        if field.name == "eic_masses":
            # Convert repeated float to repeated _EicMass.
            kwargs[field.name] = [_EicMass(value=v) for v in value]
        elif field.name == "reaction_ids":
            # Convert repeated string to repeated _ReactionId.
            kwargs[field.name] = [_ReactionId(reaction_id=v) for v in value]
        elif field.type == FieldDescriptor.TYPE_MESSAGE:
            field_mapper = getattr(mapper, field.name).mapper.class_
            if isinstance(value, Mapping):
                kwargs[field.name] = [from_proto(v, mapper=field_mapper, key=k) for k, v in value.items()]
            elif field.label == FieldDescriptor.LABEL_REPEATED:
                kwargs[field.name] = [from_proto(v, mapper=field_mapper) for v in value]
            else:
                kwargs[field.name] = from_proto(value, mapper=field_mapper)
        elif field.type == FieldDescriptor.TYPE_ENUM:
            kwargs[field.name] = field.enum_type.values_by_number[value].name
        else:
            kwargs[field.name] = value
    if isinstance(message, (reaction_pb2.Compound, reaction_pb2.ProductCompound)):
        # Add RDKit cartridge functionality.
        field_mapper = getattr(mapper, "structure").mapper.class_
        try:
            kwargs["structure"] = field_mapper(smiles=message_helpers.smiles_from_compound(message))
        except ValueError:
            pass
    return mapper(**kwargs)


def to_proto(base: Base) -> Message:
    """Converts an ORM object into a protobuf message.

    Args:
        base: ORM object.

    Returns:
        Protobuf message.
    """
    kwargs = {}
    proto = MESSAGE_TYPES[base]
    assert issubclass(proto, Message)
    for field in proto.DESCRIPTOR.fields:
        value = getattr(base, field.name)
        if field.name == "eic_masses":
            # Convert repeated FloatValue to repeated float.
            kwargs[field.name] = [v.value for v in value]
        elif field.name == "reaction_ids":
            # Convert repeated ReactionId to repeated string.
            kwargs[field.name] = [v.reaction_id for v in value]
        elif isinstance(value, list):
            if len(value) == 0:
                continue
            if hasattr(value[0], "key"):
                kwargs[field.name] = {v.key: to_proto(v) for v in value}
            else:
                kwargs[field.name] = [to_proto(v) for v in value]
        elif isinstance(value, Base):
            kwargs[field.name] = to_proto(value)
        else:
            kwargs[field.name] = value
    return proto(**kwargs)
