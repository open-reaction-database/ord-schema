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
    * We use inheritance to handle messages that appear in more than one context; see
      https://docs.sqlalchemy.org/en/14/orm/inheritance.html. The possible constraints are:
        - Some message types are used more than once in a parent (such as Time in ReactionInput), which forces
          the use of either joined or single table inheritance.
        - Some messages can be repeated, while others are unique to their parent. These uniqueness constraints force
          the use of joined table inheritance since they are specific to their polymorphic type.
      For convenience and consistency, we use single table inheritance for *all* message types, regardless of whether
      they are used in one context or more than one context. This means that we do not enforce the second constraint in
      the database.
"""
from collections import defaultdict
from hashlib import md5
from operator import attrgetter
from typing import Any, Mapping, Optional, Type

from google.protobuf.descriptor import Descriptor, FieldDescriptor
from google.protobuf.message import Message
from inflection import underscore
from sqlalchemy import Boolean, Column, Enum, Float, ForeignKey, Integer, LargeBinary, String, Text
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import relationship

import ord_schema.orm.rdkit_mappers  # pylint: disable=unused-import
from ord_schema import message_helpers
from ord_schema.logging import get_logger
from ord_schema.orm import Base
from ord_schema.proto import dataset_pb2, reaction_pb2

logger = get_logger(__name__)


def get_message_type(full_name: str) -> Any:
    """Fetches the class for a protocol buffer message type."""
    if not full_name.startswith("ord."):
        raise ValueError(f"Unexpected message type: {full_name}")
    operator = attrgetter(full_name[4:])
    try:
        return operator(reaction_pb2)
    except AttributeError:
        return operator(dataset_pb2)


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
    if descriptor is None:
        raise ValueError((descriptor, parent, field_name, unique))
    counts = {(descriptor.full_name, parent, field_name, unique)}
    for field in descriptor.fields:
        if field.type == FieldDescriptor.TYPE_MESSAGE:
            if set(field.message_type.fields_by_name.keys()) == {"key", "value"}:
                # Check for maps.
                logger.debug(f"Found map: ({descriptor.full_name}, {field.name})")
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


def build_mappers() -> dict[Type[Message], Type]:
    """Creates ORM mapper classes for protocol buffer message types.

    Returns:
        Dict mapping protocol buffer message types to mapper classes.
    """
    logger.debug("Building ORM mappers")
    mappers = {}
    parents = get_parents(dataset_pb2.Dataset)
    for message_type in sorted(parents, key=lambda x: x.DESCRIPTOR.name):
        logger.debug(f"Building mapper for {message_type}")
        mappers[message_type] = build_mapper(message_type, parents=parents)
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


def build_mapper(  # pylint: disable=too-many-branches
    message_type: Type[Message], parents: dict[Type[Message], list[tuple[Type[Message], str, bool]]]
) -> Type:
    """Creates a mapper class for a specific protocol buffer message type.

    Args:
        message_type: Protocol buffer message type.
        parents: Dict mapping message types to lists of (parent message type, field name, unique) tuples.

    Returns:
        Generated mapper class.
    """
    attrs = {
        "__tablename__": underscore(message_type.DESCRIPTOR.name),
        "id": Column(Integer, primary_key=True),
        "ord_schema_context": Column(Text, nullable=False),
        "__table_args__": ({"schema": "ord"},),
    }
    attrs["__mapper_args__"] = {
        "polymorphic_on": attrs["ord_schema_context"],
        "polymorphic_identity": message_type.DESCRIPTOR.name,
        "with_polymorphic": "*",
    }
    add_key = False
    for parent_type, field_name, _ in parents[message_type]:
        if isinstance(getattr(parent_type(), field_name), Mapping):
            add_key = True
    if add_key:
        attrs["key"] = Column(Text)  # Map key.
    for field in message_type.DESCRIPTOR.fields:
        if field.name == "eic_masses":
            attrs[field.name] = Column(ARRAY(Float))
        elif field.name == "reaction_ids":
            attrs[field.name] = Column(ARRAY(Text))
        elif field.type == FieldDescriptor.TYPE_MESSAGE:
            kwargs = {}
            if field.label != FieldDescriptor.LABEL_REPEATED:
                kwargs["uselist"] = False
            # All relationships are to polymorphic child classes.
            child_class_name = f"_{message_type.DESCRIPTOR.name}{field.name.capitalize()}"
            attrs[field.name] = relationship(child_class_name, back_populates="parent", **kwargs)
        elif field.type == FieldDescriptor.TYPE_ENUM:
            attrs[field.name] = Column(Enum(*field.enum_type.values_by_name.keys(), name=field.enum_type.name))
        else:
            attrs[field.name] = Column(_FIELD_TYPES[field.type])
    if message_type == dataset_pb2.Dataset:
        # Make dataset IDs globally unique.
        attrs["dataset_id"] = Column(Text, nullable=False, unique=True)
        # Track the MD5 hash so we can quickly identify changes.
        attrs["md5"] = Column(String(32), nullable=False)
        # Track the number of reactions for quicker browsing.
        attrs["num_reactions"] = Column(Integer, nullable=False)
    elif message_type == reaction_pb2.Reaction:
        # Make reaction IDs globally unique.
        attrs["reaction_id"] = Column(Text, nullable=False, unique=True)
        # Serialize and store the entire Reaction proto.
        attrs["proto"] = Column(LargeBinary, nullable=False)
        attrs["reaction_smiles"] = Column(Text, index=True)
        attrs["rdkit_reaction_id"] = Column(Integer, ForeignKey("rdkit.reactions.id", ondelete="CASCADE"), index=True)
        attrs["rdkit_reaction"] = relationship("RDKitReactions")
    elif message_type in {reaction_pb2.Compound, reaction_pb2.ProductCompound}:
        attrs["smiles"] = Column(Text, index=True)
        attrs["rdkit_mol_id"] = Column(Integer, ForeignKey("rdkit.mols.id", ondelete="CASCADE"), index=True)
        attrs["rdkit_mol"] = relationship("RDKitMols")
    elif message_type in {reaction_pb2.CompoundPreparation, reaction_pb2.CrudeComponent}:
        # Add foreign key to reaction.reaction_id.
        kwargs = {}
        if message_type == reaction_pb2.CrudeComponent:
            kwargs["nullable"] = False
        attrs["reaction_id"] = Column(
            Text, ForeignKey("ord.reaction.reaction_id", ondelete="CASCADE"), index=True, **kwargs
        )
    logger.debug(f"Creating mapper {message_type.DESCRIPTOR.name}: {attrs}")
    mapper_class = type(message_type.DESCRIPTOR.name, (Base,), attrs)
    # Create polymorphic child classes.
    for parent_type, field_name, _ in parents[message_type]:
        foreign_table_name = underscore(parent_type.DESCRIPTOR.name)
        foreign_key = f"ord.{foreign_table_name}.id"
        child_attrs = {
            "__mapper_args__": {"polymorphic_identity": f"{parent_type.DESCRIPTOR.name}.{field_name}"},
            # Use get() to avoid column conflicts; see
            # https://docs.sqlalchemy.org/en/14/orm/inheritance.html#resolving-column-conflicts.
            #
            # NOTE(skearnes): We are not enforcing unique constraints on this column; see the module docstring.
            f"{foreign_table_name}_id": mapper_class.__table__.c.get(
                f"{foreign_table_name}_id",
                Column(Integer, ForeignKey(foreign_key, ondelete="CASCADE"), index=True),
            ),
            "parent": relationship(parent_type.DESCRIPTOR.name, back_populates=field_name, uselist=False),
        }
        child_class_name = f"_{parent_type.DESCRIPTOR.name}{field_name.capitalize()}"
        logger.debug(f"Creating child mapper {child_class_name}: {child_attrs}")
        type(child_class_name, (mapper_class,), child_attrs)
    return mapper_class


_MESSAGE_TO_MAPPER: dict[Type[Message], Type] = build_mappers()
_MAPPER_TO_MESSAGE: dict[Type, Type[Message]] = {value: key for key, value in _MESSAGE_TO_MAPPER.items()}


class _MappersMeta(type):
    """Metaclass for Mappers; see https://stackoverflow.com/a/3155493."""

    def __getattr__(cls, item):
        return cls._MAPPERS[item]

    def __getitem__(cls, item):
        return cls._MAPPERS[item]


class Mappers(metaclass=_MappersMeta):
    """Container for generated mapper classes."""

    _MAPPERS: dict[str, Type] = {c.__name__: c for c in _MAPPER_TO_MESSAGE}


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
        mapper = _MESSAGE_TO_MAPPER[type(message)]
    kwargs = {}
    if key is not None:
        kwargs["key"] = key
    for field, value in message.ListFields():
        if field.type == FieldDescriptor.TYPE_MESSAGE:
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
    if isinstance(message, dataset_pb2.Dataset):
        kwargs["md5"] = md5(message.SerializeToString(deterministic=True)).hexdigest()
        assert hasattr(message, "reactions") and hasattr(message, "reaction_ids")  # Type hints.
        kwargs["num_reactions"] = len(message.reactions) or len(message.reaction_ids)
    elif isinstance(message, reaction_pb2.Reaction):
        kwargs["proto"] = message.SerializeToString(deterministic=True)
        try:
            reaction_smiles = message_helpers.get_reaction_smiles(
                message, generate_if_missing=True, allow_incomplete=False, validate=True
            )
        except ValueError as error:
            assert hasattr(message, "reaction_id")  # Type hint.
            logger.debug(f"Error generating reaction SMILES for {message.reaction_id}: {error}")
            reaction_smiles = None
        if reaction_smiles is not None:
            kwargs["reaction_smiles"] = reaction_smiles.split()[0]  # Handle CXSMILES.
    elif isinstance(message, (reaction_pb2.Compound, reaction_pb2.ProductCompound)):
        try:
            kwargs["smiles"] = message_helpers.smiles_from_compound(message)
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
    try:
        proto = _MAPPER_TO_MESSAGE[type(base)]
    except KeyError:
        proto = _MAPPER_TO_MESSAGE[type(base).__bases__[0]]
    assert issubclass(proto, Message)
    for field in proto.DESCRIPTOR.fields:
        value = getattr(base, field.name)
        if isinstance(value, list) and len(value) == 0:
            continue
        if isinstance(value, list) and isinstance(value[0], Base):
            if hasattr(value[0], "key"):
                kwargs[field.name] = {v.key: to_proto(v) for v in value}
            else:
                kwargs[field.name] = [to_proto(v) for v in value]
        elif isinstance(value, Base):
            kwargs[field.name] = to_proto(value)
        else:
            kwargs[field.name] = value
    return proto(**kwargs)
