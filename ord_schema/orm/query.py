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

"""ORM queries."""
from typing import Type

from sqlalchemy import inspect
from sqlalchemy.orm.attributes import InstrumentedAttribute
from sqlalchemy.orm.relationships import RelationshipProperty
from sqlalchemy.orm.util import AliasedClass

from ord_schema.logging import get_logger
from ord_schema.orm.base import Base
from ord_schema.orm.mappers import POLYMORPHIC_MAPPERS

logger = get_logger(__name__)


def get_polymorphic_mappers() -> dict[Type[Base], AliasedClass]:
    """Returns a dict mapping child mapper classes to parent (polymorphic) mapper classes."""
    mappers = {}
    for mapper in POLYMORPHIC_MAPPERS:
        mappers.update({m.class_: mapper for m in inspect(mapper).with_polymorphic_mappers})
    return mappers


_POLYMORPHIC_MAPPERS = get_polymorphic_mappers()


def get_join_path(
    source: Type[Base],
    target: InstrumentedAttribute,
    raise_on_multiple: bool = True,
    raise_on_failure: bool = True,
    include_source: bool = True,
) -> list[Type[Base] | AliasedClass] | None:
    """Uses depth-first search to find a join path between a source mapper and a target attribute.

    Args:
        source: Source mapper; the object type that you want to return from the SQL query.
        target: Target attribute; the attribute you want to filter on.
        raise_on_multiple: If True, raise if a found path may not be unique.
        raise_on_failure: If True, raise if a path cannot be found.
        include_source: If True, include the source in the list of returned classes.

    Returns:
        List of mapper classes in the order they should be joined.

    Raises:
        ValueError: If raise_on_failure is True and a path cannot be found.
    """
    for value in source.__dict__.values():
        if not isinstance(value, InstrumentedAttribute):
            continue
        if value.property == target.property:
            if isinstance(value.property, RelationshipProperty):
                mapper = value.property.mapper.class_
                mapper = _POLYMORPHIC_MAPPERS.get(mapper, mapper)
                return [mapper]
            return []
        if not isinstance(value.property, RelationshipProperty):
            continue
        mapper = value.property.mapper.class_
        child_classes = get_join_path(
            source=mapper,
            target=target,
            raise_on_multiple=raise_on_multiple,
            raise_on_failure=False,
            include_source=False,
        )
        if child_classes is not None:
            classes = [source] if include_source else []
            if mapper in _POLYMORPHIC_MAPPERS:
                if raise_on_multiple:
                    raise ValueError(
                        f"path between {source} and {target} may not be unique; "
                        "consider making the target more precise"
                    )
                logger.warning(f"path between {source} and {target} may not be unique")
                mapper = _POLYMORPHIC_MAPPERS[mapper]
            classes.append(mapper)
            classes.extend(child_classes)
            return classes
    if raise_on_failure:
        raise ValueError(f"could not find a path between {source} and {target}")
    return None
