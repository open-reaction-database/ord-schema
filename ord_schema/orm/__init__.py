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

"""Base ORM objects."""
from inspect import getmro
from inflection import underscore
from sqlalchemy import Column, Integer, ForeignKey, String
from sqlalchemy.orm import declarative_base, declarative_mixin, declared_attr


class Base:
    """See https://docs.sqlalchemy.org/en/14/orm/declarative_mixins.html#augmenting-the-base.

    * The table name is a snake_case rendering of the class name.
    * Every table has an `id` column used as the primary key. Note that _Child tables do not have this column.
    """

    id = Column(Integer, primary_key=True)

    @declared_attr
    def __tablename__(cls):  # pylint: disable=no-self-argument
        name = cls.__name__  # pylint: disable=no-member
        if name.startswith("_"):
            name = name[1:]  # Remove leading underscore.
        return underscore(name)  # pylint: disable=no-member


Base = declarative_base(cls=Base)


@declarative_mixin
class Parent:
    """Mixin class for parent classes.

    * Polymorphic types are identified in the `_type` column.
    """

    _type = Column(String(255), nullable=False)

    @declared_attr
    def __mapper_args__(cls):  # pylint: disable=no-self-argument
        name = cls.__name__  # pylint: disable=no-member
        if name.startswith("_"):
            name = name[1:]  # Remove leading underscore.
        return {
            "polymorphic_identity": underscore(name),  # pylint: disable=no-member
            "polymorphic_on": cls._type,
        }


@declarative_mixin
class Child:
    """Mixin class for child classes.

    * The parent ID is stored in the `parent_id` column (also used as the primary key).
    """

    @declared_attr
    def parent_id(cls):  # pylint: disable=no-self-argument
        """Creates a `parent_id` column referring to the parent table."""
        parent_table = None
        for parent_class in getmro(cls)[1:]:
            if issubclass(parent_class, Base):
                parent_table = parent_class.__tablename__
                break
        assert parent_table is not None, (cls, getmro(cls))
        if parent_table in ["structure"]:
            foreign_key = ForeignKey(f"rdkit.{parent_table}.id")
        else:
            foreign_key = ForeignKey(f"{parent_table}.id")
        return Column(f"{parent_table}_id", Integer, foreign_key, primary_key=True)

    @declared_attr
    def __mapper_args__(cls):  # pylint: disable=no-self-argument
        return {"polymorphic_identity": underscore(cls.__name__)}  # pylint: disable=no-member
