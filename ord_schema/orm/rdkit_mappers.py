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

"""RDKit PostgreSQL cartridge functionality.

Notes:
  * These tables live in a separate "rdkit" schema to avoid name conflicts between tables and extension types.
  * The RDKit-specific columns are populated by ord_schema.orm.database.add_rdkit; this allows the ORM to function
    normally even if if the RDKit PostgreSQL cartridge is not installed (the `smiles` column will be populated and
    the other columns will be empty).
  * Objects with this type are added to the ORM in from_proto() using the `rdkit` field.
"""
from __future__ import annotations

import os
from distutils.util import strtobool  # pylint: disable=deprecated-module
from enum import Enum

from sqlalchemy import Column, Index, Integer, Text, cast, func
from sqlalchemy.sql.expression import ColumnElement
from sqlalchemy.types import UserDefinedType

from ord_schema.orm import Base


def rdkit_cartridge() -> bool:
    """Returns whether to use RDKit PostgreSQL cartridge functionality."""
    return bool(strtobool(os.environ.get("ORD_POSTGRES_RDKIT", "1")))


class RDKitMol(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L4."""

    cache_ok = True

    @property
    def python_type(self):
        raise NotImplementedError

    def get_col_spec(self, **kwargs):
        """Returns the column type."""
        del kwargs  # Unused.
        return "mol" if rdkit_cartridge() else "bytea"


class RDKitReaction(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L129."""

    cache_ok = True

    @property
    def python_type(self):
        raise NotImplementedError

    def get_col_spec(self, **kwargs):
        """Returns the column type."""
        del kwargs  # Unused.
        return "reaction" if rdkit_cartridge() else "bytea"


class RDKitBfp(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L81."""

    cache_ok = True

    @property
    def python_type(self):
        raise NotImplementedError

    def get_col_spec(self, **kwargs):
        """Returns the column type."""
        del kwargs  # Unused.
        return "bfp" if rdkit_cartridge() else "bytea"


class RDKitSfp(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L105."""

    cache_ok = True

    @property
    def python_type(self):
        raise NotImplementedError

    def get_col_spec(self, **kwargs):
        """Returns the column type."""
        del kwargs  # Unused.
        return "sfp" if rdkit_cartridge() else "bytea"


class CString(UserDefinedType):
    """PostgreSQL cstring."""

    cache_ok = True

    @property
    def python_type(self):
        raise NotImplementedError

    def get_col_spec(self, **kwargs):
        """Returns the column type."""
        del kwargs  # Unused.
        return "cstring"


class FingerprintType(Enum):
    """RDKit PostgreSQL fingerprint types."""

    # NOTE(skearnes): Add Column and Index entries for each member to RDKitMol below.
    MORGAN_BFP = func.morganbv_fp
    MORGAN_SFP = func.morgan_fp

    def __call__(self, *args, **kwargs):
        return self.value(*args, **kwargs)


class RDKitMols(Base):
    """Table for storing compound structures and associated RDKit cartridge data."""

    __tablename__ = "mols"
    id = Column(Integer, primary_key=True)
    smiles = Column(Text, index=True, unique=True)
    mol = Column(RDKitMol)
    morgan_bfp = Column(RDKitBfp)
    morgan_sfp = Column(RDKitSfp)

    __table_args__ = (
        Index("mol_index", "mol", postgresql_using="gist"),
        Index("morgan_bfp_index", "morgan_bfp", postgresql_using="gist"),
        Index("morgan_sfp_index", "morgan_sfp", postgresql_using="gist"),
        {"schema": "rdkit"},
    )

    @classmethod
    def tanimoto(cls, other: str, fp_type: FingerprintType = FingerprintType.MORGAN_BFP) -> ColumnElement[float]:
        return func.tanimoto_sml(getattr(cls, fp_type.name.lower()), fp_type(cast(other, RDKitMol)))

    @classmethod
    def contains_substructure(cls, pattern: str) -> ColumnElement[bool]:
        return func.substruct(cls.mol, cast(pattern, RDKitMol))

    @classmethod
    def matches_smarts(cls, pattern: str) -> ColumnElement[bool]:
        return func.substruct(cls.mol, func.qmol_from_smarts(cast(pattern, CString)))


class RDKitReactions(Base):
    """Table for storing reaction objects and associated RDKit cartridge data."""

    __tablename__ = "reactions"
    id = Column(Integer, primary_key=True)
    reaction_smiles = Column(Text, index=True, unique=True)
    reaction = Column(RDKitReaction)

    __table_args__ = (
        Index("reaction_index", "reaction", postgresql_using="gist"),
        {"schema": "rdkit"},
    )

    @classmethod
    def matches_smarts(cls, pattern: str) -> ColumnElement[bool]:
        return func.substruct(cls.reaction, func.reaction_from_smarts(cast(pattern, CString)))
