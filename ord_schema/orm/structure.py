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

"""RDKit PostgreSQL cartridge functionality."""
from __future__ import annotations

import os
from distutils.util import strtobool

from sqlalchemy import Column, Index, Integer, ForeignKey, Text, func
from sqlalchemy.orm import with_polymorphic
from sqlalchemy.types import UserDefinedType

from ord_schema.orm import Base, Child, Parent


def rdkit_cartridge() -> bool:
    """Returns whether to use RDKit PostgreSQL cartridge functionality."""
    return bool(strtobool(os.environ.get("ORD_POSTGRES_RDKIT", "1")))


class _RDKitMol(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L4."""

    cache_ok = True

    @property
    def python_type(self):
        raise NotImplementedError

    def get_col_spec(self, **kwargs):
        """Returns the column type."""
        del kwargs  # Unused.
        return "rdkit.mol" if rdkit_cartridge() else "bytea"


class _RDKitBfp(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L81."""

    @property
    def python_type(self):
        raise NotImplementedError

    def get_col_spec(self, **kwargs):
        """Returns the column type."""
        del kwargs  # Unused.
        return "rdkit.bfp" if rdkit_cartridge() else "bytea"


class _MorganBinaryFingerprint(_RDKitBfp):
    """rdkit.bfp subclass for Morgan binary fingerprints."""

    cache_ok = True

    class comparator_factory(_RDKitBfp.Comparator):  # pylint: disable=invalid-name # pytype: disable=attribute-error
        def __mod__(self, other):
            return self.bool_op("operator(rdkit.%)")(func.rdkit.morganbv_fp(other))

    @property
    def python_type(self):
        raise NotImplementedError


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


class _Structure(Parent, Base):
    """Table for storing compound structures and associated RDKit cartridge data.

    Notes:
      * This table lives in a separate "rdkit" schema to avoid name conflicts between tables and extension types.
      * The RDKit-specific columns are populated by ord_schema.orm.database.add_rdkit; this allows the ORM to function
        normally even if if the RDKit PostgreSQL cartridge is not installed (the `smiles` column will be populated and
        the other columns will be empty).
      * Objects with this type are added to the ORM in from_proto() using the `structure` field.
    """

    smiles = Column(Text)
    mol = Column(_RDKitMol)
    morgan_binary_fingerprint = Column(_MorganBinaryFingerprint)

    __table_args__ = (
        Index("mol_index", "mol", postgresql_using="gist"),
        Index("morgan_binary_fingerprint_index", "morgan_binary_fingerprint", postgresql_using="gist"),
        {"schema": "rdkit"},
    )


class _CompoundStructure(Child, _Structure):
    compound_id = Column(Integer, ForeignKey("compound.id"), nullable=False)

    __table_args__ = {"schema": "rdkit"}


class _ProductCompoundStructure(Child, _Structure):
    product_compound_id = Column(Integer, ForeignKey("product_compound.id"), nullable=False)

    __table_args__ = {"schema": "rdkit"}


Structure = with_polymorphic(_Structure, "*")
