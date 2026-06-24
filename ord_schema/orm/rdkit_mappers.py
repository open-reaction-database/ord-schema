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
from enum import Enum
from typing import TYPE_CHECKING, Any

from sqlalchemy import Column, Index, Integer, Text, cast, func
from sqlalchemy.types import UserDefinedType

from ord_schema.orm import Base

if TYPE_CHECKING:
    from sqlalchemy.sql.expression import ColumnElement


def rdkit_cartridge() -> bool:
    """Returns whether to use RDKit PostgreSQL cartridge functionality."""
    return os.environ.get("ORD_POSTGRES_RDKIT", "1").lower() in (
        "1",
        "true",
        "yes",
        "y",
        "on",
    )


class RDKitMol(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L4."""

    cache_ok = True

    @property
    def python_type(self) -> type:
        """Raises NotImplementedError; this type has no Python equivalent."""
        raise NotImplementedError

    def get_col_spec(self, **kwargs: Any) -> str:
        """Returns the column type."""
        del kwargs  # Unused.
        return "mol" if rdkit_cartridge() else "bytea"


class RDKitReaction(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L129."""

    cache_ok = True

    @property
    def python_type(self) -> type:
        """Raises NotImplementedError; this type has no Python equivalent."""
        raise NotImplementedError

    def get_col_spec(self, **kwargs: Any) -> str:
        """Returns the column type."""
        del kwargs  # Unused.
        return "reaction" if rdkit_cartridge() else "bytea"


class RDKitBfp(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L81."""

    cache_ok = True

    @property
    def python_type(self) -> type:
        """Raises NotImplementedError; this type has no Python equivalent."""
        raise NotImplementedError

    def get_col_spec(self, **kwargs: Any) -> str:
        """Returns the column type."""
        del kwargs  # Unused.
        return "bfp" if rdkit_cartridge() else "bytea"


class RDKitSfp(UserDefinedType):
    """https://github.com/rdkit/rdkit/blob/master/Code/PgSQL/rdkit/rdkit.sql.in#L105."""

    cache_ok = True

    @property
    def python_type(self) -> type:
        """Raises NotImplementedError; this type has no Python equivalent."""
        raise NotImplementedError

    def get_col_spec(self, **kwargs: Any) -> str:
        """Returns the column type."""
        del kwargs  # Unused.
        return "sfp" if rdkit_cartridge() else "bytea"


class CString(UserDefinedType):
    """PostgreSQL cstring."""

    cache_ok = True

    @property
    def python_type(self) -> type:
        """Raises NotImplementedError; this type has no Python equivalent."""
        raise NotImplementedError

    def get_col_spec(self, **kwargs: Any) -> str:
        """Returns the column type."""
        del kwargs  # Unused.
        return "cstring"


class FingerprintType(Enum):
    """RDKit PostgreSQL fingerprint types."""

    # NOTE(skearnes): Add Column and Index entries for each member to RDKitMol below.
    MORGAN_BFP = func.morganbv_fp
    MORGAN_SFP = func.morgan_fp

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        """Invokes the wrapped RDKit fingerprint function."""
        return self.value(*args, **kwargs)


class RDKitMols(Base):
    """Table for storing compound structures and associated RDKit cartridge data."""

    __tablename__ = "mols"
    id = Column(Integer, primary_key=True)
    smiles = Column(Text, unique=True)
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
    def tanimoto(
        cls, other: str, fp_type: FingerprintType = FingerprintType.MORGAN_BFP
    ) -> ColumnElement[float]:
        """Returns the Tanimoto similarity value between the stored fingerprint and ``other``.

        This computes the similarity for every row (no index assistance) and is
        intended for selecting or ordering by similarity. To *filter* by similarity,
        use ``is_similar``, which uses the GiST fingerprint index.
        """
        return func.tanimoto_sml(
            getattr(cls, fp_type.name.lower()), fp_type(cast(other, RDKitMol))
        )

    @classmethod
    def is_similar(
        cls, other: str, fp_type: FingerprintType = FingerprintType.MORGAN_BFP
    ) -> ColumnElement[bool]:
        """Returns an expression testing whether the stored fingerprint is similar to ``other``.

        Uses the ``%`` operator, which is backed by the GiST fingerprint index. For
        ``MORGAN_BFP`` the cutoff is read from the ``rdkit.tanimoto_threshold`` session
        setting (default 0.5). For ``MORGAN_SFP`` (sparse fingerprints) the cartridge
        uses Dice similarity and reads ``rdkit.dice_threshold`` instead. Set the
        relevant GUC via ``set_config`` before executing.
        """
        return getattr(cls, fp_type.name.lower()).bool_op("%")(
            fp_type(cast(other, RDKitMol))
        )

    @classmethod
    def contains_substructure(cls, pattern: str) -> ColumnElement[bool]:
        """Returns an expression testing whether the stored mol contains the ``pattern`` substructure.

        Uses the ``@>`` operator, which is backed by the GiST mol index.
        """
        return cls.mol.bool_op("@>")(cast(pattern, RDKitMol))

    @classmethod
    def matches_smarts(cls, pattern: str) -> ColumnElement[bool]:
        """Returns an expression testing whether the stored mol matches the SMARTS ``pattern``.

        Uses the ``@>`` operator, which is backed by the GiST mol index.
        """
        return cls.mol.bool_op("@>")(func.qmol_from_smarts(cast(pattern, CString)))


class RDKitReactions(Base):
    """Table for storing reaction objects and associated RDKit cartridge data."""

    __tablename__ = "reactions"
    id = Column(Integer, primary_key=True)
    reaction_smiles = Column(Text, unique=True)
    reaction = Column(RDKitReaction)

    __table_args__ = (
        Index("reaction_index", "reaction", postgresql_using="gist"),
        {"schema": "rdkit"},
    )

    @classmethod
    def matches_smarts(cls, pattern: str) -> ColumnElement[bool]:
        """Returns an expression testing whether the stored reaction matches the reaction SMARTS ``pattern``.

        Uses the ``@>`` operator, which is backed by the GiST reaction index.
        """
        return cls.reaction.bool_op("@>")(
            func.reaction_from_smarts(cast(pattern, CString))
        )
