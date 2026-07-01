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

"""ORM tables for derived, best-effort data that is not part of the proto.

These live in a separate "derived" schema to keep the proto-mirroring "ord" schema
canonical. Populated by post-passes (not from_proto); the ORM works without them.
"""

from sqlalchemy import Column, ForeignKey, Index, Integer, Text, Uuid, text
from sqlalchemy.orm import relationship

from ord_schema.orm import Base


class ReactionSmiles(Base):
    """Generated reaction SMILES and its link to the deduplicated RDKit reaction.

    One row per reaction: reaction_smiles is filled by update_derived_tables, rdkit_reaction_id
    by the RDKit linking pass (rdkit.reactions is deduplicated by SMILES).
    """

    __tablename__ = "reaction_smiles"
    # FK to the ord.reaction UUIDv7 surrogate key (see mappers).
    reaction_id = Column(
        Uuid, ForeignKey("ord.reaction.id", ondelete="CASCADE"), primary_key=True
    )
    reaction_smiles = Column(Text, index=True)
    rdkit_reaction_id = Column(
        Integer, ForeignKey("rdkit.reactions.id", ondelete="CASCADE"), index=True
    )
    rdkit_reaction = relationship("RDKitReactions")

    __table_args__ = (
        # Partial index on the RDKit join key over not-yet-linked rows; the RDKit linking
        # pass joins these to rdkit.reactions by reaction_smiles and filters unlinked, so
        # this stays tiny once most rows are linked. See ord_schema/orm/README.md.
        Index(
            "reaction_smiles_unlinked_index",
            "reaction_smiles",
            postgresql_where=text("rdkit_reaction_id IS NULL"),
        ),
        {"schema": "derived"},
    )


class CompoundSmiles(Base):
    """Generated compound SMILES and its link to the deduplicated RDKit mol.

    One row per ord.compound: smiles is filled by update_derived_tables, rdkit_mol_id by
    the RDKit linking pass.
    """

    __tablename__ = "compound_smiles"
    compound_id = Column(
        Uuid, ForeignKey("ord.compound.id", ondelete="CASCADE"), primary_key=True
    )
    smiles = Column(Text, index=True)
    rdkit_mol_id = Column(
        Integer, ForeignKey("rdkit.mols.id", ondelete="CASCADE"), index=True
    )
    rdkit_mol = relationship("RDKitMols")

    __table_args__ = (
        # Partial index on the RDKit join key (smiles) over not-yet-linked rows; see README.
        Index(
            "compound_smiles_unlinked_index",
            "smiles",
            postgresql_where=text("rdkit_mol_id IS NULL"),
        ),
        {"schema": "derived"},
    )


class ProductCompoundSmiles(Base):
    """Generated product-compound SMILES and its link to the deduplicated RDKit mol.

    One row per ord.product_compound; see CompoundSmiles.
    """

    __tablename__ = "product_compound_smiles"
    product_compound_id = Column(
        Uuid,
        ForeignKey("ord.product_compound.id", ondelete="CASCADE"),
        primary_key=True,
    )
    smiles = Column(Text, index=True)
    rdkit_mol_id = Column(
        Integer, ForeignKey("rdkit.mols.id", ondelete="CASCADE"), index=True
    )
    rdkit_mol = relationship("RDKitMols")

    __table_args__ = (
        # Partial index on the RDKit join key (smiles) over not-yet-linked rows; see README.
        Index(
            "product_compound_smiles_unlinked_index",
            "smiles",
            postgresql_where=text("rdkit_mol_id IS NULL"),
        ),
        {"schema": "derived"},
    )


class ReactionClasses(Base):
    """Best-effort reaction classification, one row per classified reaction.

    A row's presence marks that classification was attempted; NULL reaction_class /
    reaction_name mean Rxn-INSIGHT could not assign a value. The post-pass classifies
    only reactions without a row, so it is idempotent.
    """

    __tablename__ = "reaction_classes"
    # FK to the ord.reaction UUIDv7 surrogate key (see mappers).
    reaction_id = Column(
        Uuid,
        ForeignKey("ord.reaction.id", ondelete="CASCADE"),
        primary_key=True,
    )
    # Coarse category (e.g. "C-C Coupling").
    reaction_class = Column(Text, index=True)
    # Specific named reaction within the class (e.g. "Suzuki coupling with boronic acids").
    reaction_name = Column(Text, index=True)

    __table_args__ = ({"schema": "derived"},)
