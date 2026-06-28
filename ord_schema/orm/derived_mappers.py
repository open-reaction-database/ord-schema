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

These live in a separate "derived" schema (keyed by the canonical reaction_id) to
keep the proto-mirroring "ord" schema canonical. They are populated by post-passes
rather than from_proto, and the ORM functions normally when they are absent.
"""

from sqlalchemy import Column, ForeignKey, Integer, Text

from ord_schema.orm import Base


class ReactionClasses(Base):
    """Best-effort reaction classification, keyed by the reaction's integer id.

    A row's presence means classification was attempted for that reaction; NULL
    reaction_class/reaction_name mean Rxn-INSIGHT could not assign a value. Keying on
    ord.reaction.id -- the ORM's convention for reaction foreign keys, not the derived
    reaction_smiles -- keeps the post-pass idempotent: it classifies only reactions
    without a row here.
    """

    __tablename__ = "reaction_classes"
    # Integer FK to ord.reaction.id, per the ORM's reaction foreign-key convention
    # (see ord_schema.orm.mappers); the text ord.reaction.reaction_id is reserved for
    # cases that specifically need the ORD reaction ID.
    reaction_id = Column(
        Integer,
        ForeignKey("ord.reaction.id", ondelete="CASCADE"),
        primary_key=True,
    )
    # Coarse category (e.g. "C-C Coupling").
    reaction_class = Column(Text, index=True)
    # Specific named reaction within the class (e.g. "Suzuki coupling with boronic acids").
    reaction_name = Column(Text, index=True)

    __table_args__ = ({"schema": "derived"},)
