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

from sqlalchemy import Column, ForeignKey, Integer, Text

from ord_schema.orm import Base


class ReactionClasses(Base):
    """Best-effort reaction classification, one row per classified reaction.

    A row's presence marks that classification was attempted; NULL reaction_class /
    reaction_name mean Rxn-INSIGHT could not assign a value. The post-pass classifies
    only reactions without a row, so it is idempotent.
    """

    __tablename__ = "reaction_classes"
    # Integer FK to ord.reaction.id, per the ORM's reaction FK convention (see mappers).
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
