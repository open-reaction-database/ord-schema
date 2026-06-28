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

"""ORM table for the served Reaction proto payload, in the public schema.

The ord.* schema is the search index over reactions; the canonical serialized Reaction
proto that the API serves lives here, keyed by the public reaction_id.
"""

from sqlalchemy import Column, ForeignKey, LargeBinary, Text

from ord_schema.orm import Base


class ReactionProto(Base):
    """Serialized Reaction proto (the served payload), keyed by reaction_id."""

    __tablename__ = "reactions"
    reaction_id = Column(
        Text,
        ForeignKey("ord.reaction.reaction_id", ondelete="CASCADE"),
        primary_key=True,
    )
    proto = Column(LargeBinary, nullable=False)

    __table_args__ = ({"schema": "public"},)
