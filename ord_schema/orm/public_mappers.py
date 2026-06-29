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

from sqlalchemy import Column, Date, ForeignKey, Integer, LargeBinary, String, Text

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


class DatasetMetadata(Base):
    """Per-dataset metadata for the API, keyed by dataset_id.

    md5 hashes the original Dataset proto (computed at ingest, not reconstructed) and is
    used to skip re-importing unchanged datasets. num_reactions and submitted_at support
    dataset browsing; submitted_at is filled by database.set_submitted_at.
    """

    __tablename__ = "datasets"
    dataset_id = Column(
        Text,
        ForeignKey("ord.dataset.dataset_id", ondelete="CASCADE"),
        primary_key=True,
    )
    md5 = Column(String(32), nullable=False)
    num_reactions = Column(Integer, nullable=False)
    submitted_at = Column(Date, index=True)

    __table_args__ = ({"schema": "public"},)
