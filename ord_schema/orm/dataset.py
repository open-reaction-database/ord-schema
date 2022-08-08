"""Table mappings for Dataset protos."""
from sqlalchemy import (
    Column,
    Integer,
    ForeignKey,
    String,
    Text,
)
from sqlalchemy.orm import relationship

from ord_schema.orm.reaction import Base

class Dataset(Base):
    name = Column(Text)
    description = Column(Text)
    reactions = relationship("Reaction")
    reaction_ids = relationship("ReactionId")
    dataset_id = Column(String(255))


class ReactionId(Base):
    dataset_id = Column(Integer, ForeignKey("dataset.id"), nullable=False)

    reaction_id = Column(String(255), ForeignKey("reaction.reaction_id"), nullable=False)


class DatasetExample(Base):
    dataset_id = Column(String(255), ForeignKey("dataset.dataset_id"), nullable=False)
    description = Column(Text)
    url = Column(Text)
    created = relationship("DatasetExampleCreated", uselist=False)
