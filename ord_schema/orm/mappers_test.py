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

"""Tests for ord_schema.orm."""
import os
import pytest

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session
from testing.postgresql import Postgresql

from ord_schema.message_helpers import load_message
from ord_schema.orm.mappers import (
    Base,
    CompoundIdentifiers,
    Reaction,
    ReactionInputs,
    ReactionInputComponents,
    from_proto,
    to_proto,
)
from ord_schema.proto.dataset_pb2 import Dataset


@pytest.mark.parametrize(
    "filename",
    (
        os.path.join(os.path.dirname(__file__), "testdata", "empty.pbtxt"),
        os.path.join(os.path.dirname(__file__), "testdata", "full.pbtxt"),
        os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"),
    ),
)
def test_round_trip(filename):
    dataset = load_message(filename, Dataset)
    assert dataset == to_proto(from_proto(dataset))


def test_orm():
    dataset = load_message(os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"), Dataset)
    with Postgresql() as postgres:
        engine = create_engine(postgres.url())
        Base.metadata.create_all(engine)
        with Session(engine) as session:
            session.add(from_proto(dataset))
            # TODO(skearnes): Use with_polymorphic to simplify joins; e.g. with_polymorphic(Compound, "*").
            # https://docs.sqlalchemy.org/en/14/orm/inheritance_loading.html#with-polymorphic.
            query = (
                select(Reaction)
                .join(ReactionInputs)
                .join(ReactionInputComponents)
                .join(CompoundIdentifiers)
                .where(CompoundIdentifiers.type == "SMILES", CompoundIdentifiers.value == "c1ccccc1CCC(O)C")
            )
            for row in session.execute(query):
                print(row[0].reaction_id)
