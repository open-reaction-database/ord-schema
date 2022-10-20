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

"""Tests for ord_schema.orm.database."""
from sqlalchemy import select

from ord_schema.orm.database import delete_dataset
from ord_schema.orm.mappers import Percentage, ProductCompound, ProductMeasurement, Reaction, ReactionOutcome
from ord_schema.proto import reaction_pb2


def test_orm(test_session):
    query = (
        select(Reaction)
        .join(ReactionOutcome)
        .join(ProductCompound)
        .join(ProductMeasurement)
        .join(Percentage)
        .where(ProductMeasurement.type == "YIELD", Percentage.value >= 70)
    )
    results = test_session.execute(query)
    reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]
    assert len(reactions) == 12


def test_delete_dataset(test_session):
    assert test_session.query(Reaction).count() == 80
    delete_dataset("test_dataset", test_session)
    assert test_session.query(Reaction).count() == 0
