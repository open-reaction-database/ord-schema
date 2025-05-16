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
import pytest
from sqlalchemy import select

from ord_schema.orm.database import delete_dataset, get_dataset_md5, get_dataset_size
from ord_schema.orm.mappers import Mappers
from ord_schema.proto import reaction_pb2


def test_orm(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionOutcome)
        .join(Mappers.ProductCompound)
        .join(Mappers.ProductMeasurement)
        .join(Mappers.Percentage)
        .where(Mappers.ProductMeasurement.type == "YIELD", Mappers.Percentage.value >= 70)
    )
    results = test_session.execute(query)
    reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]
    assert len(reactions) == 12


def test_delete_dataset(test_session):
    assert test_session.query(Mappers.Reaction).count() == 80
    delete_dataset("test_dataset", test_session)
    assert test_session.query(Mappers.Reaction).count() == 0


def test_get_dataset_md5(test_session):
    assert get_dataset_md5("test_dataset", test_session) == "0343d39a98d38eb39abd69d899af2bdf"
    assert get_dataset_md5("other_dataset", test_session) is None


def test_get_dataset_size(test_session):
    assert get_dataset_size("test_dataset", test_session) == 80
    with pytest.raises(ValueError):
        get_dataset_size("other_dataset", test_session)
