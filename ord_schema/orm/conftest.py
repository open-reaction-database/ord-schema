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

"""Pytest fixtures."""
import os
import re
from typing import Iterator

import pytest
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session
from testing.postgresql import Postgresql

from ord_schema.message_helpers import load_message
from ord_schema.orm.database import add_dataset, prepare_database
from ord_schema.proto import dataset_pb2


@pytest.fixture(name="test_engine")
def test_engine_fixture() -> Iterator[Engine]:
    with Postgresql() as postgres:
        # See https://docs.sqlalchemy.org/en/20/dialects/postgresql.html#module-sqlalchemy.dialects.postgresql.psycopg.
        url = re.sub("postgresql://", "postgresql+psycopg://", postgres.url())
        engine = create_engine(url, future=True)
        yield engine


@pytest.fixture(name="test_session")
def test_session_fixture(test_engine) -> Iterator[Session]:
    datasets = [
        load_message(
            os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"), dataset_pb2.Dataset
        )
    ]
    rdkit_cartridge = prepare_database(test_engine)
    with Session(test_engine) as session:
        for dataset in datasets:
            with session.begin():
                add_dataset(dataset, session, rdkit_cartridge=rdkit_cartridge)
    with Session(test_engine) as session:
        yield session
