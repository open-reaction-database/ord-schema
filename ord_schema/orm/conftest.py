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

import pathlib
import re
from collections.abc import Iterator

import pytest
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session
from testing.postgresql import Postgresql

from ord_schema.datasets import load_dataset
from ord_schema.orm.database import add_dataset, prepare_database, update_derived_data


@pytest.fixture(name="test_engine")
def test_engine_fixture() -> Iterator[Engine]:
    with Postgresql() as postgres:
        # See https://docs.sqlalchemy.org/en/20/dialects/postgresql.html#module-sqlalchemy.dialects.postgresql.psycopg.
        url = re.sub("postgresql://", "postgresql+psycopg://", postgres.url())
        engine = create_engine(url, future=True)
        yield engine


@pytest.fixture(name="prepared_engine")
def prepared_engine_fixture(test_engine: Engine) -> Engine:
    """``test_engine`` with the ORM schema (and RDKit cartridge) installed."""
    assert prepare_database(test_engine)
    return test_engine


@pytest.fixture(name="test_session")
def test_session_fixture(test_engine: Engine) -> Iterator[Session]:
    datasets = [
        load_dataset(
            pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt"
        )
    ]
    rdkit_cartridge = prepare_database(test_engine)
    with Session(test_engine) as session:
        for dataset in datasets:
            with session.begin():
                add_dataset(dataset, session)
                update_derived_data(
                    dataset.dataset_id, session, rdkit_cartridge=rdkit_cartridge
                )
    with Session(test_engine) as session:
        yield session
