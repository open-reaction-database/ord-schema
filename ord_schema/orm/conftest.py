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
from typing import Iterator

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from testing.postgresql import Postgresql

from ord_schema.message_helpers import load_message
from ord_schema.orm.database import add_dataset, add_rdkit, prepare_database
from ord_schema.proto.dataset_pb2 import Dataset


@pytest.fixture
def test_session() -> Iterator[Session]:
    dataset = load_message(os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"), Dataset)
    with Postgresql() as postgres:
        engine = create_engine(postgres.url(), future=True)
        rdkit_cartridge = prepare_database(engine)
        with Session(engine) as session:
            add_dataset(dataset, session)
            session.flush()
            if rdkit_cartridge:
                add_rdkit(session)
            session.commit()
        with Session(engine) as session:
            yield session
