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
from sqlalchemy.orm import Session

from ord_schema.message_helpers import load_message
from ord_schema.orm.database import add_dataset, update_rdkit_ids, update_rdkit_tables
from ord_schema.orm.testing import get_test_engine
from ord_schema.proto import dataset_pb2


@pytest.fixture(name="test_session")
def test_session_fixture() -> Iterator[Session]:
    datasets = [
        load_message(
            os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"), dataset_pb2.Dataset
        )
    ]
    with get_test_engine() as (engine, rdkit_cartridge):
        with Session(engine) as session:
            for dataset in datasets:
                with session.begin():
                    add_dataset(dataset, session)
                if rdkit_cartridge:
                    with session.begin():
                        update_rdkit_tables(dataset.dataset_id, session)
                    with session.begin():
                        update_rdkit_ids(dataset.dataset_id, session)
        with Session(engine) as session:
            yield session
