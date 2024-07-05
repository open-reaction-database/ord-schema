# Copyright 2024 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Test utilities."""
import re
from contextlib import contextmanager
from typing import Iterator

from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from testing.postgresql import Postgresql

from ord_schema.orm.database import prepare_database


@contextmanager
def get_test_engine(**kwargs) -> Iterator[tuple[Engine, bool]]:
    """Creates a test database with no initial data."""
    with Postgresql() as postgres:
        # See https://docs.sqlalchemy.org/en/20/dialects/postgresql.html#module-sqlalchemy.dialects.postgresql.psycopg.
        url = re.sub("postgresql://", "postgresql+psycopg://", postgres.url())
        engine = create_engine(url, **kwargs)
        rdkit_cartridge = prepare_database(engine)
        yield engine, rdkit_cartridge
