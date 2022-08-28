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

"""Functions for creating/managing the PostgreSQL database."""
import os
from unittest.mock import patch

from sqlalchemy import cast, func, text, update
from sqlalchemy.engine import Engine
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import Session

from ord_schema.orm import mappers
from ord_schema.proto import dataset_pb2


def get_connection_string(
    database: str, username: str, password: str, host: str = "localhost", port: int = 5432
) -> str:
    """Creates an SQLAlchemy connection string."""
    return f"postgresql://{username}:{password}@{host}:{port}/{database}"


def prepare_database(engine: Engine) -> bool:
    """Prepares the database and creates the ORM table structure.

    Args:
        engine: SQLAlchemy Engine.

    Returns:
        Whether the RDKit PostgreSQL cartridge is installed.
    """
    with engine.begin() as connection:
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS rdkit"))
        try:
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA rdkit"))
            rdkit_cartridge = True
        except OperationalError:
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS btree_gist WITH SCHEMA rdkit"))
            rdkit_cartridge = False
    with patch.dict(os.environ, {"ORD_POSTGRES_RDKIT": "1" if rdkit_cartridge else "0"}):
        mappers.Base.metadata.create_all(engine)
    return rdkit_cartridge


def add_datasets(datasets: list[dataset_pb2.Dataset], engine: Engine) -> None:
    """Adds datasets to the database."""
    with Session(engine) as session:
        for dataset in datasets:
            session.add(mappers.from_proto(dataset))
        session.commit()


def add_rdkit(engine: Engine) -> None:
    """Adds RDKit PostgreSQL cartridge data."""
    table = mappers.Structure.__table__
    with Session(engine) as session:
        session.execute(update(table).values(mol=func.rdkit.mol_from_smiles(cast(table.c.smiles, mappers.CString))))
        session.execute(update(table).values(morgan_binary_fingerprint=func.rdkit.morganbv_fp(table.c.mol)))
        session.commit()
