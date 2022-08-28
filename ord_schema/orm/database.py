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
from sqlalchemy import cast, func, text, update
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session

from ord_schema.orm.mappers import Base, CString, Structure, from_proto
from ord_schema.proto import dataset_pb2


def get_connection_string(
    database: str, username: str, password: str, host: str = "localhost", port: int = 5432
) -> str:
    """Creates an SQLAlchemy connection string."""
    return f"postgresql://{username}:{password}@{host}:{port}/{database}"


def prepare_database(engine: Engine) -> None:
    """Prepares the database and creates the ORM table structure."""
    with engine.begin() as connection:
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS rdkit"))
        connection.execute(text("CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA rdkit"))
    Base.metadata.create_all(engine)


def add_datasets(datasets: list[dataset_pb2.Dataset], engine: Engine) -> None:
    """Adds datasets to the database."""
    with Session(engine) as session:
        for dataset in datasets:
            session.add(from_proto(dataset))
        session.commit()


def add_rdkit(engine: Engine) -> None:
    """Adds RDKit PostgreSQL cartridge data."""
    table = Structure.__table__
    with Session(engine) as session:
        session.execute(update(table).values(mol=func.rdkit.mol_from_smiles(cast(Structure.smiles, CString))))
        session.execute(update(table).values(morgan_binary_fingerprint=func.rdkit.morganbv_fp(Structure.mol)))
        session.commit()
