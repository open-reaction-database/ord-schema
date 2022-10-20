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
import time
import os
from unittest.mock import patch

from sqlalchemy import cast, delete, func, text, update
from sqlalchemy.engine import Engine
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger
from ord_schema.orm.mappers import Base, Dataset, from_proto
from ord_schema.orm.structure import CString, FingerprintType, Structure
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


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
            logger.warning("RDKit PostgreSQL cartridge is not installed; structure search will be disabled")
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS btree_gist WITH SCHEMA rdkit"))
            rdkit_cartridge = False
    with patch.dict(os.environ, {"ORD_POSTGRES_RDKIT": "1" if rdkit_cartridge else "0"}):
        Base.metadata.create_all(engine)
    return rdkit_cartridge


def add_dataset(dataset: dataset_pb2.Dataset, session: Session) -> None:
    """Adds a dataset to the database."""
    logger.info(f"Adding dataset {dataset.dataset_id}")
    start = time.time()
    mapped_dataset = from_proto(dataset)
    logger.info(f"from_proto() took {time.time() - start}s")
    start = time.time()
    session.add(mapped_dataset)
    logger.info(f"session.add() took {time.time() - start}s")


def delete_dataset(dataset_id: str, session: Session) -> None:
    """Deletes a dataset from the database."""
    logger.info(f"Deleting dataset {dataset_id}")
    start = time.time()
    session.execute(delete(Dataset).where(Dataset.dataset_id == dataset_id))
    logger.info(f"delete took {time.time() - start}s")


def add_rdkit(session: Session) -> None:
    """Adds RDKit PostgreSQL cartridge data to any null-valued cells in the structure table."""
    logger.info("Populating RDKit columns")
    assert hasattr(Structure, "__table__")  # Type hint.
    table = Structure.__table__
    start = time.time()
    session.execute(
        update(table).where(table.c.mol.is_(None)).values(mol=func.rdkit.mol_from_smiles(cast(table.c.smiles, CString)))
    )
    logger.info(f"Adding mol took {time.time() - start}s")
    for fp_type in FingerprintType:
        start = time.time()
        column = fp_type.name.lower()
        session.execute(
            update(table).where(getattr(table.c, column).is_(None)).values(**{column: fp_type(table.c.mol)})
        )
        logger.info(f"Adding {fp_type} took {time.time() - start}s")
