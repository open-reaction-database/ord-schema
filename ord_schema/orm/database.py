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
import time
from unittest.mock import patch

from sqlalchemy import cast, delete, func, select, text, update
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.engine import Engine
from sqlalchemy.exc import NotSupportedError, OperationalError
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger
from ord_schema.orm.mappers import Base, Mappers, from_proto
from ord_schema.orm.rdkit_mappers import CString, FingerprintType, RDKitMol, RDKitReaction
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def get_connection_string(
    database: str, username: str, password: str, host: str = "localhost", port: int = 5432
) -> str:
    """Creates an SQLAlchemy connection string."""
    return f"postgresql://{username}:{password}@{host}:{port}/{database}?client_encoding=utf-8"


def prepare_database(engine: Engine) -> bool:
    """Prepares the database and creates the ORM table structure.

    Args:
        engine: SQLAlchemy Engine.

    Returns:
        Whether the RDKit PostgreSQL cartridge is installed.
    """
    with engine.begin() as connection:
        try:
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS tsm_system_rows"))  # For random sampling.
        except OperationalError:
            logger.warning("tsm_system_rows cartridge is not installed; random sampling will be disabled")
    with engine.begin() as connection:
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS ord"))
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS rdkit"))
    try:
        with engine.begin() as connection:
            # NOTE(skearnes): The RDKit PostgreSQL extension works best in the public schema.
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS rdkit"))
        rdkit_cartridge = True
    except (OperationalError, NotSupportedError):
        with engine.begin() as connection:
            logger.warning("RDKit PostgreSQL cartridge is not installed; structure search will be disabled")
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS btree_gist"))
        rdkit_cartridge = False
    with patch.dict(os.environ, {"ORD_POSTGRES_RDKIT": "1" if rdkit_cartridge else "0"}):
        Base.metadata.create_all(engine)
    return rdkit_cartridge


def add_dataset(dataset: dataset_pb2.Dataset, session: Session, rdkit_cartridge: bool = True) -> None:
    """Adds a dataset to the database."""
    logger.info(f"Adding dataset {dataset.dataset_id}")
    start = time.time()
    mapped_dataset = from_proto(dataset)
    logger.info(f"from_proto() took {time.time() - start:g}s")
    start = time.time()
    session.add(mapped_dataset)
    logger.info(f"session.add() took {time.time() - start:g}s")
    if rdkit_cartridge:
        session.flush()
        update_rdkit_tables(dataset.dataset_id, session)
        session.flush()
        update_rdkit_ids(dataset.dataset_id, session)


def get_dataset_md5(dataset_id: str, session: Session) -> str | None:
    """Returns the MD5 hash of the current version of a dataset, if it exists in the database."""
    result = session.execute(select(Mappers.Dataset.md5).where(Mappers.Dataset.dataset_id == dataset_id))
    row = result.first()
    return row[0] if row else None


def get_dataset_size(dataset_id: str, session: Session) -> int:
    """Returns the number of reactions in a dataset."""
    result = session.execute(select(Mappers.Dataset.num_reactions).where(Mappers.Dataset.dataset_id == dataset_id))
    row = result.first()
    if row is None:
        raise ValueError(dataset_id)
    return row[0]


def delete_dataset(dataset_id: str, session: Session) -> None:
    """Deletes a dataset from the database."""
    logger.info(f"Deleting dataset {dataset_id}")
    start = time.time()
    session.execute(delete(Mappers.Dataset).where(Mappers.Dataset.dataset_id == dataset_id))
    logger.info(f"delete took {time.time() - start}s")


def update_rdkit_tables(dataset_id: str, session: Session) -> None:
    """Updates RDKit PostgreSQL cartridge data."""
    _update_rdkit_reactions(dataset_id, session)
    _update_rdkit_mols(dataset_id, session)


def _update_rdkit_reactions(dataset_id: str, session: Session) -> None:
    """Updates the RDKit reactions table."""
    logger.info("Updating RDKit reactions")
    assert hasattr(RDKitReaction, "__table__")  # Type hint.
    table = RDKitReaction.__table__
    start = time.time()
    session.execute(
        insert(table)
        .from_select(
            ["reaction_smiles"],
            select(Mappers.Reaction.reaction_smiles)
            .join(Mappers.Dataset)
            .where(Mappers.Dataset.dataset_id == dataset_id, Mappers.Reaction.reaction_smiles.is_not(None))
            .distinct(),
        )
        .on_conflict_do_nothing(index_elements=["reaction_smiles"])
    )
    session.execute(
        update(table)
        .where(table.c.reaction.is_(None))
        .values(reaction=func.reaction_from_smiles(cast(table.c.reaction_smiles, CString)))
    )
    logger.info(f"Updating reactions took {time.time() - start:g}s")


def _update_rdkit_mols(dataset_id: str, session: Session) -> None:
    """Updates the RDKit mols table."""
    logger.info("Updating RDKit mols")
    assert hasattr(RDKitMol, "__table__")  # Type hint.
    table = RDKitMol.__table__
    start = time.time()
    # NOTE(skearnes): This join path will not include non-input compounds like workups, internal standards, etc.
    session.execute(
        insert(table)
        .from_select(
            ["smiles"],
            select(Mappers.Compound.smiles)
            .join(Mappers.ReactionInput)
            .join(Mappers.Reaction)
            .join(Mappers.Dataset)
            .where(
                Mappers.Dataset.dataset_id == dataset_id,
                Mappers.Compound.smiles.is_not(None),
                # See https://github.com/open-reaction-database/ord-schema/issues/672.
                Mappers.Compound.smiles.not_like("%[Ti+5]%"),
            )
            .distinct(),
        )
        .on_conflict_do_nothing(index_elements=["smiles"])
    )
    session.execute(
        insert(table)
        .from_select(
            ["smiles"],
            select(Mappers.ProductCompound.smiles)
            .join(Mappers.ReactionOutcome)
            .join(Mappers.Reaction)
            .join(Mappers.Dataset)
            .where(Mappers.Dataset.dataset_id == dataset_id, Mappers.ProductCompound.smiles.is_not(None))
            .distinct(),
        )
        .on_conflict_do_nothing(index_elements=["smiles"])
    )
    session.execute(
        update(table).where(table.c.mol.is_(None)).values(mol=func.mol_from_smiles(cast(table.c.smiles, CString)))
    )
    logger.info(f"Updating mols took {time.time() - start:g}s")
    logger.info("Updating fingerprints")
    for fp_type in FingerprintType:
        start = time.time()
        column = fp_type.name.lower()
        session.execute(
            update(table)
            .where(getattr(table.c, column).is_(None), table.c.mol.is_not(None))
            .values(**{column: fp_type(table.c.mol)})
        )
        logger.info(f"Updating {fp_type} took {time.time() - start:g}s")


def update_rdkit_ids(dataset_id: str, session: Session) -> None:
    """Updates RDKit reaction and mol ID associations in the ORD tables."""
    logger.info("Updating RDKit ID associations")
    start = time.time()
    # Update Reaction.
    query = session.execute(
        select(Mappers.Reaction.id, RDKitReaction.id)
        .join(RDKitReaction, Mappers.Reaction.reaction_smiles == RDKitReaction.reaction_smiles)
        .join(Mappers.Dataset)
        .where(Mappers.Dataset.dataset_id == dataset_id)
    )
    updates = []
    for ord_id, rdkit_id in query.fetchall():
        updates.append({"id": ord_id, "rdkit_reaction_id": rdkit_id})
    session.execute(update(Mappers.Reaction), updates)
    # Update Compound.
    query = session.execute(
        select(Mappers.Compound.id, RDKitMol.id)
        .join(RDKitMol, Mappers.Compound.smiles == RDKitMol.smiles)
        .join(Mappers.ReactionInput)
        .join(Mappers.Reaction)
        .join(Mappers.Dataset)
        .where(Mappers.Dataset.dataset_id == dataset_id)
    )
    updates = []
    for ord_id, rdkit_id in query.fetchall():
        updates.append({"id": ord_id, "rdkit_mol_id": rdkit_id})
    session.execute(update(Mappers.Compound), updates)
    # Update ProductCompound.
    query = session.execute(
        select(Mappers.ProductCompound.id, RDKitMol.id)
        .join(RDKitMol, Mappers.ProductCompound.smiles == RDKitMol.smiles)
        .join(Mappers.ReactionOutcome)
        .join(Mappers.Reaction)
        .join(Mappers.Dataset)
        .where(Mappers.Dataset.dataset_id == dataset_id)
    )
    updates = []
    for ord_id, rdkit_id in query.fetchall():
        updates.append({"id": ord_id, "rdkit_mol_id": rdkit_id})
    session.execute(update(Mappers.ProductCompound), updates)
    logger.info(f"Updating RDKit IDs took {time.time() - start:g}s")
