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

from sqlalchemy import delete, select, text
from sqlalchemy.engine import Engine
from sqlalchemy.exc import NotSupportedError, OperationalError
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger
from ord_schema.orm.mappers import Base, Mappers, from_proto
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def get_connection_string(
    database: str, username: str, password: str, host: str = "localhost", port: int = 5432
) -> str:
    """Creates an SQLAlchemy connection string."""
    return f"postgresql+psycopg://{username}:{password}@{host}:{port}/{database}?client_encoding=utf8"


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
    logger.debug(f"Adding dataset {dataset.dataset_id}")
    start = time.time()
    mapped_dataset = from_proto(dataset)
    logger.debug(f"from_proto() took {time.time() - start:g}s")
    session.add(mapped_dataset)
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
    logger.debug(f"Deleting dataset {dataset_id}")
    start = time.time()
    session.execute(delete(Mappers.Dataset).where(Mappers.Dataset.dataset_id == dataset_id))
    logger.debug(f"delete took {time.time() - start}s")


def update_rdkit_tables(dataset_id: str, session: Session) -> None:
    """Updates RDKit PostgreSQL cartridge data."""
    logger.debug(f"Updating RDKit tables for {dataset_id=}")
    _update_rdkit_reactions(dataset_id, session)
    _update_rdkit_mols(dataset_id, session)


def _update_rdkit_reactions(dataset_id: str, session: Session) -> None:
    """Updates the RDKit reactions table."""
    logger.debug("Updating RDKit reactions")
    start = time.time()
    result = session.execute(
        text(
            """
            INSERT INTO rdkit.reactions (reaction_smiles, reaction)
            SELECT reaction_smiles, reaction_from_smiles(reaction_smiles::cstring)
            FROM (
                SELECT reaction_smiles
                    FROM ord.reaction
                    JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                    WHERE ord.dataset.dataset_id = :dataset_id
                      AND ord.reaction.rdkit_reaction_id IS NULL
                EXCEPT 
                SELECT reaction_smiles
                    FROM rdkit.reactions
            ) subquery
            """
        ),
        {"dataset_id": dataset_id},
    )
    logger.debug(f"Updating reactions took {time.time() - start:g}s ({result.rowcount} rows)")


def _update_rdkit_mols(dataset_id: str, session: Session) -> None:
    """Updates the RDKit mols table."""
    logger.debug("Updating RDKit mols")
    start = time.time()
    result = session.execute(
        text(
            """
            INSERT INTO rdkit.mols (smiles, mol, morgan_bfp, morgan_sfp)
            SELECT smiles, mol, morgan_bfp, morgan_sfp
            FROM (
                SELECT smiles, mol, morganbv_fp(mol) AS morgan_bfp, morgan_fp(mol) AS morgan_sfp
                FROM (
                    SELECT smiles, mol_from_smiles(smiles::cstring) AS mol
                    FROM (
                        SELECT smiles
                            -- NOTE(skearnes): This join path does not include non-input compounds like workups, 
                            -- internal standards, etc.
                            FROM ord.compound
                            JOIN ord.reaction_input ON ord.compound.reaction_input_id = ord.reaction_input.id
                            JOIN ord.reaction ON ord.reaction_input.reaction_id = ord.reaction.id
                            JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                            WHERE ord.dataset.dataset_id = :dataset_id
                              AND ord.compound.rdkit_mol_id IS NULL
                        UNION
                        SELECT smiles
                            FROM ord.product_compound
                            JOIN ord.reaction_outcome
                                ON ord.product_compound.reaction_outcome_id = ord.reaction_outcome.id
                            JOIN ord.reaction ON ord.reaction_outcome.reaction_id = ord.reaction.id
                            JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                            WHERE ord.dataset.dataset_id = :dataset_id
                              AND ord.product_compound.rdkit_mol_id IS NULL
                        EXCEPT
                        SELECT smiles
                            FROM rdkit.mols
                    ) smiles_subquery
                    -- See https://github.com/open-reaction-database/ord-schema/issues/672.
                    WHERE smiles NOT LIKE '%[Ti+5]%'
                ) mol_subquery
            ) fp_subquery
            ON CONFLICT (smiles) DO NOTHING
            """
        ),
        {"dataset_id": dataset_id},
    )
    logger.debug(f"Updating mols took {time.time() - start:g}s ({result.rowcount} rows)")


def update_rdkit_ids(dataset_id: str, session: Session) -> None:
    """Updates RDKit reaction and mol ID associations in the ORD tables."""
    logger.debug("Updating RDKit ID associations")
    start = time.time()
    # Update Reaction.
    session.execute(
        text(
            """
            UPDATE ord.reaction
            SET rdkit_reaction_id = subquery.rdkit_reaction_id
            FROM (
                SELECT ord.reaction.id, rdkit.reactions.id AS rdkit_reaction_id
                    FROM ord.reaction
                    JOIN rdkit.reactions USING (reaction_smiles)
                    JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                    WHERE ord.dataset.dataset_id = :dataset_id
                      AND ord.reaction.rdkit_reaction_id IS NULL
            ) AS subquery
            WHERE ord.reaction.id = subquery.id
            """
        ),
        {"dataset_id": dataset_id},
    )
    # Update Compound.
    session.execute(
        text(
            """
            UPDATE ord.compound
            SET rdkit_mol_id = subquery.rdkit_mol_id
            FROM (
                SELECT ord.compound.id, rdkit.mols.id AS rdkit_mol_id
                    FROM ord.compound
                    JOIN rdkit.mols USING (smiles)
                    JOIN ord.reaction_input ON ord.compound.reaction_input_id = ord.reaction_input.id
                    JOIN ord.reaction ON ord.reaction_input.reaction_id = ord.reaction.id
                    JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                    WHERE ord.dataset.dataset_id = :dataset_id
                      AND ord.compound.rdkit_mol_id IS NULL
            ) AS subquery
            WHERE ord.compound.id = subquery.id
            """
        ),
        {"dataset_id": dataset_id},
    )
    session.execute(
        text(
            """
            UPDATE ord.product_compound
            SET rdkit_mol_id = subquery.rdkit_mol_id
            FROM (
                SELECT ord.product_compound.id, rdkit.mols.id AS rdkit_mol_id
                    FROM ord.product_compound
                    JOIN rdkit.mols USING (smiles)
                    JOIN ord.reaction_outcome ON ord.product_compound.reaction_outcome_id = ord.reaction_outcome.id
                    JOIN ord.reaction ON ord.reaction_outcome.reaction_id = ord.reaction.id
                    JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                    WHERE ord.dataset.dataset_id = :dataset_id
                      AND ord.product_compound.rdkit_mol_id IS NULL
            ) AS subquery
            WHERE ord.product_compound.id = subquery.id
            """
        ),
        {"dataset_id": dataset_id},
    )
    logger.debug(f"Updating RDKit IDs took {time.time() - start:g}s")
