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
from typing import Any, cast
from unittest.mock import patch

from sqlalchemy import delete, select, text
from sqlalchemy.engine import Engine
from sqlalchemy.exc import NotSupportedError, OperationalError
from sqlalchemy.orm import Session

from ord_schema import parquet_dataset
from ord_schema.logging import get_logger
from ord_schema.orm.mappers import Base, Mappers, from_proto
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def get_connection_string(
    database: str,
    username: str,
    password: str,
    host: str = "localhost",
    port: int = 5432,
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
            connection.execute(
                text("CREATE EXTENSION IF NOT EXISTS tsm_system_rows")
            )  # For random sampling.
        except OperationalError:
            logger.warning(
                "tsm_system_rows cartridge is not installed; random sampling will be disabled"
            )
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
            logger.warning(
                "RDKit PostgreSQL cartridge is not installed; structure search will be disabled"
            )
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS btree_gist"))
        rdkit_cartridge = False
    with patch.dict(
        os.environ, {"ORD_POSTGRES_RDKIT": "1" if rdkit_cartridge else "0"}
    ):
        Base.metadata.create_all(engine)
    return rdkit_cartridge


def add_dataset(
    dataset: dataset_pb2.Dataset, session: Session, rdkit_cartridge: bool = True
) -> None:
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


# Reactions are flushed and expunged in batches of this size during streaming
# ingest. Larger values reduce SQL roundtrips at the cost of more pending ORM
# objects held in memory; 200 is unmeasured and worked fine on the existing
# ord-data fixtures.
_PARQUET_FLUSH_BATCH = 200


def add_parquet_dataset(
    path: str, session: Session, rdkit_cartridge: bool = True
) -> None:
    """Streams a Parquet-serialized Dataset into the ORM tables.

    Two streaming passes over the Parquet file:

    1. ``parquet_dataset.streaming_md5`` to compute ``md5`` and
       ``num_reactions`` without holding Reactions in memory.
    2. ``iter_reactions`` to insert Reaction rows in flush/expunge batches.

    This keeps the canonical Parquet MD5 formula in one place and bounds
    peak memory to one row group plus ``_PARQUET_FLUSH_BATCH`` ORM objects,
    regardless of dataset size.
    """
    start = time.time()
    metadata = parquet_dataset.read_metadata(path)
    logger.debug(f"Streaming Parquet Dataset {metadata.dataset_id}")
    md5_hex, num_reactions = parquet_dataset.streaming_md5(path)
    reaction_child_class = Mappers.Dataset.reactions.mapper.class_
    mapped_dataset = Mappers.Dataset(
        name=metadata.name,
        description=metadata.description,
        dataset_id=metadata.dataset_id,
        md5=md5_hex,
        num_reactions=num_reactions,
    )
    session.add(mapped_dataset)
    session.flush()
    pending: list = []

    def flush_batch() -> None:
        if not pending:
            return
        session.flush()
        for obj in pending:
            session.expunge(obj)
        pending.clear()

    for _, reaction in parquet_dataset.iter_reactions(path):
        reaction_mapper = from_proto(reaction, mapper=reaction_child_class)
        reaction_mapper.parent = mapped_dataset
        session.add(reaction_mapper)
        pending.append(reaction_mapper)
        if len(pending) >= _PARQUET_FLUSH_BATCH:
            flush_batch()
    flush_batch()
    logger.debug(
        f"add_parquet_dataset() took {time.time() - start:g}s ({num_reactions} reactions)"
    )
    if rdkit_cartridge:
        update_rdkit_tables(metadata.dataset_id, session)
        session.flush()
        update_rdkit_ids(metadata.dataset_id, session)


def get_dataset_md5(dataset_id: str, session: Session) -> str | None:
    """Returns the MD5 hash of the current version of a dataset, if it exists in the database."""
    result = session.execute(
        select(Mappers.Dataset.md5).where(Mappers.Dataset.dataset_id == dataset_id)
    )
    row = result.first()
    return row[0] if row else None


def get_dataset_size(dataset_id: str, session: Session) -> int:
    """Returns the number of reactions in a dataset."""
    result = session.execute(
        select(Mappers.Dataset.num_reactions).where(
            Mappers.Dataset.dataset_id == dataset_id
        )
    )
    row = result.first()
    if row is None:
        raise ValueError(dataset_id)
    return row[0]


def delete_dataset(dataset_id: str, session: Session) -> None:
    """Deletes a dataset from the database."""
    logger.debug(f"Deleting dataset {dataset_id}")
    start = time.time()
    session.execute(
        delete(Mappers.Dataset).where(Mappers.Dataset.dataset_id == dataset_id)
    )
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
        text("""
            INSERT INTO rdkit.reactions (reaction_smiles, reaction)
            SELECT reaction_smiles, reaction
            FROM (
                SELECT reaction_smiles, reaction_from_smiles(reaction_smiles::cstring) AS reaction
                FROM (
                    -- NOTE(skearnes): NOT EXISTS probes the unique reaction_smiles index per candidate instead
                    -- of EXCEPT-scanning all of rdkit.reactions; DISTINCT dedupes within the dataset. There is no
                    -- ON CONFLICT backstop here, so this relies on the RDKit phase running serially (see
                    -- add_datasets) to avoid unique violations on concurrent inserts of the same reaction_smiles.
                    -- reaction_smiles IS NOT NULL is required: unlike EXCEPT (which treats NULLs as equal),
                    -- NOT EXISTS never matches a NULL, so without it a no-SMILES reaction would re-insert a junk
                    -- (NULL, NULL) row on every run.
                    SELECT DISTINCT reaction_smiles
                        FROM ord.reaction
                        JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                        WHERE ord.dataset.dataset_id = :dataset_id
                          AND ord.reaction.rdkit_reaction_id IS NULL
                          AND ord.reaction.reaction_smiles IS NOT NULL
                          AND NOT EXISTS (
                              SELECT 1 FROM rdkit.reactions
                              WHERE rdkit.reactions.reaction_smiles = ord.reaction.reaction_smiles
                          )
                ) candidates
            ) computed
            -- reaction_from_smiles returns NULL for an unparseable reaction SMILES; skip those so we never
            -- insert a NULL-reaction row (mirrors the mol IS NOT NULL guard in _update_rdkit_mols).
            WHERE reaction IS NOT NULL
            """),
        {"dataset_id": dataset_id},
    )
    logger.debug(
        f"Updating reactions took {time.time() - start:g}s ({cast(Any, result).rowcount} rows)"
    )


def _update_rdkit_mols(dataset_id: str, session: Session) -> None:
    """Updates the RDKit mols table."""
    logger.debug("Updating RDKit mols")
    start = time.time()
    result = session.execute(
        text("""
            WITH new_smiles AS MATERIALIZED (
                -- Materialization barrier (AS MATERIALIZED): resolve the dataset's not-yet-linked SMILES that
                -- aren't already in rdkit.mols FIRST -- a cheap per-candidate probe of the unique smiles index --
                -- so the expensive mol_from_smiles/morgan_*_fp calls below run only on the survivors. Without the
                -- barrier the planner may evaluate those functions (and especially the WHERE mol IS NOT NULL
                -- guard's mol_from_smiles) on every candidate before the anti-join; observed ~240s for a no-op
                -- dataset. NOT EXISTS replaces the old EXCEPT, which scanned all of rdkit.mols but also happened
                -- to materialize for this same reason.
                SELECT smiles
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
                ) candidates
                WHERE smiles NOT LIKE '%[Ti+5]%'  -- See https://github.com/open-reaction-database/ord-schema/issues/672.
                  AND NOT EXISTS (SELECT 1 FROM rdkit.mols WHERE rdkit.mols.smiles = candidates.smiles)
            )
            INSERT INTO rdkit.mols (smiles, mol, morgan_bfp, morgan_sfp)
            SELECT smiles, mol, morganbv_fp(mol) AS morgan_bfp, morgan_fp(mol) AS morgan_sfp
            FROM (
                SELECT smiles, mol_from_smiles(smiles::cstring) AS mol FROM new_smiles
            ) computed
            -- mol_from_smiles returns NULL for unparseable SMILES; skip those so we never insert a NULL-mol
            -- row (which ON CONFLICT would not catch for a genuinely new SMILES) or feed NULL to the fingerprints.
            WHERE mol IS NOT NULL
            ON CONFLICT (smiles) DO NOTHING
            """),
        {"dataset_id": dataset_id},
    )
    logger.debug(
        f"Updating mols took {time.time() - start:g}s ({cast(Any, result).rowcount} rows)"
    )


def _link_mol_ids(
    session: Session,
    *,
    target_table: str,
    link_table: str,
    link_column: str,
    dataset_id: str,
) -> int:
    """Links rdkit.mols ids into a compound table for one dataset's unlinked rows.

    The Compound and ProductCompound updates are identical apart from the table
    and join path, so they share this helper. Identifier arguments are trusted
    literals supplied by ``update_rdkit_ids`` (never user input) and are
    interpolated into the statement; ``dataset_id`` is passed as a bind param.

    Args:
        session: Active SQLAlchemy session.
        target_table: Compound table to update (e.g. ``"ord.compound"``).
        link_table: Join table connecting ``target_table`` to ``ord.reaction``
            (e.g. ``"ord.reaction_input"``).
        link_column: Foreign-key column on ``target_table`` that references
            ``link_table`` (e.g. ``"reaction_input_id"``).
        dataset_id: Dataset to scope the update to.

    Returns:
        The number of rows updated.
    """
    result = session.execute(
        text(f"""
            UPDATE {target_table}
            SET rdkit_mol_id = rdkit.mols.id
            FROM rdkit.mols, {link_table}, ord.reaction, ord.dataset
            WHERE rdkit.mols.smiles = {target_table}.smiles
              AND {target_table}.{link_column} = {link_table}.id
              AND {link_table}.reaction_id = ord.reaction.id
              AND ord.reaction.dataset_id = ord.dataset.id
              AND ord.dataset.dataset_id = :dataset_id
              AND {target_table}.rdkit_mol_id IS NULL
            """),  # noqa: S608  (table/column names are internal constants, not user input)
        {"dataset_id": dataset_id},
    )
    return cast(Any, result).rowcount


def update_rdkit_ids(dataset_id: str, session: Session) -> None:
    """Updates RDKit reaction and mol ID associations in the ORD tables."""
    logger.debug("Updating RDKit ID associations")
    start = time.time()
    # NOTE(skearnes): These use the flat ``UPDATE ... FROM`` form (target updated in place) rather than
    # ``FROM (SELECT id, rdkit_id ...) WHERE target.id = subquery.id``, which materialized the pairs and then
    # re-joined the target by id. The rdkit join keys (reaction_smiles/smiles) are unique-indexed, and the
    # dataset scope is reached via the indexed foreign keys.
    # Update Reaction.
    reaction_result = session.execute(
        text("""
            UPDATE ord.reaction
            SET rdkit_reaction_id = rdkit.reactions.id
            FROM rdkit.reactions, ord.dataset
            WHERE rdkit.reactions.reaction_smiles = ord.reaction.reaction_smiles
              AND ord.reaction.dataset_id = ord.dataset.id
              AND ord.dataset.dataset_id = :dataset_id
              AND ord.reaction.rdkit_reaction_id IS NULL
            """),
        {"dataset_id": dataset_id},
    )
    reaction_rows = cast(Any, reaction_result).rowcount
    compound_rows = _link_mol_ids(
        session,
        target_table="ord.compound",
        link_table="ord.reaction_input",
        link_column="reaction_input_id",
        dataset_id=dataset_id,
    )
    product_compound_rows = _link_mol_ids(
        session,
        target_table="ord.product_compound",
        link_table="ord.reaction_outcome",
        link_column="reaction_outcome_id",
        dataset_id=dataset_id,
    )
    logger.debug(
        f"Updating RDKit IDs took {time.time() - start:g}s "
        f"(reaction={reaction_rows}, compound={compound_rows}, product_compound={product_compound_rows})"
    )
