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

import datetime
import os
import time
from typing import Any, cast
from unittest.mock import patch

from dateutil import parser
from rdkit import Chem
from sqlalchemy import Uuid, bindparam, delete, inspect, select, text
from sqlalchemy.engine import Engine
from sqlalchemy.exc import NotSupportedError, OperationalError, ProgrammingError
from sqlalchemy.orm import RelationshipDirection, Session
from tqdm import tqdm
from uuid6 import uuid7

from ord_schema import message_helpers, parquet
from ord_schema.logging import get_logger
from ord_schema.orm.mappers import Base, Mappers, from_proto, to_proto
from ord_schema.orm.public_mappers import DatasetMetadata
from ord_schema.proto import dataset_pb2, reaction_pb2

# COPY the ingest tables parents-first so foreign keys resolve; sorted_tables orders by FK
# dependency. Precomputed once; only the tables populated by a given batch are written.
_COPY_TABLE_ORDER = [table.fullname for table in Base.metadata.sorted_tables]

try:
    from ord_schema.orm.reaction_class import update_reaction_classes
except Exception as error:  # noqa: BLE001
    # The optional 'reaction-class' extra pulls in rxn-insight -> rxnmapper/torch, which
    # can fail to import for many reasons (absent package, native-library OSError, init
    # error). Catch them all so the ORM stays usable when classify_reactions is False;
    # the opt-in path re-raises with the original error as the cause.
    update_reaction_classes = None  # ty: ignore[invalid-assignment]
    _reaction_class_import_error: Exception | None = error
else:
    _reaction_class_import_error = None

logger = get_logger(__name__)


def _classify_reactions(dataset_id: str, session: Session) -> None:
    """Populates reaction class/name columns, requiring the optional extra."""
    if update_reaction_classes is None:
        raise ImportError(
            "Reaction classification requires the 'reaction-class' extra: "
            "pip install ord-schema[reaction-class]"
        ) from _reaction_class_import_error
    update_reaction_classes(dataset_id, session)


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
        # Derived, best-effort data that is not part of the proto (e.g. reaction class).
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS derived"))
    with engine.begin() as connection:
        # Pin the database's default search_path to public. The connecting role is often
        # named "ord", which would make Postgres's default ("$user", public) resolve the
        # ord schema first; instead require explicit ord/derived/rdkit qualification (the
        # ORM qualifies every table) and keep only public on the path, where the RDKit
        # cartridge functions live. Best-effort: a non-owner connection cannot ALTER
        # DATABASE, but the ORM still works since it never relies on the default path.
        try:
            connection.execute(
                text(
                    "DO $$ BEGIN EXECUTE format("
                    "'ALTER DATABASE %I SET search_path TO public', "
                    "current_database()); END $$;"
                )
            )
        except (OperationalError, ProgrammingError) as error:
            logger.warning(f"Could not set the default search_path to public: {error}")
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


def add_dataset(dataset: dataset_pb2.Dataset, session: Session) -> None:
    """Ingests a dataset, writing the ``ord.*`` search index and ``public.*`` payload.

    Derived data (SMILES, RDKit links, reaction classes) is populated separately by
    ``update_derived_data`` so ingest and derivation can run independently.

    Args:
        dataset: Dataset to add.
        session: SQLAlchemy session.
    """
    logger.debug(f"Adding dataset {dataset.dataset_id}")
    start = time.time()
    mapped_dataset = from_proto(dataset)
    logger.debug(f"from_proto() took {time.time() - start:g}s")
    session.add(mapped_dataset)
    session.flush()
    set_submitted_at(dataset.dataset_id, session)


def update_derived_data(
    dataset_id: str,
    session: Session,
    *,
    rdkit_cartridge: bool = True,
    classify_reactions: bool = False,
) -> None:
    """Populates the derived tables for an already-ingested dataset.

    Idempotent (NOT EXISTS guards), so it is safe to re-run to backfill or recompute
    derived data over datasets that are already present.

    Args:
        dataset_id: Dataset to derive.
        session: SQLAlchemy session.
        rdkit_cartridge: Whether to populate RDKit cartridge tables and links.
        classify_reactions: Whether to assign reaction class/name labels; requires the
            optional ``reaction-class`` extra.
    """
    update_derived_tables(dataset_id, session)
    if rdkit_cartridge:
        session.flush()
        update_rdkit_tables(dataset_id, session)
        session.flush()
        update_rdkit_ids(dataset_id, session)
    if classify_reactions:
        session.flush()
        _classify_reactions(dataset_id, session)


def _parse_date(value: str) -> datetime.date | None:
    """Parses the date from a free-text provenance timestamp, or None if unparseable.

    Provenance times are free text in several formats (ISO, US, ctime), so they're
    parsed with dateutil. Only the calendar date is kept; recency browsing needs no
    finer granularity, and dropping the time avoids the timezones these values may
    carry.
    """
    try:
        return parser.parse(value).date()
    except (ValueError, OverflowError):
        return None


def set_submitted_at(dataset_id: str, session: Session) -> None:
    """Sets ``dataset.submitted_at`` from the dataset's latest reaction record event.

    Uses the last ``record_modified`` entry of an arbitrary reaction in the dataset:
    record events are appended in chronological order, so the last one is the most
    recent. Falls back to ``record_created`` (required) when there is no
    ``record_modified``. Reactions in a dataset share their submission-pipeline
    timestamps, so one reaction is representative and avoids scanning every reaction.
    Leaves NULL only when the timestamp is missing or unparseable.
    """
    session.flush()  # set_submitted_at reads via raw SQL; make pending rows visible.
    value = session.execute(
        text(
            """
            SELECT date_time.value
            FROM ord.reaction
            JOIN ord.reaction_provenance
                ON reaction_provenance.reaction_id = reaction.id
            JOIN ord.record_event
                ON record_event.reaction_provenance_id = reaction_provenance.id
                AND record_event.ord_schema_context IN (
                    'ReactionProvenance.record_created',
                    'ReactionProvenance.record_modified'
                )
            JOIN ord.date_time ON date_time.record_event_id = record_event.id
            WHERE reaction.id = (
                SELECT id FROM ord.reaction
                WHERE dataset_id = (
                    SELECT id FROM ord.dataset WHERE dataset_id = :dataset_id
                )
                LIMIT 1
            )
            -- Prefer record_modified over record_created, then the last entry
            -- (record events are inserted in chronological order, so the highest
            -- id is the most recent).
            ORDER BY
                (record_event.ord_schema_context
                    = 'ReactionProvenance.record_modified') DESC,
                record_event.id DESC
            LIMIT 1
            """
        ),
        {"dataset_id": dataset_id},
    ).scalar_one_or_none()
    session.execute(
        text(
            "UPDATE public.datasets SET submitted_at = :submitted_at "
            "WHERE dataset_id = :dataset_id"
        ),
        {
            "submitted_at": _parse_date(value) if value is not None else None,
            "dataset_id": dataset_id,
        },
    )


def backfill_submission_times(session: Session) -> None:
    """Populates ``submitted_at`` for datasets that don't have it yet.

    One-time/maintenance helper for databases loaded before the column existed; new
    datasets are populated at ingest by ``add_dataset``/``add_parquet_dataset``.
    """
    dataset_ids = (
        session.execute(
            select(DatasetMetadata.dataset_id).where(
                DatasetMetadata.submitted_at.is_(None)
            )
        )
        .scalars()
        .all()
    )
    for dataset_id in dataset_ids:
        set_submitted_at(dataset_id, session)


# Reactions are built into ORM trees and COPYed in batches of this size during streaming
# ingest. Larger values reduce COPY roundtrips at the cost of more trees held in memory.
_COPY_BATCH = 1000

# Number of reactions/compounds the derived passes load and process at a time. The ids needing
# derivation are fetched up front, but the heavy work (proto deserialization, SMILES generation)
# and the pending inserts are limited to this many rows at a time.
_DERIVED_BATCH = 1000


def _collect_rows(root: Any, rows_by_table: dict[str, tuple[Any, list]]) -> None:
    """Walks one ORM object tree, minting keys and wiring FKs, into per-table column tuples.

    Each surrogate (Uuid) primary key is assigned a UUIDv7 value, and every child's foreign-key
    column is set from the parent's referenced column (from the relationship metadata), so the
    rows are self-consistent without a database round trip. The resulting tuples can be streamed
    with COPY. ``rows_by_table`` maps a table's full name to ``(table, rows)``.
    """
    stack = [root]
    seen: set[int] = set()
    while stack:
        node = stack.pop()
        if id(node) in seen:
            continue
        seen.add(id(node))
        mapper = inspect(node).mapper
        table = mapper.local_table
        for key_column in table.primary_key:
            if (
                isinstance(key_column.type, Uuid)
                and getattr(node, key_column.key) is None
            ):
                setattr(node, key_column.key, uuid7())
        for relationship in mapper.relationships:
            if relationship.direction is RelationshipDirection.MANYTOONE:
                continue  # The foreign key is on this side; the owning parent sets it.
            children = getattr(node, relationship.key)
            if children is None:
                continue
            if not isinstance(children, list):
                children = [children]
            for child in children:
                for local, remote in relationship.local_remote_pairs:
                    setattr(child, remote.key, getattr(node, local.key))
                stack.append(child)
        entry = rows_by_table.get(table.fullname)
        if entry is None:
            entry = (table, [])
            rows_by_table[table.fullname] = entry
        # Under single-table polymorphic inheritance, a sibling subclass's foreign-key column
        # shares the table but is not an attribute of this instance, so it is NULL for this row.
        entry[1].append(tuple(getattr(node, c.key, None) for c in table.columns))


def _copy_rows(
    dbapi_connection: Any, rows_by_table: dict[str, tuple[Any, list]]
) -> None:
    """Streams collected rows into each table with COPY, parents before children."""
    with dbapi_connection.cursor() as cursor:
        for fullname in _COPY_TABLE_ORDER:
            entry = rows_by_table.get(fullname)
            if entry is None:
                continue
            table, rows = entry
            columns = ", ".join(f'"{column.name}"' for column in table.columns)
            # fullname/columns are internal schema constants, not user input.
            with cursor.copy(f"COPY {fullname} ({columns}) FROM STDIN") as copy:
                for row in rows:
                    copy.write_row(row)


def add_parquet_dataset(path: str, session: Session) -> None:
    """Streams a Parquet-serialized Dataset into the ORM tables with COPY.

    Two streaming passes over the Parquet file:

    1. ``parquet.streaming_md5`` computes ``md5`` and ``num_reactions`` without holding
       Reactions in memory.
    2. ``iter_reactions`` builds ORM trees a batch at a time; each tree is assigned UUIDv7
       primary keys with foreign keys wired from the relationship metadata (``_collect_rows``),
       and the rows are streamed with ``COPY`` (``_copy_rows``) rather than the ORM unit of
       work. Peak memory is bounded to one row group plus ``_COPY_BATCH`` trees.

    Derived data is populated separately by ``update_derived_data``.

    Args:
        path: Path to the Parquet-serialized Dataset.
        session: SQLAlchemy session.
    """
    start = time.time()
    metadata = parquet.load_metadata(path)
    logger.debug(f"Streaming Parquet Dataset {metadata.dataset_id}")
    md5_hex, num_reactions = parquet.streaming_md5(path)
    reaction_child_class = Mappers.Dataset.reactions.mapper.class_
    dataset = Mappers.Dataset(
        name=metadata.name,
        description=metadata.description,
        dataset_id=metadata.dataset_id,
        metadata_row=DatasetMetadata(md5=md5_hex, num_reactions=num_reactions),
    )
    dbapi_connection = session.connection().connection
    dataset_rows: dict[str, tuple[Any, list]] = {}
    _collect_rows(dataset, dataset_rows)
    _copy_rows(dbapi_connection, dataset_rows)
    batch: list = []

    def copy_batch() -> None:
        if not batch:
            return
        rows_by_table: dict[str, tuple[Any, list]] = {}
        for reaction in batch:
            reaction_mapper = from_proto(reaction, mapper=reaction_child_class)
            reaction_mapper.dataset_id = dataset.id
            _collect_rows(reaction_mapper, rows_by_table)
        _copy_rows(dbapi_connection, rows_by_table)
        batch.clear()

    for _, reaction in tqdm(
        parquet.iter_reactions(path),
        total=num_reactions,
        desc=f"ingest {metadata.dataset_id}",
        unit="rxn",
        leave=False,
    ):
        batch.append(reaction)
        if len(batch) >= _COPY_BATCH:
            copy_batch()
    copy_batch()
    set_submitted_at(metadata.dataset_id, session)
    logger.debug(
        f"add_parquet_dataset() took {time.time() - start:g}s ({num_reactions} reactions)"
    )


def get_dataset_md5(dataset_id: str, session: Session) -> str | None:
    """Returns the MD5 hash of the current version of a dataset, if it exists in the database."""
    result = session.execute(
        select(DatasetMetadata.md5).where(DatasetMetadata.dataset_id == dataset_id)
    )
    row = result.first()
    return row[0] if row else None


def get_dataset_size(dataset_id: str, session: Session) -> int:
    """Returns the number of reactions in a dataset."""
    result = session.execute(
        select(DatasetMetadata.num_reactions).where(
            DatasetMetadata.dataset_id == dataset_id
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


def update_derived_tables(dataset_id: str, session: Session) -> None:
    """Populates the derived SMILES tables from the search index.

    Reaction SMILES come from the ground-truth proto in public.reactions; compound SMILES
    from each ord.compound row's reconstructed message. Idempotent (skips rows that already
    have a derived entry); runs before the RDKit pass, which reads the SMILES. Reactions and
    compounds are processed in batches of _DERIVED_BATCH so the heavy per-row work (proto
    deserialization, SMILES generation) and pending inserts stay within one batch.
    """
    logger.debug(f"Updating derived tables for {dataset_id=}")
    start = time.time()
    # dataset_id and the compound link columns live on the polymorphic child mappers, so
    # rows are scoped to the dataset via raw SQL (like the RDKit pass).
    dataset_pk = session.execute(
        text("SELECT id FROM ord.dataset WHERE dataset_id = :dataset_id"),
        {"dataset_id": dataset_id},
    ).scalar_one()
    # Reaction SMILES from the served proto (no ORM objects loaded), keyed by ord.reaction.id.
    # Resolve the ids needing derivation in one indexed pass, then load and parse the (large)
    # protos in batches of _DERIVED_BATCH.
    reaction_ids = (
        session.execute(
            text("""
                SELECT ord.reaction.id
                FROM ord.reaction
                JOIN public.reactions
                    ON public.reactions.reaction_id = ord.reaction.reaction_id
                WHERE ord.reaction.dataset_id = :id
                  AND NOT EXISTS (
                      SELECT 1 FROM derived.reaction_smiles
                      WHERE derived.reaction_smiles.reaction_id = ord.reaction.id
                  )
                """),
            {"id": dataset_pk},
        )
        .scalars()
        .all()
    )
    select_protos = text("""
        SELECT ord.reaction.id, public.reactions.proto
        FROM ord.reaction
        JOIN public.reactions
            ON public.reactions.reaction_id = ord.reaction.reaction_id
        WHERE ord.reaction.id IN :ids
        """).bindparams(bindparam("ids", expanding=True))
    insert_reaction_smiles = text(
        "INSERT INTO derived.reaction_smiles (reaction_id, reaction_smiles) "
        "VALUES (:reaction_id, :reaction_smiles)"
    )
    for batch_start in tqdm(
        range(0, len(reaction_ids), _DERIVED_BATCH),
        desc=f"reaction SMILES {dataset_id}",
        unit="batch",
        leave=False,
    ):
        batch_ids = reaction_ids[batch_start : batch_start + _DERIVED_BATCH]
        inserts = []
        for reaction_id, proto in session.execute(select_protos, {"ids": batch_ids}):
            try:
                reaction_smiles = message_helpers.get_reaction_smiles(
                    reaction_pb2.Reaction.FromString(proto),
                    generate_if_missing=True,
                    allow_incomplete=False,
                    validate=True,
                )
            except ValueError as error:
                logger.debug(
                    f"No reaction SMILES for reaction id={reaction_id}: {error}"
                )
                continue
            tokens = (
                reaction_smiles.split() if reaction_smiles else []
            )  # Handle CXSMILES.
            if tokens:
                inserts.append(
                    {"reaction_id": reaction_id, "reaction_smiles": tokens[0]}
                )
        if inserts:
            session.execute(insert_reaction_smiles, inserts)
    _update_compound_smiles(
        session,
        dataset_pk=dataset_pk,
        compound_table="ord.compound",
        link_table="ord.reaction_input",
        link_column="reaction_input_id",
        compound_class=Mappers.Compound,
        derived_table="derived.compound_smiles",
        derived_id="compound_id",
    )
    _update_compound_smiles(
        session,
        dataset_pk=dataset_pk,
        compound_table="ord.product_compound",
        link_table="ord.reaction_outcome",
        link_column="reaction_outcome_id",
        compound_class=Mappers.ProductCompound,
        derived_table="derived.product_compound_smiles",
        derived_id="product_compound_id",
    )
    logger.debug(f"Updating derived tables took {time.time() - start:g}s")


def _update_compound_smiles(
    session: Session,
    *,
    dataset_pk: int,
    compound_table: str,
    link_table: str,
    link_column: str,
    compound_class: Any,
    derived_table: str,
    derived_id: str,
) -> None:
    """Derives SMILES for one (product) compound table's not-yet-derived rows.

    The common case -- a compound with a stored SMILES identifier -- is served set-based: the
    first SMILES identifier for each compound in the batch is fetched from ord.compound_identifier
    in a single query and canonicalized with RDKit, avoiding a per-compound ORM load. Compounds
    without a stored SMILES fall back to reconstructing the message via to_proto and computing the
    SMILES from its other identifiers. Ids are resolved up front and processed in batches of
    _DERIVED_BATCH so memory stays bounded.
    """
    compound_ids = (
        session.execute(
            text(f"""
                SELECT {compound_table}.id
                FROM {compound_table}
                JOIN {link_table} ON {compound_table}.{link_column} = {link_table}.id
                JOIN ord.reaction ON {link_table}.reaction_id = ord.reaction.id
                WHERE ord.reaction.dataset_id = :id
                  AND NOT EXISTS (
                      SELECT 1 FROM {derived_table}
                      WHERE {derived_table}.{derived_id} = {compound_table}.id
                  )
                """),  # noqa: S608  (table/column names are internal constants, not user input)
            {"id": dataset_pk},
        )
        .scalars()
        .all()
    )
    # The first SMILES identifier (lowest id == proto order) for each requested compound. The FK
    # column is indexed, so one query per batch replaces an ORM load plus lazy child fetches per
    # compound. Interpolated names are internal constants (not user input); see S608 below.
    select_smiles = text(f"""
        SELECT DISTINCT ON (ord.compound_identifier.{derived_id})
               ord.compound_identifier.{derived_id}, ord.compound_identifier.value
        FROM ord.compound_identifier
        WHERE ord.compound_identifier.{derived_id} IN :ids
          AND ord.compound_identifier.type = 'SMILES'
        ORDER BY ord.compound_identifier.{derived_id}, ord.compound_identifier.id
        """).bindparams(bindparam("ids", expanding=True))  # noqa: S608
    insert_smiles = text(
        f"INSERT INTO {derived_table} ({derived_id}, smiles) "  # noqa: S608
        f"VALUES (:{derived_id}, :smiles)"
    )
    for batch_start in tqdm(
        range(0, len(compound_ids), _DERIVED_BATCH),
        desc=derived_table.rsplit(".", maxsplit=1)[-1],
        unit="batch",
        leave=False,
    ):
        batch_ids = compound_ids[batch_start : batch_start + _DERIVED_BATCH]
        stored_smiles = {
            row[0]: row[1] for row in session.execute(select_smiles, {"ids": batch_ids})
        }
        inserts = []
        for compound_id in batch_ids:
            value = stored_smiles.get(compound_id)
            # An absent or empty SMILES identifier both fall through to the reconstruction path,
            # matching smiles_from_compound's `get_compound_smiles(...) or ...` falsiness (an
            # empty string is not a parseable structure to canonicalize).
            if value:
                # Canonicalize the stored SMILES, matching smiles_from_compound's default.
                mol = Chem.MolFromSmiles(value)
                if mol is None:
                    logger.debug(
                        f"Cannot parse SMILES for compound id={compound_id}: {value}"
                    )
                    continue
                smiles = Chem.MolToSmiles(mol)
            else:
                # No stored SMILES: reconstruct the message and derive from other identifiers.
                compound = session.get(compound_class, compound_id)
                assert compound is not None  # Selected by id above.
                try:
                    smiles = message_helpers.smiles_from_compound(
                        cast(
                            "reaction_pb2.Compound | reaction_pb2.ProductCompound",
                            to_proto(compound),
                        )
                    )
                except ValueError:
                    continue
            inserts.append({derived_id: compound_id, "smiles": smiles})
        if inserts:
            session.execute(insert_smiles, inserts)


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
                    SELECT DISTINCT derived.reaction_smiles.reaction_smiles
                        FROM derived.reaction_smiles
                        JOIN ord.reaction ON ord.reaction.id = derived.reaction_smiles.reaction_id
                        JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                        WHERE ord.dataset.dataset_id = :dataset_id
                          AND derived.reaction_smiles.rdkit_reaction_id IS NULL
                          AND derived.reaction_smiles.reaction_smiles IS NOT NULL
                          AND NOT EXISTS (
                              SELECT 1 FROM rdkit.reactions
                              WHERE rdkit.reactions.reaction_smiles = derived.reaction_smiles.reaction_smiles
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
                    SELECT derived.compound_smiles.smiles
                        -- NOTE(skearnes): This join path does not include non-input compounds like workups,
                        -- internal standards, etc.
                        FROM derived.compound_smiles
                        JOIN ord.compound ON ord.compound.id = derived.compound_smiles.compound_id
                        JOIN ord.reaction_input ON ord.compound.reaction_input_id = ord.reaction_input.id
                        JOIN ord.reaction ON ord.reaction_input.reaction_id = ord.reaction.id
                        JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                        WHERE ord.dataset.dataset_id = :dataset_id
                          AND derived.compound_smiles.rdkit_mol_id IS NULL
                    UNION
                    SELECT derived.product_compound_smiles.smiles
                        FROM derived.product_compound_smiles
                        JOIN ord.product_compound
                            ON ord.product_compound.id = derived.product_compound_smiles.product_compound_id
                        JOIN ord.reaction_outcome
                            ON ord.product_compound.reaction_outcome_id = ord.reaction_outcome.id
                        JOIN ord.reaction ON ord.reaction_outcome.reaction_id = ord.reaction.id
                        JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
                        WHERE ord.dataset.dataset_id = :dataset_id
                          AND derived.product_compound_smiles.rdkit_mol_id IS NULL
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
    derived_table: str,
    derived_id_column: str,
    compound_table: str,
    link_table: str,
    link_column: str,
    dataset_id: str,
) -> int:
    """Links rdkit.mols ids into a derived compound-SMILES table for one dataset.

    The Compound and ProductCompound updates are identical apart from the tables
    and join path, so they share this helper. Identifier arguments are trusted
    literals supplied by ``update_rdkit_ids`` (never user input) and are
    interpolated into the statement; ``dataset_id`` is passed as a bind param.

    Args:
        session: Active SQLAlchemy session.
        derived_table: Derived table to update (e.g. ``"derived.compound_smiles"``).
        derived_id_column: ``derived_table`` foreign key to ``compound_table.id``
            (e.g. ``"compound_id"``).
        compound_table: ORD compound table (e.g. ``"ord.compound"``), joined to
            reach the dataset via ``link_table``.
        link_table: Join table connecting ``compound_table`` to ``ord.reaction``
            (e.g. ``"ord.reaction_input"``).
        link_column: Foreign-key column on ``compound_table`` that references
            ``link_table`` (e.g. ``"reaction_input_id"``).
        dataset_id: Dataset to scope the update to.

    Returns:
        The number of rows updated.
    """
    result = session.execute(
        text(f"""
            UPDATE {derived_table}
            SET rdkit_mol_id = rdkit.mols.id
            FROM rdkit.mols, {compound_table}, {link_table}, ord.reaction, ord.dataset
            WHERE rdkit.mols.smiles = {derived_table}.smiles
              AND {derived_table}.{derived_id_column} = {compound_table}.id
              AND {compound_table}.{link_column} = {link_table}.id
              AND {link_table}.reaction_id = ord.reaction.id
              AND ord.reaction.dataset_id = ord.dataset.id
              AND ord.dataset.dataset_id = :dataset_id
              AND {derived_table}.rdkit_mol_id IS NULL
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
    # Update derived.reaction_smiles (reaction_smiles and rdkit_reaction_id live there).
    reaction_result = session.execute(
        text("""
            UPDATE derived.reaction_smiles
            SET rdkit_reaction_id = rdkit.reactions.id
            FROM rdkit.reactions, ord.reaction, ord.dataset
            WHERE rdkit.reactions.reaction_smiles = derived.reaction_smiles.reaction_smiles
              AND derived.reaction_smiles.reaction_id = ord.reaction.id
              AND ord.reaction.dataset_id = ord.dataset.id
              AND ord.dataset.dataset_id = :dataset_id
              AND derived.reaction_smiles.rdkit_reaction_id IS NULL
            """),
        {"dataset_id": dataset_id},
    )
    reaction_rows = cast(Any, reaction_result).rowcount
    compound_rows = _link_mol_ids(
        session,
        derived_table="derived.compound_smiles",
        derived_id_column="compound_id",
        compound_table="ord.compound",
        link_table="ord.reaction_input",
        link_column="reaction_input_id",
        dataset_id=dataset_id,
    )
    product_compound_rows = _link_mol_ids(
        session,
        derived_table="derived.product_compound_smiles",
        derived_id_column="product_compound_id",
        compound_table="ord.product_compound",
        link_table="ord.reaction_outcome",
        link_column="reaction_outcome_id",
        dataset_id=dataset_id,
    )
    logger.debug(
        f"Updating RDKit IDs took {time.time() - start:g}s "
        f"(reaction={reaction_rows}, compound={compound_rows}, product_compound={product_compound_rows})"
    )
