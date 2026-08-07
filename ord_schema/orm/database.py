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
import uuid
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
from ord_schema.orm.sharding import shard_predicate as _shard_predicate
from ord_schema.proto import dataset_pb2, reaction_pb2

# COPY the ingest tables parents-first so foreign keys resolve; sorted_tables orders by
# FK dependency. Precomputed once; only the tables populated by a given batch are
# written.
_COPY_TABLE_ORDER = [table.fullname for table in Base.metadata.sorted_tables]

try:
    from ord_schema.orm.reaction_class import update_reaction_classes
except Exception as error:  # noqa: BLE001
    # The optional 'reaction-class' extra pulls in rxn-insight -> rxnmapper/torch, which
    # can fail to import many ways (missing package, native-library OSError). Catch all
    # so the ORM stays usable without it; _classify_reactions re-raises with this error
    # as the cause.
    update_reaction_classes = None  # ty: ignore[invalid-assignment]
    _reaction_class_import_error: Exception | None = error
else:
    _reaction_class_import_error = None

logger = get_logger(__name__)


def _classify_reactions(
    dataset_id: str, session: Session, *, shard: tuple[int, int] | None = None
) -> None:
    """Populates reaction class/name columns, requiring the optional extra."""
    if update_reaction_classes is None:
        raise ImportError(
            "Reaction classification requires the 'reaction-class' extra: "
            "pip install ord-schema[reaction-class]"
        ) from _reaction_class_import_error
    update_reaction_classes(dataset_id, session, shard=shard)


def classify_dataset(
    dataset_id: str, session: Session, *, shard: tuple[int, int] | None = None
) -> None:
    """Labels reaction class/name for an already-derived dataset (classification only).

    SMILES derivation is done separately (and, in the loader, sharded); this does not
    re-derive it, so a failed SMILES shard is not silently backfilled here. ``shard``
    (index, num_shards) restricts to one hash-partition of the dataset's reaction ids.
    Requires the ``reaction-class`` extra.
    """
    _classify_reactions(dataset_id, session, shard=shard)


def get_connection_string(
    database: str,
    username: str,
    password: str,
    host: str = "localhost",
    port: int = 5432,
) -> str:
    """Creates an SQLAlchemy connection string."""
    return f"postgresql+psycopg://{username}:{password}@{host}:{port}/{database}?client_encoding=utf8"


# PostgreSQL 12 (server_version_num 120000), when AS MATERIALIZED CTEs -- used
# throughout the derived/RDKit passes -- were introduced.
_MIN_SERVER_VERSION_NUM = 120000


def _check_server_version(connection: Any) -> None:
    """Raises if the server predates the minimum supported PostgreSQL version.

    The derived/RDKit passes use ``AS MATERIALIZED`` CTEs, which an older server
    rejects with a parse error deep inside the RDKit update. Checking at the entry
    points (prepare_database and the update_rdkit_* functions) turns that into a clear
    message before any such SQL runs. The result is memoized on the connection, since
    the version is immutable and these run once per shard.

    Raises:
        RuntimeError: If the server is older than PostgreSQL 12.
    """
    info = getattr(connection, "info", None)
    if isinstance(info, dict) and info.get("ord_server_version_checked"):
        return
    version_num = int(connection.execute(text("SHOW server_version_num")).scalar_one())
    if version_num < _MIN_SERVER_VERSION_NUM:
        raise RuntimeError(
            "PostgreSQL 12+ is required (the derived/RDKit passes use "
            "AS MATERIALIZED CTEs); "
            f"server reports server_version_num={version_num}."
        )
    if isinstance(info, dict):
        info["ord_server_version_checked"] = True


def prepare_database(engine: Engine) -> bool:
    """Prepares the database and creates the ORM table structure.

    Args:
        engine: SQLAlchemy Engine (server must be PostgreSQL 12+).

    Returns:
        Whether the RDKit PostgreSQL cartridge is installed.

    Raises:
        RuntimeError: If the server is older than PostgreSQL 12.
    """
    with engine.begin() as connection:
        _check_server_version(connection)
        try:
            connection.execute(
                text("CREATE EXTENSION IF NOT EXISTS tsm_system_rows")
            )  # For random sampling.
        except OperationalError:
            logger.warning(
                "tsm_system_rows cartridge is not installed; "
                "random sampling will be disabled"
            )
    with engine.begin() as connection:
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS ord"))
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS rdkit"))
        # Derived, best-effort data that is not part of the proto (e.g. reaction class).
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS derived"))
    with engine.begin() as connection:
        # Pin the default search_path to public. The role is often named "ord", so
        # Postgres's default ("$user", public) would resolve the ord schema first; the
        # ORM qualifies every table, so keep only public -- where the RDKit cartridge
        # functions live. Best-effort: a non-owner connection cannot ALTER DATABASE, but
        # the ORM never relies on the path anyway.
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
            # NOTE(skearnes): The RDKit PostgreSQL extension works best in the public
            # schema.
            connection.execute(text("CREATE EXTENSION IF NOT EXISTS rdkit"))
        rdkit_cartridge = True
    except (OperationalError, NotSupportedError):
        with engine.begin() as connection:
            logger.warning(
                "RDKit PostgreSQL cartridge is not installed; "
                "structure search will be disabled"
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
    rederive: bool = False,
) -> None:
    """Populates the derived tables for an already-ingested dataset.

    Idempotent (NOT EXISTS guards), so it is safe to re-run to backfill derived data
    over datasets that are already present. Backfilling is not the same as rebuilding:
    the guards mean an existing row is never revisited, so use ``rederive`` when the
    derivation itself has changed.

    Args:
        dataset_id: Dataset to derive.
        session: SQLAlchemy session.
        rdkit_cartridge: Whether to populate RDKit cartridge tables and links.
        classify_reactions: Whether to assign reaction class/name labels; requires the
            optional ``reaction-class`` extra.
        rederive: Whether to delete existing derived rows first so they are recomputed.
    """
    update_derived_tables(dataset_id, session, rederive=rederive)
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

    Provenance times are free text in several formats (ISO, US, ctime), so dateutil
    parses them. Only the calendar date is kept: recency browsing needs no finer
    granularity, and dropping the time sidesteps the timezones these values may carry.
    """
    try:
        return parser.parse(value).date()
    except (ValueError, OverflowError):
        return None


def set_submitted_at(dataset_id: str, session: Session) -> None:
    """Sets ``dataset.submitted_at`` from the dataset's latest reaction record event.

    Reactions in a dataset share submission-pipeline timestamps, so one arbitrary
    reaction is representative (no full scan). Uses its last ``record_modified`` entry
    (record events are appended chronologically), falling back to the required
    ``record_created``. Leaves NULL when the timestamp is missing or unparseable.
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

# Number of reactions/compounds the derived passes load and process at a time. The ids
# needing derivation are fetched up front, but the heavy work (proto deserialization,
# SMILES generation) and the pending inserts are limited to this many rows at a time.
_DERIVED_BATCH = 1000


def _collect_rows(root: Any, rows_by_table: dict[str, tuple[Any, list]]) -> None:
    """Walks one ORM object tree into per-table column tuples, minting keys and FKs.

    Each surrogate (Uuid) primary key is assigned a UUIDv7 value, and every child's
    foreign-key column is set from the parent's referenced column (from the relationship
    metadata), so the rows are self-consistent without a database round trip. The
    resulting tuples can be streamed with COPY. ``rows_by_table`` maps a table's full
    name to ``(table, rows)``.
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
        # Every PK must be populated now: Uuid PKs are minted above; text PKs (e.g.
        # public.reactions.reaction_id) come from the proto or FK wiring from the
        # parent. A NULL here means a new table has a PK that is neither, which would
        # otherwise fail inside COPY.
        for key_column in table.primary_key:
            assert getattr(node, key_column.key, None) is not None, (
                f"unpopulated primary key {table.fullname}.{key_column.name}; "
                "a new table needs "
                "a Uuid primary key or one wired from a parent foreign key"
            )
        # Under single-table polymorphic inheritance, a sibling subclass's foreign-key
        # column shares the table but is not an attribute of this instance, so it is
        # NULL for this row.
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


def add_parquet_dataset_row(
    path: str, dataset_uuid: uuid.UUID, session: Session
) -> None:
    """Inserts scalars-only ``ord.dataset`` search-index row (no reactions/metadata).

    The surrogate ``id`` is supplied by the caller so a subsequent sharded reaction
    load can reference it without a round trip. The ``public.datasets`` metadata row is
    written last, by ``add_parquet_dataset_metadata``, so its presence marks a fully
    loaded dataset.

    Args:
        path: Path to the Parquet-serialized Dataset.
        dataset_uuid: Surrogate primary key for the ``ord.dataset`` row.
        session: SQLAlchemy session.
    """
    metadata = parquet.DatasetView(path)
    session.add(
        Mappers.Dataset(
            id=dataset_uuid,
            name=metadata.name,
            description=metadata.description,
            dataset_id=metadata.dataset_id,
        )
    )
    session.flush()


def add_parquet_reactions(
    path: str,
    dataset_uuid: uuid.UUID,
    session: Session,
    *,
    row_group: int | None = None,
    progress_desc: str | None = None,
    progress_total: int | None = None,
) -> None:
    """Streams a Parquet dataset's Reactions into the ORM tables with COPY.

    Builds ORM trees a batch at a time via ``from_proto``; each tree is assigned UUIDv7
    primary keys with foreign keys wired from the relationship metadata
    (``_collect_rows``), and the rows are streamed with ``COPY`` (``_copy_rows``) rather
    than the ORM unit of work. Peak memory is bounded to one row group plus
    ``_COPY_BATCH`` trees. The parent ``ord.dataset`` row identified by ``dataset_uuid``
    must already exist (see ``add_parquet_dataset_row``).

    Because rows are streamed with COPY rather than the unit of work, SQLAlchemy
    instrumentation does not run for this path: ``@validates`` methods and
    ``before_insert``/``after_insert`` event hooks on the mapper classes do not fire.
    ``add_dataset`` (used by ord-interface) still goes through the ORM, so any such hook
    must be reflected here as well to keep the two paths in sync.

    Args:
        path: Path to the Parquet-serialized Dataset.
        dataset_uuid: Surrogate key of the parent ``ord.dataset`` row (foreign key
            target).
        session: SQLAlchemy session.
        row_group: If set, only that Parquet row group is loaded; the unit of
            parallelism for sharded ingest. ``None`` loads every reaction in the file.
        progress_desc: If set, show a tqdm progress bar with this label (the
            single-process path passes it; shards leave it None so the pool-level bar is
            the only one).
        progress_total: Total reaction count for the progress bar's percentage/ETA.
    """
    reaction_child_class = Mappers.Dataset.reactions.mapper.class_
    dbapi_connection = session.connection().connection
    batch: list = []

    def copy_batch() -> None:
        if not batch:
            return
        rows_by_table: dict[str, tuple[Any, list]] = {}
        for reaction in batch:
            reaction_mapper = from_proto(reaction, mapper=reaction_child_class)
            reaction_mapper.dataset_id = dataset_uuid
            _collect_rows(reaction_mapper, rows_by_table)
        _copy_rows(dbapi_connection, rows_by_table)
        batch.clear()

    reactions: Any = parquet.DatasetView(path).iter_reactions(row_group=row_group)
    if progress_desc is not None:
        reactions = tqdm(
            reactions,
            total=progress_total,
            desc=progress_desc,
            unit="rxn",
            leave=False,
        )
    for _, reaction in reactions:
        batch.append(reaction)
        if len(batch) >= _COPY_BATCH:
            copy_batch()
    copy_batch()


def add_parquet_dataset_metadata(
    dataset_id: str, md5_hex: str, num_reactions: int, session: Session
) -> None:
    """Writes the ``public.datasets`` metadata row and populates ``submitted_at``.

    Written after the reactions are loaded, so the presence of this row marks a fully
    loaded dataset (``get_dataset_md5`` reads it); a crashed load leaves an
    ``ord.dataset`` row with no metadata row, which the next ingest deletes and reloads.

    Args:
        dataset_id: Dataset ID (``public.datasets`` primary key).
        md5_hex: Streaming MD5 of the Parquet file.
        num_reactions: Reaction count from the streaming pass.
        session: SQLAlchemy session.
    """
    session.add(
        DatasetMetadata(dataset_id=dataset_id, md5=md5_hex, num_reactions=num_reactions)
    )
    session.flush()
    set_submitted_at(dataset_id, session)


def add_parquet_dataset(path: str, session: Session) -> None:
    """Streams a Parquet Dataset into the ORM tables with COPY in one transaction.

    Composes ``add_parquet_dataset_row`` (search-index row), ``add_parquet_reactions``
    (the reactions), and ``add_parquet_dataset_metadata`` (the ``public.datasets``
    marker). Sharded ingest calls the same primitives across worker processes; here they
    run serially so the whole dataset lands atomically. Derived data is populated
    separately by ``update_derived_data``.

    Args:
        path: Path to the Parquet-serialized Dataset.
        session: SQLAlchemy session.
    """
    start = time.time()
    view = parquet.DatasetView(path)
    logger.debug(f"Streaming Parquet Dataset {view.dataset_id}")
    md5_hex = view.md5()
    num_reactions = len(view.reactions)
    dataset_uuid = uuid7()
    add_parquet_dataset_row(path, dataset_uuid, session)
    add_parquet_reactions(
        path,
        dataset_uuid,
        session,
        progress_desc=f"ingest {view.dataset_id}",
        progress_total=num_reactions,
    )
    add_parquet_dataset_metadata(view.dataset_id, md5_hex, num_reactions, session)
    logger.debug(
        f"add_parquet_dataset() took {time.time() - start:g}s "
        f"({num_reactions} reactions)"
    )


def get_dataset_md5(dataset_id: str, session: Session) -> str | None:
    """Returns the MD5 hash of a dataset's current version, or None if not stored."""
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


# The derived/RDKit passes scope to a dataset by walking from a compound up to
# ord.reaction, but ord.compound hangs off three parents, so no single join reaches all
# of them:
#   * a reaction input (reaction_input.reaction_id set);
#   * a workup input (that reaction_input has reaction_id NULL, reaction_workup_id set);
#   * a product measurement (compound.reaction_input_id NULL) -- authentic standards.
# Callers UNION the paths rather than joining through COALESCE of the parent keys, which
# is not sargable and would keep the planner off the foreign-key indexes. The paths are
# disjoint (a compound sets exactly one of reaction_input_id/product_measurement_id; a
# reaction_input belongs to a reaction xor a workup), so UNION ALL never double-counts.
_COMPOUND_REACTION_JOINS: tuple[str, ...] = (
    """
    JOIN ord.reaction_input ON ord.compound.reaction_input_id = ord.reaction_input.id
    JOIN ord.reaction ON ord.reaction_input.reaction_id = ord.reaction.id
    """,
    """
    JOIN ord.reaction_input ON ord.compound.reaction_input_id = ord.reaction_input.id
    JOIN ord.reaction_workup
        ON ord.reaction_input.reaction_workup_id = ord.reaction_workup.id
    JOIN ord.reaction ON ord.reaction_workup.reaction_id = ord.reaction.id
    """,
    """
    JOIN ord.product_measurement
        ON ord.compound.product_measurement_id = ord.product_measurement.id
    JOIN ord.product_compound
        ON ord.product_measurement.product_compound_id = ord.product_compound.id
    JOIN ord.reaction_outcome
        ON ord.product_compound.reaction_outcome_id = ord.reaction_outcome.id
    JOIN ord.reaction ON ord.reaction_outcome.reaction_id = ord.reaction.id
    """,
)

# ord.product_compound has a single parent, so one path suffices.
_PRODUCT_COMPOUND_REACTION_JOINS: tuple[str, ...] = (
    """
    JOIN ord.reaction_outcome
        ON ord.product_compound.reaction_outcome_id = ord.reaction_outcome.id
    JOIN ord.reaction ON ord.reaction_outcome.reaction_id = ord.reaction.id
    """,
)


def _resolve_dataset_pk(dataset_id: str, session: Session) -> Any:
    """Returns the ord.dataset surrogate key for ``dataset_id``.

    The derived and RDKit passes scope to a dataset by this key
    (``ord.reaction.dataset_id = :pk``) rather than joining ord.dataset on the string
    ``dataset_id``. A literal surrogate lets the planner use per-dataset row estimates
    (a small dataset is estimated small), so it picks a bounded nested-loop plan; a
    value reached only through a join collapses to the average dataset size and drives a
    whole-table scan instead. See update_rdkit_ids.

    Raises:
        sqlalchemy.exc.NoResultFound: If no dataset has ``dataset_id`` (it must already
            be ingested).
    """
    return session.execute(
        text("SELECT id FROM ord.dataset WHERE dataset_id = :dataset_id"),
        {"dataset_id": dataset_id},
    ).scalar_one()


def delete_derived_data(
    dataset_id: str, session: Session, *, shard: tuple[int, int] | None = None
) -> None:
    """Deletes a dataset's derived SMILES rows so the derived passes recompute them.

    Every derived pass is guarded by ``NOT EXISTS``, which makes re-running cheap but
    means an existing row is never revisited: a change to how SMILES are derived reaches
    only rows that do not exist yet. Deleting first is what turns the idempotent passes
    into a rebuild, so a database can be updated in place instead of reloaded.

    ``rdkit.mols`` and ``rdkit.reactions`` are deliberately left alone. They are shared,
    deduplicated by structure, and referenced by every dataset, so deleting this
    dataset's share of them would break others. The link columns live on the rows
    removed here, so the RDKit pass re-links from scratch, reusing structures that are
    already present and inserting the ones that are not. A rebuild that changes SMILES
    therefore leaves unreferenced structures behind; they cost space, not correctness.

    ``derived.reaction_classes`` is also left alone: classification is a separate opt-in
    pass, and discarding it here would silently make every rebuild pay for it again.

    Args:
        dataset_id: Dataset whose derived rows to delete.
        session: SQLAlchemy session.
        shard: ``(index, num_shards)`` to restrict to one disjoint hash-partition, using
            the same partitioning as the derived passes, so a sharded rebuild deletes
            exactly the rows that shard is about to recompute.
    """
    logger.debug(f"Deleting derived data for {dataset_id=}")
    start = time.time()
    dataset_pk = _resolve_dataset_pk(dataset_id, session)
    reaction_shard_sql, reaction_shard_params = _shard_predicate(
        "ord.reaction.id", shard
    )
    session.execute(
        text(f"""
        DELETE FROM derived.reaction_smiles
        WHERE derived.reaction_smiles.reaction_id IN (
            SELECT ord.reaction.id
            FROM ord.reaction
            WHERE ord.reaction.dataset_id = :id
              {reaction_shard_sql}
        )
        """),  # noqa: S608  (table/column names are internal constants, not user input)
        {"id": dataset_pk, **reaction_shard_params},
    )
    for compound_table, reaction_joins, derived_table, derived_id in (
        (
            "ord.compound",
            _COMPOUND_REACTION_JOINS,
            "derived.compound_smiles",
            "compound_id",
        ),
        (
            "ord.product_compound",
            _PRODUCT_COMPOUND_REACTION_JOINS,
            "derived.product_compound_smiles",
            "product_compound_id",
        ),
    ):
        shard_sql, shard_params = _shard_predicate(f"{compound_table}.id", shard)
        select_ids = "\nUNION ALL\n".join(
            f"""
                SELECT {compound_table}.id
                FROM {compound_table}
                {reaction_join}
                WHERE ord.reaction.dataset_id = :id
                  {shard_sql}
            """  # noqa: S608
            for reaction_join in reaction_joins
        )
        session.execute(
            text(f"""
            DELETE FROM {derived_table}
            WHERE {derived_table}.{derived_id} IN ({select_ids})
            """),  # noqa: S608
            {"id": dataset_pk, **shard_params},
        )
    logger.debug(f"Deleting derived data took {time.time() - start:g}s")


def update_derived_tables(
    dataset_id: str,
    session: Session,
    *,
    shard: tuple[int, int] | None = None,
    rederive: bool = False,
) -> None:
    """Populates the derived SMILES tables from the search index.

    Reaction SMILES come from the ground-truth proto in public.reactions; compound
    SMILES from each ord.compound row's reconstructed message. Idempotent (skips rows
    that already have a derived entry); runs before the RDKit pass, which reads the
    SMILES. Reactions and compounds are processed in batches of _DERIVED_BATCH so the
    heavy per-row work (proto deserialization, SMILES generation) and pending inserts
    stay within one batch.

    When ``shard`` is ``(index, num_shards)``, only that disjoint hash-partition of the
    dataset's reaction/compound ids is derived, so a large dataset can be split across
    worker processes (the derived-stage analog of the row-group sharding used for
    ingest). Idempotency makes the shards safe to run in any order or overlap.

    Args:
        dataset_id: Dataset to derive.
        session: SQLAlchemy session.
        shard: ``(index, num_shards)`` to derive one disjoint hash-partition.
        rederive: If True, delete this dataset's existing derived rows first, so rows
            written by an earlier version of the derivation are recomputed rather than
            skipped by the ``NOT EXISTS`` guards. Sharded runs delete only their own
            partition, so workers stay disjoint. See :func:`delete_derived_data`.
    """
    logger.debug(f"Updating derived tables for {dataset_id=}")
    start = time.time()
    if rederive:
        delete_derived_data(dataset_id, session, shard=shard)
    # dataset_id and the compound link columns live on the polymorphic child mappers, so
    # rows are scoped to the dataset via raw SQL (like the RDKit pass).
    dataset_pk = _resolve_dataset_pk(dataset_id, session)
    # Reaction SMILES from the served proto (no ORM objects loaded), keyed by
    # ord.reaction.id. Resolve the ids needing derivation in one indexed pass, then load
    # and parse the (large) protos in batches of _DERIVED_BATCH.
    reaction_shard_sql, reaction_shard_params = _shard_predicate(
        "ord.reaction.id", shard
    )
    reaction_ids = (
        session.execute(
            text(f"""
                SELECT ord.reaction.id
                FROM ord.reaction
                JOIN public.reactions
                    ON public.reactions.reaction_id = ord.reaction.reaction_id
                WHERE ord.reaction.dataset_id = :id
                  AND NOT EXISTS (
                      SELECT 1 FROM derived.reaction_smiles
                      WHERE derived.reaction_smiles.reaction_id = ord.reaction.id
                  )
                  {reaction_shard_sql}
                """),  # noqa: S608  (shard predicate is an internal constant fragment)
            {"id": dataset_pk, **reaction_shard_params},
        )
        .scalars()
        .all()
    )
    select_protos = text(
        """
        SELECT ord.reaction.id, public.reactions.proto
        FROM ord.reaction
        JOIN public.reactions
            ON public.reactions.reaction_id = ord.reaction.reaction_id
        WHERE ord.reaction.id IN :ids
        """
    ).bindparams(bindparam("ids", expanding=True))
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
            reaction_smiles = message_helpers.derived_reaction_smiles(
                reaction_pb2.Reaction.FromString(proto)
            )
            if reaction_smiles is None:
                logger.debug(f"No reaction SMILES for reaction id={reaction_id}")
                continue
            # Stored whole: an extension block is chemistry the source recorded, and
            # reaction_from_smiles reads one, so rdkit.reactions can be keyed by it.
            inserts.append(
                {"reaction_id": reaction_id, "reaction_smiles": reaction_smiles}
            )
        if inserts:
            session.execute(insert_reaction_smiles, inserts)
    _update_compound_smiles(
        session,
        dataset_pk=dataset_pk,
        compound_table="ord.compound",
        reaction_joins=_COMPOUND_REACTION_JOINS,
        compound_class=Mappers.Compound,
        derived_table="derived.compound_smiles",
        derived_id="compound_id",
        shard=shard,
    )
    _update_compound_smiles(
        session,
        dataset_pk=dataset_pk,
        compound_table="ord.product_compound",
        reaction_joins=_PRODUCT_COMPOUND_REACTION_JOINS,
        compound_class=Mappers.ProductCompound,
        derived_table="derived.product_compound_smiles",
        derived_id="product_compound_id",
        shard=shard,
    )
    logger.debug(f"Updating derived tables took {time.time() - start:g}s")


def _update_compound_smiles(
    session: Session,
    *,
    dataset_pk: int,
    compound_table: str,
    reaction_joins: tuple[str, ...],
    compound_class: Any,
    derived_table: str,
    derived_id: str,
    shard: tuple[int, int] | None = None,
) -> None:
    """Derives SMILES for one (product) compound table's not-yet-derived rows.

    The common case -- a compound with a stored SMILES identifier -- is served
    set-based: the first SMILES identifier for each compound in the batch is fetched
    from ord.compound_identifier in a single query and canonicalized with RDKit,
    avoiding a per-compound ORM load. Compounds without a stored SMILES but with some
    other structural identifier (INCHI/MOLBLOCK) fall back to reconstructing the message
    via to_proto and computing the SMILES from those identifiers; a second set-based
    query flags which batch ids carry any structural identifier at all, so a compound
    with only non-structural identifiers (e.g. NAME) -- which never yields a derived row
    and is therefore re-selected every run -- is skipped without the ORM load plus
    to_proto that is certain to raise. Ids are resolved up front and processed in
    batches of _DERIVED_BATCH so memory stays bounded. ``shard`` (index, num_shards)
    restricts to one disjoint hash-partition of this table's ids so large datasets can
    be split across workers.

    ``reaction_joins`` are the disjoint join paths from ``compound_table`` to
    ord.reaction; one id query runs per path and the results are concatenated. See
    _COMPOUND_REACTION_JOINS.
    """
    compound_shard_sql, compound_shard_params = _shard_predicate(
        f"{compound_table}.id", shard
    )
    select_ids = "\nUNION ALL\n".join(
        f"""
                SELECT {compound_table}.id
                FROM {compound_table}
                {reaction_join}
                WHERE ord.reaction.dataset_id = :id
                  AND NOT EXISTS (
                      SELECT 1 FROM {derived_table}
                      WHERE {derived_table}.{derived_id} = {compound_table}.id
                  )
                  {compound_shard_sql}
        """  # noqa: S608  (table/column names are internal constants, not user input)
        for reaction_join in reaction_joins
    )
    compound_ids = (
        session.execute(
            text(select_ids),
            {"id": dataset_pk, **compound_shard_params},
        )
        .scalars()
        .all()
    )
    # The first SMILES identifier (lowest id == proto order) for each requested
    # compound. The FK column is indexed, so one query per batch replaces an ORM load
    # plus lazy child fetches per compound. Interpolated names are internal constants
    # (not user input); see S608 below.
    select_smiles = text(f"""
        SELECT DISTINCT ON (ord.compound_identifier.{derived_id})
               ord.compound_identifier.{derived_id}, ord.compound_identifier.value
        FROM ord.compound_identifier
        WHERE ord.compound_identifier.{derived_id} IN :ids
          AND ord.compound_identifier.type = 'SMILES'
        ORDER BY ord.compound_identifier.{derived_id}, ord.compound_identifier.id
        """).bindparams(bindparam("ids", expanding=True))  # noqa: S608
    # Ids in the batch that carry a non-empty structural identifier, i.e. one of the
    # types smiles_from_compound can build a Mol from. A compound with only
    # non-structural identifiers (e.g. NAME) can never yield a SMILES, so the
    # reconstruction fallback below is skipped for it -- keeping the re-attempt of a
    # permanently underivable compound (no derived row is ever inserted, so it is
    # re-selected every run) cheap instead of paying an ORM load plus to_proto that is
    # certain to raise. The type list is bound from message_helpers, so it tracks the
    # loaders that back the reconstruction rather than a hand-maintained literal that
    # could drift from them; the names, not the enum numbers, are what this column
    # stores.
    select_structural = text(f"""
        SELECT DISTINCT ord.compound_identifier.{derived_id}
        FROM ord.compound_identifier
        WHERE ord.compound_identifier.{derived_id} IN :ids
          AND ord.compound_identifier.type IN :structural_types
          AND ord.compound_identifier.value <> ''
        """).bindparams(  # noqa: S608
        bindparam("ids", expanding=True),
        bindparam("structural_types", expanding=True),
    )
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
        derivable = {
            row[0]
            for row in session.execute(
                select_structural,
                {
                    "ids": batch_ids,
                    "structural_types": list(
                        message_helpers.STRUCTURAL_IDENTIFIER_TYPE_NAMES
                    ),
                },
            )
        }
        inserts = []
        for compound_id in batch_ids:
            value = stored_smiles.get(compound_id)
            # An absent or empty SMILES identifier both fall through to the
            # reconstruction path: an empty string is not a structure to canonicalize.
            if value:
                # Canonicalize the stored SMILES. Plain MolToSmiles, so this path
                # drops the enhanced stereochemistry smiles_from_compound keeps; see
                # #936, which is about closing that gap.
                mol = Chem.MolFromSmiles(value)
                if mol is None:
                    logger.debug(
                        f"Cannot parse SMILES for compound id={compound_id}: {value}"
                    )
                    continue
                smiles = Chem.MolToSmiles(mol)
            elif compound_id not in derivable:
                # Only non-structural identifiers, so smiles_from_compound has nothing
                # to read; skip the reconstruction (see select_structural).
                continue
            else:
                # No stored SMILES: reconstruct the message and derive from other
                # identifiers.
                compound = session.get(compound_class, compound_id)
                assert compound is not None  # Selected by id above.
                smiles = message_helpers.smiles_from_compound(
                    cast(
                        "reaction_pb2.Compound | reaction_pb2.ProductCompound",
                        to_proto(compound),
                    )
                )
                if smiles is None:
                    continue
            inserts.append({derived_id: compound_id, "smiles": smiles})
        if inserts:
            session.execute(insert_smiles, inserts)


def update_rdkit_tables(
    dataset_id: str, session: Session, *, shard: tuple[int, int] | None = None
) -> None:
    """Updates RDKit PostgreSQL cartridge data.

    ``shard`` (index, num_shards) restricts the inserts to one hash-partition of the
    SMILES values, so a dataset's RDKit work can be split across workers. Because every
    RDKit sub-step (here and in ``update_rdkit_ids``) partitions by the *same* SMILES
    hash, shard ``k`` inserts exactly the structures its links will reference --
    disjoint key sets across shards, so concurrent shards never collide on the shared
    ``rdkit.*`` unique indexes (datasets still run serially).

    Raises:
        RuntimeError: If the server is older than PostgreSQL 12 (the queries use AS
            MATERIALIZED).
    """
    logger.debug(f"Updating RDKit tables for {dataset_id=}")
    _check_server_version(session.connection())
    _update_rdkit_reactions(dataset_id, session, shard=shard)
    _update_rdkit_mols(dataset_id, session, shard=shard)


def _update_rdkit_reactions(
    dataset_id: str, session: Session, *, shard: tuple[int, int] | None = None
) -> None:
    """Updates the RDKit reactions table."""
    logger.debug("Updating RDKit reactions")
    start = time.time()
    dataset_pk = _resolve_dataset_pk(dataset_id, session)
    shard_sql, shard_params = _shard_predicate(
        "derived.reaction_smiles.reaction_smiles", shard
    )
    result = session.execute(
        text(f"""
            WITH scoped_reactions AS MATERIALIZED (
                -- Scope to the dataset via the surrogate key, carrying no
                -- rdkit_reaction_id predicate, so the planner reaches the reactions
                -- through ix_reaction_dataset_id instead of the whole-database
                -- "unlinked" partial index (which a backfill would rescan per (dataset,
                -- shard)). See _link_mol_ids and #895.
                SELECT ord.reaction.id AS id
                FROM ord.reaction
                WHERE ord.reaction.dataset_id = :dataset_pk
            )
            INSERT INTO rdkit.reactions (reaction_smiles, reaction)
            SELECT reaction_smiles, reaction
            FROM (
                SELECT reaction_smiles,
                    reaction_from_smiles(reaction_smiles::cstring) AS reaction
                FROM (
                    -- NOT EXISTS probes the unique reaction_smiles index per candidate
                    -- rather than EXCEPT-scanning all of rdkit.reactions; DISTINCT
                    -- dedupes within the dataset. reaction_smiles IS NOT NULL is
                    -- required: NOT EXISTS never matches a NULL, so without it a
                    -- no-SMILES reaction would re-insert a junk (NULL, NULL) row every
                    -- run. ON CONFLICT is a backstop; shards insert on disjoint
                    -- hash-partitions.
                    SELECT DISTINCT derived.reaction_smiles.reaction_smiles
                        FROM derived.reaction_smiles
                        JOIN scoped_reactions
                            ON scoped_reactions.id = derived.reaction_smiles.reaction_id
                        WHERE derived.reaction_smiles.rdkit_reaction_id IS NULL
                          AND derived.reaction_smiles.reaction_smiles IS NOT NULL
                          AND NOT EXISTS (
                              SELECT 1 FROM rdkit.reactions
                              WHERE rdkit.reactions.reaction_smiles
                                  = derived.reaction_smiles.reaction_smiles
                          )
                          {shard_sql}
                ) candidates
            ) computed
            -- reaction_from_smiles returns NULL for an unparseable reaction SMILES;
            -- skip those so we never insert a NULL-reaction row (mirrors the mol IS NOT
            -- NULL guard in _update_rdkit_mols).
            WHERE reaction IS NOT NULL
            ON CONFLICT (reaction_smiles) DO NOTHING
            """),  # noqa: S608  (shard predicate is an internal constant fragment)
        {"dataset_pk": dataset_pk, **shard_params},
    )
    logger.debug(
        f"Updating reactions took {time.time() - start:g}s "
        f"({cast(Any, result).rowcount} rows)"
    )


def _update_rdkit_mols(
    dataset_id: str, session: Session, *, shard: tuple[int, int] | None = None
) -> None:
    """Updates the RDKit mols table."""
    logger.debug("Updating RDKit mols")
    start = time.time()
    dataset_pk = _resolve_dataset_pk(dataset_id, session)
    shard_sql, shard_params = _shard_predicate("candidates.smiles", shard)
    # Scope each compound table to the dataset via the surrogate key, in a MATERIALIZED
    # CTE with no rdkit_mol_id predicate, so the planner reaches the compounds through
    # ix_reaction_dataset_id rather than the whole-database "unlinked rows" partial
    # index; new_smiles then joins the derived tables to these bounded id sets. See
    # _link_mol_ids and #895. UNION ALL is safe: the join paths select disjoint
    # compounds (see _COMPOUND_REACTION_JOINS).
    scoped_compound_ids = "\nUNION ALL\n".join(
        f"""
                SELECT ord.compound.id AS id
                FROM ord.compound
                {reaction_join}
                WHERE ord.reaction.dataset_id = :dataset_pk
        """  # noqa: S608  (the join path is an internal constant, not user input)
        for reaction_join in _COMPOUND_REACTION_JOINS
    )
    scoped_product_compound_ids = "\nUNION ALL\n".join(
        f"""
                SELECT ord.product_compound.id AS id
                FROM ord.product_compound
                {reaction_join}
                WHERE ord.reaction.dataset_id = :dataset_pk
        """  # noqa: S608  (the join path is an internal constant, not user input)
        for reaction_join in _PRODUCT_COMPOUND_REACTION_JOINS
    )
    result = session.execute(
        text(f"""
            WITH scoped_compounds AS MATERIALIZED (
                {scoped_compound_ids}
            ),
            scoped_product_compounds AS MATERIALIZED (
                {scoped_product_compound_ids}
            ),
            new_smiles AS MATERIALIZED (
                -- MATERIALIZED barrier: resolve the not-yet-linked SMILES absent from
                -- rdkit.mols FIRST (a cheap probe of the unique smiles index) so the
                -- expensive mol_from_smiles/morgan_*_fp calls below run only on the
                -- survivors. Without it the planner may run those functions on every
                -- candidate before the anti-join (~240s for a no-op dataset). UNION
                -- dedupes SMILES a dataset repeats across compounds.
                SELECT smiles
                FROM (
                    SELECT derived.compound_smiles.smiles
                        FROM derived.compound_smiles
                        JOIN scoped_compounds
                            ON scoped_compounds.id = derived.compound_smiles.compound_id
                        WHERE derived.compound_smiles.rdkit_mol_id IS NULL
                    UNION
                    SELECT derived.product_compound_smiles.smiles
                        FROM derived.product_compound_smiles
                        JOIN scoped_product_compounds
                            ON scoped_product_compounds.id
                                = derived.product_compound_smiles.product_compound_id
                        WHERE derived.product_compound_smiles.rdkit_mol_id IS NULL
                ) candidates
                WHERE smiles NOT LIKE '%[Ti+5]%'  -- See https://github.com/open-reaction-database/ord-schema/issues/672.
                  AND NOT EXISTS (
                      SELECT 1 FROM rdkit.mols
                      WHERE rdkit.mols.smiles = candidates.smiles
                  )
                  {shard_sql}
            )
            INSERT INTO rdkit.mols (smiles, mol, morgan_bfp, morgan_sfp)
            SELECT smiles, mol,
                morganbv_fp(mol) AS morgan_bfp, morgan_fp(mol) AS morgan_sfp
            FROM (
                SELECT smiles, mol_from_smiles(smiles::cstring) AS mol FROM new_smiles
            ) computed
            -- mol_from_smiles returns NULL for unparseable SMILES; skip those so we
            -- never insert a NULL-mol row (which ON CONFLICT would not catch for a
            -- genuinely new SMILES) or feed NULL to the fingerprints.
            WHERE mol IS NOT NULL
            ON CONFLICT (smiles) DO NOTHING
            """),  # noqa: S608  (shard predicate is an internal constant fragment)
        {"dataset_pk": dataset_pk, **shard_params},
    )
    logger.debug(
        f"Updating mols took {time.time() - start:g}s "
        f"({cast(Any, result).rowcount} rows)"
    )


def _link_mol_ids(
    session: Session,
    *,
    derived_table: str,
    derived_id_column: str,
    compound_table: str,
    reaction_joins: tuple[str, ...],
    dataset_pk: Any,
    shard: tuple[int, int] | None = None,
) -> int:
    """Links rdkit.mols ids into a derived compound-SMILES table for one dataset.

    The Compound and ProductCompound updates are identical apart from the tables and
    join paths, so they share this helper. Table/column arguments are trusted literals
    from ``update_rdkit_ids`` (never user input); ``dataset_pk`` is a bind param.

    One statement runs per join path. Each scopes to the dataset in a MATERIALIZED CTE
    keyed on the surrogate ``dataset_pk`` and carrying no rdkit_mol_id predicate, so the
    planner uses a per-dataset row estimate and reaches the compounds through
    ix_reaction_dataset_id; the outer UPDATE then touches only those compounds' derived
    rows. Without this a backfill -- every dataset's rows unlinked at once -- rescans
    the whole unlinked set once per (dataset, shard); see
    https://github.com/open-reaction-database/ord-schema/issues/895. The paths select
    disjoint compounds, so a row is linked by exactly one of them.

    Args:
        session: Active SQLAlchemy session.
        derived_table: Derived table to update (e.g. ``"derived.compound_smiles"``).
        derived_id_column: ``derived_table`` foreign key to ``compound_table.id``
            (e.g. ``"compound_id"``).
        compound_table: ORD compound table (e.g. ``"ord.compound"``).
        reaction_joins: Disjoint join paths from ``compound_table`` to ord.reaction
            (e.g. ``_COMPOUND_REACTION_JOINS``).
        dataset_pk: ord.dataset surrogate key to scope the update to (see
            _resolve_dataset_pk).
        shard: Optional ``(index, num_shards)`` partition (by the derived table's SMILES
            hash) to restrict the update to, matching the mol-insert partition so a
            shard links exactly the rows it inserted.

    Returns:
        The number of rows updated.
    """
    shard_sql, shard_params = _shard_predicate(f"{derived_table}.smiles", shard)
    rows = 0
    for reaction_join in reaction_joins:
        result = session.execute(
            text(f"""
                WITH scoped_compounds AS MATERIALIZED (
                    SELECT {compound_table}.id AS id
                    FROM {compound_table}
                    {reaction_join}
                    WHERE ord.reaction.dataset_id = :dataset_pk
                )
                UPDATE {derived_table}
                SET rdkit_mol_id = rdkit.mols.id
                FROM rdkit.mols, scoped_compounds
                WHERE {derived_table}.{derived_id_column} = scoped_compounds.id
                  AND rdkit.mols.smiles = {derived_table}.smiles
                  AND {derived_table}.rdkit_mol_id IS NULL
                  {shard_sql}
                """),  # noqa: S608
            # (table/column names are internal constants, not user input)
            {"dataset_pk": dataset_pk, **shard_params},
        )
        rows += cast(Any, result).rowcount
    return rows


def update_rdkit_ids(
    dataset_id: str, session: Session, *, shard: tuple[int, int] | None = None
) -> None:
    """Updates RDKit reaction and mol ID associations in the ORD tables.

    ``shard`` (index, num_shards) restricts the links to one SMILES-hash partition,
    matching the insert partition in ``update_rdkit_tables`` so shard ``k`` links
    exactly the structures it inserted.

    Raises:
        RuntimeError: If the server is older than PostgreSQL 12 (the queries use AS
            MATERIALIZED).
    """
    logger.debug("Updating RDKit ID associations")
    _check_server_version(session.connection())
    start = time.time()
    # Each UPDATE scopes to the dataset in a MATERIALIZED CTE keyed on the surrogate,
    # carrying no rdkit_*_id predicate, so the planner reaches the rows through
    # ix_reaction_dataset_id. The earlier flat ``UPDATE ... FROM ord.reaction WHERE
    # ... rdkit_*_id IS NULL`` let it drive from the "unlinked rows" partial index --
    # fine on an end-to-end load (the set drains dataset by dataset) but quadratic on a
    # backfill, where every dataset is unlinked at once and each (dataset, shard)
    # rescans the whole set. See #895.
    dataset_pk = _resolve_dataset_pk(dataset_id, session)
    reaction_shard_sql, reaction_shard_params = _shard_predicate(
        "derived.reaction_smiles.reaction_smiles", shard
    )
    reaction_result = session.execute(
        text(f"""
            WITH scoped_reactions AS MATERIALIZED (
                SELECT ord.reaction.id AS id
                FROM ord.reaction
                WHERE ord.reaction.dataset_id = :dataset_pk
            )
            UPDATE derived.reaction_smiles
            SET rdkit_reaction_id = rdkit.reactions.id
            FROM rdkit.reactions, scoped_reactions
            WHERE derived.reaction_smiles.reaction_id = scoped_reactions.id
              AND rdkit.reactions.reaction_smiles
                  = derived.reaction_smiles.reaction_smiles
              AND derived.reaction_smiles.rdkit_reaction_id IS NULL
              {reaction_shard_sql}
            """),  # noqa: S608  (shard predicate is an internal constant fragment)
        {"dataset_pk": dataset_pk, **reaction_shard_params},
    )
    reaction_rows = cast(Any, reaction_result).rowcount
    compound_rows = _link_mol_ids(
        session,
        derived_table="derived.compound_smiles",
        derived_id_column="compound_id",
        compound_table="ord.compound",
        reaction_joins=_COMPOUND_REACTION_JOINS,
        dataset_pk=dataset_pk,
        shard=shard,
    )
    product_compound_rows = _link_mol_ids(
        session,
        derived_table="derived.product_compound_smiles",
        derived_id_column="product_compound_id",
        compound_table="ord.product_compound",
        reaction_joins=_PRODUCT_COMPOUND_REACTION_JOINS,
        dataset_pk=dataset_pk,
        shard=shard,
    )
    logger.debug(
        f"Updating RDKit IDs took {time.time() - start:g}s "
        f"(reaction={reaction_rows}, compound={compound_rows}, "
        f"product_compound={product_compound_rows})"
    )
