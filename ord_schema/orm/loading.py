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

"""Staged loading of ORD datasets into the ORM database.

Loading is two independent stages: *ingest* writes the ``ord.*`` search index and
``public.*`` payload, and *derivation* writes the ``derived.*`` SMILES, RDKit links,
and (optionally) reaction classes. Either stage can run without the other; derivation
is idempotent, so the derived-only stage backfills or recomputes derived data over
already-ingested datasets. ``load_datasets`` orchestrates both stages over a glob of
dataset files; the per-dataset helpers (``ingest_dataset``, ``derive_dataset``,
``add_rdkit``) are exposed for callers that compose their own pipeline.
"""

import dataclasses
import math
import time
import uuid
from collections.abc import Callable, Iterable
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from contextlib import ExitStack
from functools import partial
from glob import glob
from hashlib import md5
from typing import TypeVar

from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session
from tqdm import tqdm
from uuid6 import uuid7

from ord_schema import parquet
from ord_schema.logging import get_logger
from ord_schema.message_helpers import load_message
from ord_schema.orm import database
from ord_schema.proto import dataset_pb2

_Item = TypeVar("_Item")
_Result = TypeVar("_Result")

logger = get_logger(__name__)

STAGES = ("ingest", "derived")

# SMILES derivation is sharded within a dataset so one large dataset does not pin a single worker
# (the derived-stage analog of row-group sharding for ingest). A dataset is split into one shard
# per this many reactions, capped, so small datasets stay a single shard.
_DERIVE_SHARD_SIZE = 50_000
_DERIVE_SHARD_CAP = 32


def dataset_id_for_file(filename: str) -> str:
    """Returns the dataset ID for a dataset file without ingesting it."""
    if filename.endswith(".parquet"):
        return parquet.load_metadata(filename).dataset_id
    return load_message(filename, dataset_pb2.Dataset).dataset_id


def ingest_dataset(filename: str, *, dsn: str, overwrite: bool) -> str:
    """Ingests a single dataset (``ord.*`` + ``public.*``); skips unchanged datasets.

    Derived data is populated separately by ``derive_dataset``/``add_rdkit``.

    Args:
        filename: Dataset filename. ``.parquet`` inputs are streamed via
            ``database.add_parquet_dataset``; other formats load the whole
            Dataset proto and call ``database.add_dataset``.
        dsn: Database connection string.
        overwrite: If True, update the dataset if the MD5 hash has changed.

    Returns:
        Dataset ID.

    Raises:
        ValueError: If the dataset already exists in the database and `overwrite` is not set.
    """
    if filename.endswith(".parquet"):
        logger.debug(f"Streaming {filename}")
        dataset_id = parquet.load_metadata(filename).dataset_id

        def compute_md5() -> str:
            md5_hex, _ = parquet.streaming_md5(filename)
            return md5_hex

        def insert(session: Session) -> None:
            database.add_parquet_dataset(filename, session)
    else:
        logger.debug(f"Loading {filename}")
        dataset = load_message(filename, dataset_pb2.Dataset)
        dataset_id = dataset.dataset_id

        def compute_md5() -> str:
            return md5(
                dataset.SerializeToString(deterministic=True), usedforsecurity=False
            ).hexdigest()

        def insert(session: Session) -> None:
            database.add_dataset(dataset, session)

    # NOTE(skearnes): Multiprocessing is hard to get right for shared connection pools, so we don't even try; see
    # https://docs.sqlalchemy.org/en/20/core/pooling.html#using-connection-pools-with-multiprocessing-or-os-fork.
    # Each call owns its engine and disposes it so pooled connections don't accumulate in reused pool workers.
    engine = create_engine(dsn)
    try:
        with Session(engine) as session:
            with session.begin():
                existing_md5 = database.get_dataset_md5(dataset_id, session)
            if existing_md5 is not None:
                if compute_md5() != existing_md5:
                    if not overwrite:
                        raise ValueError(
                            f"`overwrite` is required when a dataset already exists: {dataset_id}"
                        )
                    logger.debug(f"existing dataset {dataset_id} changed; updating")
                    with session.begin():
                        database.delete_dataset(dataset_id, session)
                else:
                    logger.debug(f"existing dataset {dataset_id} unchanged; skipping")
                    return dataset_id
            start = time.time()
            with session.begin():
                insert(session)
            logger.debug(f"insert took {time.time() - start:g}s")
    finally:
        engine.dispose()
    return dataset_id


def derive_dataset(dataset_id: str, *, dsn: str, classify_reactions: bool) -> str:
    """Runs the parallel-safe derived passes (SMILES + optional classification).

    The RDKit linking pass is excluded here and run serially in ``add_rdkit`` to avoid
    deadlocks. Idempotent, so it is safe to re-run over already-derived datasets.

    Args:
        dataset_id: Dataset to derive.
        dsn: Database connection string.
        classify_reactions: Whether to assign reaction class/name labels.

    Returns:
        Dataset ID.
    """
    engine = create_engine(dsn)
    try:
        with Session(engine) as session, session.begin():
            database.update_derived_data(
                dataset_id,
                session,
                rdkit_cartridge=False,  # Done serially in add_rdkit() to avoid deadlocks.
                classify_reactions=classify_reactions,
            )
    finally:
        engine.dispose()
    return dataset_id


def _derive_shard_items(
    dataset_ids: Iterable[str], dsn: str
) -> list[tuple[str, int, int]]:
    """Builds ``(dataset_id, shard_index, num_shards)`` SMILES-derivation work items.

    A dataset is split into ``ceil(num_reactions / _DERIVE_SHARD_SIZE)`` shards (capped at
    ``_DERIVE_SHARD_CAP``, at least 1), so a large dataset fans out across the pool while small
    ones stay a single shard. Datasets without a size (no metadata row) default to one shard.
    """
    items: list[tuple[str, int, int]] = []
    engine = create_engine(dsn)
    try:
        with Session(engine) as session:
            for dataset_id in dataset_ids:
                try:
                    size = database.get_dataset_size(dataset_id, session)
                except ValueError:
                    size = 0
                num_shards = max(
                    1, min(_DERIVE_SHARD_CAP, math.ceil(size / _DERIVE_SHARD_SIZE))
                )
                items.extend(
                    (dataset_id, shard_index, num_shards)
                    for shard_index in range(num_shards)
                )
    finally:
        engine.dispose()
    return items


def _derive_smiles_shard(
    item: tuple[str, int, int], *, dsn: str
) -> tuple[str, int, int]:
    """Derives one hash-partition of a dataset's SMILES (the parallel-safe, shardable pass)."""
    dataset_id, shard_index, num_shards = item
    engine = create_engine(dsn)
    try:
        with Session(engine) as session, session.begin():
            database.update_derived_tables(
                dataset_id, session, shard=(shard_index, num_shards)
            )
    finally:
        engine.dispose()
    return item


def _classify_dataset(dataset_id: str, *, dsn: str) -> str:
    """Assigns reaction class/name labels for a dataset (SMILES already derived by the shard pass).

    Classification only -- it does not re-derive SMILES, so a failed SMILES shard stays visibly
    incomplete rather than being silently backfilled here (keeping shard-failure semantics the same
    whether or not classification is enabled).
    """
    engine = create_engine(dsn)
    try:
        with Session(engine) as session, session.begin():
            database.classify_dataset(dataset_id, session)
    finally:
        engine.dispose()
    return dataset_id


def add_rdkit(engine: Engine, dataset_id: str) -> None:
    """Populates and links the RDKit cartridge tables for a dataset."""
    with Session(engine) as session:
        with session.begin():
            database.update_rdkit_tables(dataset_id, session)
        with session.begin():
            database.update_rdkit_ids(dataset_id, session)


def _run_parallel(
    func: Callable[[_Item], _Result],
    items: Iterable[_Item],
    *,
    n_jobs: int,
    desc: str,
) -> tuple[list[_Result], list[_Item]]:
    """Runs ``func`` over ``items`` with ``n_jobs`` workers.

    Returns ``(results, failures)``: the successful return values and the inputs whose
    call raised.
    """
    items = list(items)
    results: list[_Result] = []
    failures: list[_Item] = []
    with ExitStack() as stack:
        if n_jobs > 1:
            executor = stack.enter_context(ProcessPoolExecutor(n_jobs))
        else:
            executor = stack.enter_context(ThreadPoolExecutor(n_jobs))
        futures = {executor.submit(func, item): item for item in items}
        for future in tqdm(as_completed(futures), total=len(futures), desc=desc):
            try:
                results.append(future.result())
            except Exception:
                item = futures[future]
                failures.append(item)
                logger.exception(f"{desc} failed for {item}")
    return results, failures


@dataclasses.dataclass(frozen=True)
class _ParquetPlan:
    """Per-dataset outcome of the prep phase of sharded Parquet ingest.

    ``needs_load`` is False for datasets skipped as unchanged; those carry only ``dataset_id``
    (for the derived stage). Datasets that need loading also carry the surrogate ``dataset_uuid``
    the shard phase wires reactions to and the ``num_row_groups`` to shard over.
    """

    filename: str
    dataset_id: str
    md5_hex: str
    num_reactions: int
    needs_load: bool
    dataset_uuid: uuid.UUID | None = None
    num_row_groups: int = 0


def _prep_parquet_dataset(filename: str, *, dsn: str, overwrite: bool) -> _ParquetPlan:
    """Prep phase: decide skip/overwrite and insert the empty ``ord.dataset`` row.

    Skips datasets whose streaming MD5 is unchanged; otherwise deletes any existing full or
    partial rows and inserts a fresh search-index row whose surrogate id the shard phase
    references. The ``public.datasets`` marker is written only after all shards succeed
    (``_finalize_parquet_dataset``), so a crashed load leaves no marker and is cleanly redone.

    Raises:
        ValueError: If the dataset exists with changed content and ``overwrite`` is not set.
    """
    footer = parquet.load_footer(filename)
    dataset_id = footer.dataset.dataset_id
    md5_hex, num_reactions = parquet.streaming_md5(filename)
    engine = create_engine(dsn)
    try:
        with Session(engine) as session, session.begin():
            existing_md5 = database.get_dataset_md5(dataset_id, session)
            if existing_md5 is not None:
                if existing_md5 == md5_hex:
                    logger.debug(f"existing dataset {dataset_id} unchanged; skipping")
                    return _ParquetPlan(
                        filename, dataset_id, md5_hex, num_reactions, needs_load=False
                    )
                if not overwrite:
                    raise ValueError(
                        f"`overwrite` is required when a dataset already exists: {dataset_id}"
                    )
                logger.debug(f"existing dataset {dataset_id} changed; reloading")
            # Wipe any full or partial prior rows, then insert a fresh search-index row.
            database.delete_dataset(dataset_id, session)
            dataset_uuid = uuid7()
            database.add_parquet_dataset_row(filename, dataset_uuid, session)
    finally:
        engine.dispose()
    return _ParquetPlan(
        filename,
        dataset_id,
        md5_hex,
        num_reactions,
        needs_load=True,
        dataset_uuid=dataset_uuid,
        num_row_groups=footer.num_row_groups,
    )


def _ingest_parquet_shard(item: tuple[str, uuid.UUID, int], *, dsn: str) -> None:
    """Shard phase: COPY-load one Parquet row group's reactions in its own transaction."""
    filename, dataset_uuid, row_group = item
    engine = create_engine(dsn)
    try:
        with Session(engine) as session, session.begin():
            database.add_parquet_reactions(
                filename, dataset_uuid, session, row_group=row_group
            )
    finally:
        engine.dispose()


def _finalize_parquet_dataset(plan: _ParquetPlan, *, dsn: str) -> str:
    """Finalize phase: write the ``public.datasets`` completeness marker and ``submitted_at``."""
    engine = create_engine(dsn)
    try:
        with Session(engine) as session, session.begin():
            database.add_parquet_dataset_metadata(
                plan.dataset_id, plan.md5_hex, plan.num_reactions, session
            )
    finally:
        engine.dispose()
    return plan.dataset_id


def _ingest_parquet_sharded(
    filenames: list[str], *, dsn: str, overwrite: bool, n_jobs: int
) -> tuple[list[str], list[str]]:
    """Ingests Parquet datasets by sharding their row groups across a worker pool.

    Three phases, each a barrier so the next reads the previous one's committed rows: prep inserts
    each dataset's search-index row, shards COPY the reactions row-group by row-group, and finalize
    writes the ``public.datasets`` marker for datasets whose shards all succeeded. The row group is
    the unit of parallelism, so one large dataset's shards fill the pool instead of pinning a
    single worker. A dataset with any failed shard is left unmarked (skipped in finalize) and
    reported, so a re-run redoes it from scratch.

    Returns ``(dataset_ids, failures)``: ids of skipped and fully loaded datasets (for the derived
    stage), and the filenames that failed in any phase.
    """
    failures: list[str] = []
    plans, prep_failures = _run_parallel(
        partial(_prep_parquet_dataset, dsn=dsn, overwrite=overwrite),
        filenames,
        n_jobs=n_jobs,
        desc="Prep",
    )
    failures.extend(prep_failures)
    skipped_ids = [plan.dataset_id for plan in plans if not plan.needs_load]
    to_load = [plan for plan in plans if plan.needs_load]
    shard_items: list[tuple[str, uuid.UUID, int]] = []
    for plan in to_load:
        if plan.dataset_uuid is None:
            # Prep sets dataset_uuid for every needs_load plan; enforce the contract explicitly
            # (a bare assert is dropped under -O and would surface later as an obscure error).
            raise ValueError(
                f"prep returned needs_load with no dataset_uuid: {plan.filename}"
            )
        shard_items.extend(
            (plan.filename, plan.dataset_uuid, row_group)
            for row_group in range(plan.num_row_groups)
        )
    _, shard_failures = _run_parallel(
        partial(_ingest_parquet_shard, dsn=dsn),
        shard_items,
        n_jobs=n_jobs,
        desc="Ingest",
    )
    failed_uuids = {dataset_uuid for _, dataset_uuid, _ in shard_failures}
    failures.extend(sorted({filename for filename, _, _ in shard_failures}))
    finalize_plans = [plan for plan in to_load if plan.dataset_uuid not in failed_uuids]
    finalized_ids, finalize_failures = _run_parallel(
        partial(_finalize_parquet_dataset, dsn=dsn),
        finalize_plans,
        n_jobs=n_jobs,
        desc="Finalize",
    )
    failures.extend(plan.filename for plan in finalize_failures)
    return skipped_ids + finalized_ids, failures


def load_datasets(
    pattern: str,
    dsn: str,
    *,
    stages: Iterable[str] = STAGES,
    overwrite: bool = False,
    classify_reactions: bool = False,
    n_jobs: int = 1,
) -> None:
    """Runs the selected stages over the datasets matching ``pattern``.

    The parallel-safe work (ingest, SMILES/classification) runs in a worker pool; the
    RDKit linking pass runs serially to avoid deadlocks. When ingest is skipped, dataset
    IDs are resolved from ``pattern`` so derivation runs over the already-ingested rows.

    Args:
        pattern: Glob for dataset filenames.
        dsn: Database connection string.
        stages: Stages to run (any of ``STAGES``); defaults to all.
        overwrite: If True, update changed datasets during ingest.
        classify_reactions: If True, assign reaction class/name labels during derivation.
        n_jobs: Number of parallel workers.

    Raises:
        ValueError: If ``stages`` is empty or names an unknown stage.
        ImportError: If ``classify_reactions`` is set without the ``reaction-class`` extra.
        RuntimeError: If any dataset failed in any stage; the inputs are the message.
    """
    stages = tuple(stages)
    unknown = set(stages) - set(STAGES)
    if unknown:
        raise ValueError(
            f"unknown stages {sorted(unknown)}; choose from {list(STAGES)}"
        )
    if not stages:
        raise ValueError(f"stages must include at least one of {list(STAGES)}")
    if classify_reactions and database.update_reaction_classes is None:
        raise ImportError(
            "reaction classification requires the 'reaction-class' extra: "
            "pip install ord-schema[reaction-class]"
        )
    filenames = sorted(glob(pattern))
    failures: list[str] = []

    if "ingest" in stages:
        logger.info("Ingesting datasets")
        dataset_ids: list[str] = []
        # Parquet datasets are sharded by row group across the pool so one large dataset does not
        # pin a single worker; the atomic single-process path handles n_jobs == 1 and non-parquet
        # inputs (which cannot be sharded).
        if n_jobs > 1:
            parquet_files = [f for f in filenames if f.endswith(".parquet")]
            other_files = [f for f in filenames if not f.endswith(".parquet")]
        else:
            parquet_files, other_files = [], filenames
        if parquet_files:
            ids, stage_failures = _ingest_parquet_sharded(
                parquet_files, dsn=dsn, overwrite=overwrite, n_jobs=n_jobs
            )
            dataset_ids.extend(ids)
            failures.extend(stage_failures)
        if other_files:
            ids, stage_failures = _run_parallel(
                partial(ingest_dataset, dsn=dsn, overwrite=overwrite),
                other_files,
                n_jobs=n_jobs,
                desc="Ingest",
            )
            dataset_ids.extend(ids)
            failures.extend(stage_failures)
    else:
        dataset_ids = [dataset_id_for_file(filename) for filename in filenames]

    if "derived" in stages:
        logger.info("Deriving SMILES")
        # Shard SMILES derivation within datasets so one large dataset does not pin a single
        # worker; the pool sees a flat (dataset, shard) work list across all datasets.
        shard_items = _derive_shard_items(dataset_ids, dsn)
        _, shard_failures = _run_parallel(
            partial(_derive_smiles_shard, dsn=dsn),
            shard_items,
            n_jobs=n_jobs,
            desc="Derive",
        )
        failures.extend(sorted({dataset_id for dataset_id, _, _ in shard_failures}))
        if classify_reactions:
            logger.info("Classifying reactions")
            _, classify_failures = _run_parallel(
                partial(_classify_dataset, dsn=dsn),
                dataset_ids,
                n_jobs=n_jobs,
                desc="Classify",
            )
            failures.extend(classify_failures)
        logger.info("Adding RDKit functionality")
        engine = create_engine(dsn)
        try:
            for dataset_id in tqdm(dataset_ids, desc="RDKit"):
                try:
                    add_rdkit(
                        engine, dataset_id
                    )  # NOTE(skearnes): Do this serially to avoid deadlocks.
                except Exception:
                    failures.append(dataset_id)
                    logger.exception(
                        f"Adding RDKit functionality for {dataset_id} failed"
                    )
        finally:
            engine.dispose()

    if failures:
        raise RuntimeError(failures)
