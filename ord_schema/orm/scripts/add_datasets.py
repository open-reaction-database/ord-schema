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

"""Adds datasets to the ORM database.

Ingest (``ord.*`` search index + ``public.*`` payload) and derivation (``derived.*``
SMILES, RDKit links, and reaction classes) are separate stages selected with
``--stages``; either can run without the other. Derivation is idempotent, so the
derived-only stage backfills or recomputes derived data over already-ingested datasets.
"""

import argparse
import logging
import os
import time
from collections.abc import Callable, Iterable
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from contextlib import ExitStack
from functools import partial
from glob import glob
from hashlib import md5

from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session
from tqdm import tqdm

from ord_schema import parquet
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.message_helpers import load_message
from ord_schema.orm import database
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)

_STAGES = ("ingest", "derived")


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
    # NOTE(skearnes): Multiprocessing is hard to get right for shared connection pools, so we don't even try; see
    # https://docs.sqlalchemy.org/en/20/core/pooling.html#using-connection-pools-with-multiprocessing-or-os-fork.
    engine = create_engine(dsn)
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
    with Session(engine) as session, session.begin():
        database.update_derived_data(
            dataset_id,
            session,
            rdkit_cartridge=False,  # Done serially in add_rdkit() to avoid deadlocks.
            classify_reactions=classify_reactions,
        )
    return dataset_id


def add_rdkit(engine: Engine, dataset_id: str) -> None:
    """Updates RDKit tables."""
    with Session(engine) as session:
        with session.begin():
            database.update_rdkit_tables(dataset_id, session)
        with session.begin():
            database.update_rdkit_ids(dataset_id, session)


def _run_parallel(
    func: Callable[[str], str],
    items: Iterable[str],
    *,
    n_jobs: int,
    desc: str,
    failures: list[str],
) -> list[str]:
    """Runs ``func`` over ``items`` with ``n_jobs`` workers.

    Returns the successful results; appends each failed input to ``failures``.
    """
    items = list(items)
    results = []
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
    return results


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Add datasets to the ORM database")
    parser.add_argument(
        "--pattern", required=True, help="Pattern for dataset filenames"
    )
    parser.add_argument(
        "--stages",
        default=",".join(_STAGES),
        help="Comma-separated stages to run: 'ingest' (ord.*/public.*) and/or 'derived' "
        "(derived.* SMILES, RDKit links, reaction classes). Derived-only runs over "
        "already-ingested datasets matching --pattern.",
    )
    parser.add_argument(
        "--classify_reactions",
        action="store_true",
        help="Assign reaction class/name labels in the derived stage "
        "(requires the 'reaction-class' extra)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Update changed datasets"
    )
    parser.add_argument("--dsn", default=None, help="Postgres connection string")
    parser.add_argument("--database", default="orm", help="Database")
    parser.add_argument("--username", default="postgres", help="Database username")
    parser.add_argument("--password", default=None, help="Database password")
    parser.add_argument("--host", default="localhost", help="Database host")
    parser.add_argument("--port", type=int, default=5432, help="Database port")
    parser.add_argument(
        "--n_jobs", type=int, default=1, help="Number of parallel workers"
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Runs the selected ingest and/or derived stages over the matching datasets."""
    silence_rdkit_logs()
    if args.debug:
        get_logger(database.__name__, level=logging.DEBUG)
    stages = [stage.strip() for stage in args.stages.split(",") if stage.strip()]
    unknown = set(stages) - set(_STAGES)
    if unknown:
        raise SystemExit(
            f"Unknown --stages values {sorted(unknown)}; choose from {list(_STAGES)}"
        )
    if not stages:
        raise SystemExit(f"--stages must select at least one of {list(_STAGES)}")
    if args.classify_reactions and database.update_reaction_classes is None:
        raise SystemExit(
            "--classify_reactions requires the 'reaction-class' extra: "
            "pip install ord-schema[reaction-class]"
        )
    if args.dsn:
        dsn = args.dsn
    else:
        dsn = database.get_connection_string(
            database=args.database,
            username=args.username,
            password=args.password or os.environ["PGPASSWORD"],
            host=args.host,
            port=args.port,
        )
    filenames = sorted(glob(args.pattern))
    failures: list[str] = []

    if "ingest" in stages:
        logger.info("Ingesting datasets")
        dataset_ids = _run_parallel(
            partial(ingest_dataset, dsn=dsn, overwrite=args.overwrite),
            filenames,
            n_jobs=args.n_jobs,
            desc="Ingest",
            failures=failures,
        )
    else:
        dataset_ids = [dataset_id_for_file(filename) for filename in filenames]

    if "derived" in stages:
        logger.info("Deriving SMILES and reaction classes")
        dataset_ids = _run_parallel(
            partial(
                derive_dataset, dsn=dsn, classify_reactions=args.classify_reactions
            ),
            dataset_ids,
            n_jobs=args.n_jobs,
            desc="Derive",
            failures=failures,
        )
        logger.info("Adding RDKit functionality")
        engine = create_engine(dsn)
        for dataset_id in tqdm(dataset_ids, desc="RDKit"):
            try:
                add_rdkit(
                    engine, dataset_id
                )  # NOTE(skearnes): Do this serially to avoid deadlocks.
            except Exception:
                failures.append(dataset_id)
                logger.exception(f"Adding RDKit functionality for {dataset_id} failed")

    if failures:
        raise RuntimeError(failures)


if __name__ == "__main__":
    main(parse_args())
