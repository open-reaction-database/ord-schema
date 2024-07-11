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

Usage:
    add_datasets.py --pattern=<str> [options]
    add_datasets.py -h | --help

Options:
    --pattern=<str>         Pattern for dataset filenames
    --overwrite             Update changed datasets
    --url=<str>             Postgres connection string
    --database=<str>        Database [default: orm]
    --username=<str>        Database username [default: postgres]
    --password=<str>        Database password
    --host=<str>            Database host [default: localhost]
    --port=<int>            Database port [default: 5432]
    --n_jobs=<int>          Number of parallel workers [default: 1]
    --debug                 Enable debug logging.
"""
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob
from hashlib import md5

from docopt import docopt
from rdkit import RDLogger
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session
from tqdm import tqdm

from ord_schema.logging import get_logger
from ord_schema.message_helpers import load_message
from ord_schema.orm import database
from ord_schema.proto import dataset_pb2

engine: Engine = None
logger = get_logger(__name__)


def initializer():
    """Initializer for child processes."""
    # See https://docs.sqlalchemy.org/en/20/core/pooling.html#using-connection-pools-with-multiprocessing-or-os-fork.
    engine.dispose(close=False)


def add_dataset(filename: str, overwrite: bool) -> str:
    """Adds a single dataset to the database.

    Args:
        filename: Dataset filename.
        overwrite: If True, update the dataset if the MD5 hash has changed.

    Returns:
        Dataset ID.

    Raises:
        ValueError: If the dataset already exists in the database and `overwrite` is not set.
    """
    logger.debug(f"Loading {filename}")
    dataset = load_message(filename, dataset_pb2.Dataset)
    with Session(engine) as session:
        with session.begin():
            dataset_md5 = database.get_dataset_md5(dataset.dataset_id, session)
        if dataset_md5 is not None:
            this_md5 = md5(dataset.SerializeToString(deterministic=True)).hexdigest()
            if this_md5 != dataset_md5:
                if not overwrite:
                    raise ValueError(f"`overwrite` is required when a dataset already exists: {dataset.dataset_id}")
                logger.debug(f"existing dataset {dataset.dataset_id} changed; updating")
                with session.begin():
                    database.delete_dataset(dataset.dataset_id, session)
            else:
                logger.debug(f"existing dataset {dataset.dataset_id} unchanged; skipping")
                return dataset.dataset_id
        with session.begin():
            database.add_dataset(dataset, session, rdkit_cartridge=False)  # Do this separately in add_rdkit().
    return dataset.dataset_id


def add_rdkit(dataset_id: str) -> None:
    """Updates RDKit tables."""
    with Session(engine) as session:
        with session.begin():
            database.update_rdkit_tables(dataset_id, session)
        with session.begin():
            database.update_rdkit_ids(dataset_id, session)


def main(**kwargs):
    RDLogger.DisableLog("rdApp.*")
    if kwargs["--debug"]:
        get_logger(database.__name__, level=logging.DEBUG)
    if kwargs["--url"]:
        url = kwargs["--url"]
    else:
        url = database.get_connection_string(
            database=kwargs["--database"],
            username=kwargs["--username"],
            password=kwargs["--password"] or os.environ["PGPASSWORD"],
            host=kwargs["--host"],
            port=int(kwargs["--port"]),
        )
    global engine  # pylint: disable=global-statement
    engine = create_engine(url)
    filenames = sorted(glob(kwargs["--pattern"]))
    with ProcessPoolExecutor(initializer=initializer, max_workers=int(kwargs["--n_jobs"])) as executor:
        logger.info("Adding datasets")
        futures = {}
        for filename in filenames:
            future = executor.submit(add_dataset, filename=filename, overwrite=kwargs["--overwrite"])
            futures[future] = filename
        dataset_ids = []
        failures = []
        for future in tqdm(as_completed(futures), total=len(futures)):
            try:
                dataset_ids.append(future.result())
            except Exception as error:  # pylint: disable=broad-exception-caught
                filename = futures[future]
                failures.append(filename)
                logger.error(f"Adding dataset {filename} failed: {error}")
        logger.info("Adding RDKit functionality")
        for dataset_id in tqdm(dataset_ids):
            try:
                add_rdkit(dataset_id)  # NOTE(skearnes): Do this serially to avoid deadlocks.
            except Exception as error:  # pylint: disable=broad-exception-caught
                failures.append(dataset_id)
                logger.error(f"Adding RDKit functionality for {dataset_id} failed: {error}")
        if failures:
            raise RuntimeError(failures)


if __name__ == "__main__":
    main(**docopt(__doc__))
