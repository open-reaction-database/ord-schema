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
"""
import os
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from glob import glob
from hashlib import md5

from docopt import docopt
from rdkit import RDLogger
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger
from ord_schema.message_helpers import load_message
from ord_schema.orm.database import (
    add_dataset,
    delete_dataset,
    get_connection_string,
    get_dataset_md5,
    update_rdkit_ids,
    update_rdkit_tables,
)
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def _add_dataset(filename: str, url: str, overwrite: bool) -> None:
    """Adds a single dataset to the database.

    Args:
        filename: Dataset filename.
        url: Database connection string.
        overwrite: If True, update the dataset if the MD5 hash has changed.

    Raises:
        ValueError: If the dataset already exists in the database and `overwrite` is not set.
    """
    logger.info(f"Loading {filename}")
    start = time.time()
    dataset = load_message(filename, dataset_pb2.Dataset)
    logger.info(f"load_message() took {time.time() - start:g}s")
    engine = create_engine(url, future=True)
    with Session(engine) as session:
        dataset_md5 = get_dataset_md5(dataset.dataset_id, session)
        if dataset_md5 is not None:
            this_md5 = md5(dataset.SerializeToString(deterministic=True)).hexdigest()
            if this_md5 != dataset_md5:
                if not overwrite:
                    raise ValueError(f"`overwrite` is required when a dataset already exists: {dataset.dataset_id}")
                logger.info(f"existing dataset {dataset.dataset_id} changed; updating")
                delete_dataset(dataset.dataset_id, session)
            else:
                logger.info(f"existing dataset {dataset.dataset_id} unchanged; skipping")
                return
        add_dataset(dataset, session)
        session.flush()
        update_rdkit_tables(dataset.dataset_id, session=session)
        session.flush()
        update_rdkit_ids(dataset.dataset_id, session=session)
        start = time.time()
        session.commit()
        logger.info(f"session.commit() took {time.time() - start:g}s")


def main(**kwargs):
    RDLogger.DisableLog("rdApp.*")
    if kwargs.get("--url"):
        url = kwargs["--url"]
    else:
        url = get_connection_string(
            database=kwargs["--database"],
            username=kwargs["--username"],
            password=kwargs["--password"] or os.environ["PGPASSWORD"],
            host=kwargs["--host"],
            port=int(kwargs["--port"]),
        )
    function = partial(_add_dataset, url=url, overwrite=kwargs["--overwrite"])
    filenames = glob(kwargs["--pattern"])
    with ProcessPoolExecutor(max_workers=int(kwargs["--n_jobs"])) as executor:
        for _ in executor.map(function, filenames):
            pass  # Must iterate over results to raise exceptions.


if __name__ == "__main__":
    main(**docopt(__doc__))
