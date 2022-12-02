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
    --overwrite             Overwrite existing datasets
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
from functools import partial
from glob import glob
from concurrent.futures import ProcessPoolExecutor

from docopt import docopt
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger
from ord_schema.message_helpers import load_message
from ord_schema.orm.database import add_dataset, add_rdkit, delete_dataset, get_connection_string
from ord_schema.proto.dataset_pb2 import Dataset

logger = get_logger(__name__)


def _add_dataset(filename: str, url: str, overwrite: bool) -> None:
    """Adds a single dataset to the database.

    Args:
        filename: Dataset filename.
        url: Database connection string.
        overwrite: If True, overwrite an existing dataset.
    """
    logger.info(f"Loading {filename}")
    start = time.time()
    dataset = load_message(filename, Dataset)
    logger.info(f"load_message() took {time.time() - start}s")
    engine = create_engine(url, future=True)
    with Session(engine) as session:
        if overwrite:
            delete_dataset(dataset.dataset_id, session)
        add_dataset(dataset, session)
        start = time.time()
        session.commit()
        logger.info(f"session.commit() took {time.time() - start}s")


def main(**kwargs):
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
    engine = create_engine(url, future=True)
    with Session(engine) as session:
        add_rdkit(session)
        start = time.time()
        session.commit()
        logger.info(f"session.commit() took {time.time() - start}s")


if __name__ == "__main__":
    main(**docopt(__doc__))
