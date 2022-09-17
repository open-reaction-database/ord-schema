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

"""Adds a dataset to the ORM database.

Usage:
    add_dataset.py --pattern=<str> [options]
    add_dataset.py -h | --help

    --pattern=<str>         Pattern for dataset filenames

Options:
    --database=<str>        Database [default: orm]
    --username=<str>        Database username [default: postgres]
    --password=<str>        Database password
    --host=<str>            Database host [default: localhost]
    --port=<int>            Database port [default: 5432]
"""
import os
import time
from docopt import docopt
from glob import glob

from sqlalchemy import create_engine

from ord_schema.logging import get_logger
from ord_schema.message_helpers import load_message
from ord_schema.orm.database import add_datasets, add_rdkit, get_connection_string
from ord_schema.proto.dataset_pb2 import Dataset

logger = get_logger(__name__)


def main(**kwargs):
    engine = create_engine(
        get_connection_string(
            database=kwargs["--database"],
            username=kwargs["--username"],
            password=kwargs["--password"] or os.environ["PGPASSWORD"],
            host=kwargs["--host"],
            port=int(kwargs["--port"]),
        ),
        future=True,
    )
    for filename in sorted(glob(kwargs["--pattern"])):
        logger.info(f"Loading {filename}")
        t0 = time.time()
        dataset = load_message(filename, Dataset)
        logger.info(f"load_message() took {time.time() - t0}s")
        add_datasets([dataset], engine=engine)
    logger.info("Updating RDKit functionality")
    add_rdkit(engine)


if __name__ == "__main__":
    main(**docopt(__doc__))
