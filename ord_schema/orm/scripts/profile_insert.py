# Copyright 2024 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Profiles the ORM.

Usage:
    profile.py --pattern=<str> --output=<str>
    profile.py -h | --help

Options:
    --pattern=<str>         Pattern for dataset filenames
    --output=<str>          Output filename for profile information.
"""
import cProfile
import time
from glob import glob

from docopt import docopt
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger
from ord_schema.message_helpers import load_message
from ord_schema.orm.database import add_dataset, update_rdkit_ids, update_rdkit_tables
from ord_schema.orm.testing import get_test_engine
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def main(**kwargs):
    with get_test_engine(echo=False) as (engine, rdkit_cartridge):
        assert rdkit_cartridge
        for filename in glob(kwargs["--pattern"]):
            with Session(engine) as session:
                start = time.time()
                dataset = load_message(filename, dataset_pb2.Dataset)
                logger.debug(f"load_message() took {time.time() - start:g}s")
                with session.begin():
                    add_dataset(dataset, session)
                with session.begin():
                    update_rdkit_tables(dataset.dataset_id, session)
                with session.begin():
                    update_rdkit_ids(dataset.dataset_id, session)


if __name__ == "__main__":
    docopt_kwargs = docopt(__doc__)
    cProfile.run("main(**docopt_kwargs)", docopt_kwargs["--output"])
