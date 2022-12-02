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

"""Tests for ord_schema.orm.add_datasets."""
import os

import docopt
import pytest
from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError
from testing.postgresql import Postgresql

from ord_schema.orm import add_datasets
from ord_schema.orm.database import prepare_database


def test_main():
    with Postgresql() as postgres:
        engine = create_engine(postgres.url(), future=True)
        assert prepare_database(engine)  # Requires RDKit.
        argv = [
            "--url",
            postgres.url(),
            "--pattern",
            os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"),
        ]
        add_datasets.main(**docopt.docopt(add_datasets.__doc__, argv))
        with pytest.raises(IntegrityError, match="violates unique constraint"):
            add_datasets.main(**docopt.docopt(add_datasets.__doc__, argv))
        argv.append("--overwrite")
        add_datasets.main(**docopt.docopt(add_datasets.__doc__, argv))
