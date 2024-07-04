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

"""Tests for ord_schema.orm.scripts.add_datasets."""
import os

import docopt
import pytest

from ord_schema.orm.database import prepare_database
from ord_schema.orm.scripts import add_datasets


def test_main(test_engine):
    if not prepare_database(test_engine):
        pytest.skip("RDKit cartridge is required")
    argv = [
        "--url",
        test_engine.url,
        "--pattern",
        os.path.join(os.path.dirname(__file__), "..", "testdata", "ord-nielsen-example.pbtxt"),
    ]
    add_datasets.main(**docopt.docopt(add_datasets.__doc__, argv))
