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

"""Tests for ord_schema.orm.scripts.profile_insert."""
import os

import docopt
import pytest

from ord_schema.orm.scripts import profile_insert
from ord_schema.orm.testing import get_test_engine


def test_main(tmp_path):
    with get_test_engine() as (_, rdkit_cartridge):
        if not rdkit_cartridge:
            pytest.skip("RDKit cartridge is required")
    filename = (tmp_path / "test.profile").as_posix()
    argv = [
        "--pattern",
        os.path.join(os.path.dirname(__file__), "..", "testdata", "ord-nielsen-example.pbtxt"),
        "--output",
        filename,
    ]
    profile_insert.main(**docopt.docopt(profile_insert.__doc__, argv))
