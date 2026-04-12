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

"""Tests for ord_schema.orm.mappers."""

import os
import subprocess
import sys

import pytest

from ord_schema.message_helpers import load_message
from ord_schema.orm.mappers import from_proto, to_proto
from ord_schema.proto.dataset_pb2 import Dataset


def test_mappers_registers_rdkit_relationship_targets_in_clean_interpreter():
    """Ensure mappers.py stays self-contained for string relationship() names.

    Pytest collection loads ord_schema.orm.rdkit_mappers via rdkit_mappers_test,
    which registers RDKitMols / RDKitReactions and hides a missing side-effect
    import in mappers.py. A fresh interpreter only imports mappers.
    """
    script = "from sqlalchemy.orm import configure_mappers\nimport ord_schema.orm.mappers\nconfigure_mappers()\n"
    proc = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert proc.returncode == 0, f"stderr:\n{proc.stderr}\nstdout:\n{proc.stdout}"


@pytest.mark.parametrize(
    "filename",
    (
        os.path.join(os.path.dirname(__file__), "testdata", "empty.pbtxt"),
        os.path.join(os.path.dirname(__file__), "testdata", "full.pbtxt"),
        os.path.join(os.path.dirname(__file__), "testdata", "ord-nielsen-example.pbtxt"),
    ),
)
def test_round_trip(filename):
    dataset = load_message(filename, Dataset)
    assert dataset == to_proto(from_proto(dataset))
