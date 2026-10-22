# Copyright 2026 Open Reaction Database Project Authors
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

"""Tests for ord_schema.orm."""

import os
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    "statement",
    ["import ord_schema.orm", "from ord_schema.orm.database import prepare_database"],
)
def test_import_warns_importer_under_default_filters(statement):
    """Ensure importing the ORM shows a DeprecationWarning naming the importer.

    The package is already imported in this process, so a fresh interpreter imports
    it. PYTHONWARNINGS is dropped so the interpreter runs Python's default filters,
    which show a DeprecationWarning only when it is attributed to __main__.
    """
    env = {key: value for key, value in os.environ.items() if key != "PYTHONWARNINGS"}
    proc = subprocess.run(
        [sys.executable, "-c", statement],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
        env=env,
    )
    assert proc.returncode == 0, f"stderr:\n{proc.stderr}\nstdout:\n{proc.stdout}"
    assert "<string>:1: DeprecationWarning: ord_schema.orm is deprecated" in proc.stderr
