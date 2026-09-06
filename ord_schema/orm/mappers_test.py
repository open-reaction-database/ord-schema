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

import pathlib
import subprocess
import sys

import pytest

from ord_schema.message_helpers import load_message
from ord_schema.orm.mappers import from_proto, to_proto
from ord_schema.proto.dataset_pb2 import Dataset


def test_mappers_registers_rdkit_relationship_targets_in_clean_interpreter():
    """Ensure mappers.py stays self-contained for string relationship() names.

    Pytest collection loads ord_schema.orm.rdkit_mappers via rdkit_mappers_test, which
    registers RDKitMols / RDKitReactions and hides a missing side-effect import in
    mappers.py. A fresh interpreter only imports mappers.
    """
    script = (
        "from sqlalchemy.orm import configure_mappers\n"
        "import ord_schema.orm.mappers\n"
        "configure_mappers()\n"
    )
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
    [
        pathlib.Path(__file__).parent / "testdata" / "empty.pbtxt",
        pathlib.Path(__file__).parent / "testdata" / "full.pbtxt",
        pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
    ],
)
def test_round_trip(filename):
    dataset = load_message(filename, Dataset)
    assert dataset == to_proto(from_proto(dataset))


def test_full_fixture_preserves_section_metadata_paths():
    """Distinct metadata keys under nested contexts must not swap or drop."""
    dataset = load_message(
        pathlib.Path(__file__).parent / "testdata" / "full.pbtxt", Dataset
    )
    reaction = to_proto(from_proto(dataset)).reactions[0]
    assert (
        reaction.inputs["4"].metadata["input_meta"].string_value == "input_meta_value"
    )
    assert reaction.setup.metadata["setup_meta"].string_value == "setup_meta_value"
    assert (
        reaction.conditions.metadata["conditions_meta"].string_value
        == "conditions_meta_value"
    )
    assert (
        reaction.conditions.temperature.metadata["temp_meta"].string_value
        == "temp_meta_value"
    )
    assert (
        reaction.conditions.pressure.metadata["pressure_meta"].string_value
        == "pressure_meta_value"
    )
    assert (
        reaction.conditions.stirring.metadata["stirring_meta"].string_value
        == "stirring_meta_value"
    )
    assert (
        reaction.conditions.illumination.metadata["illum_meta"].string_value
        == "illum_meta_value"
    )
    assert (
        reaction.conditions.electrochemistry.metadata["electro_meta"].string_value
        == "electro_meta_value"
    )
    assert (
        reaction.conditions.flow.metadata["flow_meta"].string_value == "flow_meta_value"
    )
    assert reaction.notes.metadata["notes_meta"].string_value == "notes_meta_value"
    workup = reaction.workups[0]
    assert workup.metadata["workup_meta"].string_value == "workup_meta_value"
    assert (
        workup.temperature.metadata["workup_temp_meta"].string_value
        == "workup_temp_meta_value"
    )
    assert (
        workup.stirring.metadata["workup_stirring_meta"].string_value
        == "workup_stirring_meta_value"
    )
