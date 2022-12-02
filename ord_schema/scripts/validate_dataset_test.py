# Copyright 2020 Open Reaction Database Project Authors
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
"""Tests for ord_schema.scripts.validate_dataset."""
import os

import docopt
import pytest

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import validate_dataset


@pytest.fixture
def setup(tmp_path) -> str:
    test_subdirectory = tmp_path.as_posix()
    reaction1 = reaction_pb2.Reaction()
    dummy_input = reaction1.inputs["dummy_input"]
    dummy_component = dummy_input.components.add()
    dummy_component.identifiers.add(type="CUSTOM")
    dummy_component.identifiers[0].details = "custom_identifier"
    dummy_component.identifiers[0].value = "custom_value"
    dummy_component.is_limiting = True
    dummy_component.amount.mass.value = 1
    dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
    reaction1.outcomes.add().conversion.value = 75
    dataset1 = dataset_pb2.Dataset(reactions=[reaction1])
    dataset1_filename = os.path.join(test_subdirectory, "dataset1.pbtxt")
    message_helpers.write_message(dataset1, dataset1_filename)
    # reaction2 is empty.
    reaction2 = reaction_pb2.Reaction()
    dataset2 = dataset_pb2.Dataset(reactions=[reaction1, reaction2])
    dataset2_filename = os.path.join(test_subdirectory, "dataset2.pbtxt")
    message_helpers.write_message(dataset2, dataset2_filename)
    yield test_subdirectory


def test_simple(setup):
    test_subdirectory = setup
    argv = ["--input", os.path.join(test_subdirectory, "dataset1.pbtxt")]
    validate_dataset.main(docopt.docopt(validate_dataset.__doc__, argv))


def test_filter(setup):
    test_subdirectory = setup
    argv = ["--input", os.path.join(test_subdirectory, "dataset1.pbtxt"), "--filter", "dataset"]
    validate_dataset.main(docopt.docopt(validate_dataset.__doc__, argv))


def test_validation_errors(setup):
    test_subdirectory = setup
    argv = ["--input", os.path.join(test_subdirectory, "dataset*.pbtxt")]
    with pytest.raises(
        validations.ValidationError,
        match="Reactions should have at least 1 reaction input",
    ):
        validate_dataset.main(docopt.docopt(validate_dataset.__doc__, argv))


@pytest.mark.parametrize(
    "pattern,expected",
    (
        (r"data/\d{2}", ["data/11/foo.pb"]),
        (r"data/\d[a-z]", ["data/1a/foo.pb"]),
        (r"data/[a-z]\d", ["data/a1/foo.pb"]),
        (r"data/[a-z]{2}", ["data/aa/foo.pb"]),
    ),
)
def test_filter_filenames(pattern, expected):
    filenames = [
        "dataset.pb",
        "data/a1/foo.pb",
        "data/aa/foo.pb",
        "data/11/foo.pb",
        "data/1a/foo.pb",
    ]
    assert expected == validate_dataset.filter_filenames(filenames, pattern)
