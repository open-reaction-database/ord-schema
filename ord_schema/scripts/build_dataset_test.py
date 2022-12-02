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
"""Tests for ord_schema.scripts.build_dataset."""
import os

import docopt
import pytest

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import build_dataset


@pytest.fixture
def dirname(tmp_path) -> str:
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
    dirname = tmp_path.as_posix()
    message_helpers.write_message(reaction1, os.path.join(dirname, "reaction-1.pbtxt"))
    # reaction2 is empty.
    reaction2 = reaction_pb2.Reaction()
    message_helpers.write_message(reaction2, os.path.join(dirname, "reaction-2.pbtxt"))
    yield dirname


def test_simple(dirname):
    input_pattern = os.path.join(dirname, "reaction-1.pbtxt")
    output_filename = os.path.join(dirname, "dataset.pbtxt")
    argv = [
        "--input",
        input_pattern,
        "--name",
        "test dataset",
        "--description",
        "this is a test dataset",
        "--output",
        output_filename,
    ]
    build_dataset.main(docopt.docopt(build_dataset.__doc__, argv))
    assert os.path.exists(output_filename)
    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    assert dataset.name == "test dataset"
    assert dataset.description == "this is a test dataset"
    assert len(dataset.reactions) == 1


def test_validation(dirname):
    input_pattern = os.path.join(dirname, "reaction-?.pbtxt")
    output_filename = os.path.join(dirname, "dataset.pbtxt")
    argv = [
        "--input",
        input_pattern,
        "--name",
        "test dataset",
        "--description",
        "this is a test dataset",
        "--output",
        output_filename,
    ]
    with pytest.raises(
        validations.ValidationError,
        match="Reactions should have at least 1 reaction input",
    ):
        build_dataset.main(docopt.docopt(build_dataset.__doc__, argv))
    # Make sure disabling validation works.
    argv = [
        "--input",
        input_pattern,
        "--name",
        "test dataset",
        "--description",
        "this is a test dataset",
        "--output",
        output_filename,
        "--no-validate",
    ]
    build_dataset.main(docopt.docopt(build_dataset.__doc__, argv))
