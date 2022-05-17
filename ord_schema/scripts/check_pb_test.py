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
"""Tests for ord_schema.scripts.check_pb."""
import docopt
import pytest

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.scripts import check_pb


def test_main_pass(tmp_path):
    pb_filename = (tmp_path / "test.pb").as_posix()
    pbtxt_filename = (tmp_path / "test.pbtxt").as_posix()
    dataset = dataset_pb2.Dataset()
    reaction = dataset.reactions.add()
    component = reaction.inputs["test"].components.add()
    component.identifiers.add(value="c1ccccc1", type="SMILES")
    message_helpers.write_message(dataset, pb_filename)
    message_helpers.write_message(dataset, pbtxt_filename)
    check_pb.main(docopt.docopt(check_pb.__doc__, ["--pb", pb_filename, "--pbtxt", pbtxt_filename]))


def test_main_fail(tmp_path):
    pb_filename = (tmp_path / "test.pb").as_posix()
    pbtxt_filename = (tmp_path / "test.pbtxt").as_posix()
    dataset = dataset_pb2.Dataset()
    reaction = dataset.reactions.add()
    component = reaction.inputs["test"].components.add()
    component.identifiers.add(value="c1ccccc1", type="SMILES")
    message_helpers.write_message(dataset, pb_filename)
    component.identifiers.add(value="benzene", type="NAME")
    message_helpers.write_message(dataset, pbtxt_filename)
    with pytest.raises(ValueError, match="Datasets differ"):
        check_pb.main(docopt.docopt(check_pb.__doc__, ["--pb", pb_filename, "--pbtxt", pbtxt_filename]))
