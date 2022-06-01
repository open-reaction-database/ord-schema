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
"""Tests for ord_schema.proto.dataset_pb2."""

import pytest

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


@pytest.fixture
def serialized_dataset() -> bytes:
    dataset = dataset_pb2.Dataset()
    dataset.name = "test"
    dataset.description = "test dataset"
    # Add a reaction directly to the dataset.
    reaction1 = dataset.reactions.add()
    reaction1.identifiers.add(value="C(C)Cl.Br>>C(C)Br.Cl", type="REACTION_SMILES")
    # Copy a reaction created elsewhere.
    reaction2 = reaction_pb2.Reaction()
    reaction2.identifiers.add(value="amide coupling", type="NAME")
    dataset.reactions.add().CopyFrom(reaction2)
    yield dataset.SerializeToString()


def test_dataset(serialized_dataset):
    dataset = dataset_pb2.Dataset.FromString(serialized_dataset)
    assert dataset.name == "test"
    assert dataset.description == "test dataset"
    assert len(dataset.reactions) == 2
    assert dataset.reactions[0].identifiers[0].type == reaction_pb2.ReactionIdentifier.REACTION_SMILES
    assert dataset.reactions[1].identifiers[0].type == reaction_pb2.ReactionIdentifier.NAME
