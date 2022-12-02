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
"""Tests for ord_schema.proto.reaction_pb2."""

import pytest

from ord_schema.proto import reaction_pb2


def test_simple():
    reaction = reaction_pb2.Reaction()
    reaction.identifiers.add(value="C(C)Cl.Br>>C(C)Br.Cl", type="REACTION_SMILES")
    assert reaction.IsInitialized()
    assert len(reaction.identifiers) == 1
    assert not reaction.HasField("setup")
    with pytest.raises(ValueError, match="not_a_field"):
        reaction.HasField("not_a_field")
