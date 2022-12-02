# Copyright 2020 The Open Reaction Database Authors
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
"""Tests for ord_schema.frozen_storage."""

import dataclasses
import os
import tempfile

import pytest

from ord_schema import frozen_message
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.proto import test_pb2


def _freeze(message) -> frozen_message.FrozenMessage:
    """Runs a round-trip to disk as an extra check for FrozenMessage."""
    with tempfile.TemporaryDirectory() as tempdir:
        filename = os.path.join(tempdir, "message.pbtxt")
        message_helpers.write_message(message, filename)
        loaded_message = message_helpers.load_message(filename, type(message))
        return frozen_message.FrozenMessage(loaded_message)


def test_scalar_without_optional():
    """Illustrates the problem with non-optional scalar fields.

    This is what happens if scalar fields are not labeled `optional`:
      * Unset scalar fields return True with hasattr().
      * It's not possible to know whether a value matching the default was
        set explicitly.
    """
    frozen = _freeze(reaction_pb2.StirringConditions.StirringRate())
    assert hasattr(frozen, "rpm")  # Wrong!
    assert round(abs(frozen.rpm - 0.0), 7) == 0


def test_access_optional_scalar():
    frozen = _freeze(reaction_pb2.Percentage(value=12.3))
    assert hasattr(frozen, "value")
    assert round(abs(frozen.value - 12.3), 6) == 0
    assert not hasattr(frozen, "precision")
    with pytest.raises(AttributeError):
        _ = frozen.precision


def test_access_optional_bool():
    frozen = _freeze(test_pb2.Scalar(string_value="test"))
    assert not hasattr(frozen, "bool_value")
    with pytest.raises(AttributeError):
        _ = frozen.bool_value
    frozen = _freeze(test_pb2.Scalar(bool_value=False))
    assert hasattr(frozen, "bool_value")
    assert not frozen.bool_value
    frozen = _freeze(test_pb2.Scalar(bool_value=True))
    assert hasattr(frozen, "bool_value")
    assert frozen.bool_value


def test_access_submessage():
    frozen = _freeze(reaction_pb2.Reaction())
    assert not hasattr(frozen, "setup")
    with pytest.raises(AttributeError):
        _ = frozen.setup


def test_access_repeated_scalar():
    message = reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails()
    frozen_empty = _freeze(message)
    assert len(frozen_empty.eic_masses) == 0
    with pytest.raises(IndexError):
        _ = frozen_empty.eic_masses[0]

    message.eic_masses.append(1.5)
    frozen = _freeze(message)
    assert len(frozen.eic_masses) == 1
    assert frozen.eic_masses[0] == 1.5
    with pytest.raises(IndexError):
        _ = frozen.eic_masses[1]


def test_access_repeated_submessage():
    message = reaction_pb2.Reaction()
    message.observations.add(comment="test")
    frozen = _freeze(message)
    assert len(frozen.observations) == 1
    assert frozen.observations[0].comment == "test"
    with pytest.raises(IndexError):
        _ = frozen.observations[1]
    assert len(frozen.outcomes) == 0
    with pytest.raises(IndexError):
        _ = frozen.outcomes[0]


def test_access_map_value():
    message = reaction_pb2.Reaction()
    message.inputs["test"].addition_order = 1
    frozen = _freeze(message)
    assert hasattr(frozen, "inputs")
    assert "test" in frozen.inputs
    assert hasattr(frozen.inputs["test"], "addition_order")
    assert frozen.inputs["test"].addition_order == 1
    assert "missing" not in frozen.inputs
    with pytest.raises(KeyError):
        _ = frozen.inputs["missing"]


def test_set_scalar_value():
    frozen = _freeze(reaction_pb2.StirringConditions.StirringRate())
    with pytest.raises(dataclasses.FrozenInstanceError):
        frozen.rpm = 123


def test_set_bool_value():
    frozen = _freeze(test_pb2.Scalar())
    with pytest.raises(dataclasses.FrozenInstanceError):
        frozen.bool_value = False


def test_set_submessage():
    frozen = _freeze(reaction_pb2.Reaction())
    with pytest.raises(AttributeError):
        frozen.setup.automation_platform = "test"
    with pytest.raises(AttributeError):
        frozen.setup.CopyFrom(reaction_pb2.ReactionSetup(automation_platform="test"))


def test_add_repeated_scalar():
    frozen = _freeze(reaction_pb2.ProductCompound())
    with pytest.raises(AttributeError):
        frozen.analysis_identity.append("test")


def test_add_repeated_submessage():
    frozen = _freeze(reaction_pb2.Reaction())
    with pytest.raises(AttributeError):
        frozen.observations.add(comment="test")


def test_set_map_value():
    frozen = _freeze(reaction_pb2.Reaction())
    with pytest.raises(KeyError):
        frozen.inputs["test"].addition_order = 1
    with pytest.raises(KeyError):
        frozen.inputs["test"].CopyFrom(reaction_pb2.ReactionInput(addition_order=1))


def test_modify_repeated_submessage():
    """See https://git.io/JfPf9."""
    message = reaction_pb2.Reaction()
    message.workups.add(type="ADDITION")
    frozen = _freeze(message)
    with pytest.raises(dataclasses.FrozenInstanceError):
        frozen.workups[0].type = reaction_pb2.ReactionWorkup.TEMPERATURE
