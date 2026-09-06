# Copyright 2024 Open Reaction Database Project Authors
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

from ord_schema import validations
from ord_schema.proto import reaction_pb2

# Section-level metadata field numbers must stay stable for wire compatibility.
_SECTION_METADATA_FIELDS = (
    (reaction_pb2.ReactionNotes, 10),
    (reaction_pb2.ReactionInput, 11),
    (reaction_pb2.ReactionSetup, 6),
    (reaction_pb2.ReactionConditions, 11),
    (reaction_pb2.TemperatureConditions, 4),
    (reaction_pb2.PressureConditions, 5),
    (reaction_pb2.StirringConditions, 4),
    (reaction_pb2.IlluminationConditions, 6),
    (reaction_pb2.ElectrochemistryConditions, 10),
    (reaction_pb2.FlowConditions, 5),
    (reaction_pb2.ReactionWorkup, 11),
)


def test_simple():
    reaction = reaction_pb2.Reaction()
    reaction.identifiers.add(value="C(C)Cl.Br>>C(C)Br.Cl", type="REACTION_SMILES")
    assert reaction.IsInitialized()
    assert len(reaction.identifiers) == 1
    assert not reaction.HasField("setup")
    with pytest.raises(ValueError, match="not_a_field"):
        reaction.HasField("not_a_field")


@pytest.mark.parametrize(("message_cls", "number"), _SECTION_METADATA_FIELDS)
def test_section_metadata_field_numbers(message_cls, number):
    field = message_cls.DESCRIPTOR.fields_by_name["metadata"]
    assert field.number == number
    assert field.message_type.GetOptions().map_entry
    value_type = field.message_type.fields_by_name["value"].message_type
    assert value_type == reaction_pb2.Data.DESCRIPTOR


def test_section_metadata_serialize_round_trip():
    reaction = reaction_pb2.Reaction()
    reaction.inputs["in"].metadata["k"].string_value = "input"
    reaction.setup.metadata["k"].string_value = "setup"
    reaction.setup.automation_code["script"].string_value = "print(1)"
    reaction.conditions.metadata["k"].string_value = "conditions"
    reaction.conditions.temperature.metadata["k"].string_value = "temp"
    reaction.conditions.pressure.metadata["k"].string_value = "pressure"
    reaction.conditions.stirring.metadata["k"].string_value = "stirring"
    reaction.conditions.illumination.metadata["k"].string_value = "illum"
    reaction.conditions.electrochemistry.metadata["k"].string_value = "electro"
    reaction.conditions.flow.metadata["k"].string_value = "flow"
    reaction.notes.metadata["k"].string_value = "notes"
    workup = reaction.workups.add()
    workup.metadata["k"].string_value = "workup"
    workup.temperature.metadata["k"].string_value = "workup_temp"
    workup.stirring.metadata["k"].string_value = "workup_stirring"
    workup.stirring.type = reaction_pb2.StirringConditions.STIR_BAR

    restored = reaction_pb2.Reaction()
    restored.ParseFromString(reaction.SerializeToString())
    assert restored.inputs["in"].metadata["k"].string_value == "input"
    assert restored.setup.metadata["k"].string_value == "setup"
    assert restored.setup.automation_code["script"].string_value == "print(1)"
    assert restored.conditions.metadata["k"].string_value == "conditions"
    assert restored.conditions.temperature.metadata["k"].string_value == "temp"
    assert restored.conditions.pressure.metadata["k"].string_value == "pressure"
    assert restored.conditions.stirring.metadata["k"].string_value == "stirring"
    assert restored.conditions.illumination.metadata["k"].string_value == "illum"
    assert restored.conditions.electrochemistry.metadata["k"].string_value == "electro"
    assert restored.conditions.flow.metadata["k"].string_value == "flow"
    assert restored.notes.metadata["k"].string_value == "notes"
    assert restored.workups[0].metadata["k"].string_value == "workup"
    assert restored.workups[0].temperature.metadata["k"].string_value == "workup_temp"
    assert restored.workups[0].stirring.metadata["k"].string_value == "workup_stirring"


def test_section_metadata_data_validation():
    hosts = [
        reaction_pb2.ReactionNotes(),
        reaction_pb2.ReactionInput(
            components=[
                reaction_pb2.Compound(
                    identifiers=[
                        reaction_pb2.CompoundIdentifier(type="NAME", value="solvent")
                    ]
                )
            ]
        ),
        reaction_pb2.ReactionSetup(),
        reaction_pb2.ReactionConditions(),
        reaction_pb2.TemperatureConditions(),
        reaction_pb2.PressureConditions(),
        reaction_pb2.StirringConditions(type=reaction_pb2.StirringConditions.STIR_BAR),
        reaction_pb2.IlluminationConditions(
            type=reaction_pb2.IlluminationConditions.LED
        ),
        reaction_pb2.ElectrochemistryConditions(
            type=reaction_pb2.ElectrochemistryConditions.CONSTANT_CURRENT
        ),
        reaction_pb2.FlowConditions(type=reaction_pb2.FlowConditions.PLUG_FLOW_REACTOR),
        reaction_pb2.ReactionWorkup(type=reaction_pb2.ReactionWorkup.FILTRATION),
    ]
    for host in hosts:
        host.ClearField("metadata")
        host.metadata["empty"]  # Creates an empty Data entry.
        with pytest.raises(validations.ValidationError, match="requires one of"):
            validations.validate_message(host)

        host.ClearField("metadata")
        host.metadata["bytes"].bytes_value = b"payload"
        with pytest.raises(validations.ValidationError, match="format is required"):
            validations.validate_message(host)

        host.ClearField("metadata")
        host.metadata["ok"].string_value = "scalar"
        host.metadata["attach"].bytes_value = b"payload"
        host.metadata["attach"].format = "bin"
        output = validations.validate_message(host)
        assert not output.errors

        host.ClearField("metadata")
        host.metadata["bad_url"].url = "not-a-url"
        output = validations.validate_message(host)
        assert not output.errors
        assert any(
            "does not look like a valid URL" in warning for warning in output.warnings
        )
