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
"""Tests for ord_schema.message_helpers."""
import pytest

from ord_schema import units
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import generate_text


@pytest.fixture
def reaction() -> reaction_pb2.Reaction:
    resolver = units.UnitResolver()
    reaction = reaction_pb2.Reaction()
    reaction.setup.is_automated = True
    reaction.inputs['dummy_input'].components.add().CopyFrom(
        message_helpers.build_compound(
            name='n-hexane',
            smiles='CCCCCC',
            role='reactant',
            amount='1 milliliters',
        ))
    reaction.inputs['dummy_input'].components.add().CopyFrom(
        message_helpers.build_compound(
            name='THF',
            smiles='C1OCCC1',
            role='solvent',
            amount='40 liters',
        ))
    reaction.inputs['dummy_input2'].components.add().CopyFrom(
        message_helpers.build_compound(
            name='Pd',
            smiles='[Pd]',
            role='catalyst',
            amount='catalytic',
        ))
    reaction.conditions.pressure.atmosphere.type = (
        reaction_pb2.PressureConditions.Atmosphere.OXYGEN)
    reaction.conditions.stirring.rate.rpm = 100
    reaction.conditions.temperature.control.type = (
        reaction_pb2.TemperatureConditions.TemperatureControl.OIL_BATH)
    reaction.conditions.temperature.setpoint.CopyFrom(
        reaction_pb2.Temperature(value=100,
                                 units=reaction_pb2.Temperature.CELSIUS))
    outcome = reaction.outcomes.add()
    outcome.reaction_time.CopyFrom(resolver.resolve('40 minutes'))
    outcome.products.add().identifiers.extend(
        message_helpers.build_compound(
            name='hexanone',
            smiles='CCCCC(=O)C',
        ).identifiers)
    yield reaction


def test_text(reaction):
    text = generate_text.generate_text(reaction)
    expected = ['vessel', 'oil bath', 'after 40 min', 'as a solvent', 'hexanone', 'automatically', 'mL', 'catalytic']
    for value in expected:
        assert value in text


def test_html(reaction):
    html = generate_text.generate_html(reaction)
    expected = ['<table', 'hexanone', 'under oxygen', '100 rpm', '40 min', 'solvent', '100 °C', 'catalytic']
    for value in expected:
        assert value in html


def test_compact_html(reaction):
    html = generate_text.generate_html(reaction, compact=True)
    assert '<table' in html
    expected = ['hexanone', 'under_oxygen', '100 rpm', '40 min', 'solvent', '100 °C']
    for value in expected:
        assert value not in html


def test_reaction_smiles_html():
    reaction = reaction_pb2.Reaction()
    reaction.identifiers.add(value='C>C>C', type='REACTION_SMILES')
    assert generate_text.generate_html(reaction) is not None
