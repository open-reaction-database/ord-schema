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

from absl.testing import absltest

from ord_schema import units
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import generate_text


class GenerateTextTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self._resolver = units.UnitResolver()

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
        reaction.conditions.pressure.atmosphere.type = (
            reaction_pb2.PressureConditions.Atmosphere.OXYGEN)
        reaction.conditions.stirring.rate.rpm = 100
        reaction.conditions.temperature.control.type = (
            reaction_pb2.TemperatureConditions.TemperatureControl.OIL_BATH)
        reaction.conditions.temperature.setpoint.CopyFrom(
            reaction_pb2.Temperature(value=100,
                                     units=reaction_pb2.Temperature.CELSIUS))
        outcome = reaction.outcomes.add()
        outcome.reaction_time.CopyFrom(self._resolver.resolve('40 minutes'))
        outcome.products.add().identifiers.extend(
            message_helpers.build_compound(
                name='hexanone',
                smiles='CCCCC(=O)C',
            ).identifiers)
        self._reaction = reaction

    def test_text(self):
        text = generate_text.generate_text(self._reaction)
        self.assertRegex(text, 'vessel')
        self.assertRegex(text, 'oil bath')
        self.assertRegex(text, 'after 40 min')
        self.assertRegex(text, 'as a solvent')
        self.assertRegex(text, 'hexanone')
        self.assertRegex(text, 'automatically')
        self.assertRegex(text, 'mL')

    def test_html(self):
        html = generate_text.generate_html(self._reaction)
        self.assertRegex(html, '<table')
        self.assertRegex(html, 'hexanone')
        self.assertRegex(html, 'under oxygen')
        self.assertRegex(html, '100 rpm')
        self.assertRegex(html, '40 min')
        self.assertRegex(html, 'solvent')
        self.assertRegex(html, '100 °C')

    def test_compact_html(self):
        html = generate_text.generate_html(self._reaction, compact=True)
        self.assertRegex(html, '<table')
        self.assertNotRegex(html, 'hexanone')
        self.assertNotRegex(html, 'under oxygen')
        self.assertNotRegex(html, '100 rpm')
        self.assertNotRegex(html, '40 min')
        self.assertNotRegex(html, 'solvent')
        self.assertNotRegex(html, '100 °C')


if __name__ == '__main__':
    absltest.main()
