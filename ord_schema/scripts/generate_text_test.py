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

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from ord_schema import units
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import generate_text


class GenerateTextTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
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
            ).identifers)
        reaction.reaction_id = 'dummy_reaction_id'
        self._reaction = reaction
        self._input = os.path.join(self.test_subdirectory, 'reaction.pbtxt')
        message_helpers.write_message(self._reaction, self._input)

    @parameterized.parameters(
        ('text', None),
        ('text', 'test.txt'),
        ('html', None),
        ('html', 'test.html'),
    )
    def test_main(self, output_type, output):
        with flagsaver.flagsaver(input=self._input,
                                 output_type=output_type,
                                 output=output):
            generate_text.main(())
        if output:
            self.assertTrue(os.path.exists(output))


if __name__ == '__main__':
    absltest.main()
