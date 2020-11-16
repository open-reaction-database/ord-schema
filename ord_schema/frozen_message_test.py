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

from absl import flags
from absl.testing import absltest

from ord_schema import frozen_message
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.proto import test_pb2


class FrozenMessageTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)

    def _freeze(self, message):
        """Runs a round-trip to disk as an extra check for FrozenMessage."""
        filename = os.path.join(self.test_subdirectory, 'message.pbtxt')
        message_helpers.write_message(message, filename)
        loaded_message = message_helpers.load_message(filename, type(message))
        return frozen_message.FrozenMessage(loaded_message)

    def test_scalar_without_optional(self):
        """Illustrates the problem with non-optional scalar fields.

        This is what happens if scalar fields are not labeled `optional`:
          * Unset scalar fields return True with hasattr().
          * It's not possible to know whether a value matching the default was
            set explicitly.
        """
        frozen = self._freeze(reaction_pb2.StirringConditions.StirringRate())
        self.assertTrue(hasattr(frozen, 'rpm'))  # Wrong!
        self.assertAlmostEqual(frozen.rpm, 0.0)

    def test_access_optional_scalar(self):
        frozen = self._freeze(reaction_pb2.Percentage(value=12.3))
        self.assertTrue(hasattr(frozen, 'value'))
        self.assertAlmostEqual(frozen.value, 12.3, places=6)
        self.assertFalse(hasattr(frozen, 'precision'))
        with self.assertRaises(AttributeError):
            _ = frozen.precision

    def test_access_optional_bool(self):
        frozen = self._freeze(test_pb2.Scalar(string_value='test'))
        self.assertFalse(hasattr(frozen, 'bool_value'))
        with self.assertRaises(AttributeError):
            _ = frozen.bool_value
        frozen = self._freeze(test_pb2.Scalar(bool_value=False))
        self.assertTrue(hasattr(frozen, 'bool_value'))
        self.assertFalse(frozen.bool_value)
        frozen = self._freeze(test_pb2.Scalar(bool_value=True))
        self.assertTrue(hasattr(frozen, 'bool_value'))
        self.assertTrue(frozen.bool_value)

    def test_access_submessage(self):
        frozen = self._freeze(reaction_pb2.Reaction())
        self.assertFalse(hasattr(frozen, 'setup'))
        with self.assertRaises(AttributeError):
            _ = frozen.setup

    def test_access_repeated_scalar(self):
        message = reaction_pb2.ProductMeasurement.MassSpecMeasurementType()
        frozen_empty = self._freeze(message)
        self.assertEmpty(frozen_empty.eic_masses)
        with self.assertRaises(IndexError):
            _ = frozen_empty.eic_masses[0]

        message.eic_masses.append(1.5)
        frozen = self._freeze(message)
        self.assertLen(frozen.eic_masses, 1)
        self.assertEqual(frozen.eic_masses[0], 1.5)
        with self.assertRaises(IndexError):
            _ = frozen.eic_masses[1]

    def test_access_repeated_submessage(self):
        message = reaction_pb2.Reaction()
        message.observations.add(comment='test')
        frozen = self._freeze(message)
        self.assertLen(frozen.observations, 1)
        self.assertEqual(frozen.observations[0].comment, 'test')
        with self.assertRaises(IndexError):
            _ = frozen.observations[1]
        self.assertEmpty(frozen.outcomes)
        with self.assertRaises(IndexError):
            _ = frozen.outcomes[0]

    def test_access_map_value(self):
        message = reaction_pb2.Reaction()
        message.inputs['test'].addition_order = 1
        frozen = self._freeze(message)
        self.assertTrue(hasattr(frozen, 'inputs'))
        self.assertIn('test', frozen.inputs)
        self.assertTrue(hasattr(frozen.inputs['test'], 'addition_order'))
        self.assertEqual(frozen.inputs['test'].addition_order, 1)
        self.assertNotIn('missing', frozen.inputs)
        with self.assertRaises(KeyError):
            _ = frozen.inputs['missing']

    def test_set_scalar_value(self):
        frozen = self._freeze(reaction_pb2.StirringConditions.StirringRate())
        with self.assertRaises(dataclasses.FrozenInstanceError):
            frozen.rpm = 123

    def test_set_bool_value(self):
        frozen = self._freeze(test_pb2.Scalar())
        with self.assertRaises(dataclasses.FrozenInstanceError):
            frozen.bool_value = False

    def test_set_submessage(self):
        frozen = self._freeze(reaction_pb2.Reaction())
        with self.assertRaises(AttributeError):
            frozen.setup.automation_platform = 'test'
        with self.assertRaises(AttributeError):
            frozen.setup.CopyFrom(
                reaction_pb2.ReactionSetup(automation_platform='test'))

    def test_add_repeated_scalar(self):
        frozen = self._freeze(reaction_pb2.ReactionProduct())
        with self.assertRaises(AttributeError):
            frozen.analysis_identity.append('test')

    def test_add_repeated_submessage(self):
        frozen = self._freeze(reaction_pb2.Reaction())
        with self.assertRaises(AttributeError):
            frozen.observations.add(comment='test')

    def test_set_map_value(self):
        frozen = self._freeze(reaction_pb2.Reaction())
        with self.assertRaises(KeyError):
            frozen.inputs['test'].addition_order = 1
        with self.assertRaises(KeyError):
            frozen.inputs['test'].CopyFrom(
                reaction_pb2.ReactionInput(addition_order=1))

    def test_modify_repeated_submessage(self):
        """See https://git.io/JfPf9."""
        message = reaction_pb2.Reaction()
        message.workup.add(type='ADDITION')
        frozen = self._freeze(message)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            frozen.workup[0].type = reaction_pb2.ReactionWorkup.TEMPERATURE


if __name__ == '__main__':
    absltest.main()
