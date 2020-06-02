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

from absl.testing import absltest

from ord_schema import frozen_message
from ord_schema.proto import reaction_pb2


class FrozenMessageTest(absltest.TestCase):
    def test_scalar_without_optional(self):
        """Illustrates the problem with non-optional scalar fields.

        This is what happens if scalar fields are not labeled `optional`:
          * Unset scalar fields return True with hasattr().
          * It's not possible to know whether a value matching the default was
            set explicitly.
        """
        frozen = frozen_message.FrozenMessage(
            reaction_pb2.StirringConditions.StirringRate())
        self.assertTrue(hasattr(frozen, 'rpm'))  # Wrong!
        self.assertAlmostEqual(frozen.rpm, 0.0)

    def test_access_optional_scalar(self):
        frozen = frozen_message.FrozenMessage(
            reaction_pb2.Percentage(value=12.3))
        self.assertTrue(hasattr(frozen, 'value'))
        self.assertAlmostEqual(frozen.value, 12.3, places=6)
        self.assertFalse(hasattr(frozen, 'precision'))
        with self.assertRaises(AttributeError):
            _ = frozen.precision

    def test_access_submessage(self):
        frozen = frozen_message.FrozenMessage(reaction_pb2.Reaction())
        self.assertFalse(hasattr(frozen, 'setup'))
        with self.assertRaises(AttributeError):
            _ = frozen.setup

    def test_access_repeated_scalar(self):
        message = reaction_pb2.ReactionProduct()
        message.analysis_identity.append('test')
        frozen = frozen_message.FrozenMessage(message)
        self.assertLen(frozen.analysis_identity, 1)
        self.assertEqual(frozen.analysis_identity[0], 'test')
        with self.assertRaises(IndexError):
            _ = frozen.analysis_identity[1]
        self.assertEmpty(frozen.analysis_yield)
        with self.assertRaises(IndexError):
            _ = frozen.analysis_yield[0]

    def test_access_repeated_submessage(self):
        message = reaction_pb2.Reaction()
        message.observations.add(comment='test')
        frozen = frozen_message.FrozenMessage(message)
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
        frozen = frozen_message.FrozenMessage(message)
        self.assertTrue(hasattr(frozen, 'inputs'))
        self.assertIn('test', frozen.inputs)
        self.assertTrue(hasattr(frozen.inputs['test'], 'addition_order'))
        self.assertEqual(frozen.inputs['test'].addition_order, 1)
        self.assertNotIn('missing', frozen.inputs)
        with self.assertRaises(KeyError):
            _ = frozen.inputs['missing']

    def test_set_scalar_value(self):
        frozen = frozen_message.FrozenMessage(
            reaction_pb2.StirringConditions.StirringRate())
        with self.assertRaises(dataclasses.FrozenInstanceError):
            frozen.rpm = 123

    def test_set_submessage(self):
        frozen = frozen_message.FrozenMessage(reaction_pb2.Reaction())
        with self.assertRaises(AttributeError):
            frozen.setup.automation_platform = 'test'
        with self.assertRaises(AttributeError):
            frozen.setup.CopyFrom(
                reaction_pb2.ReactionSetup(automation_platform='test'))

    def test_add_repeated_scalar(self):
        frozen = frozen_message.FrozenMessage(reaction_pb2.ReactionProduct())
        with self.assertRaises(AttributeError):
            frozen.analysis_identity.append('test')

    def test_add_repeated_submessage(self):
        frozen = frozen_message.FrozenMessage(reaction_pb2.Reaction())
        with self.assertRaises(AttributeError):
            frozen.observations.add(comment='test')

    def test_set_map_value(self):
        frozen = frozen_message.FrozenMessage(reaction_pb2.Reaction())
        with self.assertRaises(KeyError):
            frozen.inputs['test'].addition_order = 1
        with self.assertRaises(KeyError):
            frozen.inputs['test'].CopyFrom(
                reaction_pb2.ReactionInput(addition_order=1))

    def test_modify_repeated_submessage(self):
        """See https://git.io/JfPf9."""
        message = reaction_pb2.Reaction()
        message.workup.add(type='ADDITION')
        frozen = frozen_message.FrozenMessage(message)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            frozen.workup[0].type = reaction_pb2.ReactionWorkup.TEMPERATURE


if __name__ == '__main__':
    absltest.main()
