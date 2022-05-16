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
"""Tests for ord_schema.updates."""

from absl.testing import absltest

from ord_schema import updates
from ord_schema.proto import reaction_pb2
from ord_schema.proto import dataset_pb2


class UpdateReactionTest(absltest.TestCase):

    def test_with_updates_simple(self):
        message = reaction_pb2.Reaction()
        updates.update_reaction(message)
        assert message != reaction_pb2.Reaction()
        self.assertLen(message.provenance.record_modified, 1)

    def test_with_no_updates(self):
        message = reaction_pb2.Reaction()
        message.provenance.record_created.time.value = '2020-05-08'
        message.reaction_id = 'ord-c0bbd41f095a44a78b6221135961d809'
        copied = reaction_pb2.Reaction()
        copied.CopyFrom(message)
        updates.update_reaction(copied)
        assert copied == message

    def test_add_reaction_id(self):
        message = reaction_pb2.Reaction()
        updates.update_reaction(message)
        self.assertNotEmpty(message.reaction_id)
        self.assertLen(message.provenance.record_modified, 1)

    def test_keep_existing_reaction_id(self):
        message = reaction_pb2.Reaction()
        message.reaction_id = 'ord-c0bbd41f095a44a78b6221135961d809'
        message.provenance.record_created.time.value = '2020-01-01'
        updates.update_reaction(message)
        assert message.reaction_id == \
                         'ord-c0bbd41f095a44a78b6221135961d809'
        self.assertLen(message.provenance.record_modified, 0)


class UpdateDatasetTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.dataset = dataset_pb2.Dataset()
        reaction = self.dataset.reactions.add()
        # Minimal reaction.
        dummy_input = reaction.inputs['dummy_input']
        reaction.outcomes.add()
        dummy_component = dummy_input.components.add()
        dummy_component.identifiers.add(type='CUSTOM')
        dummy_component.identifiers[0].details = 'custom_identifier'
        dummy_component.identifiers[0].value = 'custom_value'
        dummy_component.amount.mass.value = 1
        dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
        # Placeholders for referenced reactions.
        reaction2 = self.dataset.reactions.add()
        reaction3 = self.dataset.reactions.add()
        reaction2.CopyFrom(self.dataset.reactions[0])
        reaction3.CopyFrom(self.dataset.reactions[0])

    def test_crossferences(self):
        dummy_input = self.dataset.reactions[0].inputs['dummy_input']
        dummy_input.crude_components.add(reaction_id='crude',
                                         has_derived_amount=True)
        self.dataset.reactions[1].reaction_id = 'crude'
        dummy_input.components[0].preparations.add(type='SYNTHESIZED',
                                                   reaction_id='synthesized')
        self.dataset.reactions[2].reaction_id = 'synthesized'
        # Check updated values.
        updates.update_dataset(self.dataset)
        assert dummy_input.crude_components[0].reaction_id == \
                         self.dataset.reactions[1].reaction_id
        assert dummy_input.crude_components[0].reaction_id != \
                            'crude'
        assert dummy_input.components[0].preparations[0].reaction_id == \
                         self.dataset.reactions[2].reaction_id
        assert dummy_input.components[0].preparations[0].reaction_id != \
            'synthesized'

    def test_crossferences_are_proper_ids(self):
        dummy_input = self.dataset.reactions[0].inputs['dummy_input']
        dummy_input.crude_components.add(
            reaction_id='ord-c0bbd41f095a44a78b6221135961d809',
            has_derived_amount=True)
        self.dataset.reactions[
            1].reaction_id = 'ord-c0bbd41f095a44a78b6221135961d809'
        dummy_input.components[0].preparations.add(
            type='SYNTHESIZED',
            reaction_id='ord-d0bbd41f095a44a78b6221135961d809')
        self.dataset.reactions[
            2].reaction_id = 'ord-d0bbd41f095a44a78b6221135961d809'
        # Check updated values.
        updates.update_dataset(self.dataset)
        assert dummy_input.crude_components[0].reaction_id == \
                         self.dataset.reactions[1].reaction_id
        assert dummy_input.crude_components[0].reaction_id == \
                         'ord-c0bbd41f095a44a78b6221135961d809'
        assert dummy_input.components[0].preparations[0].reaction_id == \
                         self.dataset.reactions[2].reaction_id
        assert dummy_input.components[0].preparations[0].reaction_id == \
                         'ord-d0bbd41f095a44a78b6221135961d809'


if __name__ == '__main__':
    absltest.main()
