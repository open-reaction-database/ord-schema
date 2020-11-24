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
        self.assertNotEqual(message, reaction_pb2.Reaction())
        self.assertLen(message.provenance.record_modified, 1)

    def test_with_no_updates(self):
        message = reaction_pb2.Reaction()
        message.provenance.record_created.time.value = '2020-05-08'
        message.reaction_id = 'ord-c0bbd41f095a44a78b6221135961d809'
        copied = reaction_pb2.Reaction()
        copied.CopyFrom(message)
        updates.update_reaction(copied)
        self.assertEqual(copied, message)

    def test_with_resolve_names(self):
        reaction = reaction_pb2.Reaction()
        component = reaction.inputs['ethylamine'].components.add()
        component.identifiers.add(type='NAME', value='ethylamine')
        updates.update_reaction(reaction)
        self.assertLen(component.identifiers, 2)
        self.assertEqual(component.identifiers[1].value, 'CCN')
        self.assertEqual(component.identifiers[1].type,
                         reaction_pb2.CompoundIdentifier.IdentifierType.SMILES)
        self.assertRegex(component.identifiers[1].details, 'NAME resolved')

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
        self.assertEqual(message.reaction_id,
                         'ord-c0bbd41f095a44a78b6221135961d809')
        self.assertLen(message.provenance.record_modified, 0)


class UpdateDatasetTest(absltest.TestCase):

    def test_crossferences(self):
        message = dataset_pb2.Dataset()
        reaction1 = message.reactions.add()
        reaction2 = message.reactions.add()
        reaction3 = message.reactions.add()
        # Minimal reaction 1
        dummy_input = reaction1.inputs['dummy_input']
        reaction1.outcomes.add()
        dummy_component = dummy_input.components.add()
        dummy_component.identifiers.add(type='CUSTOM')
        dummy_component.identifiers[0].details = 'custom_identifier'
        dummy_component.identifiers[0].value = 'custom_value'
        dummy_component.amount.mass.value = 1
        dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
        reaction2.CopyFrom(reaction1)
        reaction3.CopyFrom(reaction1)
        dummy_component.preparations.add(type='SYNTHESIZED')
        dummy_component.preparations[0].reaction_id = 'placeholder_id'
        reaction2.reaction_id = 'placeholder_id'
        dummy_input.crude_components.add(reaction_id='crude-making step',
                                         has_derived_amount=True)
        reaction3.reaction_id = 'crude-making step'

        updates.update_dataset(message)
        self.assertEqual(dummy_component.preparations[0].reaction_id,
                         reaction2.reaction_id)
        self.assertEqual(dummy_input.crude_components[0].reaction_id,
                         reaction3.reaction_id)
        self.assertNotEqual(dummy_component.preparations[0].reaction_id,
                            'placeholder_id')
        self.assertNotEqual(dummy_input.crude_components[0].reaction_id,
                            'crude-making step')


if __name__ == '__main__':
    absltest.main()
