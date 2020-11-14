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
"""Tests for ord_schema.interface.ord_client."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema.interface import ord_client


class OrdClientTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.client = ord_client.OrdClient()

    @parameterized.parameters(
        ('ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c', 264))
    def test_fetch_dataset(self, dataset_id, expected_num_reactions):
        dataset = self.client.fetch_dataset(dataset_id)
        self.assertLen(dataset.reactions, expected_num_reactions)

    @parameterized.parameters(
        (['ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c'], [264]))
    def test_fetch_datasets(self, dataset_ids, expected_num_reactions):
        datasets = self.client.fetch_datasets(dataset_ids)
        self.assertLen(dataset_ids, len(expected_num_reactions))
        self.assertLen(datasets, len(expected_num_reactions))
        for dataset, expected in zip(datasets, expected_num_reactions):
            self.assertLen(dataset.reactions, expected)

    @parameterized.parameters(
        ('ord-f0621fa47ac74fd59f9da027f6d13fc4', 'Jun Li'),
        ('ord-c6fbf2aab30841d198a27068a65a9a98', 'Steven Kearnes'))
    def test_fetch_reaction(self, reaction_id, created_by):
        reaction = self.client.fetch_reaction(reaction_id)
        self.assertEqual(reaction.provenance.record_created.person.name,
                         created_by)

    @parameterized.parameters(([
        'ord-f0621fa47ac74fd59f9da027f6d13fc4',
        'ord-c6fbf2aab30841d198a27068a65a9a98'
    ], ['Jun Li', 'Steven Kearnes']))
    def test_fetch_reactions(self, reaction_ids, created_by):
        reactions = self.client.fetch_reactions(reaction_ids)
        self.assertLen(reaction_ids, len(created_by))
        self.assertLen(reactions, len(created_by))
        for reaction, expected in zip(reactions, created_by):
            self.assertEqual(reaction.provenance.record_created.person.name,
                             expected)

    def test_query_reaction_ids(self):
        dataset = self.client.query(reaction_ids=[
            'ord-f0621fa47ac74fd59f9da027f6d13fc4',
            'ord-c6fbf2aab30841d198a27068a65a9a98'
        ])
        self.assertLen(dataset.reactions, 2)

    def test_query_reaction_smarts(self):
        dataset = self.client.query(
            reaction_smarts='FC(F)(F)c1ccc([F,Cl,Br,I])cc1>CS(=O)C>')
        self.assertLen(dataset.reactions, 862)

    def test_query_single_component(self):
        component = ord_client.ComponentQuery('FC(F)(F)c1ccc(Br)cc1',
                                              source='input',
                                              mode='exact')
        dataset = self.client.query(components=[component])
        self.assertLen(dataset.reactions, 287)

    def test_query_multiple_components(self):
        component1 = ord_client.ComponentQuery('FC(F)(F)c1ccc(Br)cc1',
                                               source='input',
                                               mode='exact')
        component2 = ord_client.ComponentQuery('Cc1cc(C)on1',
                                               source='input',
                                               mode='exact')
        dataset = self.client.query(components=[component1, component2])
        self.assertLen(dataset.reactions, 12)

    def test_query_stereochemistry(self):
        # TODO(kearnes): Add a test once we have more chiral stuff.
        pass

    def test_query_similarity(self):
        component = ord_client.ComponentQuery('FC(F)(F)c1ccc(Br)cc1',
                                              source='input',
                                              mode='similar')
        dataset = self.client.query(components=[component], similarity=0.6)
        self.assertLen(dataset.reactions, 575)


if __name__ == '__main__':
    absltest.main()
