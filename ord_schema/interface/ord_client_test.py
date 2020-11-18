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

    @parameterized.parameters(
        ('ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c', 264))
    def test_fetch_dataset(self, dataset_id, expected_num_reactions):
        client = ord_client.OrdClient()
        dataset = client.fetch_dataset(dataset_id)
        self.assertLen(dataset.reactions, expected_num_reactions)

    @parameterized.parameters(
        (['ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c'], [264]))
    def test_fetch_datasets(self, dataset_ids, expected_num_reactions):
        client = ord_client.OrdClient()
        datasets = client.fetch_datasets(dataset_ids)
        self.assertLen(dataset_ids, len(expected_num_reactions))
        self.assertLen(datasets, len(expected_num_reactions))
        for dataset, expected in zip(datasets, expected_num_reactions):
            self.assertLen(dataset.reactions, expected)

    @absltest.skip('Temporarily disabled for Reaction schema migration.')
    @parameterized.parameters(
        ('ord-f0621fa47ac74fd59f9da027f6d13fc4', 'Jun Li'),
        ('ord-c6fbf2aab30841d198a27068a65a9a98', 'Steven Kearnes'))
    def test_fetch_reaction(self, reaction_id, created_by):
        client = ord_client.OrdClient()
        reaction = client.fetch_reaction(reaction_id)
        self.assertEqual(reaction.provenance.record_created.person.name,
                         created_by)

    @absltest.skip('Temporarily disabled for Reaction schema migration.')
    @parameterized.parameters(([
        'ord-f0621fa47ac74fd59f9da027f6d13fc4',
        'ord-c6fbf2aab30841d198a27068a65a9a98'
    ], ['Jun Li', 'Steven Kearnes']))
    def test_fetch_reactions(self, reaction_ids, created_by):
        client = ord_client.OrdClient()
        reactions = client.fetch_reactions(reaction_ids)
        self.assertLen(reaction_ids, len(created_by))
        self.assertLen(reactions, len(created_by))
        for reaction, expected in zip(reactions, created_by):
            self.assertEqual(reaction.provenance.record_created.person.name,
                             expected)


if __name__ == '__main__':
    absltest.main()
