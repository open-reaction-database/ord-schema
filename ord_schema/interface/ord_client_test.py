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


if __name__ == '__main__':
    absltest.main()
