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
"""Tests for ord_schema.interface.query."""

import time

from absl import logging
from absl.testing import absltest
from absl.testing import parameterized
import docker
import numpy as np
import psycopg2

from ord_schema.interface import query

# Avoid conflicts with any running PostgreSQL server(s).
_POSTGRES_PORT = 5430


class QueryTest(parameterized.TestCase, absltest.TestCase):
    @classmethod
    def setUpClass(cls):
        client = docker.from_env()
        # TODO(kearnes): Use a smaller test image if needed.
        client.images.pull('openreactiondatabase/ord-postgres')
        cls._container = client.containers.run(
            'openreactiondatabase/ord-postgres',
            ports={'5432/tcp': _POSTGRES_PORT},
            detach=True)

    @classmethod
    def tearDownClass(cls):
        cls._container.stop()

    def setUp(self):
        super().setUp()
        self.postgres = query.OrdPostgres(host='localhost',
                                          port=_POSTGRES_PORT)
        while True:
            try:
                with self.postgres:
                    break
            except psycopg2.OperationalError:
                logging.info('waiting for database to be ready')
                time.sleep(1)
                continue

    @parameterized.named_parameters(('smiles', 'C', False),
                                    ('smarts', '[#6]', True))
    def test_simple(self, pattern, use_smarts):
        with self.postgres as postgres:
            results = postgres.substructure_search(pattern=pattern,
                                                   table='rdk.inputs',
                                                   limit=100,
                                                   use_smarts=use_smarts)
        self.assertLen(results.reactions, 100)
        reaction_ids = [reaction.reaction_id for reaction in results.reactions]
        # Check that we remove redundant reaction IDs.
        self.assertCountEqual(reaction_ids, np.unique(reaction_ids))

    def test_bad_smiles(self):
        with self.postgres as postgres:
            with self.assertRaisesRegex(psycopg2.errors.DataException,
                                        'could not create molecule'):
                postgres.substructure_search('invalid', 'rdk.inputs')


if __name__ == '__main__':
    absltest.main()
