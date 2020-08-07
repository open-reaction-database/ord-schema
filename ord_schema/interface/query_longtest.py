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

from ord_schema import interface
from ord_schema.interface import query


class QueryTest(parameterized.TestCase, absltest.TestCase):

    @classmethod
    def setUpClass(cls):
        client = docker.from_env()
        client.images.pull('openreactiondatabase/ord-postgres')
        cls._container = client.containers.run(
            'openreactiondatabase/ord-postgres',
            ports={'5432/tcp': interface.POSTGRES_PORT},
            detach=True,
            remove=True)

    @classmethod
    def tearDownClass(cls):
        cls._container.stop()

    def setUp(self):
        super().setUp()
        num_attempts = 0
        while True:
            num_attempts += 1
            if num_attempts > 30:
                raise RuntimeError('failed to connect to the database')
            try:
                self.postgres = query.OrdPostgres(
                    dbname=interface.POSTGRES_DB,
                    user=interface.POSTGRES_USER,
                    password=interface.POSTGRES_PASSWORD,
                    host='localhost',
                    port=interface.POSTGRES_PORT)
                break
            except psycopg2.OperationalError as error:
                logging.info('waiting for database to be ready: %s', error)
                time.sleep(1)
                continue

    def test_predicate_search(self):
        # No constraint.
        predicate = query.Predicate()
        results = self.postgres.predicate_search(predicate)
        self.assertLen(len(results.reactions) > 100)

        # Look up some specific reactions.
        predicate = query.Predicate()
        predicate.add_reaction_id('ord-00386496302144278c262d284c3bc9f0')
        predicate.add_reaction_id('ord-003de5cd29a541bdb3279081ebdefd06')
        predicate.add_reaction_id('ord-00469722514e4e148db9fae89586f4a6')
        results = self.postgres.predicate_search(predicate)
        self.assertLen(results.reactions, 3)

        # Match on input substructure.
        predicate = query.Predicate()
        predicate.add_input('CCC1=CC(=CC=C1)CC',
                            query.Predicate.MatchMode.SUBSTRUCTURE)
        results = self.postgres.predicate_search(predicate)
        self.assertTrue(len(results.reactions) > 0)

        # Match on output similarity with a specified threshold.
        predicate = query.Predicate()
        predicate.add_output(
            'CC1=CC=C2C(C=NN2C3OCCCC3)=C1C4=CC=C(N=CC=C5)C5=C4',
            query.Predicate.MatchMode.SIMILAR)
        predicate.set_tanimoto_threshold(0.6)
        results = self.postgres.predicate_search(predicate)
        self.assertTrue(len(results.reactions) > 0)

    @parameterized.named_parameters(('smiles', 'C', False),
                                    ('smarts', '[#6]', True))
    def test_substructure_search(self, pattern, use_smarts):
        results = self.postgres.substructure_search(pattern=pattern,
                                                    table='rdk.inputs',
                                                    limit=100,
                                                    use_smarts=use_smarts)
        self.assertLen(results.reactions, 100)
        reaction_ids = [reaction.reaction_id for reaction in results.reactions]
        # Check that we remove redundant reaction IDs.
        self.assertCountEqual(reaction_ids, np.unique(reaction_ids))

    def test_similarity_search(self):
        results = self.postgres.similarity_search(smiles='CC=O',
                                                  table='rdk.inputs',
                                                  limit=100,
                                                  threshold=0.5)
        self.assertEmpty(results.reactions)
        results = self.postgres.similarity_search(smiles='CC=O',
                                                  table='rdk.inputs',
                                                  limit=100,
                                                  threshold=0.05)
        self.assertLen(results.reactions, 100)
        reaction_ids = [reaction.reaction_id for reaction in results.reactions]
        # Check that we remove redundant reaction IDs.
        self.assertCountEqual(reaction_ids, np.unique(reaction_ids))

    def test_bad_smiles(self):
        with self.assertRaisesRegex(psycopg2.errors.DataException,
                                    'could not create molecule'):
            self.postgres.substructure_search('invalid', 'rdk.inputs')
        with self.assertRaisesRegex(psycopg2.errors.DataException,
                                    'could not create molecule'):
            self.postgres.similarity_search('invalid', 'rdk.inputs')


if __name__ == '__main__':
    absltest.main()
