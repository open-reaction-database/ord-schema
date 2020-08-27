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
        num_attempts = 0
        while True:
            num_attempts += 1
            if num_attempts > 30:
                raise RuntimeError('failed to connect to the database')
            try:
                cls.postgres = query.OrdPostgres(
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

    @classmethod
    def tearDownClass(cls):
        cls._container.stop()

    def test_reaction_id_query(self):
        reaction_ids = [
            'ord-00386496302144278c262d284c3bc9f0',
            'ord-003de5cd29a541bdb3279081ebdefd06',
            'ord-00469722514e4e148db9fae89586f4a6'
        ]
        command = query.ReactionIdQuery(reaction_ids)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertCountEqual(results.reaction_ids, reaction_ids)

    def test_reaction_smarts_query(self):
        pattern = '[#6]>>[#7]'
        command = query.ReactionSmartsQuery(pattern)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results.reaction_ids, 10)

    def test_substructure_query(self):
        pattern = 'C'
        mode = query.ReactionComponentPredicate.MatchMode.SUBSTRUCTURE
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results.reaction_ids, 10)
        # Check that we remove redundant reaction IDs.
        self.assertCountEqual(results.reaction_ids,
                              np.unique(results.reaction_ids))

    def test_smarts_query(self):
        pattern = '[#6]'
        mode = query.ReactionComponentPredicate.MatchMode.SMARTS
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results.reaction_ids, 10)
        # Check that we remove redundant reaction IDs.
        self.assertCountEqual(results.reaction_ids,
                              np.unique(results.reaction_ids))

    def test_similarity_query(self):
        pattern = 'CC=O'
        mode = query.ReactionComponentPredicate.MatchMode.SIMILAR
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates,
                                               tanimoto_threshold=0.5)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertEmpty(results.reaction_ids)
        command = query.ReactionComponentQuery(predicates,
                                               tanimoto_threshold=0.05)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results.reaction_ids, 10)
        # Check that we remove redundant reaction IDs.
        self.assertCountEqual(results.reaction_ids,
                              np.unique(results.reaction_ids))

    def test_bad_smiles(self):
        pattern = 'invalid_smiles'
        mode = query.ReactionComponentPredicate.MatchMode.SUBSTRUCTURE
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates)
        with self.assertRaisesRegex(psycopg2.errors.DataException,
                                    'could not create molecule'):
            self.postgres.run_query(command)


if __name__ == '__main__':
    absltest.main()
