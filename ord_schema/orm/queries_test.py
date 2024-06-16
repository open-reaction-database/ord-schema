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
"""Tests for ord_schema.orm.queries."""
import numpy as np
import pytest

import ord_interface
from ord_interface.client.queries import (
    OrdPostgres,
    ReactionIdQuery,
    ReactionComponentQuery,
    ReactionYieldQuery,
    ReactionSmartsQuery,
    ReactionConversionQuery,
    RandomSampleQuery,
    DoiQuery,
    DatasetIdQuery,
)


@pytest.fixture(scope="module")
def connection() -> OrdPostgres:
    yield OrdPostgres(
        dbname=ord_interface.client.POSTGRES_DB,
        user=ord_interface.client.POSTGRES_USER,
        password=ord_interface.client.POSTGRES_PASSWORD,
        host="localhost",
        port=ord_interface.client.POSTGRES_PORT,
    )


def test_random_sample_query(connection):
    command = RandomSampleQuery(16)
    results = connection.run_query(command, return_ids=True)
    assert len(results) == 16


def test_dataset_id_query(connection):
    dataset_ids = ["ord_dataset-89b083710e2d441aa0040c361d63359f"]
    command = DatasetIdQuery(dataset_ids)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_id_query(connection):
    reaction_ids = ["ord-cf0d04017ede4c8aab8a15119c53e57b"]
    command = ReactionIdQuery(reaction_ids)
    results = connection.run_query(command, limit=10, return_ids=False)
    assert [result.reaction_id for result in results] == reaction_ids


def test_reaction_smarts_query(connection):
    pattern = "[#6]>>[#7]"
    command = ReactionSmartsQuery(pattern)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_conversion_query(connection):
    command = ReactionConversionQuery(min_conversion=50, max_conversion=90)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_yield_query(connection):
    command = ReactionYieldQuery(min_yield=50, max_yield=90)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_doi_query(connection):
    dois = ["10.1126/science.1255525"]
    command = DoiQuery(dois)
    results = connection.run_query(command, limit=10)
    assert len(results) == 10
    for result in results:
        assert result.reaction.provenance.doi in dois


def test_exact_query(connection):
    pattern = "[Br]C1=CC=C(C(C)=O)C=C1"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            target=ReactionComponentPredicate.Target.INPUT,
            mode=ReactionComponentPredicate.MatchMode.EXACT,
        )
    ]
    command = ReactionComponentQuery(predicates)
    results = connection.run_query(command, limit=5, return_ids=True)
    assert len(results) == 5
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_substructure_query(connection):
    pattern = "C"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            target=ReactionComponentPredicate.Target.INPUT,
            mode=ReactionComponentPredicate.MatchMode.SUBSTRUCTURE,
        )
    ]
    command = ReactionComponentQuery(predicates)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_smarts_query(connection):
    pattern = "[#6]"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            target=ReactionComponentPredicate.Target.INPUT,
            mode=ReactionComponentPredicate.MatchMode.SMARTS,
        )
    ]
    command = ReactionComponentQuery(predicates)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_similarity_query(connection):
    pattern = "CC=O"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            ReactionComponentPredicate.MatchMode.SMARTS,
            mode=ReactionComponentPredicate.MatchMode.SIMILAR,
        )
    ]
    command = ReactionComponentQuery(predicates, tanimoto_threshold=0.5)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert not results
    command = ReactionComponentQuery(predicates, tanimoto_threshold=0.05)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_bad_smiles(connection):
    pattern = "invalid_smiles"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            ReactionComponentPredicate.MatchMode.SIMILAR,
            mode=ReactionComponentPredicate.MatchMode.SUBSTRUCTURE,
        )
    ]
    command = ReactionComponentQuery(predicates)
    with pytest.raises(QueryException, match="cannot parse pattern: invalid_smiles"):
        connection.run_query(command)
