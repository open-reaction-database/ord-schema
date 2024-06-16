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

from ord_schema.orm.queries import (
    DatasetIdQuery,
    DoiQuery,
    RandomSampleQuery,
    ReactionComponentQuery,
    ReactionConversionQuery,
    ReactionIdQuery,
    ReactionSmartsQuery,
    ReactionYieldQuery,
    run,
)


@pytest.mark.skip("tsm_system_rows is not part of testing.postgresql")
def test_random_sample_query(test_cursor):
    query = RandomSampleQuery(16)
    results = run(test_cursor, query, return_ids=True)
    assert len(results) == 16


def test_dataset_id_query(test_cursor):
    dataset_ids = ["test_dataset"]
    with pytest.raises(ValueError, match="Invalid dataset ID"):
        DatasetIdQuery(dataset_ids)
    query = DatasetIdQuery(dataset_ids, validate=False)
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_id_query(test_cursor):
    reaction_ids = ["test_reaction-79"]
    with pytest.raises(ValueError, match="Invalid reaction ID"):
        ReactionIdQuery(reaction_ids)
    query = ReactionIdQuery(reaction_ids, validate=False)
    results = run(test_cursor, query, limit=10, return_ids=False)
    assert [result.reaction_id for result in results] == reaction_ids


def test_reaction_smarts_query(test_cursor):
    query = ReactionSmartsQuery("[#6]>>[#7]")
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_conversion_query(test_cursor):
    query = ReactionConversionQuery(min_conversion=50, max_conversion=90)
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_yield_query(test_cursor):
    query = ReactionYieldQuery(min_yield=50, max_yield=90)
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_doi_query(test_cursor):
    dois = ["10.1126/science.1255525"]
    query = DoiQuery(dois)
    results = run(test_cursor, query, limit=10)
    assert len(results) == 10
    for result in results:
        assert result.reaction.provenance.doi in dois


def test_exact_query(test_cursor):
    pattern = "[Br]C1=CC=C(C(C)=O)C=C1"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            target=ReactionComponentPredicate.Target.INPUT,
            mode=ReactionComponentPredicate.MatchMode.EXACT,
        )
    ]
    query = ReactionComponentQuery(predicates)
    results = run(test_cursor, query, limit=5, return_ids=True)
    assert len(results) == 5
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_substructure_query(test_cursor):
    pattern = "C"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            target=ReactionComponentPredicate.Target.INPUT,
            mode=ReactionComponentPredicate.MatchMode.SUBSTRUCTURE,
        )
    ]
    query = ReactionComponentQuery(predicates)
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_smarts_query(test_cursor):
    pattern = "[#6]"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            target=ReactionComponentPredicate.Target.INPUT,
            mode=ReactionComponentPredicate.MatchMode.SMARTS,
        )
    ]
    query = ReactionComponentQuery(predicates)
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_similarity_query(test_cursor):
    pattern = "CC=O"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            ReactionComponentPredicate.MatchMode.SMARTS,
            mode=ReactionComponentPredicate.MatchMode.SIMILAR,
        )
    ]
    query = ReactionComponentQuery(predicates, tanimoto_threshold=0.5)
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert not results
    query = ReactionComponentQuery(predicates, tanimoto_threshold=0.05)
    results = run(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_bad_smiles(test_cursor):
    pattern = "invalid_smiles"
    predicates = [
        ReactionComponentPredicate(
            pattern,
            ReactionComponentPredicate.MatchMode.SIMILAR,
            mode=ReactionComponentPredicate.MatchMode.SUBSTRUCTURE,
        )
    ]
    query = ReactionComponentQuery(predicates)
    with pytest.raises(QueryException, match="cannot parse pattern: invalid_smiles"):
        run(test_cursor, query)
