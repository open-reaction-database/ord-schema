# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for ord_schema.agent.sql."""

import pytest

from ord_schema.agent import query, sql

_LAMBDA_QUERY = (
    "SELECT reaction_id FROM reactions WHERE len(list_filter(flatten(list_transform("
    "map_values(inputs), i -> i.components)), c -> c.smiles = $thf)) > 0"
)
_UNNEST_QUERY = (
    "SELECT DISTINCT reaction_id FROM reactions, UNNEST(map_values(inputs)) t(i), "
    "UNNEST(i.components) u(c) WHERE c.smiles = $thf"
)


@pytest.mark.parametrize(
    ("query", "parameters"),
    [
        ("SELECT reaction_id FROM reactions", None),
        (
            "SELECT reaction_id FROM reactions "
            "WHERE conditions.temperature.setpoint_kelvin > 300",
            None,
        ),
        (
            "SELECT provenance.record_created.person.organization, count(*) "
            "FROM reactions GROUP BY 1",
            None,
        ),
        (_LAMBDA_QUERY, {"thf": "C1CCOC1"}),
        # UNNEST in the select list projects the elements; it is the FROM clause form
        # that materializes the exploded intermediate.
        ("SELECT unnest(map_keys(inputs)) FROM reactions", None),
    ],
)
def test_valid(query, parameters):
    sql.validate(query, parameters=parameters)


@pytest.mark.parametrize(
    ("query", "parameters", "match"),
    [
        ("SELECT FROM WHERE", None, "does not parse"),
        (
            "SELECT 1 FROM reactions; SELECT 2 FROM reactions",
            None,
            "single statement",
        ),
        ("DROP TABLE reactions", None, "expected a SELECT"),
        ("SELECT nope FROM reactions", None, "cannot be planned"),
        # A placeholder the model referenced but did not declare: caught before the
        # query reaches data, so the translation can be retried.
        (
            "SELECT reaction_id FROM reactions WHERE smiles = $missing",
            None,
            "cannot be planned",
        ),
        (_UNNEST_QUERY, {"thf": "C1CCOC1"}, "UNNEST in a FROM clause"),
    ],
)
def test_invalid(query, parameters, match):
    with pytest.raises(sql.InvalidQueryError, match=match):
        sql.validate(query, parameters=parameters)


@pytest.mark.parametrize(
    "query",
    [
        "SELECT * FROM read_parquet('/etc/hosts')",
        "SELECT * FROM read_csv('/etc/hosts')",
    ],
)
def test_no_file_access(query):
    with pytest.raises(sql.InvalidQueryError, match="cannot be planned"):
        sql.validate(query)


def test_placeholders_are_bound_not_interpolated():
    # A value that would end the string literal and append a statement if it were
    # interpolated. Bound, it is just a string that matches nothing.
    sql.validate(
        "SELECT reaction_id FROM reactions WHERE smiles = $name",
        parameters={"name": "'; DROP TABLE reactions; --"},
    )


def test_validation_needs_no_data():
    # The whole point: the schema is generated, so an empty table has the real shape
    # and planning resolves every column without a corpus.
    sql.validate(
        "SELECT reaction_id FROM reactions "
        "WHERE list_contains(list_transform(coalesce(outcomes, []), "
        "o -> o.reaction_time_seconds), 3600.0)"
    )


@pytest.mark.parametrize(
    "query",
    [
        "SELECT a.reaction_id FROM reactions a, reactions b WHERE a.smiles = b.smiles",
        "SELECT reaction_id FROM reactions, range(1000000000)",
        "WITH RECURSIVE t(n) AS ("
        "  SELECT 1 UNION ALL SELECT n + 1 FROM t WHERE n < 1000000000"
        ") SELECT count(*) FROM t",
        "SELECT * FROM reactions",
    ],
)
def test_cost_is_not_bounded(query):
    # Pinning the documented limit rather than a behavior worth having: refusing UNNEST
    # is one entry in a blacklist, and the set of expensive programs expressible in SQL
    # cannot be listed. A caller needs its own timeout or row cap; if a later change
    # does start bounding cost, this test should fail and be rewritten, not deleted.
    sql.validate(query)


@pytest.mark.parametrize("padding", [0, 11, 30])
def test_unnest_is_refused_however_wide_the_plan(padding):
    # The default EXPLAIN draws the plan into a fixed-width canvas and drops whole
    # subtrees once a query has more parallel pipelines than fit -- no ellipsis, no
    # error. Eleven UNION ALL branches were enough to make the UNNEST vanish from it
    # and the query pass, so the check reads the JSON plan instead.
    filler = "SELECT reaction_id FROM reactions WHERE reaction_id = 'nomatch'"
    branches = [filler] * padding + [_UNNEST_QUERY]
    with pytest.raises(sql.InvalidQueryError, match="UNNEST in a FROM clause"):
        sql.validate(" UNION ALL ".join(branches), parameters={"thf": "C1CCOC1"})


def test_an_operator_name_in_a_literal_is_not_an_operator():
    # The plan carries filter expressions alongside operator names, so a substring
    # search over the JSON would refuse this.
    sql.validate("SELECT reaction_id FROM reactions WHERE smiles = 'INOUT_FUNCTION'")


def test_a_structure_query_validates_against_the_executable_schema():
    # A compiled structure predicate references the executor's offset column, which
    # the bare projection schema cannot bind.
    compiled = query.compile_query(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "inputs.components",
                    "where": {
                        "op": "substructure",
                        "path": "smiles",
                        "smarts": "c1ccncc1",
                    },
                }
            }
        )
    )
    parameters = {compiled.structures[0].name: "0"}
    with pytest.raises(sql.InvalidQueryError, match="structure_offset"):
        sql.validate(compiled.sql, parameters=parameters)
    sql.validate(compiled.sql, parameters=parameters, schema=query.executable_schema())
