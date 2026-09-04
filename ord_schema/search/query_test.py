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

"""Tests for ord_schema.search.query."""

import warnings

import duckdb
import pytest
from pydantic import ValidationError
from rdkit import Chem

from ord_schema.artifacts import pivot, projection
from ord_schema.search import query, sql


def _compile(payload, **kwargs):
    return query.compile_query(query.Query.model_validate(payload), **kwargs)


def _run(compiled, parameters=None):
    """Executes a compiled query against an empty projection, proving it is runnable."""
    connection = duckdb.connect()
    connection.register(query.TABLE, projection.SCHEMA.empty_table())
    try:
        return connection.execute(compiled.sql, parameters or {}).fetchall()
    finally:
        connection.close()


# Path resolution


def test_a_scalar_path_reaches_a_leaf():
    resolved = query.resolve("conditions.temperature.setpoint_kelvin")
    assert resolved.expression == "conditions.temperature.setpoint_kelvin"
    assert not resolved.repeated


def test_a_map_becomes_a_list_of_its_values():
    resolved = query.resolve("inputs.components")
    assert resolved.expression == (
        "flatten(list_transform(map_values(inputs), x -> x.components))"
    )
    assert resolved.repeated


def test_a_path_through_two_repeated_levels_flattens_once_each():
    resolved = query.resolve("outcomes.products.measurements")
    assert resolved.repeated
    assert resolved.expression.count("flatten") == 2


def test_an_unknown_field_suggests_a_near_one():
    with pytest.raises(query.QueryError, match="did you mean 'temperature'"):
        query.resolve("conditions.temprature.setpoint_kelvin")


def test_an_unknown_field_with_no_near_match_still_names_the_path():
    with pytest.raises(query.QueryError, match="no field named 'zzz'"):
        query.resolve("conditions.zzz")


def test_descending_into_a_leaf_is_refused():
    with pytest.raises(query.QueryError, match="no fields to descend into"):
        query.resolve("reaction_id.nope")


# Compilation


def test_a_comparison_compiles_to_a_scalar_test():
    compiled = _compile(
        {
            "where": {
                "op": "gt",
                "path": "conditions.temperature.setpoint_kelvin",
                "value": {"literal": 300.0},
            }
        }
    )
    assert compiled.sql == (
        "SELECT reaction_id FROM reactions "
        "WHERE conditions.temperature.setpoint_kelvin > 300.0"
    )
    assert compiled.compounds == ()


def test_exists_compiles_to_a_list_lambda_never_to_unnest():
    # The whole reason the model does not write SQL: the fast idiom is the only one the
    # compiler can emit, so the 27-200x cliff is unreachable rather than discouraged.
    compiled = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {"op": "eq", "path": "smiles", "value": {"literal": "CCO"}},
            }
        }
    )
    assert "list_filter" in compiled.sql
    assert "UNNEST" not in compiled.sql.upper()
    sql.validate(compiled.sql)


def test_forall_is_the_negation_of_a_counterexample():
    compiled = _compile(
        {
            "where": {
                "op": "forall",
                "path": "outcomes.products",
                "where": {"op": "not_null", "path": "smiles"},
            }
        }
    )
    assert "NOT (" in compiled.sql
    # Vacuously true of a level with no elements, which the projection writes as NULL,
    # so the count comparison is folded to what an absent level means rather than left
    # to answer neither way.
    assert compiled.sql.endswith(")) = 0, true)")


def test_a_nested_quantifier_binds_its_own_element():
    # The bug this pins: an inner path resolved against the row instead of the bound
    # element silently reads the reaction's own column and answers the wrong question.
    compiled = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "exists",
                    "path": "identifiers",
                    "where": {
                        "op": "eq",
                        "path": "value",
                        "value": {"literal": "THF"},
                    },
                },
            }
        }
    )
    assert "e0.identifiers" in compiled.sql
    assert "e1.value" in compiled.sql


def test_co_membership_and_mere_co_occurrence_compile_differently():
    both = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "and",
                    "clauses": [
                        {"op": "eq", "path": "smiles", "value": {"literal": "CCO"}},
                        {
                            "op": "eq",
                            "path": "reaction_role",
                            "value": {"literal": "SOLVENT"},
                        },
                    ],
                },
            }
        }
    )
    either = _compile(
        {
            "where": {
                "op": "and",
                "clauses": [
                    {
                        "op": "exists",
                        "path": "inputs.components",
                        "where": {
                            "op": "eq",
                            "path": "smiles",
                            "value": {"literal": "CCO"},
                        },
                    },
                    {
                        "op": "exists",
                        "path": "inputs.components",
                        "where": {
                            "op": "eq",
                            "path": "reaction_role",
                            "value": {"literal": "SOLVENT"},
                        },
                    },
                ],
            }
        }
    )
    assert both.sql != either.sql
    assert both.sql.count("list_filter") == 1
    assert either.sql.count("list_filter") == 2


def test_compounds_are_bound_not_spelled():
    compiled = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {"op": "eq", "path": "smiles", "value": {"compound": "thf"}},
            }
        }
    )
    assert "$thf" in compiled.sql
    assert compiled.compounds == ("thf",)
    assert _run(compiled, {"thf": "C1CCOC1"}) == []


def test_a_string_literal_carrying_a_quote_cannot_close_it():
    compiled = _compile(
        {
            "where": {
                "op": "eq",
                "path": "reaction_id",
                "value": {"literal": "'; DROP TABLE reactions; --"},
            }
        }
    )
    assert _run(compiled) == []


# Aggregation


def test_an_aggregate_groups_and_measures():
    compiled = _compile(
        {
            "aggregate": {
                "group_by": ["conditions.stirring.type"],
                "measures": [
                    {
                        "fn": "avg",
                        "path": "conditions.temperature.setpoint_kelvin",
                        "name": "mean_k",
                    },
                    {"fn": "count", "name": "n"},
                ],
            },
            "order_by": [{"key": "mean_k", "descending": True}],
            "limit": 5,
        }
    )
    assert compiled.sql == (
        "SELECT conditions.stirring.type, "
        "avg(conditions.temperature.setpoint_kelvin) AS mean_k, count(*) AS n "
        "FROM reactions GROUP BY 1 ORDER BY mean_k DESC LIMIT 5"
    )
    assert _run(compiled) == []


def test_count_distinct_counts_distinctly():
    compiled = _compile(
        {
            "aggregate": {
                "measures": [
                    {"fn": "count_distinct", "path": "provenance.doi", "name": "dois"}
                ]
            }
        }
    )
    assert "count(DISTINCT provenance.doi)" in compiled.sql


# Refusals


@pytest.mark.parametrize(
    ("payload", "match"),
    [
        # The refusal that has no SQL equivalent: UNNEST would silently mean "any".
        (
            {
                "where": {
                    "op": "eq",
                    "path": "inputs.components.smiles",
                    "value": {"literal": "CCO"},
                }
            },
            "crosses a repeated level",
        ),
        (
            {
                "where": {
                    "op": "exists",
                    "path": "reaction_id",
                    "where": {"op": "not_null", "path": "reaction_id"},
                }
            },
            "needs a repeated level",
        ),
        (
            {
                "where": {
                    "op": "gt",
                    "path": "conditions.temperature.setpoint_kelvin",
                    "value": {"literal": "hot"},
                }
            },
            "compared against a string",
        ),
        (
            {
                "where": {
                    "op": "contains",
                    "path": "conditions.temperature.setpoint_kelvin",
                    "value": {"literal": "3"},
                }
            },
            "needs a text column",
        ),
        (
            {"where": {"op": "gt", "path": "reaction_id", "value": {"literal": 1}}},
            "needs a numeric column",
        ),
        # A compound resolves to a SMILES, so it can only compare against text.
        (
            {
                "where": {
                    "op": "eq",
                    "path": "conditions.temperature.setpoint_kelvin",
                    "value": {"compound": "thf"},
                }
            },
            "compares against text",
        ),
        (
            {
                "aggregate": {
                    "group_by": ["inputs.components.smiles"],
                    "measures": [{"fn": "count", "name": "n"}],
                }
            },
            "needs a scalar column",
        ),
        (
            {
                "aggregate": {"measures": [{"fn": "count", "name": "n"}]},
                "order_by": [{"key": "nope"}],
            },
            "measure name or a group_by path",
        ),
    ],
)
def test_refused(payload, match):
    with pytest.raises(query.QueryError, match=match):
        _compile(payload)


@pytest.mark.parametrize(
    "payload",
    [
        # A measure name reaches the SQL as text, so it is held to an identifier shape.
        {"aggregate": {"measures": [{"fn": "count", "name": "n; DROP TABLE x"}]}},
        {
            "where": {
                "op": "eq",
                "path": "reaction_id",
                "value": {"compound": "x'; DROP TABLE y; --"},
            }
        },
        # Neither a literal nor a compound, or both at once.
        {"where": {"op": "eq", "path": "reaction_id", "value": {}}},
        {
            "where": {
                "op": "eq",
                "path": "reaction_id",
                "value": {"literal": "a", "compound": "b"},
            }
        },
        {"limit": 0},
        {"aggregate": {"measures": [{"fn": "sum", "name": "n"}]}},
        {"where": {"op": "and", "clauses": []}},
    ],
)
def test_rejected_before_compilation(payload):
    with pytest.raises(ValidationError):
        query.Query.model_validate(payload)


@pytest.mark.parametrize(
    "payload",
    [
        # A misspelled top-level key: dropped, this is Query() and matches everything.
        {"structure": {"path": "inputs.components.smiles", "smarts": "c1ccccc1"}},
        {"wehre": {"op": "is_null", "path": "reaction_id"}},
        # Misspelled inside a predicate, where the surviving model is still valid: the
        # quantifier keeps its body but loses the pattern, and the substructure node
        # then has neither smarts nor compound.
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {"op": "substructure", "path": "smiles", "smart": "c1ccccc1"},
            }
        },
        {"where": {"op": "eq", "path": "reaction_id", "value": {"literl": "a"}}},
        {"limit": 10, "ordr_by": []},
    ],
)
def test_an_unknown_key_is_refused_rather_than_dropped(payload):
    # Dropping one does not fail here, it widens: every narrowing field is optional or
    # sits under a discriminated union, so the survivor is a valid query missing the
    # clause the caller wrote -- and a Query with no ``where`` is every reaction in the
    # corpus, returned without a warning.
    with pytest.raises(ValidationError, match=r"[Ee]xtra"):
        query.Query.model_validate(payload)


def test_a_query_with_no_where_matches_everything():
    # What the payloads above validate to once their unknown key is dropped, and why
    # refusing them is more than typo-catching: no predicate compiles to no WHERE.
    assert "WHERE" not in query.compile_query(query.Query()).sql


# Reductions over a repeated level


def test_ordering_by_a_reduction_over_a_repeated_path():
    # A yield lives under outcomes, products, and measurements, so ordering by it needs
    # the list reduced to the one number the reaction is judged by.
    compiled = _compile(
        {
            "order_by": [
                {
                    "key": {
                        "reduce": "max",
                        "path": "outcomes.products.measurements.percentage.value",
                    },
                    "descending": True,
                }
            ],
            "limit": 10,
        }
    )
    assert "list_max(" in compiled.sql
    assert compiled.sql.endswith("DESC LIMIT 10")


def test_max_rows_bounds_a_query_that_asked_for_no_limit():
    # The grammar leaves limit optional, so the unbounded query is the ordinary one: a
    # predicate matching much of the corpus returns millions of rows to whoever ran it.
    compiled = _compile({}, max_rows=100)
    assert compiled.sql.endswith(" LIMIT 100")
    assert compiled.limit == 100


def test_max_rows_clamps_a_query_asking_for_more():
    # A bound an explicit limit can exceed bounds nothing.
    compiled = _compile({"limit": 5000}, max_rows=100)
    assert compiled.sql.endswith(" LIMIT 100")
    assert compiled.limit == 100


def test_max_rows_leaves_a_smaller_limit_alone():
    compiled = _compile({"limit": 10}, max_rows=100)
    assert compiled.sql.endswith(" LIMIT 10")
    assert compiled.limit == 10


def test_a_limit_of_zero_is_not_read_as_no_limit():
    # The trap in `query.limit or max_rows`: a falsy limit reads as absent, so a query
    # asking for no rows would be served max_rows of them. Held off by the model's
    # gt=0 today, which is not where the compiler should rest.
    compiled = query.compile_query(query.Query.model_construct(limit=0), max_rows=100)
    assert compiled.limit == 0


def test_a_max_rows_no_query_can_satisfy_is_refused():
    # Zero compiles to a LIMIT nothing satisfies and a negative one to SQL DuckDB
    # refuses at execution, both far from whoever passed it.
    for value in (0, -5):
        with pytest.raises(ValueError, match="which no query can return"):
            _compile({}, max_rows=value)


def test_without_max_rows_a_query_is_as_stated():
    assert _compile({}).limit is None
    assert " LIMIT " not in _compile({}).sql
    assert _compile({"limit": 5000}).limit == 5000


def test_a_reduction_reaches_the_same_rows_the_elements_do():
    # The reduction is over the reaction's own elements, so it agrees with the list the
    # projection holds rather than with a corpus-wide aggregate.
    compiled = _compile(
        {
            "order_by": [
                {
                    "key": {
                        "reduce": "max",
                        "path": "outcomes.products.measurements.percentage.value",
                    }
                }
            ]
        }
    )
    resolved = query.resolve("outcomes.products.measurements.percentage.value")
    assert resolved.expression in compiled.sql


def test_a_measure_may_reduce_a_repeated_path():
    compiled = _compile(
        {
            "aggregate": {
                "group_by": ["conditions.temperature.setpoint_kelvin"],
                "measures": [
                    {
                        "fn": "avg",
                        "path": {
                            "reduce": "max",
                            "path": "outcomes.products.measurements.percentage.value",
                        },
                        "name": "best_yield",
                    }
                ],
            }
        }
    )
    assert "avg(list_max(" in compiled.sql


@pytest.mark.parametrize(
    ("reducer", "expected"),
    [
        ("min", "list_min("),
        ("max", "list_max("),
        ("avg", "list_avg("),
        ("sum", "list_sum("),
        # Counting what is there means filtering the nulls a list may hold, since len()
        # would count them.
        ("count", "len(list_filter("),
    ],
)
def test_each_reducer_compiles_to_its_list_aggregate(reducer, expected):
    compiled = _compile(
        {
            "order_by": [
                {
                    "key": {
                        "reduce": reducer,
                        "path": "outcomes.products.measurements.percentage.value",
                    }
                }
            ]
        }
    )
    assert expected in compiled.sql


def test_a_reduction_over_a_scalar_path_is_refused():
    # A scalar needs no reducing, and accepting one would give the same query two
    # spellings, one of which wraps a value in a single-element list.
    with pytest.raises(query.QueryError, match="already scalar"):
        _compile(
            {
                "order_by": [
                    {
                        "key": {
                            "reduce": "max",
                            "path": "conditions.temperature.setpoint_kelvin",
                        }
                    }
                ]
            }
        )


@pytest.mark.parametrize("reducer", ["min", "max", "avg", "sum"])
def test_an_arithmetic_reduction_over_text_is_refused(reducer):
    # Summing a list of strings is a DuckDB error rather than an answer, and a query
    # that cannot mean anything is worth refusing where it is compiled.
    with pytest.raises(query.QueryError, match="needs a numeric column"):
        _compile(
            {
                "order_by": [
                    {
                        "key": {
                            "reduce": reducer,
                            "path": "outcomes.products.identifiers.value",
                        }
                    }
                ]
            }
        )


def test_counting_a_repeated_text_path_is_allowed():
    # How many identifiers a reaction's products carry is a number, whatever the
    # identifiers themselves hold.
    compiled = _compile(
        {
            "order_by": [
                {
                    "key": {
                        "reduce": "count",
                        "path": "outcomes.products.identifiers.value",
                    }
                }
            ]
        }
    )
    assert "len(list_filter(" in compiled.sql


def test_an_arithmetic_measure_over_text_is_refused():
    # The same rule reaches a measure that names a scalar column directly.
    with pytest.raises(query.QueryError, match="needs a numeric column"):
        _compile(
            {
                "aggregate": {
                    "measures": [
                        {"fn": "avg", "path": "reaction_id", "name": "mean_id"}
                    ]
                }
            }
        )


def test_an_aggregated_query_cannot_order_by_a_reduction():
    # After grouping there is no reaction left to reduce over; the reduction belongs
    # inside a measure, where it is one input to the aggregate.
    with pytest.raises(query.QueryError, match="reduce inside a measure"):
        _compile(
            {
                "aggregate": {"measures": [{"fn": "count", "name": "n"}]},
                "order_by": [
                    {
                        "key": {
                            "reduce": "max",
                            "path": "outcomes.products.measurements.percentage.value",
                        }
                    }
                ],
            }
        )


def test_a_reduction_runs():
    # The expression has to be DuckDB the planner accepts, not merely a string.
    compiled = _compile(
        {
            "order_by": [
                {
                    "key": {
                        "reduce": "max",
                        "path": "outcomes.products.measurements.percentage.value",
                    },
                    "descending": True,
                }
            ],
            "limit": 5,
        }
    )
    assert _run(compiled) == []


# The whole point


def test_no_expressible_query_needs_a_join_or_recursion():
    # The cost bound is structural: the grammar has one relation and no way to name
    # another, no recursion, and no set-returning function. This asserts the surface
    # rather than any one query -- a node type that broke it would fail here.
    fields = set(query.Query.model_fields)
    assert fields == {"where", "aggregate", "order_by", "limit"}
    operators = set()
    for member in (
        query.And,
        query.Or,
        query.Not,
        query.Quantifier,
        query.Comparison,
        query.NullCheck,
    ):
        operators.update(member.model_fields["op"].annotation.__args__)
    assert operators == {
        "and",
        "or",
        "not",
        "exists",
        "forall",
        "eq",
        "ne",
        "lt",
        "le",
        "gt",
        "ge",
        "contains",
        "starts_with",
        "ends_with",
        "is_null",
        "not_null",
    }


@pytest.mark.parametrize(
    "table",
    [
        # The cost bound is the point of the grammar, so the one string a caller can
        # put into the FROM clause must not be able to add a relation to it.
        "reactions, range(1000000000)",
        "reactions AS a, reactions AS b",
        "read_parquet('/etc/hosts')",
    ],
)
def test_a_table_that_is_not_an_identifier_is_refused(table):
    with pytest.raises(query.QueryError, match="not an identifier"):
        query.compile_query(query.Query(), table=table)


def test_a_plain_relation_name_is_still_allowed():
    assert query.compile_query(query.Query(), table="projections").sql.endswith(
        "FROM projections"
    )


def test_an_internal_column_is_refused():
    # structure_id is real in the schema but internal to the artifacts: its values are
    # not stable across builds, so a comparison against one is refused at resolution.
    with pytest.raises(query.QueryError, match="internal"):
        query.resolve("inputs.components.structure_id")


def test_suggestions_do_not_name_internal_columns():
    # A near-miss like 'structure' is exactly what a model gropes for when asked a
    # structure question; the suggestion must not teach it the withheld name.
    with pytest.raises(query.QueryError, match="no field named") as excinfo:
        query.resolve("inputs.components.structure")
    assert "structure_id" not in str(excinfo.value)


# Structure predicates


def _substructure(smarts="c1ccncc1", path="smiles"):
    return {"op": "substructure", "path": path, "smarts": smarts}


def test_substructure_compiles_to_a_bitmap_test():
    compiled = query.compile_query(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "inputs.components",
                    "where": _substructure(),
                }
            }
        )
    )
    assert "get_bit(CAST($structure_0 AS BITSTRING)" in compiled.sql
    assert "e0.structure_id IS NOT NULL" in compiled.sql
    assert f"+ {query.STRUCTURE_OFFSET})::INTEGER" in compiled.sql
    (parameter,) = compiled.structures
    assert parameter.op == "substructure"
    assert parameter.pattern == "c1ccncc1"
    assert parameter.compound is None


def test_equal_structure_predicates_share_one_parameter():
    # The executor screens and verifies each distinct predicate once, however many
    # times the query states it.
    compiled = query.compile_query(
        query.Query.model_validate(
            {
                "where": {
                    "op": "or",
                    "clauses": [
                        {
                            "op": "exists",
                            "path": "inputs.components",
                            "where": _substructure(),
                        },
                        {
                            "op": "exists",
                            "path": "outcomes.products",
                            "where": _substructure(),
                        },
                    ],
                }
            }
        )
    )
    assert len(compiled.structures) == 1
    assert compiled.sql.count("$structure_0") == 2


def test_similarity_compiles_with_its_threshold():
    compiled = query.compile_query(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "inputs.components",
                    "where": {
                        "op": "similarity",
                        "path": "smiles",
                        "smiles": "CCO",
                        "threshold": 0.7,
                    },
                }
            }
        )
    )
    (parameter,) = compiled.structures
    assert parameter.op == "similarity"
    assert parameter.pattern == "CCO"
    assert parameter.threshold == 0.7


def test_a_smarts_with_an_explicit_hydrogen_is_merged_with_a_warning():
    # Measured: [H]OC fails the fingerprint screen against methanol even though the
    # merged query matches it, so the hydrogens are folded into heavy-atom H counts
    # rather than silently missing true hits.
    with pytest.warns(UserWarning, match="rewritten"):
        model = query.Substructure.model_validate(
            {"op": "substructure", "path": "smiles", "smarts": "[H]OC"}
        )
    merged = Chem.MolFromSmarts(model.smarts)
    assert merged is not None
    assert not any(atom.GetAtomicNum() == 1 for atom in merged.GetAtoms())
    assert Chem.MolFromSmiles("CO").HasSubstructMatch(merged)


@pytest.mark.parametrize("smarts", ["[H][H]", "[2H]OC", "[2H]"])
def test_an_unfoldable_explicit_hydrogen_passes_through(smarts):
    # These hydrogens cannot be implicit in a stored molecule either -- an isotope and
    # H2 have no heavy atom to hide in -- so they are real graph atoms the query
    # already matches. 28,297 ORD reactions have an [H][H] component, so rewriting or
    # refusing them would lose exactly the queries that work.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        model = query.Substructure.model_validate(
            {"op": "substructure", "path": "smiles", "smarts": smarts}
        )
    assert model.smarts == smarts


@pytest.mark.parametrize("smarts", ["", "[]"])
def test_a_smarts_with_no_atoms_is_refused(smarts):
    # An empty query fingerprints to no bits, so the screen admits the whole corpus
    # and verification then rejects all of it: an empty answer at full cost.
    with pytest.raises(ValueError, match=r"no atoms|does not parse"):
        query.Substructure.model_validate(
            {"op": "substructure", "path": "smiles", "smarts": smarts}
        )


def test_a_similarity_smiles_with_no_atoms_is_refused():
    with pytest.raises(ValueError, match="no atoms"):
        query.Similarity.model_validate(
            {"op": "similarity", "path": "smiles", "smiles": "", "threshold": 0.5}
        )


def test_an_h_count_smarts_is_accepted():
    query.Substructure.model_validate(
        {"op": "substructure", "path": "smiles", "smarts": "[OX2H]C"}
    )


def test_an_unparseable_smarts_is_refused():
    with pytest.raises(ValueError, match="does not parse"):
        query.Substructure.model_validate(
            {"op": "substructure", "path": "smiles", "smarts": "not-smarts"}
        )


def test_a_substructure_on_a_non_smiles_column_is_refused():
    with pytest.raises(query.QueryError, match="smiles"):
        query.compile_query(
            query.Query.model_validate(
                {
                    "where": {
                        "op": "exists",
                        "path": "inputs.components",
                        "where": _substructure(path="reaction_role"),
                    }
                }
            )
        )


def test_a_substructure_on_the_reaction_smiles_is_refused():
    # The reaction-level smiles is a reaction, not a molecule; it has no structure ID
    # and reaction substructure search is a different operation.
    with pytest.raises(query.QueryError, match="no structure ID"):
        query.compile_query(
            query.Query.model_validate({"where": _substructure(path="smiles")})
        )


def test_an_unquantified_substructure_is_refused():
    with pytest.raises(query.QueryError, match="exists or forall"):
        query.compile_query(
            query.Query.model_validate(
                {"where": _substructure(path="inputs.components.smiles")}
            )
        )


def test_a_similarity_threshold_is_required_and_bounded():
    with pytest.raises(ValueError, match="threshold"):
        query.Similarity.model_validate(
            {"op": "similarity", "path": "smiles", "smiles": "CCO"}
        )
    with pytest.raises(ValueError, match="less than or equal to 1"):
        query.Similarity.model_validate(
            {"op": "similarity", "path": "smiles", "smiles": "CCO", "threshold": 1.5}
        )


def test_a_compound_named_like_a_structure_parameter_is_refused():
    with pytest.raises(query.QueryError, match="collide"):
        query.compile_query(
            query.Query.model_validate(
                {
                    "where": {
                        "op": "exists",
                        "path": "inputs.components",
                        "where": {
                            "op": "and",
                            "clauses": [
                                _substructure(),
                                {
                                    "op": "eq",
                                    "path": "smiles",
                                    "value": {"compound": "structure_0"},
                                },
                            ],
                        },
                    }
                }
            )
        )


def _pivot_sql(body, pivot=None):
    return query.compile_query(query.Query.model_validate(body), pivot=pivot).sql


def _every_level(path):
    return pivot.table_name(path) if path in pivot.LEVELS else None


def test_exists_over_a_pivoted_level_becomes_a_semi_join():
    sql = _pivot_sql(
        {
            "where": {
                "op": "exists",
                "path": "workups",
                "where": {
                    "op": "eq",
                    "path": "type",
                    "value": {"literal": "EXTRACTION"},
                },
            }
        },
        pivot=_every_level,
    )
    assert "EXISTS (SELECT 1 FROM pivot_workups AS x0" in sql
    # Qualified on both sides. A pivot carries a reaction_id of its own, so an
    # unqualified outer reference binds to the inner one and the correlation compares a
    # column to itself -- true of every reaction that has any element at all.
    assert "x0.reaction_id = reactions.reaction_id" in sql
    assert "x0.element.type" in sql
    assert "list_filter" not in sql


def test_forall_over_a_pivoted_level_becomes_a_negated_semi_join():
    sql = _pivot_sql(
        {
            "where": {
                "op": "forall",
                "path": "workups",
                "where": {
                    "op": "eq",
                    "path": "type",
                    "value": {"literal": "EXTRACTION"},
                },
            }
        },
        pivot=_every_level,
    )
    assert sql.count("NOT EXISTS (SELECT 1 FROM pivot_workups AS x0") == 1
    assert "NOT IN" not in sql
    assert "list_filter" not in sql


def test_a_nested_quantifier_joins_on_the_whole_ordinal_prefix():
    # The correctness point the ordinals exist for. Correlating a measurement to its
    # product on the reaction alone returned 23,608 reactions over ORD where the answer
    # is 22,666 -- and dropping only the outcome ordinal changed nothing at all, since
    # the corpus is effectively single-outcome, which is what made it look right.
    sql = _pivot_sql(
        {
            "where": {
                "op": "exists",
                "path": "outcomes.products",
                "where": {
                    "op": "and",
                    "clauses": [
                        {
                            "op": "eq",
                            "path": "is_desired_product",
                            "value": {"literal": True},
                        },
                        {
                            "op": "exists",
                            "path": "measurements",
                            "where": {
                                "op": "eq",
                                "path": "type",
                                "value": {"literal": "YIELD"},
                            },
                        },
                    ],
                },
            }
        },
        pivot=_every_level,
    )
    assert "FROM pivot_outcomes_products_measurements AS x1" in sql
    assert "x1.reaction_id = x0.reaction_id" in sql
    assert "x1.outcome_index = x0.outcome_index" in sql
    assert "x1.product_index = x0.product_index" in sql
    assert "list_filter" not in sql


def test_a_quantifier_inside_a_list_lambda_is_left_to_the_elements():
    # With the outer level unavailable there is no alias to correlate to, so the inner
    # quantifier compiles over the elements rather than against a pivot keyed by
    # nothing.
    sql = _pivot_sql(
        {
            "where": {
                "op": "exists",
                "path": "outcomes.products",
                "where": {
                    "op": "exists",
                    "path": "measurements",
                    "where": {
                        "op": "eq",
                        "path": "type",
                        "value": {"literal": "YIELD"},
                    },
                },
            }
        },
        pivot=lambda path: None,
    )
    assert "pivot_" not in sql
    assert sql.count("list_filter") == 2


def test_a_declined_body_leaves_no_parameters_behind():
    # The structure clause compiles -- and names its parameter -- before the nested
    # quantifier reaches a field the pruned type dropped. That name must not survive
    # the decline, or the query would carry a parameter nothing binds.
    compiled = query.compile_query(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "outcomes.products",
                    "where": {
                        "op": "and",
                        "clauses": [
                            {
                                "op": "substructure",
                                "path": "smiles",
                                "smarts": "c1ccncc1",
                            },
                            {
                                "op": "exists",
                                "path": "measurements",
                                "where": {
                                    "op": "eq",
                                    "path": "type",
                                    "value": {"literal": "YIELD"},
                                },
                            },
                        ],
                    },
                }
            }
        ),
        # The measurements level is refused -- a budget may do exactly this -- so the
        # inner quantifier declines, the outer body fails to resolve against the pruned
        # type, and the outer declines with it.
        pivot=lambda path: (
            None if path == "outcomes.products.measurements" else _every_level(path)
        ),
    )
    assert "pivot_" not in compiled.sql
    assert len(compiled.structures) == 1


def test_without_a_pivot_the_sql_is_unchanged():
    body = {
        "where": {
            "op": "exists",
            "path": "workups",
            "where": {
                "op": "eq",
                "path": "type",
                "value": {"literal": "EXTRACTION"},
            },
        }
    }
    assert _pivot_sql(body) == _pivot_sql(body, pivot=lambda path: None)


def test_a_structure_predicate_routes_through_the_pivot_with_the_outer_offset():
    # The occurrence index declines a body that is not one structure predicate plus
    # string equalities, so this lands on the pivot. Its rows carry structure_id but no
    # structure_offset, and the bitmap is indexed by the two added: the offset resolves
    # outward to the reaction the element belongs to, which is the one its ID is meant
    # to be read against. Pinned here because it holds by name resolution -- a pivot
    # that ever carried an offset column of its own would bind inward instead.
    compiled = query.compile_query(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "inputs.components",
                    "where": {
                        "op": "and",
                        "clauses": [
                            {
                                "op": "substructure",
                                "path": "smiles",
                                "smarts": "c1ccncc1",
                            },
                            {"op": "not_null", "path": "smiles"},
                        ],
                    },
                }
            }
        ),
        pivot=_every_level,
    )
    assert "EXISTS (SELECT 1 FROM pivot_inputs_components AS x0" in compiled.sql
    assert "x0.element.structure_id" in compiled.sql
    # Unqualified, so it is the enclosing relation's column rather than one of the
    # pivot's; the pivot has none.
    assert "+ structure_offset)" in compiled.sql
    assert "x0.structure_offset" not in compiled.sql
    assert len(compiled.structures) == 1


def test_a_declining_body_is_never_charged_for_a_pivot():
    # Supplying a pivot is what builds it, and over ORD that is minutes per level. A
    # body reaching a repeated field declines, so it must not ask for one first.
    asked = []

    def watched(path):
        asked.append(path)
        return _every_level(path)

    query.compile_query(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "workups",
                    "where": {
                        "op": "exists",
                        "path": "input.components",
                        "where": {
                            "op": "eq",
                            "path": "reaction_role",
                            "value": {"literal": "SOLVENT"},
                        },
                    },
                }
            }
        ),
        pivot=lambda path: (
            None if path == "workups.input.components" else watched(path)
        ),
    )
    # The inner level was refused, so the outer body could not resolve and the outer
    # level was never asked for -- which is what stops a declining query from paying
    # for a table it will not use.
    assert asked == []


def test_a_singular_struct_routes_to_its_ancestor_level_s_pivot():
    sql = _pivot_sql(
        {
            "where": {
                "op": "exists",
                "path": "outcomes.products.measurements.authentic_standard",
                "where": {"op": "not_null", "path": "smiles"},
            }
        },
        pivot=_every_level,
    )
    assert "FROM pivot_outcomes_products_measurements AS x0" in sql
    assert "x0.element.authentic_standard.smiles" in sql


# Comparing compounds rather than spellings


def test_same_compound_compiles_to_a_bitmap_test():
    compiled = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "same_compound",
                    "path": "smiles",
                    "smiles": "CC(=O)O",
                },
            }
        }
    )
    assert "get_bit" in compiled.sql
    assert [(p.op, p.pattern) for p in compiled.structures] == [
        ("same_compound", "CC(=O)O")
    ]


def test_same_compound_by_name_is_bound_rather_than_spelled():
    compiled = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "same_compound",
                    "path": "smiles",
                    "compound": "pyridine",
                },
            }
        }
    )
    assert [(p.op, p.compound) for p in compiled.structures] == [
        ("same_compound", "pyridine")
    ]


@pytest.mark.parametrize(
    "payload",
    [
        # Neither a smiles nor a compound, or both at once.
        {"op": "same_compound", "path": "smiles"},
        {
            "op": "same_compound",
            "path": "smiles",
            "smiles": "CC(=O)O",
            "compound": "pyridine",
        },
        {"op": "same_compound", "path": "smiles", "smiles": "not a molecule"},
        # An empty molecule matches nothing and still costs a pass over the corpus.
        {"op": "same_compound", "path": "smiles", "smiles": ""},
        {"op": "same_compound", "path": "smiles", "compound": "not an identifier"},
    ],
)
def test_a_malformed_same_compound_is_refused_before_compilation(payload):
    with pytest.raises(ValidationError):
        query.Query.model_validate(
            {"where": {"op": "exists", "path": "inputs.components", "where": payload}}
        )


def test_same_compound_needs_a_compound_smiles_column():
    # The hash lives beside a structure, and reaction_id has none.
    with pytest.raises(query.QueryError, match="applies to a compound"):
        _compile(
            {
                "where": {
                    "op": "same_compound",
                    "path": "reaction_id",
                    "smiles": "CC(=O)O",
                }
            }
        )


def test_same_parent_compiles_to_a_bitmap_test():
    compiled = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "same_parent",
                    "path": "smiles",
                    "smiles": "CC(=O)O",
                },
            }
        }
    )
    assert "get_bit" in compiled.sql
    assert [(p.op, p.pattern) for p in compiled.structures] == [
        ("same_parent", "CC(=O)O")
    ]


def test_the_two_compound_operators_bind_separate_parameters():
    # They ask different questions of the same molecule, so one bitmap cannot serve
    # both however alike the predicates look.
    compiled = _compile(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "and",
                    "clauses": [
                        {"op": "same_compound", "path": "smiles", "smiles": "CC(=O)O"},
                        {"op": "same_parent", "path": "smiles", "smiles": "CC(=O)O"},
                    ],
                },
            }
        }
    )
    assert [p.op for p in compiled.structures] == ["same_compound", "same_parent"]


@pytest.mark.parametrize(
    "payload",
    [
        {"op": "same_parent", "path": "smiles"},
        {
            "op": "same_parent",
            "path": "smiles",
            "smiles": "CC(=O)O",
            "compound": "pyridine",
        },
        {"op": "same_parent", "path": "smiles", "smiles": "not a molecule"},
        {"op": "same_parent", "path": "smiles", "smiles": ""},
        {"op": "same_parent", "path": "smiles", "compound": "not an identifier"},
    ],
)
def test_a_malformed_same_parent_is_refused_before_compilation(payload):
    with pytest.raises(ValidationError):
        query.Query.model_validate(
            {"where": {"op": "exists", "path": "inputs.components", "where": payload}}
        )


def test_same_parent_needs_a_compound_smiles_column():
    with pytest.raises(query.QueryError, match="applies to a compound"):
        _compile(
            {
                "where": {
                    "op": "same_parent",
                    "path": "reaction_id",
                    "smiles": "CC(=O)O",
                }
            }
        )


# Aggregating over a repeated level


def _over(body, pivot=_every_level):
    return _pivot_sql(body, pivot=pivot)


_SOLVENT_COUNTS = {
    "aggregate": {
        "over": "inputs.components",
        "where": {"op": "eq", "path": "reaction_role", "value": {"literal": "SOLVENT"}},
        "group_by": ["smiles"],
        "measures": [{"fn": "count_distinct", "path": "reaction_id", "name": "n"}],
    },
    "order_by": [{"key": "n", "descending": True}],
    "limit": 3,
}


def test_an_aggregate_over_a_level_reads_the_pivot_rather_than_the_reactions():
    # The rows are elements, so the relation is the pivot and the group_by path is
    # relative to the element -- the shape "the most common solvent" needs, which no
    # grouping of reactions can express.
    sql = _over(_SOLVENT_COUNTS)
    assert "FROM pivot_inputs_components AS x0" in sql
    assert "x0.element.smiles" in sql
    assert "count(DISTINCT x0.reaction_id) AS n" in sql
    assert "FROM reactions" not in sql


def test_the_aggregate_s_where_filters_the_elements_not_the_reactions():
    sql = _over(_SOLVENT_COUNTS)
    assert "WHERE x0.element.reaction_role = 'SOLVENT'" in sql
    assert "EXISTS" not in sql


def test_the_query_s_where_still_selects_reactions():
    # Two filters because there are two things to filter, and neither can say what the
    # other says: the elements grouped, and whose elements are counted.
    sql = _over(
        _SOLVENT_COUNTS
        | {
            "where": {
                "op": "gt",
                "path": "conditions.temperature.setpoint_kelvin",
                "value": {"literal": 350},
            }
        }
    )
    assert "x0.element.reaction_role = 'SOLVENT'" in sql
    assert (
        "EXISTS (SELECT 1 FROM reactions WHERE reactions.reaction_id = x0.reaction_id"
        in sql
    )
    assert "setpoint_kelvin > 350" in sql


def test_a_quantifier_in_the_reaction_filter_does_not_take_the_grouped_alias():
    # Both relations are the same pivot, so a shared alias would leave the inner
    # correlation binding to the rows being grouped.
    sql = _over(
        _SOLVENT_COUNTS
        | {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "eq",
                    "path": "reaction_role",
                    "value": {"literal": "CATALYST"},
                },
            }
        }
    )
    assert "AS x1 WHERE x1.reaction_id = reactions.reaction_id" in sql
    assert "AS x0" in sql


def test_grouping_elements_needs_a_pivot_for_the_level():
    # No fallback: reaching the elements without one is UNNEST in a FROM clause, which
    # is the spelling this grammar exists to make unwritable.
    with pytest.raises(query.QueryError, match="does not hold"):
        _over(_SOLVENT_COUNTS, pivot=lambda path: None)


def test_over_must_name_a_repeated_level():
    with pytest.raises(query.QueryError, match="no such level"):
        _over(
            {
                "aggregate": {
                    "over": "conditions.temperature.setpoint_kelvin",
                    "measures": [{"fn": "count", "name": "n"}],
                }
            }
        )


def test_over_names_the_level_itself_not_a_field_below_it():
    # reach() descends singular fields past a level, which a quantifier wants and this
    # cannot use: the rows would be that struct rather than the level named.
    with pytest.raises(query.QueryError, match="descends past"):
        _over(
            {
                "aggregate": {
                    "over": "outcomes.products.measurements.authentic_standard",
                    "measures": [{"fn": "count", "name": "n"}],
                }
            }
        )


def test_an_element_filter_without_a_level_is_refused():
    with pytest.raises(ValidationError, match="needs an over"):
        query.Query.model_validate(
            {
                "aggregate": {
                    "where": {"op": "not_null", "path": "reaction_id"},
                    "measures": [{"fn": "count", "name": "n"}],
                }
            }
        )


def test_a_reduction_is_refused_where_the_rows_are_already_elements():
    # A reduction folds a reaction's own list to one value; a row that is already one
    # element has no list to fold.
    with pytest.raises(query.QueryError, match="already elements"):
        _over(
            {
                "aggregate": {
                    "over": "outcomes.products",
                    "measures": [
                        {
                            "fn": "max",
                            "path": {
                                "reduce": "max",
                                "path": "outcomes.products.measurements."
                                "percentage.value",
                            },
                            "name": "n",
                        }
                    ],
                }
            }
        )


def test_the_elements_being_grouped_cannot_be_quantified_over():
    # Load-bearing, not incidental: a pivot prunes the repeated fields from its element
    # type, so an aggregate's own filter can never bind a variable. That is what lets it
    # and the reaction filter compile at the same depth without their aliases colliding.
    with pytest.raises(query.QueryError, match="no field named"):
        _over(
            {
                "aggregate": {
                    "over": "outcomes.products",
                    "where": {
                        "op": "exists",
                        "path": "measurements",
                        "where": {
                            "op": "eq",
                            "path": "type",
                            "value": {"literal": "YIELD"},
                        },
                    },
                    "measures": [{"fn": "count", "name": "n"}],
                }
            }
        )


def test_both_filters_allocate_structure_parameters_apart():
    # They append to one list, so a name taken by the element filter is not offered to
    # the reaction filter -- two bitmaps bound under one name would silently be one.
    compiled = query.compile_query(
        query.Query.model_validate(
            {
                "aggregate": {
                    "over": "inputs.components",
                    "where": {
                        "op": "substructure",
                        "path": "smiles",
                        "smarts": "c1ccccc1",
                    },
                    "group_by": ["smiles"],
                    "measures": [{"fn": "count", "name": "n"}],
                },
                "where": {
                    "op": "exists",
                    "path": "inputs.components",
                    "where": {
                        "op": "substructure",
                        "path": "smiles",
                        "smarts": "[Pd]",
                    },
                },
            }
        ),
        pivot=_every_level,
    )
    names = [parameter.name for parameter in compiled.structures]
    assert len(names) == 2
    assert len(set(names)) == 2
