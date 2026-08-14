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

"""Tests for ord_schema.agent.query."""

import warnings

import duckdb
import pytest
from pydantic import ValidationError
from rdkit import Chem

from ord_schema import projection
from ord_schema.agent import query, sql


def _compile(payload):
    return query.compile_query(query.Query.model_validate(payload))


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
