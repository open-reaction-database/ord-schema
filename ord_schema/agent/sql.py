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

"""Checking the SQL a model wrote before anything runs it.

Validation needs no corpus. ``projection.SCHEMA`` is generated from the proto
descriptors, so an empty Arrow table carrying it has the real 442-leaf shape, and
planning a query against that resolves every column reference and type without reading a
byte of data. The prompt in :mod:`ord_schema.agent.schema` renders the same schema
object, so the columns a model is told about and the columns its query is checked
against come from one call.

Two classes of defect are caught, and they are different in kind:

* **Not a query.** Anything but a single ``SELECT`` is refused, which is a whitelist
  rather than a search for known-bad text. The connection also has no filesystem access,
  so ``read_parquet`` and ``COPY`` fail to plan -- but that is a property of *this*
  connection, not a certificate about the SQL: a reader that binds its path lazily, such
  as ``sniff_csv``, plans cleanly here and is refused only when something executes it.
* **Not answerable.** A column that does not exist, a type that will not compare, a
  placeholder the model referenced but did not declare. These are ordinary binder
  errors, surfaced early enough to retry the translation instead of failing a user.

Placeholders are bound, never interpolated. A compound name reaches the query as a
parameter value, so it cannot carry SQL with it, and the same values used here are the
ones the query runs with.

.. warning::
   **Nothing here bounds what a query costs.** ``UNNEST`` in a ``FROM`` clause is
   refused, because it materializes the exploded intermediate -- measured at 27-200x the
   cost of the equivalent list lambdas, and over the full corpus the difference between
   0.9 seconds and no result in four minutes -- and because it is the spelling every SQL
   tutorial teaches, so a model reaches for it by default.

   That is one entry in a blacklist, not a guard. A self-join over 2.4M reactions, a
   cross join against ``range``, a recursive CTE counting to a billion, and a
   predicate-free ``SELECT *`` all plan cleanly here, and none of them is reachable by
   enumerating spellings: accepting SQL means accepting that the set of expensive
   programs cannot be listed. A caller running validated SQL against a real corpus needs
   its own bound -- a statement timeout, a row cap, a memory limit -- and must not read
   this module as having provided one.
"""

import json
from collections.abc import Iterator
from typing import Any

import duckdb
import pyarrow as pa

from ord_schema import projection

# The name the model is told to query. The only relation in scope, so no other table is
# reachable -- which is not the same as no join being expressible, since a query can
# still join this one to itself.
TABLE = "reactions"

# DuckDB compiles UNNEST in a FROM clause to a table in-out function; list lambdas
# compile to an ordinary projection over the child arrays. The operator name is
# therefore the signal, and it survives the query being spelled in ways a text search
# would miss.
_EXPLODING_OPERATOR = "INOUT_FUNCTION"


def _operators(node: Any) -> Iterator[str]:
    """Yields the operator name of every node in a DuckDB JSON plan.

    Only the node's own name is read, and only ``children`` is descended into, so a
    string literal that happens to spell an operator name cannot be mistaken for one: a
    query filtering on ``'INOUT_FUNCTION'`` carries that text in ``extra_info``, which
    is never visited. Both spellings of the key are accepted because DuckDB has used
    each.

    Args:
        node: A plan node, or the whole plan.

    Yields:
        Each operator name, outermost first.
    """
    if isinstance(node, list):
        for item in node:
            yield from _operators(item)
        return
    if not isinstance(node, dict):
        return
    name = node.get("name") or node.get("operator_type")
    if isinstance(name, str):
        yield name
    yield from _operators(node.get("children", []))


class InvalidQueryError(ValueError):
    """The generated SQL is not a query this surface will run."""


def _connect(schema: pa.Schema) -> duckdb.DuckDBPyConnection:
    """Returns a connection holding an empty ``TABLE`` and no access to anything else.

    External access is disabled after the table is registered, which leaves the
    in-memory relation queryable while ``read_parquet``, ``COPY``, and ``ATTACH`` can
    reach no path. The setting cannot be turned back on, so a query cannot lift it.

    Args:
        schema: Schema the registered table carries.

    Returns:
        A connection scoped to a single empty relation.
    """
    connection = duckdb.connect()
    connection.register(TABLE, schema.empty_table())
    connection.execute("SET enable_external_access=false")
    return connection


def validate(
    sql: str,
    *,
    parameters: dict[str, str] | None = None,
    schema: pa.Schema = projection.SCHEMA,
) -> None:
    """Checks that ``sql`` is a single read-only query this surface will run.

    Says nothing about what the query costs; see the module warning.

    Args:
        sql: The generated SQL, querying ``TABLE``.
        parameters: Values for the query's named placeholders, exactly as they will be
            bound at execution. Planning needs them, so a query referencing a
            placeholder the caller did not supply fails here rather than at execution.
            Binding different values later invalidates the result, since constant
            folding can prune a branch at plan time; nothing detects that, so validating
            once and re-binding is the caller's mistake to avoid.
        schema: Schema to plan against; the projection schema by default.

    Raises:
        InvalidQueryError: If the SQL does not parse, is more than one statement, is not
            a ``SELECT``, references something the schema does not hold, or plans to an
            exploded intermediate.
    """
    try:
        statements = duckdb.extract_statements(sql)
    except duckdb.ParserException as error:
        raise InvalidQueryError(f"query does not parse: {error}") from error
    if len(statements) != 1:
        raise InvalidQueryError(
            f"expected a single statement, got {len(statements)}; "
            "a query answers one question"
        )
    if statements[0].type != duckdb.StatementType.SELECT:
        raise InvalidQueryError(
            f"expected a SELECT, got {statements[0].type.name}; this surface reads"
        )
    connection = _connect(schema)
    try:
        plan = connection.execute(
            f"EXPLAIN (FORMAT json) {sql}", parameters or {}
        ).fetchone()
    except duckdb.Error as error:
        raise InvalidQueryError(f"query cannot be planned: {error}") from error
    finally:
        connection.close()
    # JSON rather than the default rendering, which draws the plan into a fixed-width
    # canvas and silently drops whole subtrees once a query has more parallel pipelines
    # than fit. A UNNEST padded with enough UNION ALL branches disappeared from it.
    try:
        operators = set(_operators(json.loads(plan[1]) if plan else None))
    except (TypeError, IndexError, ValueError) as error:
        # Fail closed. This check is the only affordability signal the module has, so a
        # plan it cannot read is a query it cannot vouch for.
        raise InvalidQueryError(f"query produced no readable plan: {error}") from error
    if _EXPLODING_OPERATOR in operators:
        raise InvalidQueryError(
            "query explodes the nested columns with UNNEST in a FROM clause, which "
            "materializes one row per element; rewrite it with list_filter / "
            "list_transform / flatten lambdas, which read the same values in place"
        )
