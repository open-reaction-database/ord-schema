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

"""A query a model can emit, and its compilation to DuckDB SQL.

A model emits one of these rather than SQL. The reason is not safety alone -- bound
parameters already gave that -- but **cost**. SQL's expensive programs cannot be
enumerated, so every guard on generated SQL is a blacklist someone has to keep adding
to. This grammar has one relation, no join, no recursion, and no set-returning
function, so the worst program expressible in it is one pass over the corpus and a
sort. Expensive queries are not discouraged here; they are unwritable.

It reaches every column the projection has, because a column is a *parameter* rather
than part of the grammar. That is the difference between this and the per-predicate
enumeration it replaces: adding a field upstream adds a queryable column and costs
nothing here.

Two things the compiler settles that SQL leaves to whoever writes it:

* **Quantifiers are stated, never assumed.** A path crossing a repeated level is
  refused unless an ``exists`` or ``forall`` binds it. The same intent in SQL is spelled
  ``UNNEST``, which silently means "any" *and* multiplies the row count.
* **Repeated levels compile to list lambdas**, never to ``UNNEST`` in a ``FROM`` clause
  -- measured at 27-200x the cost for identical answers, and the idiom every SQL
  tutorial teaches. A compiler cannot reach for the wrong one.

Co-membership survives, which a reaction-keyed flat table cannot express:
``exists(components, A and B)`` is one component that is both, while
``exists(components, A) and exists(components, B)`` is a reaction containing each.
They are different nodes rather than a subtlety of join keys.

Compounds are named, never spelled. A ``{"compound": "thf"}`` value compiles to a bound
parameter, so the model never writes a structure and the caller resolves the name
through :mod:`ord_schema.resolvers` before binding.

Structure predicates -- ``substructure`` and ``similarity`` -- compile to a bitmap test
rather than to chemistry. The chemistry runs *outside* the query, against the
:mod:`ord_schema.structures` artifact (:mod:`ord_schema.agent.execute` is the driver),
and the match set re-enters as a bitmap parameter indexed by the corpus-wide ID --
the projection's ``structure_id`` plus its file's offset: DuckDB permits no subquery
inside a lambda expression, so the element a quantifier binds cannot semi-join a match
table, but it can test one integer against a bitmap. The compiled SQL therefore
references ``structure_offset``, a column only the executor's relation carries -- a raw
projection cannot answer a structure query, which is the point: the IDs are
dataset-local and only the executor knows the offsets that make them one corpus-wide
space.
"""

import dataclasses
import difflib
import re
import warnings
from typing import Annotated, Any, Literal

import pyarrow as pa
from pydantic import BaseModel, Field, model_validator
from rdkit import Chem

from ord_schema import projection

# The relation a compiled query reads. The only one in scope, so nothing else is
# nameable and a join has nothing to join to.
TABLE = "reactions"

# The per-row column mapping a dataset-local structure_id into the corpus-wide ID
# space a bitmap parameter is indexed by. Supplied by the executor's relation and
# absent from the projection schema resolve() defaults to, so a model-supplied path
# does not reach it.
STRUCTURE_OFFSET = "structure_offset"


def executable_schema(schema: pa.Schema | None = None) -> pa.Schema:
    """Returns the schema of the relation a compiled query actually runs against.

    The executor's relation is the projection plus ``STRUCTURE_OFFSET``, so validating
    compiled SQL (:func:`ord_schema.agent.sql.validate`) needs this schema whenever the
    query carries a structure predicate; the projection schema alone cannot bind it.

    Args:
        schema: Base schema; the projection schema by default.

    Returns:
        The base schema with the offset column appended.
    """
    base = schema if schema is not None else projection.SCHEMA
    return base.append(pa.field(STRUCTURE_OFFSET, pa.int64()))


_COMPARISONS = {"eq": "=", "ne": "<>", "lt": "<", "le": "<=", "gt": ">", "ge": ">="}
_ORDERED = frozenset({"lt", "le", "gt", "ge"})
_TEXT = {"contains": "contains", "starts_with": "starts_with", "ends_with": "ends_with"}
_AGGREGATES = {
    "count": "count",
    "count_distinct": "count",
    "sum": "sum",
    "avg": "avg",
    "min": "min",
    "max": "max",
}

# Measure names and the bound-parameter names for compounds are the only model-supplied
# strings that reach the SQL as text rather than as parameters, so both are held to an
# identifier shape a quoting mistake cannot escape.
_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class QueryError(ValueError):
    """A query that cannot be compiled against the projection schema."""


def _explicit_hydrogens(pattern: Chem.Mol) -> int:
    """Returns how many of a query's atoms are hydrogens in their own right."""
    return sum(1 for atom in pattern.GetAtoms() if atom.GetAtomicNum() == 1)


@dataclasses.dataclass(frozen=True)
class StructureParameter:
    """A structure predicate the executor evaluates and binds as a bitmap.

    Exactly one of ``pattern`` and ``compound`` is set. ``pattern`` is a SMARTS for a
    substructure predicate and a SMILES for a similarity one, already validated;
    ``compound`` is a name still to be resolved at execution.
    """

    name: str
    op: str
    pattern: str | None
    compound: str | None
    threshold: float | None


@dataclasses.dataclass(frozen=True)
class Compiled:
    """A compiled query, and the parameters its execution still needs.

    ``compounds`` are names whose resolved SMILES the caller binds. ``structures`` are
    predicates the caller evaluates against the structures artifact, each bound as a
    bitmap over corpus-wide structure IDs; a query carrying any also references the
    ``STRUCTURE_OFFSET`` column, so it runs only against the executor's relation.
    """

    sql: str
    compounds: tuple[str, ...]
    structures: tuple[StructureParameter, ...] = ()


@dataclasses.dataclass(frozen=True)
class _Resolved:
    """Where a path landed: the expression, whether it is a list, and its type."""

    expression: str
    repeated: bool
    type: pa.DataType


def _members(current: pa.Schema | pa.DataType) -> list[str]:
    """Returns the field names available on a schema or struct type."""
    if isinstance(current, pa.Schema):
        return list(current.names)
    if pa.types.is_struct(current):
        return [field.name for field in current]
    return []


def _lookup(
    current: pa.Schema | pa.DataType, name: str, path: str, allow_internal: bool
) -> pa.Field:
    """Returns the field ``name`` within ``current``, or raises a helpful QueryError."""
    members = _members(current)
    if name not in members:
        if not members:
            raise QueryError(f"{path}: {name!r} has no fields to descend into")
        # Internal fields are excluded from the suggestions, not just refused: a
        # near-miss like 'structure' must not teach the model the one name the schema
        # description withholds.
        candidates = [
            member
            for member in members
            if not projection.is_internal(current.field(member))
        ]
        close = difflib.get_close_matches(name, candidates, n=3)
        suggestion = f"; did you mean {', '.join(map(repr, close))}?" if close else ""
        raise QueryError(f"{path}: no field named {name!r}{suggestion}")
    field = current.field(name)
    if projection.is_internal(field) and not allow_internal:
        # Artifact-internal machinery (structure_id): its values are not stable across
        # builds, so a comparison against one is a wrong answer waiting to move.
        raise QueryError(f"{path}: {name!r} is internal to the artifacts")
    return field


def resolve(
    path: str,
    *,
    schema: pa.Schema = projection.SCHEMA,
    root: str | None = None,
    allow_internal: bool = False,
) -> _Resolved:
    """Resolves a dotted path against the projection schema.

    Descending through a repeated level turns the expression into a list of the
    elements beneath it, so the caller can tell a scalar it may compare from a level it
    must quantify over.

    Args:
        path: Dotted column path, e.g. ``conditions.temperature.setpoint_kelvin``.
        schema: Schema or struct type to resolve within.
        root: Bound variable the path is relative to, inside a quantifier.
        allow_internal: Permit internal columns. Library-internal: the compiler reaches
            ``structure_id`` this way; a model-supplied path never sets it.

    Returns:
        The DuckDB expression, whether it evaluates to a list, and the type it reaches.

    Raises:
        QueryError: If the path is empty or names a field the schema does not hold.
    """
    if not path:
        raise QueryError("empty path")
    expression: str | None = None
    repeated = False
    current: Any = schema
    for part in path.split("."):
        field = _lookup(current, part, path, allow_internal)
        if expression is None:
            expression = part if root is None else f"{root}.{part}"
        elif repeated:
            expression = f"list_transform({expression}, x -> x.{part})"
        else:
            expression = f"{expression}.{part}"
        inner = field.type
        if pa.types.is_map(inner):
            expression = (
                f"flatten({expression})" if repeated else f"map_values({expression})"
            )
            repeated, current = True, inner.item_type
        elif pa.types.is_list(inner):
            if repeated:
                expression = f"flatten({expression})"
            repeated, current = True, inner.value_type
        else:
            current = inner
    return _Resolved(expression or "", repeated, current)


class Value(BaseModel):
    """A comparison operand: a literal, or a compound to resolve and bind."""

    literal: bool | int | float | str | None = None
    compound: str | None = None

    @model_validator(mode="after")
    def _exactly_one(self) -> "Value":
        if (self.literal is None) == (self.compound is None):
            raise ValueError("a value is either a literal or a compound, not both")
        if self.compound is not None and not _NAME.match(self.compound):
            raise ValueError(f"compound name is not an identifier: {self.compound!r}")
        return self


class And(BaseModel):
    """Every clause holds."""

    op: Literal["and"]
    clauses: list["Predicate"] = Field(min_length=1)


class Or(BaseModel):
    """At least one clause holds."""

    op: Literal["or"]
    clauses: list["Predicate"] = Field(min_length=1)


class Not(BaseModel):
    """The clause does not hold."""

    op: Literal["not"]
    clause: "Predicate"


class Quantifier(BaseModel):
    """Some, or every, element at a repeated level satisfies ``where``.

    Paths inside ``where`` are relative to the element bound here, which is what lets a
    single component be required to satisfy several conditions at once.
    """

    op: Literal["exists", "forall"]
    path: str
    where: "Predicate"


class Comparison(BaseModel):
    """A leaf compared against a literal or a resolved compound."""

    op: Literal[
        "eq", "ne", "lt", "le", "gt", "ge", "contains", "starts_with", "ends_with"
    ]
    path: str
    value: Value


class NullCheck(BaseModel):
    """Whether the source recorded the leaf at all."""

    op: Literal["is_null", "not_null"]
    path: str


class Substructure(BaseModel):
    """The element's structure contains the query as a subgraph.

    ``path`` names a compound's ``smiles``. The query is a SMARTS pattern, or a
    compound name resolved to a molecule at execution; exactly one is given.
    """

    op: Literal["substructure"]
    path: str
    smarts: str | None = None
    compound: str | None = None

    @model_validator(mode="after")
    def _check(self) -> "Substructure":
        if (self.smarts is None) == (self.compound is None):
            raise ValueError("a substructure query is a smarts or a compound, not both")
        if self.compound is not None and not _NAME.match(self.compound):
            raise ValueError(f"compound name is not an identifier: {self.compound!r}")
        if self.smarts is not None:
            pattern = Chem.MolFromSmarts(self.smarts)
            if pattern is None:
                raise ValueError(f"SMARTS does not parse: {self.smarts!r}")
            if not pattern.GetNumAtoms():
                # An empty SMARTS parses to a zero-atom query whose fingerprint has
                # no bits, so the screen admits every structure and verification then
                # rejects all of them: a guaranteed-empty answer at full corpus cost.
                raise ValueError(f"SMARTS has no atoms: {self.smarts!r}")
            merged = Chem.MergeQueryHs(pattern)
            if _explicit_hydrogens(merged) < _explicit_hydrogens(pattern):
                # A stored molecule is built from SMILES, so its hydrogens are
                # implicit wherever they can be, and a query naming one as its own
                # atom matches nothing -- [H]OC finds no methanol, and says nothing
                # about why. Folding them into heavy-atom H counts asks the question
                # the corpus can answer.
                rewritten = Chem.MolToSmarts(merged)
                warnings.warn(
                    f"SMARTS {self.smarts!r} names hydrogens that a stored molecule "
                    f"holds implicitly, and would match nothing; rewritten as "
                    f"{rewritten!r}, with those hydrogens folded into heavy-atom H "
                    "counts",
                    # No fixed stacklevel reaches the caller: pydantic-core invokes
                    # this validator, so anything but the default points into pydantic.
                    stacklevel=1,
                )
                self.smarts = rewritten
            # Hydrogens that do not fold -- isotopic ([2H]), or with no heavy atom to
            # fold into ([H][H]) -- are left exactly as written. They cannot be
            # implicit in a stored molecule either, so they are real graph atoms that
            # the query already matches and screens for; 28,297 ORD reactions have an
            # [H][H] component, and rewriting or refusing those would lose them.
        return self


class Similarity(BaseModel):
    """The element's structure is Tanimoto-similar to the query molecule.

    Similarity is defined on the Morgan fingerprints in the structures artifact, so
    the screen is the whole answer and no verification step exists. ``path`` names a
    compound's ``smiles``. The query is a SMILES, or a compound name resolved at
    execution; exactly one is given.
    """

    op: Literal["similarity"]
    path: str
    smiles: str | None = None
    compound: str | None = None
    threshold: float = Field(gt=0, le=1)

    @model_validator(mode="after")
    def _check(self) -> "Similarity":
        if (self.smiles is None) == (self.compound is None):
            raise ValueError("a similarity query is a smiles or a compound, not both")
        if self.compound is not None and not _NAME.match(self.compound):
            raise ValueError(f"compound name is not an identifier: {self.compound!r}")
        if self.smiles is not None:
            molecule = Chem.MolFromSmiles(self.smiles)
            if molecule is None:
                raise ValueError(f"SMILES does not parse: {self.smiles!r}")
            if not molecule.GetNumAtoms():
                # Its fingerprint has no bits, so Tanimoto against it is 0 for every
                # structure: a guaranteed-empty answer that still scans the corpus.
                raise ValueError(f"SMILES has no atoms: {self.smiles!r}")
        return self


Predicate = Annotated[
    And | Or | Not | Quantifier | Comparison | NullCheck | Substructure | Similarity,
    Field(discriminator="op"),
]

And.model_rebuild()
Or.model_rebuild()
Not.model_rebuild()
Quantifier.model_rebuild()


class Measure(BaseModel):
    """One aggregate over the matching rows."""

    fn: Literal["count", "count_distinct", "sum", "avg", "min", "max"]
    path: str | None = None
    name: str

    @model_validator(mode="after")
    def _check(self) -> "Measure":
        if not _NAME.match(self.name):
            raise ValueError(f"measure name is not an identifier: {self.name!r}")
        if self.fn != "count" and self.path is None:
            raise ValueError(f"{self.fn} needs a path")
        return self


class Aggregate(BaseModel):
    """Group the matching rows and measure each group.

    ``group_by`` paths must be scalar, so the number of groups is bounded by the values
    a column holds rather than by an explosion over a repeated level.
    """

    group_by: list[str] = Field(default_factory=list)
    measures: list[Measure] = Field(min_length=1)


class Order(BaseModel):
    """How to sort the result."""

    key: str
    descending: bool = False


class Query(BaseModel):
    """A whole query: which reactions, optionally grouped and measured."""

    where: Predicate | None = None
    aggregate: Aggregate | None = None
    order_by: list[Order] = Field(default_factory=list)
    limit: int | None = Field(default=None, gt=0)


def _literal(value: Value, compounds: list[str]) -> str:
    """Returns the SQL for a value, recording a compound that needs binding."""
    if value.compound is not None:
        if value.compound not in compounds:
            compounds.append(value.compound)
        return f"${value.compound}"
    if isinstance(value.literal, str):
        escaped = value.literal.replace("'", "''")
        return f"'{escaped}'"
    if isinstance(value.literal, bool):
        return "TRUE" if value.literal else "FALSE"
    return repr(value.literal)


def _check_operand(node: Comparison, resolved: _Resolved) -> None:
    """Raises if the operator or the value does not suit the leaf's type.

    Args:
        node: The comparison being compiled.
        resolved: Where its path landed.

    Raises:
        QueryError: If the comparison cannot be meant.
    """
    leaf = resolved.type
    if node.op in _TEXT and not pa.types.is_string(leaf):
        raise QueryError(f"{node.path}: {node.op} needs a text column, not {leaf}")
    if node.op in _ORDERED and not (
        pa.types.is_integer(leaf) or pa.types.is_floating(leaf)
    ):
        raise QueryError(f"{node.path}: {node.op} needs a numeric column, not {leaf}")
    if node.value.compound is not None and not pa.types.is_string(leaf):
        raise QueryError(f"{node.path}: a compound compares against text, not {leaf}")
    literal = node.value.literal
    if literal is None:
        return
    if isinstance(literal, str) and not pa.types.is_string(leaf):
        raise QueryError(f"{node.path}: holds {leaf}, compared against a string")
    if isinstance(literal, bool) and not pa.types.is_boolean(leaf):
        raise QueryError(f"{node.path}: holds {leaf}, compared against a boolean")
    if (
        isinstance(literal, int | float)
        and not isinstance(literal, bool)
        and not (pa.types.is_integer(leaf) or pa.types.is_floating(leaf))
    ):
        raise QueryError(f"{node.path}: holds {leaf}, compared against a number")


def _leaf(node: Any, resolved: _Resolved, compounds: list[str]) -> str:
    """Compiles a comparison or null check against an already-resolved scalar."""
    if isinstance(node, NullCheck):
        keyword = "NULL" if node.op == "is_null" else "NOT NULL"
        return f"{resolved.expression} IS {keyword}"
    _check_operand(node, resolved)
    operand = _literal(node.value, compounds)
    if node.op in _TEXT:
        return f"{_TEXT[node.op]}({resolved.expression}, {operand})"
    return f"{resolved.expression} {_COMPARISONS[node.op]} {operand}"


def _structure_parameter(
    node: "Substructure | Similarity", structures: list[StructureParameter]
) -> str:
    """Returns the parameter name for a structure predicate, reusing an equal one.

    Deduplicated by content so the executor screens and verifies each distinct
    predicate once, however many times the query states it.
    """
    pattern = node.smarts if isinstance(node, Substructure) else node.smiles
    threshold = node.threshold if isinstance(node, Similarity) else None
    for existing in structures:
        if (existing.op, existing.pattern, existing.compound, existing.threshold) == (
            node.op,
            pattern,
            node.compound,
            threshold,
        ):
            return existing.name
    parameter = StructureParameter(
        name=f"structure_{len(structures)}",
        op=node.op,
        pattern=pattern,
        compound=node.compound,
        threshold=threshold,
    )
    structures.append(parameter)
    return parameter.name


def _structure(
    node: "Substructure | Similarity",
    scope: str | None,
    schema: Any,
    structures: list[StructureParameter],
) -> str:
    """Compiles a structure predicate to a bitmap test on the element's structure ID.

    The chemistry happens in the executor; what compiles here is only the re-entry of
    its match set. The null guard keeps a compound with no recorded structure from
    reading as a match under negation-free semantics: no structure, no match.
    """
    parts = node.path.split(".")
    if parts[-1] != "smiles":
        raise QueryError(
            f"{node.path}: {node.op} applies to a compound ``smiles`` column"
        )
    resolved = resolve(node.path, schema=schema, root=scope)
    if resolved.repeated:
        raise QueryError(
            f"{node.path}: crosses a repeated level, so whether it means any or every "
            "element is unstated; wrap it in an exists or forall"
        )
    identifier_path = ".".join([*parts[:-1], "structure_id"])
    try:
        identifier = resolve(
            identifier_path, schema=schema, root=scope, allow_internal=True
        )
    except QueryError as error:
        raise QueryError(
            f"{node.path}: {node.op} needs a compound smiles, and this smiles has no "
            f"structure ID beside it ({error})"
        ) from error
    parameter = _structure_parameter(node, structures)
    return (
        f"({identifier.expression} IS NOT NULL AND get_bit(CAST(${parameter} AS "
        f"BITSTRING), ({identifier.expression} + {STRUCTURE_OFFSET})::INTEGER) = 1)"
    )


def _quantifier(
    node: "Quantifier",
    scope: str | None,
    schema: Any,
    compounds: list[str],
    structures: list[StructureParameter],
    depth: int,
) -> str:
    """Compiles a quantifier to a filter over the elements at a repeated level.

    The element is bound to a fresh variable, and the body compiles against *it* rather
    than against the row, which is what lets one component be required to satisfy
    several conditions at once.
    """
    resolved = resolve(node.path, schema=schema, root=scope)
    if not resolved.repeated:
        raise QueryError(
            f"{node.path}: {node.op} needs a repeated level, and this path reaches "
            f"{resolved.type}; compare it directly instead"
        )
    variable = f"e{depth}"
    body = _predicate(
        node.where, variable, resolved.type, compounds, structures, depth + 1
    )
    if node.op == "exists":
        return f"len(list_filter({resolved.expression}, {variable} -> {body})) > 0"
    # No counterexample, which is vacuously true of an empty level.
    return f"len(list_filter({resolved.expression}, {variable} -> NOT ({body}))) = 0"


def _predicate(
    node: Any,
    scope: str | None,
    schema: Any,
    compounds: list[str],
    structures: list[StructureParameter],
    depth: int,
) -> str:
    """Compiles one predicate node to a boolean DuckDB expression."""
    if isinstance(node, And | Or):
        keyword = " AND " if isinstance(node, And) else " OR "
        return (
            "("
            + keyword.join(
                _predicate(clause, scope, schema, compounds, structures, depth)
                for clause in node.clauses
            )
            + ")"
        )
    if isinstance(node, Not):
        inner = _predicate(node.clause, scope, schema, compounds, structures, depth)
        return f"(NOT {inner})"
    if isinstance(node, Quantifier):
        return _quantifier(node, scope, schema, compounds, structures, depth)
    if isinstance(node, Substructure | Similarity):
        return _structure(node, scope, schema, structures)
    resolved = resolve(node.path, schema=schema, root=scope)
    if resolved.repeated:
        raise QueryError(
            f"{node.path}: crosses a repeated level, so whether it means any or every "
            "element is unstated; wrap it in an exists or forall"
        )
    return _leaf(node, resolved, compounds)


def _scalar(path: str, schema: pa.Schema, what: str) -> str:
    """Returns the expression for a path that has to be scalar."""
    resolved = resolve(path, schema=schema)
    if resolved.repeated:
        raise QueryError(f"{path}: {what} needs a scalar column, not a repeated level")
    return resolved.expression


def compile_query(
    query: Query, *, schema: pa.Schema = projection.SCHEMA, table: str = TABLE
) -> Compiled:
    """Compiles a query to DuckDB SQL over the projection.

    Args:
        query: The query to compile.
        schema: Schema to resolve paths against; the projection schema by default.
        table: Relation name to read. Held to an identifier, because it reaches the SQL
            as text: a caller passing ``"reactions, range(1000000000)"`` would otherwise
            get a cross join the ``Query`` never asked for, and the single-relation cost
            bound this grammar exists to provide would hold only by convention.

    Returns:
        The SQL, the compound names whose resolved SMILES the caller binds, and the
        structure predicates the caller evaluates and binds as bitmaps.

    Raises:
        QueryError: If ``table`` is not an identifier, if any path, operator, or
            ordering key cannot be meant against the schema, or if a compound name
            collides with a generated structure parameter.
    """
    if not _NAME.match(table):
        raise QueryError(f"table is not an identifier: {table!r}")
    compounds: list[str] = []
    structures: list[StructureParameter] = []
    if query.aggregate:
        groups = [
            _scalar(path, schema, "group_by") for path in query.aggregate.group_by
        ]
        selected = list(groups)
        for measure in query.aggregate.measures:
            if measure.path is None:
                argument = "*"
            else:
                argument = _scalar(measure.path, schema, measure.fn)
                if measure.fn == "count_distinct":
                    argument = f"DISTINCT {argument}"
            selected.append(f"{_AGGREGATES[measure.fn]}({argument}) AS {measure.name}")
        names = {measure.name for measure in query.aggregate.measures}
        orderable = names | set(query.aggregate.group_by)
    else:
        selected = ["reaction_id"]
        orderable = None
    # S608: assembling SQL is this function's purpose. Every fragment in `selected` is
    # either an expression resolved against the schema or a measure name already held to
    # an identifier shape, `table` is this module's constant, and every model-supplied
    # value reached the string as a bound parameter or a quoted literal.
    sql = f"SELECT {', '.join(selected)} FROM {table}"  # noqa: S608
    if query.where is not None:
        sql += (
            f" WHERE {_predicate(query.where, None, schema, compounds, structures, 0)}"
        )
    taken = {parameter.name for parameter in structures}
    collisions = sorted(taken & set(compounds))
    if collisions:
        # Both reach the SQL as $-parameters, so a compound named like a structure
        # parameter would receive the wrong binding silently.
        raise QueryError(
            f"compound names collide with structure parameters: {collisions}"
        )
    if query.aggregate and groups:
        sql += " GROUP BY " + ", ".join(str(index + 1) for index in range(len(groups)))
    if query.order_by:
        keys = []
        for order in query.order_by:
            if orderable is None:
                key = _scalar(order.key, schema, "order_by")
            elif order.key in orderable:
                key = (
                    order.key if order.key in names else _scalar(order.key, schema, "")
                )
            else:
                raise QueryError(
                    f"{order.key}: an aggregated query orders by a measure name or a "
                    f"group_by path, and this is neither"
                )
            keys.append(f"{key} DESC" if order.descending else key)
        sql += " ORDER BY " + ", ".join(keys)
    if query.limit is not None:
        sql += f" LIMIT {query.limit}"
    return Compiled(sql, tuple(compounds), tuple(structures))
