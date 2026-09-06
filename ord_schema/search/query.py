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
than part of the grammar: adding a field upstream adds a queryable column and costs
nothing here.

Two things the compiler settles that SQL leaves to whoever writes it:

* **Quantifiers are stated, never assumed.** A path crossing a repeated level is
  refused in a predicate unless an ``exists`` or ``forall`` binds it. The same intent in
  SQL is spelled ``UNNEST``, which silently means "any" *and* multiplies the row count.
  A ``Reduction`` is the one reading that needs no quantifier, because it names what to
  do with the elements instead of asking which of them match.
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
:mod:`ord_schema.artifacts.structures` artifact (:mod:`ord_schema.search.execute` is the
driver), and the match set re-enters as a bitmap parameter indexed by the corpus-wide ID
-- the projection's ``structure_id`` plus its file's offset: DuckDB permits no subquery
inside a lambda expression, so the element a quantifier binds cannot semi-join a match
table, but it can test one integer against a bitmap. The compiled SQL therefore
references ``structure_offset``, a column only the executor's relation carries -- a raw
projection cannot answer a structure query, which is the point: the IDs are dataset-
local and only the executor knows the offsets that make them one corpus-wide space.
"""

import dataclasses
import datetime
import difflib
import re
import warnings
from collections.abc import Callable
from typing import Annotated, Any, Literal

import pyarrow as pa
from pydantic import BaseModel, ConfigDict, Field, model_validator
from rdkit import Chem

# Aliased because ``pivot`` is the parameter name a caller passes an index under, and a
# module shadowed inside the function that needs it is a bug waiting for an edit.
from ord_schema.artifacts import pivot as pivot_levels
from ord_schema.artifacts import projection

# The relation a compiled query reads. The only one in scope, so nothing else is
# nameable and a join has nothing to join to.
TABLE = "reactions"

# The per-row column mapping a dataset-local structure_id into the corpus-wide ID
# space a bitmap parameter is indexed by. Supplied by the executor's relation and
# absent from the projection schema resolve() defaults to, so a model-supplied path
# does not reach it.
STRUCTURE_OFFSET = "structure_offset"


# Consulted for an ``exists`` at the row, given the path it quantifies over, the field
# equalities its body asks for, and a thunk naming the structure parameter the match set
# will bind under. Returns a row condition standing for the whole quantifier, or None to
# leave it compiled as a filter over the elements. The thunk is called only by an index
# that takes the clause, because a parameter named and left unbound is an error rather
# than a wasted name.
ElementIndex = Callable[[str, dict[str, str], Callable[[], str]], str | None]


# Consulted for a quantifier at the row, given the path it ranges over. Returns the
# relation holding one row per element of that level, or None to leave the quantifier
# compiled as a filter over the elements. What the rows carry is ``pivot.LEVELS``, which
# this module reads from the schema rather than from the caller; the caller decides only
# whether the relation exists.
PivotIndex = Callable[[str], str | None]


@dataclasses.dataclass(frozen=True)
class _Routing:
    """Where a quantifier may be answered, besides by filtering the elements.

    Attributes:
        table: The relation the query reads. Carried because a semi-join has to
            qualify its correlation with it: the pivot holds a ``reaction_id`` of its
            own, so an unqualified outer reference binds to the *inner* one and the
            correlation compares a column to itself, which is true of every row.
        index: An occurrence index, consulted for an ``exists`` at the row.
        pivot: Pivoted levels, consulted after ``index`` and for ``forall`` too.
        enclosing: The alias and level of the pivot this clause is compiling inside,
            if it is inside one. A nested quantifier correlates to *that* rather than
            to the row: its level's ordinals extend the enclosing level's, and joining
            on the whole prefix is what keeps a measurement bound to the product it
            was reached through.
    """

    table: str
    index: ElementIndex | None = None
    pivot: PivotIndex | None = None
    enclosing: "tuple[str, pivot_levels.RepeatedLevel] | None" = None


def executable_schema(schema: pa.Schema | None = None) -> pa.Schema:
    """Returns the schema of the relation a compiled query runs against.

    The executor's relation is the projection plus ``STRUCTURE_OFFSET``, so validating
    compiled SQL (:func:`ord_schema.search.sql.validate`) needs this schema whenever the
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


def _is_temporal(leaf: pa.DataType) -> bool:
    """Returns whether ``leaf`` holds an instant or a day rather than a number."""
    return pa.types.is_timestamp(leaf) or pa.types.is_date(leaf)


def _parsed_instant(literal: str) -> datetime.datetime | None:
    """Returns the instant an ISO literal names, or None where it is not one.

    Args:
        literal: The string the query compares a temporal column against.

    Returns:
        The parsed value, with a bare date reading as that day's midnight so a query
        for a day compares against the first instant of it. Aware where the literal
        carries an offset, which the caller refuses rather than compares against a
        column holding no zone. None where the string is not ISO 8601, refused the same
        way.
    """
    try:
        return datetime.datetime.fromisoformat(literal)
    except ValueError:
        return None


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

    ``limit`` is the row bound the SQL carries, which is the query's own where it asked
    for one within ``max_rows`` and ``max_rows`` where it did not. A caller comparing it
    against ``Query.limit`` learns whether the answer was cut short.
    """

    sql: str
    compounds: tuple[str, ...]
    structures: tuple[StructureParameter, ...] = ()
    limit: int | None = None


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
    """Returns the field ``name`` within ``current``, or raises QueryError."""
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


class _Node(BaseModel):
    """Base of every model below, refusing keys the model does not declare.

    Dropping an unknown key is the wrong default here because these models are how a
    query arrives from outside -- a dict, a JSON body, a tool call. Every field that
    narrows a query is optional or lives under a discriminated union, so a misspelled
    key does not fail: it validates to a query missing the clause the caller wrote,
    which for a top-level ``where`` is every reaction in the corpus. The caller sees a
    plausible answer to a question they did not ask, with nothing logged.
    """

    model_config = ConfigDict(extra="forbid")


class Value(_Node):
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


class And(_Node):
    """Every clause holds."""

    op: Literal["and"]
    clauses: list["Predicate"] = Field(min_length=1)


class Or(_Node):
    """At least one clause holds."""

    op: Literal["or"]
    clauses: list["Predicate"] = Field(min_length=1)


class Not(_Node):
    """The clause does not hold."""

    op: Literal["not"]
    clause: "Predicate"


class Quantifier(_Node):
    """Some, or every, element at a repeated level satisfies ``where``.

    Paths inside ``where`` are relative to the element bound here, which is what lets a
    single component be required to satisfy several conditions at once.
    """

    op: Literal["exists", "forall"]
    path: str
    where: "Predicate"


class Comparison(_Node):
    """A leaf compared against a literal or a resolved compound."""

    op: Literal[
        "eq", "ne", "lt", "le", "gt", "ge", "contains", "starts_with", "ends_with"
    ]
    path: str
    value: Value


class NullCheck(_Node):
    """Whether the source recorded the leaf at all."""

    op: Literal["is_null", "not_null"]
    path: str


class Substructure(_Node):
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


class Similarity(_Node):
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


def _check_molecule_query(node: "SameCompound | SameParent") -> None:
    """Raises unless a node names exactly one query molecule RDKit can read.

    Args:
        node: The predicate to check.

    Raises:
        ValueError: If it names both a SMILES and a compound or neither, if the
            compound name is not an identifier, or if the SMILES does not parse or
            holds no atoms. An empty molecule hashes to something no real structure
            carries: a guaranteed-empty answer that still costs a pass over the corpus.
    """
    if (node.smiles is None) == (node.compound is None):
        raise ValueError("a compound query is a smiles or a compound, not both")
    if node.compound is not None and not _NAME.match(node.compound):
        raise ValueError(f"compound name is not an identifier: {node.compound!r}")
    if node.smiles is not None:
        molecule = Chem.MolFromSmiles(node.smiles)
        if molecule is None:
            raise ValueError(f"SMILES does not parse: {node.smiles!r}")
        if not molecule.GetNumAtoms():
            raise ValueError(f"SMILES has no atoms: {node.smiles!r}")


class SameCompound(_Node):
    """The element's structure is the same compound, however either was drawn.

    An ``eq`` on a ``smiles`` compares spellings, and two drawings of one reagent are
    unequal: acetic acid and acetate, an amine and its ammonium, a 2-pyridone and its
    2-hydroxypyridine tautomer. This compares the ``mol_hash`` the structures artifact
    derives, which ignores protonation state, tautomer, and atom-map labels, so those
    match. Fragments and stereochemistry are not ignored: sodium acetate stays distinct
    from acetic acid, and enantiomers from each other.

    ``path`` names a compound's ``smiles``. The query is a SMILES, or a compound name
    resolved at execution; exactly one is given.
    """

    op: Literal["same_compound"]
    path: str
    smiles: str | None = None
    compound: str | None = None

    @model_validator(mode="after")
    def _check(self) -> "SameCompound":
        _check_molecule_query(self)
        return self


class SameParent(_Node):
    """The element's structure is the same compound once counterions are set aside.

    What ``same_compound`` ignores, this ignores too, and the salt a reagent was sold
    as besides: sodium acetate matches acetic acid, and triethylamine hydrochloride
    matches triethylamine. Reagents whose parents differ still differ, so sodium
    carbonate does not match sodium acetate; reagents that *are* the salt are left
    whole, so sodium hydride does not match hydrogen and palladium acetate does not
    match acetic acid.

    It is the looser question of the two, and the looseness is the point: sodium and
    potassium carbonate share a parent and match here, which is right for "which
    reactions used carbonate" and wrong for a question about the sodium. Ask
    ``same_compound`` for that one.

    Only the counterions RDKit recognizes are set aside -- the halides, lithium,
    sodium, potassium, calcium, magnesium, and the usual acid counterions. A salt of
    something else, cesium carbonate among them, is left whole and matches only itself.

    ``path`` names a compound's ``smiles``. The query is a SMILES, or a compound name
    resolved at execution; exactly one is given.
    """

    op: Literal["same_parent"]
    path: str
    smiles: str | None = None
    compound: str | None = None

    @model_validator(mode="after")
    def _check(self) -> "SameParent":
        _check_molecule_query(self)
        return self


# The predicates the executor answers by evaluating chemistry outside the query and
# binding the match set as a bitmap over structure IDs; see ``_structure``.
StructurePredicate = Substructure | Similarity | SameCompound | SameParent

Predicate = Annotated[
    And
    | Or
    | Not
    | Quantifier
    | Comparison
    | NullCheck
    | Substructure
    | Similarity
    | SameCompound
    | SameParent,
    Field(discriminator="op"),
]

And.model_rebuild()
Or.model_rebuild()
Not.model_rebuild()
Quantifier.model_rebuild()


# DuckDB's list aggregates, which ignore the nulls a list may hold. `count` filters
# rather than taking len(), which would count them, and coalesces because the
# projection spells an absent repeated level NULL rather than empty -- 788,334
# reactions record no workups this way -- and len(NULL) is NULL, which would answer
# "how many does this reaction have" with "unknown" where the answer is none.
_REDUCERS = {
    "min": "list_min({expression})",
    "max": "list_max({expression})",
    "avg": "list_avg({expression})",
    "sum": "list_sum({expression})",
    "count": "coalesce(len(list_filter({expression}, value -> value IS NOT NULL)), 0)",
}

# The relation a pivoted reduction reads its elements from. A reduction is only ever
# compiled where the rows are reactions, and each subquery scopes its own alias, so one
# name serves every reduction in a query without shadowing the reactions it correlates
# to.
_REDUCTION_ALIAS = "r0"

# Reducers and aggregates that arithmetic has to carry. The grammar already holds
# ordered comparisons to numbers, and a sum over text is a DuckDB error rather than an
# answer, so both are refused where the query is compiled rather than where it runs.
_ARITHMETIC = frozenset({"min", "max", "avg", "sum"})


class Reduction(_Node):
    """One value per reaction, reduced from a path that crosses a repeated level.

    An ordering key and an aggregate's argument both have to be scalar, which leaves
    "the highest-yielding reactions" unwritable: a yield lives under outcomes, products,
    and measurements, so the path resolves to a list rather than a number. This reduces
    that list to the one value the reaction is judged by, leaving the aggregate to
    combine those across reactions.

    ``where`` selects which elements are reduced, and is what separates "the highest
    yield" from "the highest percentage anything was measured at": a percentage is
    recorded for selectivity and purity too, so a reduction over the bare path answers
    a broader question than the one asked. Its paths are relative to the element, as
    they are inside a quantifier, and the element is the deepest repeated level the
    path crosses -- for
    ``outcomes.products.measurements.percentage.value`` that is a measurement, so
    ``type`` there means the measurement's own.

    Attributes:
        reduce: How to reduce the list. ``count`` counts the values that are present.
        path: A dotted path crossing at least one repeated level.
        where: Selects the elements to reduce; all of them when absent.
    """

    reduce: Literal["min", "max", "avg", "sum", "count"]
    path: str
    where: "Predicate | None" = None


class Measure(_Node):
    """One aggregate over the matching rows."""

    fn: Literal["count", "count_distinct", "sum", "avg", "min", "max"]
    path: str | Reduction | None = None
    name: str

    @model_validator(mode="after")
    def _check(self) -> "Measure":
        if not _NAME.match(self.name):
            raise ValueError(f"measure name is not an identifier: {self.name!r}")
        if self.fn != "count" and self.path is None:
            raise ValueError(f"{self.fn} needs a path")
        return self


class Aggregate(_Node):
    """Group the matching rows and measure each group.

    ``group_by`` paths must be scalar, so the number of groups is bounded by the values
    a column holds rather than by an explosion over a repeated level.

    ``over`` changes what a row is. Without it the rows are reactions, and "how many
    reactions are automated" is a count of them. With it the rows are the elements of
    one repeated level, which is what "the most common solvent" asks for: solvents are
    not a column of a reaction, and no grouping of reactions can name one. Paths inside
    are then relative to the element, as they are inside a quantifier.

    Two filters, because there are two things to filter. This ``where`` selects the
    elements grouped -- the components that are solvents. The *query's* ``where``
    selects reactions, and says which reactions' elements are counted at all. "The most
    common solvents in reactions run above 350 K" needs both, and neither can say what
    the other says.

    One level, never nested, so the rows are that level's own cardinality rather than a
    product: the query stays one pass and a sort, over the pivot instead of over the
    projection. ``reaction_id`` is reachable as a measure path and is how a question
    about reactions is answered from a relation of elements -- a reaction using THF
    twice is one reaction and two rows.
    """

    over: str | None = None
    where: "Predicate | None" = None
    group_by: list[str] = Field(default_factory=list)
    measures: list[Measure] = Field(min_length=1)

    @model_validator(mode="after")
    def _element_filter_needs_elements(self) -> "Aggregate":
        if self.where is not None and self.over is None:
            raise ValueError(
                "aggregate.where selects the elements being grouped, so it needs an "
                "over; a condition on reactions is the query's own where"
            )
        return self


class Order(_Node):
    """How to sort the result."""

    key: str | Reduction
    descending: bool = False


class Query(_Node):
    """A whole query: which reactions, optionally grouped and measured."""

    where: Predicate | None = None
    aggregate: Aggregate | None = None
    order_by: list[Order] = Field(default_factory=list)
    limit: int | None = Field(default=None, gt=0)


def _literal(
    value: Value, compounds: list[str], leaf: pa.DataType | None = None
) -> str:
    """Returns the SQL for a value, recording a compound that needs binding.

    Args:
        value: The literal or compound name to compile.
        compounds: Compound names collected so far, appended to when this is one.
        leaf: The type it is compared against, which types a temporal literal. Without
            it a string compiles as text, which is what every non-temporal column wants.

    Returns:
        The SQL operand.
    """
    if value.compound is not None:
        if value.compound not in compounds:
            compounds.append(value.compound)
        return f"${value.compound}"
    if isinstance(value.literal, str):
        if leaf is not None and _is_temporal(leaf):
            # Typed rather than quoted, so the comparison is between instants. A string
            # beside a TIMESTAMP casts in DuckDB, but the same string beside a DATE
            # compares as text and answers by spelling.
            instant = _parsed_instant(value.literal)
            if instant is None:
                raise QueryError(
                    f"{value.literal!r} is not an ISO 8601 date or timestamp"
                )
            if pa.types.is_date(leaf):
                return f"DATE '{instant.date().isoformat()}'"
            return f"TIMESTAMP '{instant.isoformat(sep=' ')}'"
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
        pa.types.is_integer(leaf) or pa.types.is_floating(leaf) or _is_temporal(leaf)
    ):
        raise QueryError(
            f"{node.path}: {node.op} needs a numeric or temporal column, not {leaf}"
        )
    if node.value.compound is not None and not pa.types.is_string(leaf):
        raise QueryError(f"{node.path}: a compound compares against text, not {leaf}")
    literal = node.value.literal
    if literal is None:
        return
    if isinstance(literal, str) and _is_temporal(leaf):
        # Parsed here rather than left to _literal, so the path is in the message.
        instant = _parsed_instant(literal)
        if instant is None:
            raise QueryError(
                f"{node.path}: holds {leaf}, compared against {literal!r}, which is "
                "not an ISO 8601 date or timestamp"
            )
        if instant.tzinfo is not None:
            # The column holds no zone, because no shape the projection reads carries
            # one: a value that does is left unparsed rather than read as a wall clock.
            # DuckDB drops the offset from a TIMESTAMP literal rather than shifting by
            # it, so accepting one here would compare the caller's wall clock while
            # reading as though it compared their instant.
            raise QueryError(
                f"{node.path}: holds {leaf}, which records no time zone, and "
                f"{literal!r} carries an offset; give the instant without one"
            )
    elif isinstance(literal, str) and not pa.types.is_string(leaf):
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
    operand = _literal(node.value, compounds, resolved.type)
    if node.op in _TEXT:
        return f"{_TEXT[node.op]}({resolved.expression}, {operand})"
    return f"{resolved.expression} {_COMPARISONS[node.op]} {operand}"


def _structure_parameter(
    node: "StructurePredicate",
    structures: list[StructureParameter],
) -> str:
    """Returns the parameter name for a structure predicate, reusing an equal one.

    Deduplicated by content so the executor screens and verifies each distinct
    predicate once, however many times the query states it.

    Args:
        node: The structure predicate to name.
        structures: Parameters named so far, appended to when this one is new.

    Returns:
        The parameter name to bind the predicate's bitmap under.
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
    node: "StructurePredicate",
    scope: str | None,
    schema: Any,
    structures: list[StructureParameter],
) -> str:
    """Compiles a structure predicate to a bitmap test on the element's structure ID.

    The chemistry happens in the executor; what compiles here is only the re-entry of
    its match set. The null guard keeps a compound with no recorded structure from
    reading as a match under negation-free semantics: no structure, no match.

    Args:
        node: The structure predicate to compile.
        scope: Bound variable the path is relative to, or None at the row.
        schema: Schema or struct type the path resolves within.
        structures: Parameters named so far, appended to when this one is new.

    Returns:
        The DuckDB expression testing the element's ID against the predicate's bitmap.

    Raises:
        QueryError: If the path does not name a compound's ``smiles``, or crosses a
            repeated level without a quantifier to say which element is meant.
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


def _element_terms(
    node: "Quantifier",
) -> "tuple[StructurePredicate, dict[str, str]] | None":
    """Returns what an ``exists`` body asks of one element, or None if it asks more.

    The shape an element index can answer, stated without knowing what any index holds:
    one structure predicate on the element's own ``smiles``, and equalities on the
    element's own fields. A caller decides whether those fields are ones it carries.

    Args:
        node: The quantifier to read.

    Returns:
        ``(structure, fields)``, where ``fields`` maps a field name to the literal it
        must equal, or None if the body is anything else -- a ``forall``, a negation, a
        disjunction, a nested quantifier, a comparison that is not an equality against a
        string literal, a path that descends past the element, a field equated twice,
        or a count of structure predicates that is not exactly one.
    """
    if node.op != "exists":
        return None
    clauses = node.where.clauses if isinstance(node.where, And) else [node.where]
    structure: StructurePredicate | None = None
    fields: dict[str, str] = {}
    for clause in clauses:
        if isinstance(clause, StructurePredicate):
            if structure is not None or clause.path != "smiles":
                return None
            structure = clause
        elif (
            isinstance(clause, Comparison)
            and clause.op == "eq"
            and "." not in clause.path
            and clause.path not in fields
            and isinstance(clause.value.literal, str)
        ):
            fields[clause.path] = clause.value.literal
        else:
            return None
    if structure is None:
        return None
    return structure, fields


def _pivot_target(
    path: str, routing: _Routing
) -> "tuple[pivot_levels.RepeatedLevel, tuple[str, ...], pa.DataType] | None":
    """Returns the pivot a quantifier over ``path`` would range within, or None.

    Args:
        path: The quantifier's path, relative to the enclosing element if there is one.
        routing: Where this quantifier may be answered; its ``enclosing`` says which
            element the path is relative to.

    Returns:
        The level, the remaining segments from its element down to the value, and the
        type they reach -- or None if no pivot covers the path, or the level found is
        not a descendant of the enclosing one, which leaves no ordinal prefix to say
        which outer element an inner one belongs to.
    """
    if routing.enclosing is None:
        return pivot_levels.reach(path)
    alias_level = routing.enclosing[1]
    # A nested quantifier's path is relative to the element it sits inside, so the
    # level it names is that element's level extended by it.
    reached = pivot_levels.reach(f"{alias_level.path}.{path}")
    if reached is None:
        return None
    outer = alias_level.ordinals
    return reached if reached[0].ordinals[: len(outer)] == outer else None


def _pivoted(
    node: "Quantifier",
    compounds: list[str],
    structures: list[StructureParameter],
    depth: int,
    routing: _Routing,
) -> str | None:
    """Returns a semi-join standing for a quantifier, or None to leave it to the level.

    A pivot holds one row per element, so ``exists`` is a semi-join and ``forall`` is
    the absence of a counterexample. Both are written as ``EXISTS`` rather than ``IN``:
    a NULL ``reaction_id`` would make ``NOT IN`` answer NULL for every reaction and
    return nothing, where ``NOT EXISTS`` answers the same on a corpus that has one and
    on a corpus that does not.

    The body compiles against the element's *pruned* type, so a path reaching a field
    the pivot dropped -- anything repeated -- raises and this declines, leaving the
    quantifier to the elements. That is what keeps the covered set honest without a
    second list of what a pivot can answer.

    A level the source never recorded has no rows here, which is what the nested form
    means by folding NULL to an empty level: no element satisfies an ``exists``, and a
    ``forall`` has no counterexample.

    Args:
        node: The quantifier to compile.
        compounds: Compound names collected so far, appended to by the body.
        structures: Structure parameters collected so far, appended to by the body.
        depth: How many quantifiers enclose this one, which names its variable.
        routing: Where this quantifier may be answered; its ``pivot`` names the
            relation, and its ``table`` qualifies the correlation.

    Returns:
        The condition, or None if no pivot holds the level or the body asks for
        something its rows do not carry.
    """
    if routing.pivot is None:
        return None
    reached = _pivot_target(node.path, routing)
    if reached is None:
        return None
    level, remainder, element_type = reached
    variable = f"x{depth}"
    # Compiled into throwaway lists: a body that raises partway would otherwise leave
    # this quantifier's parameters behind for a compilation that no longer wants them.
    taken_compounds, taken_structures = list(compounds), list(structures)
    try:
        body = _predicate(
            node.where,
            ".".join([variable, pivot_levels.ELEMENT, *remainder]),
            element_type,
            taken_compounds,
            taken_structures,
            depth + 1,
            dataclasses.replace(routing, enclosing=(variable, level)),
        )
    except QueryError:
        return None
    # Asked for only once the body is known to resolve, because supplying a pivot is
    # what builds it: a level unnested over ORD is minutes, and a body reaching a
    # repeated field would otherwise pay them and then decline anyway.
    table = routing.pivot(level.path)
    if table is None:
        return None
    compounds[:], structures[:] = taken_compounds, taken_structures
    # S608: the relation is named by this module's own walk of the schema, and every
    # fragment of the body is an expression resolved against that walk's element type.
    if routing.enclosing is None:
        correlation = [
            f"{variable}.{pivot_levels.REACTION_ID} = "
            f"{routing.table}.{pivot_levels.REACTION_ID}"
        ]
    else:
        alias, outer_level = routing.enclosing
        # Joined on the whole ordinal prefix, not merely the reaction: a measurement
        # reached through a product belongs to *that* product, and correlating on
        # anything less returns reactions where the two clauses hold of different
        # elements.
        correlation = [
            f"{variable}.{pivot_levels.REACTION_ID} = "
            f"{alias}.{pivot_levels.REACTION_ID}",
            *(
                f"{variable}.{ordinal} = {alias}.{ordinal}"
                for ordinal in outer_level.ordinals
            ),
        ]
    inner = (
        f"SELECT 1 FROM {table} AS {variable} "  # noqa: S608
        f"WHERE {' AND '.join(correlation)} AND "
    )
    if node.op == "exists":
        return f"EXISTS ({inner}{body})"
    return f"NOT EXISTS ({inner}NOT ({body}))"


def _quantifier(
    node: "Quantifier",
    scope: str | None,
    schema: Any,
    compounds: list[str],
    structures: list[StructureParameter],
    depth: int,
    routing: _Routing,
) -> str:
    """Compiles a quantifier to a filter over the elements at a repeated level.

    The element is bound to a fresh variable, and the body compiles against *it* rather
    than against the row, which is what lets one component be required to satisfy
    several conditions at once.

    An ``index`` that can answer the body replaces the filter with whatever condition
    it returns, which is how a caller holding a precomputed index over elements spends
    it without the rest of the query changing. Offered only at the row, since a
    condition on the row is not a condition on an enclosing element, and only for
    ``exists``: an index shows the elements that match, never that every element does.

    Args:
        node: The quantifier to compile.
        scope: Bound variable the path is relative to, or None at the row.
        schema: Schema or struct type the path resolves within.
        compounds: Compound names collected so far, appended to by the body.
        structures: Structure parameters collected so far, appended to by the body.
        depth: How many quantifiers enclose this one, which names its variable.
        routing: Where this quantifier may be answered besides the elements. Its
            ``index`` is consulted for an ``exists`` at the row; its ``pivot`` is
            consulted after that, and for ``forall`` too, because a pivot holds every
            element rather than only the ones carrying a structure.

    Returns:
        The DuckDB expression filtering the level, true when the quantifier holds. A
        level the source never recorded is a level with no elements: the projection
        writes NULL there rather than an empty list, and ``len`` of NULL is NULL, so
        both forms are folded to what an empty level means -- no element satisfies an
        ``exists``, and a ``forall`` has no counterexample. Left as NULL, a reaction
        without workups would answer neither yes nor no to a question about them, and
        so would vanish from ``NOT exists`` -- where an indexed quantifier, which knows
        only whether the index holds a matching occurrence, says true.

    Raises:
        QueryError: If the path does not reach a repeated level.
    """
    # Offered before the path is resolved, because inside a pivot it cannot be: the
    # enclosing element's type is pruned of exactly the repeated fields a nested
    # quantifier ranges over, so resolving first would raise and decline every one.
    if routing.pivot is not None and routing.enclosing is not None:
        condition = _pivoted(node, compounds, structures, depth, routing)
        if condition is not None:
            return condition
    resolved = resolve(node.path, schema=schema, root=scope)
    if not resolved.repeated:
        raise QueryError(
            f"{node.path}: {node.op} needs a repeated level, and this path reaches "
            f"{resolved.type}; compare it directly instead"
        )
    # Offered at the row only. A nested quantifier's path is element-relative and would
    # not match an index's absolute paths anyway, so the scope check is defense in
    # depth rather than the working guard.
    if routing.index is not None and scope is None:
        terms = _element_terms(node)
        if terms is not None:
            structure, fields = terms
            # The thunk defers naming the structure parameter until the index accepts.
            # A declined clause allocates through its element compilation instead, and
            # _structure_parameter dedupes by content, so early or late it is one name.
            condition = routing.index(
                node.path, fields, lambda: _structure_parameter(structure, structures)
            )
            if condition is not None:
                return condition
    # A nested quantifier was offered the pivot above, before its path was resolved.
    # Inside a list lambda there is nothing to correlate to, so the elements answer it.
    if routing.pivot is not None and scope is None:
        condition = _pivoted(node, compounds, structures, depth, routing)
        if condition is not None:
            return condition
    variable = f"e{depth}"
    body = _predicate(
        node.where,
        variable,
        resolved.type,
        compounds,
        structures,
        depth + 1,
        routing,
    )
    if node.op == "exists":
        return (
            f"coalesce(len(list_filter({resolved.expression}, "
            f"{variable} -> {body})) > 0, false)"
        )
    # No counterexample, which is vacuously true of a level with no elements.
    return (
        f"coalesce(len(list_filter({resolved.expression}, "
        f"{variable} -> NOT ({body}))) = 0, true)"
    )


def _predicate(
    node: Any,
    scope: str | None,
    schema: Any,
    compounds: list[str],
    structures: list[StructureParameter],
    depth: int,
    routing: _Routing,
) -> str:
    """Compiles one predicate node to a boolean DuckDB expression.

    An ``index`` is carried down unchanged. A condition it returns stands for one
    quantifier and composes like any other boolean: negated, it says no element matches;
    beside another clause, it narrows the same rows the rest of the query is about.

    Args:
        node: The clause to compile: a connective, a quantifier, a structure predicate,
            or a comparison against a leaf.
        scope: Bound variable the paths within are relative to, or None at the row.
        schema: Schema or struct type the paths resolve within.
        compounds: Compound names collected so far, appended to as they are met.
        structures: Structure parameters collected so far, appended to as they are met.
        depth: How many quantifiers enclose this clause, which names their variables.
        routing: Where a quantifier may be answered besides the elements; carried
            down unchanged. See ``_Routing``.

    Returns:
        The DuckDB expression, true when the clause holds of the row or element.

    Raises:
        QueryError: If a path does not resolve, crosses a repeated level without a
            quantifier, or is compared in a way its type does not allow.
    """
    if isinstance(node, And | Or):
        keyword = " AND " if isinstance(node, And) else " OR "
        return (
            "("
            + keyword.join(
                _predicate(clause, scope, schema, compounds, structures, depth, routing)
                for clause in node.clauses
            )
            + ")"
        )
    if isinstance(node, Not):
        inner = _predicate(
            node.clause, scope, schema, compounds, structures, depth, routing
        )
        return f"(NOT {inner})"
    if isinstance(node, Quantifier):
        return _quantifier(node, scope, schema, compounds, structures, depth, routing)
    if isinstance(node, StructurePredicate):
        return _structure(node, scope, schema, structures)
    resolved = resolve(node.path, schema=schema, root=scope)
    if resolved.repeated:
        raise QueryError(
            f"{node.path}: crosses a repeated level, so whether it means any or every "
            "element is unstated; wrap it in an exists or forall"
        )
    return _leaf(node, resolved, compounds)


def _scalar(
    path: str, schema: pa.Schema, what: str, root: str | None = None
) -> _Resolved:
    """Returns where a path that has to be scalar landed."""
    resolved = resolve(path, schema=schema, root=root)
    if resolved.repeated:
        raise QueryError(f"{path}: {what} needs a scalar column, not a repeated level")
    return resolved


def _check_numeric(path: str, leaf: pa.DataType, what: str) -> None:
    """Raises unless a leaf holds numbers.

    Args:
        path: The path, named in the message.
        leaf: The type that path reaches.
        what: The reducer or aggregate asking, named in the message.

    Raises:
        QueryError: If the leaf is neither an integer nor a floating-point type.
    """
    if not (pa.types.is_integer(leaf) or pa.types.is_floating(leaf)):
        raise QueryError(f"{path}: {what} needs a numeric column, not {leaf}")


def _reduced_element(reduction: Reduction) -> tuple[Any, tuple[str, ...], Any]:
    """Returns the level a reduction's elements come from, and the descent past it.

    Args:
        reduction: What to reduce, and how.

    Returns:
        The level, the segments from its element down to the reduced value, and that
        element's type.

    Raises:
        QueryError: If no repeated level covers the path, which leaves a ``where``
            nothing to bind its paths to.
    """
    reached = pivot_levels.reach(reduction.path)
    if reached is None:
        raise QueryError(
            f"{reduction.path}: a reduction's where selects among the elements of a "
            "repeated level, and this path crosses none this grammar names"
        )
    level, remainder, _ = reached
    return level, remainder, level.element_type


def _pivoted_reduction(
    reduction: Reduction,
    table: str,
    compounds: list[str],
    structures: list["StructureParameter"],
    routing: "_Routing",
) -> str | None:
    """Returns the reduction as an aggregate over a pivot, where one covers the path.

    Reading the elements as rows costs a scan of the narrow pivot rather than a walk of
    the nested lists the projection holds, measured at 2.5s against 0.76s for a corpus
    ranking -- which is the cost of reading ``reaction_id`` alone, so the reduction
    itself becomes free.

    Args:
        reduction: What to reduce, and how.
        table: The relation of reactions the subquery correlates to.
        compounds: Compound names collected so far, appended to as they are met.
        structures: Structure parameters collected so far, appended to as they are met.
        routing: Where a filter's own clauses may be answered. See ``_Routing``.

    Returns:
        A correlated scalar subquery over the pivot, or None where no pivot covers the
        path and the projection has to answer it.
    """
    if routing.pivot is None:
        return None
    reached = pivot_levels.reach(reduction.path)
    if reached is None:
        return None
    level, remainder, _ = reached
    relation = routing.pivot(level.path)
    if relation is None:
        return None
    column = ".".join([_REDUCTION_ALIAS, pivot_levels.ELEMENT, *remainder])
    condition = (
        f"{_REDUCTION_ALIAS}.{pivot_levels.REACTION_ID} = "
        f"{table}.{pivot_levels.REACTION_ID}"
    )
    if reduction.where is not None:
        condition += " AND " + _predicate(
            reduction.where,
            f"{_REDUCTION_ALIAS}.{pivot_levels.ELEMENT}",
            level.element_type,
            compounds,
            structures,
            1,
            routing,
        )
    return (
        f"(SELECT {_AGGREGATES[reduction.reduce]}({column}) "  # noqa: S608
        f"FROM {relation} AS {_REDUCTION_ALIAS} WHERE {condition})"
    )


def _filtered_elements(
    reduction: Reduction,
    schema: pa.Schema,
    compounds: list[str],
    structures: list["StructureParameter"],
    routing: "_Routing",
) -> str:
    """Returns the projection's list of reduced values, narrowed to matching elements.

    Args:
        reduction: What to reduce, and how. Its ``where`` is compiled here.
        schema: Schema the level's own path resolves against.
        compounds: Compound names collected so far, appended to as they are met.
        structures: Structure parameters collected so far, appended to as they are met.
        routing: Where a filter's own clauses may be answered. See ``_Routing``.

    Returns:
        A DuckDB list expression holding the values the reducer combines.
    """
    level, remainder, element_type = _reduced_element(reduction)
    elements = resolve(level.path, schema=schema)
    body = _predicate(
        reduction.where,
        _REDUCTION_ALIAS,
        element_type,
        compounds,
        structures,
        1,
        routing,
    )
    expression = f"list_filter({elements.expression}, {_REDUCTION_ALIAS} -> {body})"
    for segment in remainder:
        expression = f"list_transform({expression}, x -> x.{segment})"
    return expression


def _reduced(
    reduction: Reduction,
    schema: pa.Schema,
    *,
    compounds: list[str] | None = None,
    structures: list["StructureParameter"] | None = None,
    routing: "_Routing | None" = None,
) -> str:
    """Returns the expression reducing a repeated path to one value per reaction.

    Args:
        reduction: What to reduce, and how.
        schema: Schema the path resolves against.
        compounds: Compound names collected so far, appended to as a ``where`` meets
            them; a fresh list where the caller compiles no filter.
        structures: Structure parameters, collected the same way.
        routing: Where the elements may be read besides the projection's lists, and
            the relation a pivoted reduction correlates to. See ``_Routing``.

    Returns:
        A DuckDB expression yielding one scalar per reaction. An arithmetic reducer
        yields NULL for a reaction holding no elements under that path; ``count``
        yields zero there.

    Raises:
        QueryError: If the path is already scalar, which needs no reduction (accepting
            one would give the same query two spellings), if an arithmetic reducer
            reaches a leaf that does not hold numbers, or if a ``where`` quantifies
            over a level of its own.
    """
    collected_compounds = [] if compounds is None else compounds
    collected_structures = [] if structures is None else structures
    where_to_look = _Routing(TABLE) if routing is None else routing
    resolved = resolve(reduction.path, schema=schema)
    if not resolved.repeated:
        raise QueryError(
            f"{reduction.path}: {reduction.reduce} reduces a repeated level, and this "
            f"path is already scalar; order by the path itself"
        )
    if reduction.reduce in _ARITHMETIC:
        _check_numeric(reduction.path, resolved.type, reduction.reduce)
    if reduction.where is not None:
        level, _, _ = _reduced_element(reduction)
        _refuse_quantified(reduction.where, level.path, "a reduction's where")
    pivoted = _pivoted_reduction(
        reduction,
        where_to_look.table,
        collected_compounds,
        collected_structures,
        where_to_look,
    )
    if pivoted is not None:
        return pivoted
    expression = (
        resolved.expression
        if reduction.where is None
        else _filtered_elements(
            reduction,
            schema,
            collected_compounds,
            collected_structures,
            where_to_look,
        )
    )
    return _REDUCERS[reduction.reduce].format(expression=expression)


def _measure_argument(
    measure: Measure,
    schema: pa.Schema,
    element: str | None = None,
    *,
    compounds: list[str] | None = None,
    structures: list["StructureParameter"] | None = None,
    routing: "_Routing | None" = None,
) -> str:
    """Returns the expression a measure aggregates over.

    Args:
        measure: The measure being compiled.
        schema: Schema its path resolves against.
        element: The bound element paths are relative to, for an aggregate over a
            repeated level; None where the rows are reactions. A ``Reduction`` needs
            none, since a level's elements hold no list to reduce.
        compounds: Compound names collected so far, appended to as a reduction's
            ``where`` meets them.
        structures: Structure parameters, collected the same way.
        routing: Where a reduction's elements may be read besides the projection's
            lists. See ``_Routing``.

    Returns:
        ``*`` for a bare count, a reduced expression where the path crosses a repeated
        level, and the resolved column otherwise, prefixed with ``DISTINCT`` for
        ``count_distinct``.

    Raises:
        QueryError: If the path cannot be meant as this measure's argument.
    """
    if measure.path is None:
        return "*"
    if isinstance(measure.path, Reduction):
        if element is not None:
            raise QueryError(
                f"{measure.path.path}: a reduction reads a reaction's own elements, "
                "and this aggregate's rows are already elements"
            )
        argument = _reduced(
            measure.path,
            schema,
            compounds=compounds,
            structures=structures,
            routing=routing,
        )
    elif element is not None and measure.path == pivot_levels.REACTION_ID:
        # The one path a pivot row carries outside its element, and what makes a count
        # of elements answerable as a count of reactions.
        argument = f"{element}.{pivot_levels.REACTION_ID}"
    else:
        resolved = _scalar(
            measure.path,
            schema,
            measure.fn,
            None if element is None else f"{element}.{pivot_levels.ELEMENT}",
        )
        if measure.fn in _ARITHMETIC:
            _check_numeric(measure.path, resolved.type, measure.fn)
        argument = resolved.expression
    return f"DISTINCT {argument}" if measure.fn == "count_distinct" else argument


def _refuse_quantified(
    node: "Predicate", over: str | None, what: str = "aggregate.where"
) -> None:
    """Raises where an element filter tries to bind elements of its own.

    A pivot carries an element's non-repeated fields and drops the rest, so a quantifier
    here resolves against a type missing the level it names and fails saying the field
    does not exist -- which is confusing where the field plainly does, one level up in
    the projection. Said once, here, in terms of what the rows are. Refused whichever
    route reads the elements, so that a filter compiling at all does not depend on
    whether the corpus holds the pivot.

    Args:
        node: The element filter, or a clause within it.
        over: The level whose elements are filtered, named in the message.
        what: The clause doing the filtering, named in the message.

    Raises:
        QueryError: If the filter quantifies over anything.
    """
    if isinstance(node, Quantifier):
        raise QueryError(
            f"{node.path}: {what} selects among the elements of {over}, and "
            "a pivot carries their non-repeated fields only, so there is no level "
            "here to quantify over; a condition on the reaction is the query's own "
            "where"
        )
    for clause in getattr(node, "clauses", []):
        _refuse_quantified(clause, over, what)
    inner = getattr(node, "clause", None)
    if inner is not None:
        _refuse_quantified(inner, over, what)


def _element_relation(
    query: Query, table: str, pivot: PivotIndex | None
) -> tuple[str, str, pa.DataType] | None:
    """Returns the relation an aggregate's rows come from, where they are not reactions.

    Args:
        query: The query being compiled.
        table: The relation of reactions, named in the message where the level is not
            one this grammar can group over.
        pivot: Pivoted levels; the only relation of elements there is.

    Returns:
        ``(alias, relation, element type)`` for an ``aggregate.over``, or None where the
        rows are reactions.

    Raises:
        QueryError: If the path names no repeated level, or no pivot covers it. There
            is no fallback: reaching the elements without one means ``UNNEST`` in a
            ``FROM`` clause, measured at 27-200x for the same answer and the one
            spelling this grammar exists to make unwritable.
    """
    if query.aggregate is None or query.aggregate.over is None:
        return None
    over = query.aggregate.over
    reached = pivot_levels.reach(over)
    if reached is None:
        raise QueryError(
            f"{over}: aggregate.over names a repeated level to group the elements of, "
            f"and this reaches no such level of {table}"
        )
    level, remainder, element_type = reached
    if remainder:
        # reach() descends singular fields past a level, which a quantifier wants and
        # this cannot use: the rows would be that struct, and group_by paths relative to
        # it would silently mean something narrower than the level named.
        raise QueryError(
            f"{over}: aggregate.over names the level itself, and this descends past "
            f"{level.path} into {'.'.join(remainder)}"
        )
    relation = None if pivot is None else pivot(level.path)
    if relation is None:
        raise QueryError(
            f"{over}: grouping elements reads the pivot artifact for {level.path}, "
            "which this corpus does not hold"
        )
    return "x0", relation, element_type


def compile_query(
    query: Query,
    *,
    schema: pa.Schema = projection.SCHEMA,
    table: str = TABLE,
    index: ElementIndex | None = None,
    pivot: PivotIndex | None = None,
    max_rows: int | None = None,
) -> Compiled:
    """Compiles a query to DuckDB SQL over the projection.

    Args:
        query: The query to compile.
        schema: Schema to resolve paths against; the projection schema by default.
        table: Relation name to read. Held to an identifier, because it reaches the SQL
            as text: a caller passing ``"reactions, range(1000000000)"`` would otherwise
            get a cross join the ``Query`` never asked for, and the single-relation cost
            bound this grammar exists to provide would hold only by convention.
        index: An element index to spend where it can answer a quantifier; see
            ``ElementIndex``. Without one, every quantifier compiles to a filter over
            the elements, which is what the projection alone can answer.
        pivot: Pivoted levels to spend where they can answer a quantifier; see
            ``PivotIndex``. Consulted after ``index``, and for ``forall`` as well as
            ``exists``.
        max_rows: Most rows the SQL may return, applied whether or not the query asked
            for a limit and clamping one that asks for more. The grammar leaves
            ``limit`` optional, and a query without one returns every matching
            reaction -- millions of rows at corpus scale, materialized in whatever
            process is executing. A caller serving queries it did not write wants this
            set; None leaves the query's own limit, or none, exactly as stated. It
            bounds an aggregated query's *groups*, where a truncated answer is an
            arbitrary part of a distribution rather than a sample of the matching
            reactions; ``Corpus.search`` says so when one comes back full.

    Returns:
        The SQL, the compound names whose resolved SMILES the caller binds, and the
        structure predicates the caller evaluates and binds as bitmaps.

    Raises:
        QueryError: If ``table`` is not an identifier, if any path, operator, or
            ordering key cannot be meant against the schema, or if a compound name
            collides with a generated structure parameter.
        ValueError: If ``max_rows`` is less than one, which is a bound no query can
            satisfy rather than a small one.
    """
    if not _NAME.match(table):
        raise QueryError(f"table is not an identifier: {table!r}")
    if max_rows is not None and max_rows < 1:
        # Zero compiles to a LIMIT no query can satisfy and a negative one to SQL
        # DuckDB refuses at execution, both of them far from whoever passed it.
        raise ValueError(
            f"max_rows is {max_rows}, which no query can return; leave it unset to "
            "bound nothing"
        )
    compounds: list[str] = []
    structures: list[StructureParameter] = []
    element = _element_relation(query, table, pivot)
    if query.aggregate:
        alias = None if element is None else element[0]
        element_schema = schema if element is None else element[2]
        groups = [
            _scalar(
                path,
                element_schema,
                "group_by",
                None if alias is None else f"{alias}.{pivot_levels.ELEMENT}",
            ).expression
            for path in query.aggregate.group_by
        ]
        selected = list(groups)
        for measure in query.aggregate.measures:
            argument = _measure_argument(
                measure,
                element_schema,
                alias,
                compounds=compounds,
                structures=structures,
                routing=_Routing(table, index, pivot),
            )
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
    source = table if element is None else f"{element[1]} AS {element[0]}"
    sql = f"SELECT {', '.join(selected)} FROM {source}"  # noqa: S608
    conditions = []
    if element is not None and query.aggregate is not None:
        if query.aggregate.where is not None:
            _refuse_quantified(query.aggregate.where, query.aggregate.over)
            conditions.append(
                _predicate(
                    query.aggregate.where,
                    f"{element[0]}.{pivot_levels.ELEMENT}",
                    element[2],
                    compounds,
                    structures,
                    1,
                    _Routing(table, index, pivot),
                )
            )
    if query.where is not None:
        condition = _predicate(
            query.where,
            None,
            schema,
            compounds,
            structures,
            # Past the alias this aggregate's own relation took, so a quantifier inside
            # does not bind its variable to the elements being grouped.
            0 if element is None else 1,
            _Routing(table, index, pivot),
        )
        if element is not None:
            # ``where`` selects reactions whatever the rows are, so over a relation of
            # elements it is a semi-join back to them rather than a filter here: the
            # predicate is written against reaction columns, which a pivot row does not
            # carry.
            condition = (
                f"EXISTS (SELECT 1 FROM {table} WHERE "  # noqa: S608
                f"{table}.{pivot_levels.REACTION_ID} = "
                f"{element[0]}.{pivot_levels.REACTION_ID} AND {condition})"
            )
        conditions.append(condition)
    if conditions:
        sql += " WHERE " + " AND ".join(conditions)
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
            if isinstance(order.key, Reduction):
                if orderable is not None:
                    # After grouping there is no reaction left to reduce over: the
                    # reduction is one input to a measure, not a key beside it.
                    raise QueryError(
                        "an aggregated query orders by a measure name or a group_by "
                        "path; reduce inside a measure instead"
                    )
                key = _reduced(
                    order.key,
                    schema,
                    compounds=compounds,
                    structures=structures,
                    routing=_Routing(table, index, pivot),
                )
            elif orderable is None:
                key = _scalar(order.key, schema, "order_by").expression
            elif order.key in orderable:
                key = (
                    order.key
                    if order.key in names
                    else _scalar(order.key, schema, "").expression
                )
            else:
                raise QueryError(
                    f"{order.key}: an aggregated query orders by a measure name or a "
                    f"group_by path, and this is neither"
                )
            keys.append(f"{key} DESC" if order.descending else key)
        sql += " ORDER BY " + ", ".join(keys)
    if max_rows is None:
        limit = query.limit
    elif query.limit is None:
        limit = max_rows
    else:
        limit = min(query.limit, max_rows)
    if limit is not None:
        sql += f" LIMIT {limit}"
    return Compiled(sql, tuple(compounds), tuple(structures), limit)
