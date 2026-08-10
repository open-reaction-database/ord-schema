# Reaching ORD from an agent

The serialized protocol buffers in [ord-data](https://github.com/open-reaction-database/ord-data)
are the [source of truth](https://en.wikipedia.org/wiki/Single_source_of_truth), and the
[projection](../projection.py) restates them as nested Parquet so a query engine can read one
leaf without deserializing the rest. This package is how a language model reaches that
projection: what the model is told it may query, what it is allowed to emit, and how that
becomes SQL.

Requires the `agent` extra (`pip install "ord-schema[agent]"`).

## Why the model does not write SQL

The obvious design is to let the model write DuckDB SQL and check it. That was tried, and the
check cannot be made to hold: **the expensive programs expressible in SQL cannot be
enumerated.** A validator that refuses `UNNEST` in a `FROM` clause still accepts a self-join
over 2.4M reactions, a cross join against `range`, a recursive CTE counting to a billion, and a
predicate-free `SELECT *`. Entries can be added to that blacklist, but it never finishes.

So the model emits a `Query` instead, and the library compiles it. The grammar has one
relation, no join, no recursion, and no set-returning function, which means the worst program
expressible in it is **one pass over the corpus plus a sort**. Expensive queries are not
discouraged; they are unwritable.

The same move settles two other things SQL leaves to whoever writes it:

- **Quantifiers are stated, never assumed.** A path crossing a repeated level is refused unless
  an `exists` or `forall` binds it. The equivalent SQL is `UNNEST`, which silently means "any"
  *and* multiplies the row count.
- **Repeated levels compile to list lambdas**, never to `UNNEST` in a `FROM` clause — measured
  at 27–200× the cost for identical answers, and the idiom every SQL tutorial teaches. A
  compiler cannot reach for the wrong one.

It gives up nothing in reach. A column is a *parameter* rather than part of the grammar, so all
442 leaves of the projection are queryable, and a field added upstream becomes queryable at no
cost here.

## The grammar

```text
Query      = { where?: Predicate, aggregate?: Aggregate,
               order_by?: [Order], limit?: int }

Predicate  = { op: "and" | "or", clauses: [Predicate] }
           | { op: "not", clause: Predicate }
           | { op: "exists" | "forall", path: Path, where: Predicate }
           | { op: "eq"|"ne"|"lt"|"le"|"gt"|"ge", path: Path, value: Value }
           | { op: "contains"|"starts_with"|"ends_with", path: Path, value: Value }
           | { op: "is_null" | "not_null", path: Path }
           | { op: "substructure", path: Path, smarts?: string, compound?: <name> }
           | { op: "similarity", path: Path, smiles?: string, compound?: <name>,
               threshold: float }

Value      = { literal: <scalar> } | { compound: <name> }

Aggregate  = { group_by: [Path],
               measures: [{ fn: "count"|"count_distinct"|"sum"|"avg"|"min"|"max",
                            path?: Path, name: string }] }

Order      = { key: string, descending?: bool }
```

A `Path` is a dotted column path such as `conditions.temperature.setpoint_kelvin`. Inside an
`exists` or `forall` it is **relative to the bound element**, which is what lets one component
be required to satisfy several conditions at once.

Every rule below is enforced against `projection.SCHEMA` before any SQL exists, so a bad query
is a compile error rather than a wrong answer:

- A comparison path must end at a scalar and must not cross a repeated level.
- An `exists`/`forall` path must end *at* a repeated level.
- A path naming an artifact-internal column (`structure_id`, the join key into the
  structures artifact) is refused, and internal columns are withheld from the rendered
  schema and from "did you mean" suggestions: their values are not stable across
  builds.
- Operators must suit the leaf type — `contains` is text-only, ordering is numeric-only.
- `group_by` paths must be scalar, so the number of groups is bounded by the values a column
  holds rather than by an explosion over a repeated level.
- A `{"compound": ...}` value is resolved through [`ord_schema.resolvers`](../resolvers.py) and
  **bound as a parameter**, so the model names compounds and never spells structures.
- A `substructure`/`similarity` path must name a compound's `smiles`, inside a
  quantifier like any other element predicate.
- A SMARTS naming a hydrogen the corpus stores implicitly is rewritten, with a warning,
  rather than run as written: stored molecules come from SMILES, so `[H]OC` matches no
  methanol and would return empty without saying why. `MergeQueryHs` folds it to
  `[O&!H0]C`, which matches. Hydrogens that cannot be folded — isotopic (`[2H]`), or
  with no heavy atom to fold into (`[H][H]`) — are left exactly as written, because
  those are real graph atoms in a stored molecule and the query already works;
  28,297 ORD reactions have an `[H][H]` component.

## Structure search

A structure predicate compiles to a bitmap test, not to chemistry. The chemistry runs
in [`execute.Corpus`](execute.py) against the [structures artifact](../structures.py).
Substructure runs through an RDKit `SubstructLibrary` built over the corpus: a
fingerprint **screen** — complete but not exact — over every structure row
(deduplicated per dataset, so a molecule in two datasets is screened twice), then exact
subgraph **verification** of the survivors. It runs on RDKit's own threads with the GIL
released, so one `Corpus` serves concurrent requests without forking; the library is
built on the first substructure query and costs about 8s and 2 GB for the full corpus.
Similarity needs no verification — Tanimoto is defined on the Morgan fingerprint, so
the screen is the answer — and stays in SQL. The verified match set re-enters the query
as a `BITSTRING` parameter indexed by the corpus-wide ID (`structure_id` plus the
file's offset), which is what preserves element binding: `exists(components,
substructure(pyridine) and role = SOLVENT)` means one component that is both. The
alternative — intersecting reaction-ID sets — over-returns by 94% on exactly that
query; see [ord-logbook#28](https://github.com/open-reaction-database/ord-logbook/pull/28),
finding 4.

```python
from ord_schema.agent import execute, query

corpus = execute.Corpus("projections/*/*.parquet", "structures/*/*.parquet")
table = corpus.search(
    query.Query.model_validate(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "and",
                    "clauses": [
                        {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
                        {"op": "eq", "path": "reaction_role",
                         "value": {"literal": "SOLVENT"}},
                    ],
                },
            },
            "limit": 100,
        }
    ),
    timeout_seconds=60,
)
```

`Corpus` refuses artifacts that do not pair: every projection must have a structures
artifact derived from the same source dataset, and long enough to hold every
`structure_id` the projection carries, because the IDs joining them are meaningful only
as a pair. Pairing is by source dataset rather than by filename, so the two trees need
not share a layout. Structure IDs are dataset-local; the executor's relation
carries a per-file offset column (`structure_offset`), which is why compiled SQL with a
structure predicate runs only there — validate it against `query.executable_schema()`
rather than the bare projection schema.

## Usage

### Tell the model what it may query

```python
from ord_schema.agent import schema

print(schema.describe())
```

That renders the projection as an indented type tree in DuckDB's type names — 442 leaves in 537
lines, small enough to sit in a system prompt whole, which is what lets translation stay a
single tool call rather than a retrieval loop over column metadata. Units ride along in the
column names (`setpoint_kelvin`, `mass_grams`), so nothing has to explain them, and each enum
column carries the values it may hold:

```text
      reaction_role: VARCHAR  (UNSPECIFIED | REACTANT | REAGENT | SOLVENT | CATALYST | ...)
```

An enum projects as its value *name*, so without this a model would have the column and no way
to learn the spelling to compare against. The members ride in Arrow field metadata on the
projection itself and survive the Parquet round trip at any nesting depth, so a published file
says what its own enum columns may hold — no type can do that job, since a dictionary column
carries only the values a batch happened to contain and DuckDB reads it back as `VARCHAR`
regardless.

`projection_schema.txt` alongside this file is a checked-in snapshot of that output. It is a
review artifact, never a source: `describe()` generates what actually ships, and a test keeps
the two identical so a proto field added upstream shows up in a diff rather than in a model's
behavior.

### Compile a query

```python
from ord_schema.agent import query

compiled = query.compile_query(
    query.Query.model_validate(
        {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {"op": "eq", "path": "smiles", "value": {"compound": "thf"}},
            },
            "limit": 100,
        }
    )
)
print(compiled.sql)
print(compiled.compounds)  # ('thf',) -- resolve and bind these
```

```sql
SELECT reaction_id FROM reactions
WHERE len(list_filter(flatten(list_transform(map_values(inputs), x -> x.components)),
                      e0 -> e0.smiles = $thf)) > 0
LIMIT 100
```

`compiled.compounds` names the parameters still to be bound. Resolve each through
`ord_schema.resolvers` and pass the SMILES at execution — the name never reaches the SQL as
text.

### Run it

```python
import duckdb

from ord_schema.resolvers import resolve_name

connection = duckdb.connect()
connection.execute(
    "CREATE VIEW reactions AS SELECT * FROM read_parquet('projections/*/*.parquet')"
)
parameters = {name: resolve_name("name", name)[0] for name in compiled.compounds}
reaction_ids = [row[0] for row in connection.execute(compiled.sql, parameters).fetchall()]
```

### Ask for the same component, not merely the same reaction

The distinction a reaction-keyed flat table cannot express, and the reason quantifiers bind an
element rather than a row:

```python
# One component that is BOTH tetrahydrofuran AND used in quantity.
{"op": "exists", "path": "inputs.components", "where": {"op": "and", "clauses": [
    {"op": "eq", "path": "smiles", "value": {"compound": "thf"}},
    {"op": "gt", "path": "amount.volume_liters", "value": {"literal": 0.005}}]}}

# A reaction containing each, possibly as different components.
{"op": "and", "clauses": [
    {"op": "exists", "path": "inputs.components",
     "where": {"op": "eq", "path": "smiles", "value": {"compound": "thf"}}},
    {"op": "exists", "path": "inputs.components",
     "where": {"op": "gt", "path": "amount.volume_liters", "value": {"literal": 0.005}}}]}
```

### Group and measure

```python
{
    "aggregate": {
        "group_by": ["conditions.stirring.type"],
        "measures": [
            {"fn": "avg", "path": "conditions.temperature.setpoint_kelvin", "name": "mean_k"},
            {"fn": "count", "name": "n"},
        ],
    },
    "order_by": [{"key": "mean_k", "descending": True}],
    "limit": 10,
}
```

An aggregated query orders by a measure name or a `group_by` path, and by nothing else.

## What it deliberately cannot express

Arbitrary expressions (`yield / temperature`), window functions, and joins against anything
else. That is the line drawn in the design log: lookup and search stay mediated, analysis goes
direct. An analyst wanting a window function opens the Parquet in DuckDB, where the user is a
human who is not the injection risk and whose slow query costs only their own session.

## The backstop

[`sql.validate`](sql.py) plans a query against an *empty* Arrow table carrying the projection
schema, so it needs no corpus. It refuses anything that is not a single `SELECT`, refuses
`UNNEST` in a `FROM` clause, and runs on a connection with no filesystem access.

It is a check on the compiler's own output, not the guard it once was — the guard is the
grammar. It still **does not bound cost**, and says so: a caller running arbitrary SQL against
a real corpus needs its own statement timeout, row cap, or memory limit.

## Not yet solved

Execution has no sandbox. `enable_external_access=false` cannot be combined with a lazy Parquet
view — only a materialized table survives it — so running against the full corpus needs
DuckDB's `allowed_directories` or a separate process. Compiling from an IR removes the reasons
to fear the *query*; it does not remove the reasons to contain the *process*.
