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
fingerprint **screen** — complete but not exact — over every distinct molecule in the
corpus, then exact subgraph **verification** of the survivors. It runs on RDKit's own
threads with the GIL released, so one `Corpus` serves concurrent requests without
forking; the library is built on the first substructure query and costs seconds and
about 1.5 GB for the full corpus.
Each search takes its own DuckDB cursor, since a connection holds the pending result of
its last `execute` and concurrent searches sharing one would read each other's rows.
`threads` is per search rather than per corpus, so a server expecting several at once
should set it to the core count divided by that concurrency instead of leaving it at
`-1`.

Finding which reactions hold a match is a separate cost from the chemistry, and the
larger one. An **occurrence index** — one row per structure occurrence, carrying the
corpus-wide ID, the path, and the element's own `reaction_role` — makes that a semi-join
against a narrow table rather than a scan of every reaction. Keeping the role beside the
structure is what keeps a bound query a condition on one row. At corpus scale the index
is 14.1M rows and about 130 MB, resident for the life of the `Corpus`. A search for
pyridine among the inputs returns in 1.46s **end to end**, pyridine as the solvent in
1.74s, and a boronic acid in 0.15s — and for a common pattern nearly all of that is the
RDKit screen and verify, which the index does not touch: pyridine's match alone is about
1.4s of its 1.46s. What the index cut is the reaction lookup, from roughly 3.5s to 0.2s.

The index answers one *quantifier* at a time, not one query: an `exists` whose body asks
for one structure and at most one role becomes `reaction_id IN (...)`, and the rest of
the query compiles as if the index did not exist. So "reactions with yield > 50% where
pyridine is the solvent" spends the index on its pyridine clause and answers the yield
clause from the projection, in one query — and aggregates, orderings, limits, negations,
disjunctions, and second structure predicates all compose the same way. A quantifier the
index cannot carry — one binding another element field, holding no structure predicate
or two, or any `forall`, which needs every element rather than some — compiles over the
elements, and the projection answers it. Either way it is one compiled query, screened
and verified identically; the log line says whether the index took a clause. A level the
source never recorded — most reactions have no workups and no authentic standards — is a
level with no elements: nothing satisfies an `exists` there and nothing contradicts a
`forall`, so "reactions **without** pyridine in the workup" includes every reaction that
has no workup at all, whichever way the clause was answered. The index
is built on the first query that spends it, one pass over the projections per indexed
path, so a server that wants its first real query to be fast should issue a throwaway
structure query at startup.

What remains is the chemistry itself, and two things cut it. Structures are deduplicated
per dataset, so the library holds one entry per **distinct** molecule — 1,435,426 of the
corpus's 2,016,224 rows, so 29% of the matching disappears — and maps each entry back to
every structure ID sharing it. And a match set depends on the query molecule, the
operation, and the threshold, so recent ones are **cached**: pyridine costs 1.05s the
first time and 0.02s the next. A compound is keyed by what it resolved to rather than by
its name, and by which parser reads it — the same text is one molecule as SMILES and
another as SMARTS. A predicate asked
again while the first pass is still running waits for that pass, so a burst of identical
requests costs one match; unrelated searches still overlap.

Queries the index declines read the projection, and read it faster from a **materialized
column set** holding only the top-level columns they name. Reading a handful of columns
out of a 442-leaf projection spread over 53 files costs mostly per-file overhead, which
is the same whatever the query asks for; paying it once and answering from memory
afterwards takes a temperature filter from 1.24s to 0.21s and a group-by on stirring
type from 0.75s to 0.003s. The columns are read back off the compiled SQL and the query
is then compiled again against the table, so a column it names and the table lacks is a
catalog error rather than a wrong answer, and no SQL text is edited. Sets are held to a
memory budget and evicted least-recently-used, skipping any a search is still reading;
one too large to keep is not kept, the projection answers directly, and that is
remembered rather than rediscovered by building it again. Only builds wait on builds: a
search whose columns are already materialized is answered while another is being built. The rows are the
same rows either way, in an order neither relation promises — a query wanting one has to
say so.

**The budget is the cliff, and it says so.** A refused column set logs a **warning** on
every query that names it, giving what the set costs, what the budget is, and the
argument that changes them — so a query that suddenly takes seconds explains itself in
the log rather than by profiling:

```text
WARNING materializing ['outcomes', 'reaction_id'] takes 3.3 GB, over this corpus's
2.0 GB budget, so every query naming those columns reads the Parquet files instead --
seconds rather than tenths of a second at ORD's scale. Open the Corpus with a larger
narrow_budget_bytes, or give the process more memory: the default budget is a quarter
of what it may hold.
```

Materialized, ORD's `inputs` is 4.07 GB and its `outcomes`
3.26 GB, against 1.71 GB for `conditions` and 1.46 GB for `notes`. A column set the
budget refuses is read from the 53 Parquet files on every query that names it, and
scattered rows do not let a reader skip row groups, so the same question costs seconds
instead of tenths — "pyridine as the solvent with a yield above 50%" is 3.34s where
`outcomes` is refused and 0.41s where it is held. The default is therefore a **quarter of
the machine's memory** rather than a fixed figure, and `Corpus(narrow_budget_bytes=...)`
states it outright where the machine is shared or more can be spent. This is also the
second thing the index buys: a query whose structure clause becomes a semi-join never
names `inputs` at all, so what has to be resident is far smaller than the same question
answered from the elements.

Over the whole corpus (2,428,291 reactions, 53 files, a 6 GB budget, warm), a mixed set
of queries pairing an indexed clause with one only the projection can answer:

| query | index | projection |
| --- | --- | --- |
| pyridine solvent, above 350 K | 0.02s | 0.20s |
| pyridine solvent, yield > 50% | 0.41s | 4.82s |
| pyridine solvent, white product | 0.24s | 4.30s |
| pyridine solvent, "reflux" in the procedure | 0.05s | 0.20s |
| yields by product color (grouped) | 0.40s | 1.53s |
| hottest with a yield (ordered, limited) | 0.41s | 0.99s |
| boronic acid, pyridine solvent, yield > 50% | 0.43s | 12.75s |
| any aromatic carbon, yield > 50% | 0.45s | 14.79s |
| **not** pyridine anywhere, with a yield | 0.55s | 14.86s |

The screen itself is not a lever. It is the same `PatternFingerprint` the RDKit
PostgreSQL cartridge screens with (`rdkit.sss_fp_size`, via `makeMolSignature`), and for
a small common query it admits most of the corpus — 80% for pyridine, of which 73% then
fail verification. Measured and rejected: more bits does nothing (2048 → 8192 moves
pyridine 80% → 79%), and a circular fingerprint is not usable at all, since it cannot be
computed from a SMARTS and its bits are not preserved under subgraph extraction — a
radius-2 Morgan screen dropped 4 of 5 true matches. Verification of the survivors is the
real cost, and it is intrinsic.

Similarity needs no verification — Tanimoto is defined on the Morgan fingerprint, so
the screen is the answer — and stays in SQL. The verified match set re-enters the query
as a `BITSTRING` parameter indexed by the corpus-wide ID (`structure_id` plus the file's
offset), which is how the projection path tests a single element: `exists(components,
substructure(pyridine) and role = SOLVENT)` means one component that is both. The
alternative — intersecting reaction-ID sets — over-returns by 94% on exactly that
query; see [ord-logbook#28](https://github.com/open-reaction-database/ord-logbook/pull/28),
finding 4. Every corpus-scale figure on this page was measured against that logbook
entry's snapshot of ORD.

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
