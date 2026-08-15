# Pivoted element index

- **Date:** 2026-08-15
- **Status:** approved, phase 1 and phase 2 in implementation

## Problem

The agent executor answers a quantifier by filtering the elements of a repeated level
with a `list_filter` lambda over the projection. To keep that from re-reading Parquet on
every query it materializes whole top-level columns as in-memory narrow tables, which
costs 18.46 GiB resident for the projection and sets how large a serving container must
be.

That cost is decoded *nesting*, not the information a predicate reads. Measured over the
full corpus (2,428,291 reactions), pivoting the repeated levels to one row per element
gives the same query surface in 2.68 GiB resident and 381 MiB as Parquet, answering the
benchmark queries in 6–48 ms where the nested form needs 1.3–3.4 s, with identical
answers. See the logbook entry
[Where the agent search cache can live](https://github.com/open-reaction-database/ord-logbook/blob/main/entries/2026-08-15-where-the-search-cache-lives.md).

## What a pivot is

For a repeated level at dotted path `P`, the pivot is a table holding one row per
element of `P`, carrying:

1. `reaction_id`.
2. The **ordinal of every repeated level from the root down to and including `P`**.
3. The element's own fields, with every **repeated** field — list or map — removed
   recursively. Struct nesting is kept.

It is **unfiltered**: every element gets a row, including elements whose fields are all
NULL.

### Why the full ordinal prefix

A spike compared the nested and pivoted spellings of "a desired product with a yield
above 50%" over the whole corpus, as reaction-ID sets:

| spelling | reactions | false positives |
| --- | --- | --- |
| nested (the answer) | 22,666 | — |
| prefix join on `(reaction_id, outcome_index, product_index)` | 22,666 | 0 |
| outcome-only join | 23,608 | 942 |
| reaction-only join | 23,608 | 942 |

Correlating on anything short of the full prefix over-returns by 4.2%. The corpus is
effectively single-outcome, so dropping the *outcome* ordinal changes nothing — a test
that checked only outcome-level correlation would pass while the product-level defect
shipped. The ordinals are what the nested form gets structurally by binding an element,
and they must be reconstructed explicitly.

### Why struct nesting is kept

The pivot's win comes from removing repeated levels, not structs. Building a table that
pulls `conditions.temperature.setpoint_kelvin` out of a deep struct across 2,428,291
rows took 2 seconds; DuckDB stores a `STRUCT` as separate child vectors, so reading one
field does not touch its siblings. Keeping the nesting means paths resolve identically
on both routes — `percentage.value` is `x.percentage.value` against a pivot and
`e0.percentage.value` against the projection — which needs no flattening, no name
collisions, and no change to `resolve`.

It also makes declining structural. A body reaching `measurements` finds no such field
on the pruned element type, so `resolve` refuses and the quantifier falls back to the
projection. The phase 1 boundary enforces itself rather than being a flag to keep
honest.

### Why completeness matters

`_quantifier`'s docstring says an index is offered "only for `exists`: an index shows
the elements that match, never that every element does." That holds for the occurrence
index because it filters to `structure_id IS NOT NULL`. A pivot is unfiltered and
therefore complete for its path, which is what makes `forall` indexable. Completeness,
not shape, is the deciding property.

## Phase 1 — the query path

### New module: `ord_schema/agent/pivot.py`

Pure schema functions, importing `pyarrow` and `ord_schema.projection` only, so that
`query.py` may import it without a cycle:

- `repeated_levels(schema) -> dict[str, RepeatedLevel]` — walks the schema for every
  list or map level, in the spirit of `structures.structure_levels`: derived from one
  walk rather than a hand-kept list, so a repeated field added upstream is covered
  without anyone updating a list. Each entry carries the dotted path, the ordinals its
  rows need, and the pruned element type.
- `element_type(path) -> pa.DataType` — the element's struct type with repeated fields
  removed recursively.
- `table_name(path) -> str` — a deterministic identifier for the pivot of `path`.

Nothing here reads the corpus or emits SQL.

### `ord_schema/agent/query.py`

- `PivotIndex = Callable[[str], str | None]` — a path to the pivot table holding it, or
  None. The compiler derives the pruned element type from the schema itself, so this is
  the only thing it needs from the executor.
- `compile_query(..., pivot: PivotIndex | None = None)`, threaded to `_predicate` and
  `_quantifier` exactly as `index` already is.
- `_quantifier` routing order, at the row only:
  1. The occurrence index, unchanged.
  2. The pivot: if `pivot(path)` returns a table and the body resolves against the
     pruned element type with `root="x"`, emit the semi-join below.
  3. The existing lambda over the elements.

```sql
-- exists
EXISTS (SELECT 1 FROM <table> AS x
        WHERE x.reaction_id = reaction_id AND <body>)
-- forall
NOT EXISTS (SELECT 1 FROM <table> AS x
            WHERE x.reaction_id = reaction_id AND NOT (<body>))
```

`NOT EXISTS`, never `NOT IN`. The spike found `NOT IN` agrees only because no
`reaction_id` in the corpus is NULL; a single NULL would make it answer NULL for every
reaction and return nothing.

A body that fails to resolve against the pruned type declines rather than raising: the
projection can still answer it.

### `ord_schema/agent/execute.py`

- `_pivot(path) -> str | None` — builds the pivot for one path on first use and returns
  its table name, or None if the budget refuses it.
- Built under `_narrow_build_lock`, cost measured as a `duckdb_memory()` delta, charged
  against the same budget as `_Narrow`, evicted by the same LRU, and not evicted while
  it has readers. A pivot is a narrow table of a different shape and reuses that
  machinery rather than duplicating it.
- A refusal falls back to the projection route and warns as `_warn_refused` does,
  naming the path and `narrow_budget_bytes`.
- The build is minutes per path on the full corpus (`outcomes.products` 469 s,
  `inputs.components` 692 s), so it logs before starting, the way the executor already
  announces that the occurrence index answers part of a query.

### Not in phase 1

- **Nested correlation.** A quantifier inside a quantifier keeps compiling to the
  lambda form. The prefix-join machinery is real work, and the single-level route
  already covers yield, color, temperature, role, and amounts.
- **Subsuming the occurrence index.** A pivot carries `structure_id`, so
  `get_bit($p, x.structure_id + offset) = 1` works against it, and it carries every
  other scalar field besides — so it answers strictly more than `_index_condition`,
  including bodies that decline today for binding a second field. Collapsing them is
  the likely end state, but it should follow a measurement showing the occurrence index
  no longer earns its 130 MiB, not happen as a side effect of this change.

## Phase 2 — the derived artifact

The phase 1 build cost belongs offline. Phase 2 derives the pivots as stamped artifacts
beside the projections and structures, so a container loads instead of builds.

- A build script writing one Parquet file per repeated level, sorted on the columns
  predicates use so row-group statistics can skip rather than merely decompress.
- Stamped through `ord_schema.artifacts` with the same version and source-content
  metadata the other derived artifacts carry, so `is_current` decides staleness the
  same way.
- The executor prefers a present artifact and falls back to building in process, which
  keeps phase 1's path alive as the degraded mode rather than deleting it.

## Testing

The corpus barely exercises the edge cases — 0 reactions with NULL or empty `outcomes`,
119 with NULL products under an outcome — so corpus agreement is weak evidence and
fixtures do the real work.

- **Synthetic cases** in the test corpus: level absent (NULL), level present but empty,
  NULL leaf inside a present element, several elements where only one matches, and
  **multiple outcomes**, since the real corpus is effectively single-outcome and that is
  exactly why a missing ordinal was invisible.
- **A differential test**: a table of query bodies, each asserted to return identical
  reaction-ID sets with the pivot route on and off. This is what catches the 942.
- **Sabotage checks**: dropping an ordinal from the prefix join must fail the suite.
  That precise mutation shipped silently in the probe that motivated this design.

## Risks

- **Size.** Keeping every non-repeated field makes a pivot larger than the 2.68 GiB
  measured, since those probes hand-picked leaves. Query time should not move, since
  unread struct children are not scanned. If it is too large the answer is pruning to
  referenced subtrees, which is `_Narrow`'s existing job.
- **Build latency in phase 1.** Minutes per path, paid by whichever query touches the
  path first. Phase 2 is the fix; until then it is logged, not hidden. Whether this is
  a regression against `_Narrow` materializing the same column is unmeasured and should
  be measured rather than assumed.
- **Two routes to keep in agreement.** The differential test is the guard, and it is
  the reason that test is a table rather than a handful of cases.
