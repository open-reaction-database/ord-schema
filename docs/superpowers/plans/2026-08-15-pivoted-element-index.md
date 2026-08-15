# Pivoted Element Index Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Answer quantifiers over repeated levels with a semi-join against a flat
one-row-per-element table instead of a `list_filter` lambda over the nested projection.

**Architecture:** A schema walk derives, for every repeated level, a pivot whose rows
carry `reaction_id`, the ordinals of every repeated level from the root down to and
including that level, and the element's non-repeated fields with struct nesting intact.
The compiler routes a quantifier to `EXISTS` / `NOT EXISTS` against that table when the
body resolves against the pruned element type; the executor builds the table lazily,
charged against the existing narrow-table budget. Phase 2 ships the tables as stamped
artifacts so a container loads instead of builds.

**Tech Stack:** Python 3.11+, pyarrow, DuckDB 1.5.5, pytest, ruff, ty.

**Spec:** `docs/superpowers/specs/2026-08-15-pivoted-element-index-design.md`

## Global Constraints

- Google-style docstrings; every parameter gets an `Args:` entry. No lazy imports.
- Comments describe the present, never the change that produced it.
- American English. Ruff line length 88.
- Free functions over `@staticmethod`.
- `pytest -n auto` for the suite.
- Emitted SQL uses `NOT EXISTS`, never `NOT IN`: a NULL `reaction_id` would make
  `NOT IN` answer NULL for every reaction and return nothing.
- A pivot is unfiltered — every element gets a row — because completeness is what makes
  `forall` indexable.

## File Structure

| file | responsibility |
| --- | --- |
| `ord_schema/agent/pivot.py` (create) | pure schema functions: the repeated-level walk, the pruned element type, table naming. Imports pyarrow and `ord_schema.projection` only, so `query.py` may import it without a cycle. |
| `ord_schema/agent/pivot_test.py` (create) | tests for the walk |
| `ord_schema/agent/query.py` (modify) | `PivotIndex`, `compile_query(pivot=...)`, `_quantifier` routing |
| `ord_schema/agent/execute.py` (modify) | generalized cache key, `_pivot`, wiring |
| `ord_schema/agent/execute_test.py` (modify) | build, budget, and differential tests |

---

### Task 1: The repeated-level walk

**Files:**

- Create: `ord_schema/agent/pivot.py`
- Test: `ord_schema/agent/pivot_test.py`

**Interfaces:**

- Consumes: `ord_schema.projection.SCHEMA`, `projection.is_internal`
- Produces:

```python
@dataclasses.dataclass(frozen=True)
class RepeatedLevel:
    path: str                    # dotted query path, e.g. "outcomes.products"
    ordinals: tuple[str, ...]    # ("outcome_index", "product_index")
    element_type: pa.DataType    # struct, repeated fields removed recursively

def repeated_levels(schema: pa.Schema = projection.SCHEMA) -> dict[str, RepeatedLevel]
def table_name(path: str) -> str          # "outcomes.products" -> "pivot_outcomes_products"
LEVELS: dict[str, RepeatedLevel]          # module-level, built at import
```

- [ ] **Step 1: Write the failing tests**

```python
def test_covers_known_levels():
    levels = pivot.repeated_levels()
    assert "inputs.components" in levels
    assert "outcomes.products" in levels
    assert "outcomes.products.measurements" in levels
    assert "workups" in levels


def test_ordinals_accumulate_down_the_path():
    levels = pivot.repeated_levels()
    assert levels["workups"].ordinals == ("workup_index",)
    assert levels["outcomes.products"].ordinals == ("outcome_index", "product_index")
    assert levels["outcomes.products.measurements"].ordinals == (
        "outcome_index",
        "product_index",
        "measurement_index",
    )


def test_element_type_drops_repeated_fields_but_keeps_structs():
    element = pivot.repeated_levels()["outcomes.products"].element_type
    names = [field.name for field in element]
    assert "isolated_color" in names          # scalar kept
    assert "texture" in names                 # struct kept
    assert "measurements" not in names        # list dropped
    assert "identifiers" not in names         # list dropped
    assert "analyses" not in names            # map dropped
    texture = element.field("texture").type
    assert [field.name for field in texture] == ["type", "details"]


def test_nested_repeated_inside_a_struct_is_dropped_recursively():
    element = pivot.repeated_levels()["outcomes.products.measurements"].element_type
    standard = element.field("authentic_standard").type
    names = [field.name for field in standard]
    assert "smiles" in names
    assert "identifiers" not in names          # list inside a kept struct


def test_a_map_level_is_covered():
    # inputs is MAP<VARCHAR, STRUCT>, so its components list sits under a map.
    assert "inputs.components" in pivot.repeated_levels()


def test_table_name_is_an_identifier():
    assert pivot.table_name("outcomes.products") == "pivot_outcomes_products"
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest ord_schema/agent/pivot_test.py -x -q`
Expected: FAIL, no module named `pivot`.

- [ ] **Step 3: Implement the walk**

Recurse the schema. On a list or map, the level's path is the path so far; recurse into
the value type carrying an ordinal named `<singular>_index` derived from the last path
segment. On a struct, recurse into children. `element_type` prunes by rebuilding the
struct with list/map fields omitted, recursing into kept struct children. Skip
`projection.is_internal` fields nowhere — `structure_id` must be kept, since a pivot can
answer structure predicates.

- [ ] **Step 4: Run to verify they pass**

Run: `pytest ord_schema/agent/pivot_test.py -q`

- [ ] **Step 5: Lint and commit**

```bash
ruff format ord_schema/agent/pivot.py ord_schema/agent/pivot_test.py
ruff check ord_schema/agent/pivot.py ord_schema/agent/pivot_test.py
git add ord_schema/agent/pivot.py ord_schema/agent/pivot_test.py
git commit -m "Walk the schema for every repeated level"
```

---

### Task 2: Compiler routing

**Files:**

- Modify: `ord_schema/agent/query.py`
- Test: `ord_schema/agent/query_test.py`

**Interfaces:**

- Consumes: `pivot.LEVELS`, `pivot.table_name`, `resolve`
- Produces:

```python
PivotIndex = Callable[[str], str | None]
def compile_query(query, *, schema=..., table=TABLE, index=None, pivot=None) -> Compiled
```

- [ ] **Step 1: Write the failing tests**

```python
def _sql(body, pivot=None):
    return query.compile_query(query.Query.model_validate(body), pivot=pivot).sql


def _everything(path):
    return pivot_module.table_name(path)


def test_exists_over_a_covered_path_becomes_a_semi_join():
    sql = _sql(
        {"where": {"op": "exists", "path": "workups",
                   "where": {"op": "eq", "path": "type",
                             "value": {"literal": "EXTRACTION"}}}},
        pivot=_everything,
    )
    assert "EXISTS (SELECT 1 FROM pivot_workups AS x" in sql
    assert "list_filter" not in sql


def test_forall_becomes_a_negated_semi_join():
    sql = _sql(
        {"where": {"op": "forall", "path": "workups",
                   "where": {"op": "eq", "path": "type",
                             "value": {"literal": "EXTRACTION"}}}},
        pivot=_everything,
    )
    assert "NOT EXISTS (SELECT 1 FROM pivot_workups AS x" in sql
    assert "NOT IN" not in sql


def test_a_body_reaching_a_repeated_field_declines():
    # workups.input.components is a list on the workup element, so the pruned type
    # has no such field and the quantifier stays a lambda over the projection.
    sql = _sql(
        {"where": {"op": "exists", "path": "workups",
                   "where": {"op": "exists", "path": "input.components",
                             "where": {"op": "eq", "path": "reaction_role",
                                       "value": {"literal": "SOLVENT"}}}}},
        pivot=_everything,
    )
    assert "list_filter" in sql


def test_without_a_pivot_nothing_changes():
    body = {"where": {"op": "exists", "path": "workups",
                      "where": {"op": "eq", "path": "type",
                                "value": {"literal": "EXTRACTION"}}}}
    assert _sql(body) == _sql(body, pivot=lambda path: None)
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest ord_schema/agent/query_test.py -x -q -k pivot`

- [ ] **Step 3: Implement**

Add `PivotIndex`, thread `pivot` through `compile_query` → `_predicate` → `_quantifier`
beside `index`. In `_quantifier`, after the occurrence route declines and while
`scope is None`, look up `pivot(node.path)`; if it returns a table and
`pivot.LEVELS` has the path, compile the body with `_predicate(node.where, "x",
level.element_type, ...)` inside a `try/except QueryError` — a body the pruned type
cannot resolve declines to the lambda form. Emit the `EXISTS` / `NOT EXISTS` forms.

- [ ] **Step 4: Run to verify they pass, then the whole query suite**

Run: `pytest ord_schema/agent/query_test.py -q`

- [ ] **Step 5: Lint and commit**

```bash
ruff format ord_schema/agent/query.py ord_schema/agent/query_test.py
ruff check ord_schema/agent/query.py ord_schema/agent/query_test.py
git commit -am "Compile a covered quantifier to a semi-join over its pivot"
```

---

### Task 3: Executor build and budget

**Files:**

- Modify: `ord_schema/agent/execute.py` (`_Narrow`, `_materialize`, `_cached`,
  `_build`, `_evict`, `_warn_refused`, `search`)
- Test: `ord_schema/agent/execute_test.py`

**Interfaces:**

- Produces:

```python
@dataclasses.dataclass(frozen=True)
class _Held:
    """What a cache entry holds: a column set, or a pivot over one repeated level."""
    kind: str                 # "columns" | "pivot"
    parts: tuple[str, ...]    # sorted column names, or (path,)

def _pivot(self, path: str) -> str | None
```

- [ ] **Step 1: Write the failing tests**

```python
def test_a_pivot_is_built_once_and_reused(deep_corpus, caplog):
    where = {"op": "exists", "path": "workups",
             "value": None, "where": {"op": "eq", "path": "type",
                                      "value": {"literal": "EXTRACTION"}}}
    first = _search(deep_corpus, where)
    second = _search(deep_corpus, where)
    assert first == second == {"ord-cc01"}
    assert sum("pivot" in record.message for record in caplog.records) == 1


def test_a_zero_budget_builds_no_pivot(corpus_dir):
    with execute.Corpus(..., narrow_budget_bytes=0) as value:
        assert value._pivot("workups") is None
```

- [ ] **Step 2: Run to verify they fail**

- [ ] **Step 3: Implement**

Replace the `frozenset[str]` cache key with `_Held` throughout, so one budget and one
LRU cover both kinds. `_build` takes the `CREATE TABLE ... AS <select>` body from its
caller. `_pivot(path)` assembles the select with `unnest(... ) WITH ORDINALITY` for each
repeated level on the path, projecting the pruned element's fields plus the ordinals and
`reaction_id`. Wire `pivot=self._pivot` into `compile_query` in `search`, and hold a
read on the entry for the query's lifetime the way `_narrowed_table` does.

- [ ] **Step 4: Run**

Run: `pytest ord_schema/agent/execute_test.py -q -n auto`

- [ ] **Step 5: Commit**

---

### Task 4: Differential test and synthetic fixtures

**Files:**

- Modify: `ord_schema/agent/execute_test.py`

- [ ] **Step 1: Add a fixture the real corpus does not provide**

`wide_corpus`: one reaction with **two outcomes**, the first carrying a non-desired
product with a 90% yield and the second a desired product with a 10% yield — so an
outcome-blind or product-blind correlation answers differently from the nested form.
Plus a reaction with an empty outcomes list, one with a NULL leaf inside a present
product, and one with several products where only the last matches.

- [ ] **Step 2: Write the differential test**

```python
@pytest.mark.parametrize("where", _DIFFERENTIAL_CASES)
def test_pivot_and_projection_agree(wide_corpus, where):
    with_pivot = _search(wide_corpus, where)
    with mock.patch.object(execute.Corpus, "_pivot", lambda self, path: None):
        without = _search(wide_corpus, where)
    assert with_pivot == without
```

`_DIFFERENTIAL_CASES` covers `exists`, `forall`, `not exists`, a body with two
conditions on one element, and a body whose leaf is NULL for some elements.

- [ ] **Step 3: Run, then sabotage**

Drop an ordinal from the emitted prefix join and confirm the suite fails. Restore.

- [ ] **Step 4: Commit**

---

### Task 5: Pivot artifacts

**Files:**

- Modify: `ord_schema/agent/pivot.py` (add the writer)
- Create: `ord_schema/scripts/build_pivots.py`
- Test: `ord_schema/agent/pivot_test.py`

**Interfaces:**

```python
def write_pivot(source, output, path, *, source_md5=None, compression="zstd") -> int
```

One Parquet file per (projection file, repeated level), laid out as
`pivots/<level>/<source>.parquet` so a level globs to one relation. Stamped through
`ord_schema.artifacts.to_metadata` with the same fields `write_structures` uses, so
`is_current` decides staleness identically.

- [ ] **Step 1:** test that a written pivot round-trips its rows and carries stamps
- [ ] **Step 2:** implement over one projection row group at a time
- [ ] **Step 3:** script looping levels × sources
- [ ] **Step 4:** commit

---

### Task 6: Load artifacts in preference to building

**Files:**

- Modify: `ord_schema/agent/execute.py`

- [ ] **Step 1:** test that a corpus opened with a pivots pattern registers views over
  the artifacts and never builds in process
- [ ] **Step 2:** `Corpus(..., pivots_pattern=None)`; when given, `_pivot(path)` returns
  a view over `pivots/<level>/*.parquet` rather than building, and pairs them the way
  projections and structures are paired
- [ ] **Step 3:** keep the in-process build as the fallback when no pattern is given
- [ ] **Step 4:** commit

---

## Self-Review

**Spec coverage:** walk → Task 1; compiler routing and `NOT EXISTS` → Task 2; lazy build,
budget reuse, refusal warning → Task 3; synthetic fixtures, differential test, sabotage
→ Task 4; artifact and loader → Tasks 5–6. Non-goals (nested correlation, occurrence
subsumption) are not implemented anywhere, as intended.

**Types:** `RepeatedLevel.element_type` is consumed by Task 2's `_predicate` call and
Task 3's select assembly; `table_name` by Tasks 2, 3, and 6; `_Held` by Task 3 only.

**Known gap:** Task 3's `_pivot` select assembly is the one step whose SQL is not written
out here, because it depends on the `resolve` traversal for each level and is best read
off `_occurrences`, which already assembles the analogous `unnest` per indexed path.
