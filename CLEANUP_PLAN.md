# Plan: ord-schema Cleanup

A deep-dive cleanup pass across the repo. Findings are grouped by confidence and
leverage. **The headline is that the codebase is already in good shape** — ruff
reports zero unused imports/vars, there is no dead code, no skipped/xfail tests,
no tracked build artifacts, and dependency/CI config is modern. So this is a
short, high-signal list, not a rewrite.

## Recommendation

Do **Tier 1** — small, mechanical, high-confidence wins with near-zero risk.
**Tier 2** is worthwhile but more opinionated (touches structure or test infra);
do it if we're already in these files. **Tier 3** is "noted, but stay put."

This is independent of the in-flight `speedup-add-rdkit` branch. The separate
validation effort (`VALIDATION_PLAN.md`) is now **complete** — every Tier 1 and
Tier 2 check landed in PR #817 (merged to `main`), Tier 3 was intentionally
skipped, and that planning doc has been deleted. Line references below were
re-verified against the post-merge tree.

---

## Tier 1 — mechanical, high-confidence, do these

1. **Drop the 11 copy-pasted unit validators into a factory.**
   [validations.py:1226-1301](ord_schema/validations.py#L1226) — `validate_time`,
   `validate_mass`, `validate_moles`, `validate_volume`, `validate_concentration`,
   `validate_pressure`, `validate_current`, `validate_voltage`, `validate_length`,
   `validate_wavelength`, `validate_flow_rate` are byte-for-byte identical
   (`check_value_and_units` + two `ensure_float_nonnegative`). `validate_temperature`
   is the one real exception (range per unit) and stays as-is.
   *Fix:* one `_validate_simple_unit(message)` helper registered directly in
   `_VALIDATOR_SWITCH` for those 11 types. ~75 lines → ~5. No behavior change.

2. **Collapse the ~20 `check_type_and_details`-only wrappers.**
   [validations.py:680-927](ord_schema/validations.py#L680) — many validators are
   nothing but `check_type_and_details(message)` (addition_device, vessel,
   vessel_material/attachment/preparation, temperature/pressure control,
   atmosphere, pressure_measurement, stirring, etc.).
   *Fix:* register `check_type_and_details` (and similarly the bare-`pass` no-ops)
   **directly as the value in `_VALIDATOR_SWITCH`** instead of via named wrappers.
   This *preserves* the deliberate "every message type must be listed" invariant
   documented at [validations.py:195-198](ord_schema/validations.py#L195)
   (`NOTE(ccoley)`) while deleting ~20 trivial functions. **Do not** simply delete
   the entries — the switch intentionally forces a decision per message type.

3. **Fix typo.** [validations.py:674](ord_schema/validations.py#L674) —
   `"the ReationInput"` → `"the ReactionInput"`.

4. **Remove stale working-tree cruft (local only, not tracked).**
   `build/lib/ord_schema/` (732K stale duplicate of the source tree, regenerated
   by `python -m build`) and the root `package-lock.json` (108K leftover; no root
   `package.json` — the real JS lives in `js/*/`). Both are already gitignored, so
   this is just `rm -rf build package-lock.json` to declutter the checkout. No
   commit needed.

---

## Tier 2 — worthwhile, slightly opinionated

6. **Turn on parallel tests in CI.** [run_tests.yml:56](.github/workflows/run_tests.yml#L56)
   runs `uv run pytest -vv --cov=ord_schema --durations=20` — serial — even though
   `pytest-xdist` is already a declared dev dep and the global default is
   `pytest -n auto`. Add `-n auto`. **Caveat:** the `orm/` postgres tests share a
   database fixture; verify they're either function-isolated or marked to run
   serially before flipping the switch (run locally with `-n auto` first). If they
   conflict, scope `-n auto` to the non-ORM suite or use `--dist loadgroup`.

7. **De-duplicate the two near-identical UPDATEs in `update_rdkit_ids`.**
   [orm/database.py:270-321](ord_schema/orm/database.py#L270) — the `Compound` and
   `ProductCompound` UPDATEs differ only in table/join path
   (`reaction_input` vs `reaction_outcome`). Extract a small helper that takes the
   table + join path. Also log per-statement rowcounts (the sibling
   `_update_rdkit_mols`/`_update_rdkit_reactions` already do) for parity and
   easier debugging of incremental runs.

8. **Document the partial indexes for existing databases.**
   [orm/mappers.py:240-260](ord_schema/orm/mappers.py#L240) adds partial indexes
   over unlinked rows (the recent speedup) and the comment says to apply them with
   `CREATE INDEX CONCURRENTLY` on existing DBs — but there's no migration
   snippet/SQL or README note. Add a short "applying to an existing database"
   section so deployed instances actually get the indexes.

9. **Replace `time.sleep(1)` in a test.**
   [message_helpers_test.py:467](ord_schema/message_helpers_test.py#L467) sleeps a
   full second in `test_gzip_reproducibility`. Control the mtime directly (or mock
   it) instead — removes a fixed 1s from every suite run.

---

## Tier 3 — noted, stay put

- **Pre-compiling validation regexes** (CAS `\d{2,7}-\d{2}-\d` at
  [validations.py:778](ord_schema/validations.py#L778), InChIKey at
  [:784](ord_schema/validations.py#L784), email at
  [:1222](ord_schema/validations.py#L1222)). Originally flagged as a perf win — it
  isn't. `re` already caches compiled patterns (`re._cache`, up to `re._MAXCACHE`
  = 512 entries), so a literal pattern passed to `re.fullmatch` is compiled once
  and reused; hoisting to a module constant saves only a dict lookup, not a
  recompile. With a handful of distinct patterns there is no measurable win. Skip
  unless touching these lines anyway for readability.
- **Compound get/set identifier wrappers** ([message_helpers.py:559-721](ord_schema/message_helpers.py#L559),
  `get/set_compound_smiles|name|molblock`). A factory would cut boilerplate, but
  these are public API with distinct docstrings/type signatures and good
  discoverability; a factory would *hurt* IDE/readability. Leave them.
- **CTE/`NOT EXISTS`/`MATERIALIZED` nesting in `_update_rdkit_mols`**
  ([orm/database.py:223-259](ord_schema/orm/database.py#L223)). Complex but
  correct, recently tuned, and well-commented. Don't touch working perf-critical
  SQL for aesthetics; just trim the now-historical "we used to use EXCEPT"
  asides if they distract.
- **Adding docstrings to every one-line validator.** Real but low-value churn;
  fold into Tier 1 #2 (most disappear) rather than a standalone pass.
- **Python 3.14 in the CI matrix** — verify rdkit/tensorflow/psycopg wheels
  exist before treating it as supported; this is a watch item, not a cleanup.

---

## Mechanics

- Tier 1 #1-#4 are independent; each is a self-contained commit.
- After the `validations.py` changes, run
  `pytest -n auto ord_schema/validations_test.py` — the switch refactor is
  behavior-preserving, so the existing suite is the regression guard.
- For Tier 2 #6, validate locally with `pytest -n auto ord_schema/orm/` before
  changing CI.
- ORM changes (#7-#8) should ride with or after the `speedup-add-rdkit` branch to
  avoid conflicts in `database.py`/`mappers.py`.
