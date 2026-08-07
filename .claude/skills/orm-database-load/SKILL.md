---
name: orm-database-load
description: Load or backfill the ORD ORM Postgres database and verify it before cutover. Use when running ord_schema.orm.scripts.add_datasets, populating a new ord_<date> search database, backfilling derived SMILES or RDKit links after a pipeline change, or checking a candidate database for parity with the live one.
---

# Loading and backfilling the ORM database

The ORM database is populated in two stages. **ingest** writes the `ord.*` search index and the
`public.*` payload; **derived** writes `derived.*` SMILES, the `rdkit.*` structures and their
links, and (opt-in) reaction classes. Every pass is guarded by `NOT EXISTS`, so all of it is
idempotent: safe to re-run, safe to kill and resume.

## Where it runs

Not on a laptop. Loads run on a dev VM inside the database's VPC, because the Aurora cluster has
no public endpoint and the job is chatty. The VM is normally stopped — start it, use it, stop it.

```sh
aws ec2 describe-instances --filters "Name=tag:Name,Values=dev-vm" \
  --query 'Reservations[].Instances[].{Id:InstanceId,State:State.Name}'
```

Reach it with SSM (`aws ssm send-command`, document `AWS-RunShellScript`); there is no SSH. Two
traps: SSM runs commands as **root under `dash`**, so wrap real work in
`sudo -u ubuntu -H bash -l`; and `fs.protected_regular` stops even root from reopening a
`ubuntu`-owned file in `/tmp`, so give each uploaded script a unique name.

### The RDKit version must match the database

`derived.*.smiles` and `rdkit.mols.smiles` are written by **Python RDKit**, not by the Postgres
cartridge, and the two do not agree on every molecule. The database stores `[Br][Ag]`; the
cartridge canonicalizes that same structure to `Br[Ag]`. Load from a machine with a different
RDKit and a few thousand organometallics will canonicalize differently, inserting duplicate
structures under different strings.

Any SMILES already in `rdkit.mols` must be a fixed point of the writer you are about to use.
Check before writing anything:

```python
from rdkit import Chem
for smiles in ("[Br][Ag]", "[B]#[Ni]", "[Cl][Ti]([Cl])([Cl])[Cl]"):
    assert Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == smiles, smiles
```

## Full load

Put `PGHOST` / `PGPORT` / `PGUSER` / `PGPASSWORD` / `PGSSLMODE` / `PGDATABASE` in the environment
(the master password lives in Secrets Manager) and let the DSN inherit them:

```sh
python -u -m ord_schema.orm.scripts.add_datasets \
  --pattern "$ORD_DATA/data/*/*.parquet" \
  --dsn "postgresql+psycopg://" \
  --n_jobs 16
```

`--classify_reactions` is opt-in, needs the optional reaction-class extra, and is slow. Omit it
unless you actually want to (re)classify.

## Backfill

`--stages derived` fills in derived rows over already-ingested datasets. It runs the same
per-dataset, per-shard RDKit passes as an end-to-end load, so nothing special is required:

```sh
python -u -m ord_schema.orm.scripts.add_datasets \
  --pattern "$ORD_DATA/data/*/*.parquet" --stages derived \
  --dsn "postgresql+psycopg://" --n_jobs 16
```

### Backfilling is not rebuilding

Every derived pass is guarded by `NOT EXISTS`. That is what makes the load idempotent, and it
also means an existing row is **never revisited**: after a change to how SMILES are derived,
the command above writes nothing and exits clean, leaving every stale row in place. Silence
here is not success.

`--rederive` deletes a dataset's derived SMILES rows before deriving them, which is what turns
the idempotent passes into a rebuild:

```sh
python -u -m ord_schema.orm.scripts.add_datasets \
  --pattern "$ORD_DATA/data/*/*.parquet" --stages derived --rederive \
  --dsn "postgresql+psycopg://" --n_jobs 16
```

Do not confuse it with `--overwrite`, which is an *ingest* flag meaning "re-ingest a dataset
whose MD5 changed."

Two things it deliberately leaves alone. `rdkit.mols` and `rdkit.reactions` are shared and
deduplicated by structure, so deleting one dataset's share would break every other dataset;
the link columns live on the deleted rows, so the RDKit pass re-links from scratch, reusing
structures already present. A rebuild that changes SMILES therefore leaves unreferenced
structures behind — space, not correctness. And `derived.reaction_classes` survives, because
classification is a slow opt-in pass that a rebuild should not silently redo.

Verify with counts before and after: a rebuild should end with the row counts it started with,
plus a spot-check of a value you expect to have changed.

### If you are on an ord-schema without the #895 fix

Before that fix, the RDKit passes scoped `rdkit_*_id IS NULL` per `(dataset, shard)` but the
planner drove each from the whole-database unlinked-rows partial index, applying the dataset and
shard predicates as filters. On an end-to-end load that never bites (the unlinked set drains
dataset by dataset), but a backfill leaves every dataset's rows unlinked at once, so each
`(dataset, shard)` rescans the entire unlinked set. Measured 2026-07-09 (~4.7M compounds, 53
datasets × 32 shards): **~1,084 s per shard, >9 h projected**, versus **372 s** for one global
pass. The #895 fix scopes each pass to the dataset via the resolved surrogate key, so this no
longer happens.

To recover a run that is exhibiting this — or on any version, to link a large standing unlinked
set in one pass:

1. Watch the log. When the `compound_smiles` bars finish and the `RDKit` bar appears, the derive
   is done and only the link remains.
2. Kill the job. Killing the client leaves server backends still executing their `UPDATE`s —
   terminate them, scoped to the target database and never the live one:

   ```sql
   SELECT pg_terminate_backend(pid) FROM pg_stat_activity
   WHERE datname = current_database() AND pid <> pg_backend_pid();
   ```

3. Run `scripts/global_link.py`, which inserts the missing `rdkit.mols` (keeping the `[Ti+5]`
   guard from issue #672) and links every unlinked row in one pass, chunked to bound WAL, then
   `ANALYZE`s. It is idempotent, so it is also a safe way to finish any partially-linked database.

## Should you scale up the writer?

Decide from metrics, not instinct.

- The writer is a **single instance with no reader**. Changing its class restarts it: measured
  ~3.5 minutes of downtime, and you pay that twice (up, then back down). The outage hits the
  *live* search database, which shares the instance.
- Check `CPUCreditBalance` first on burstable classes. Observed 2026-07-09: a 4.7M-row backfill
  ran at ~93% CPU and burned ~65–80 credits/hour against a full balance of 864 — hours of runway.
  The instance class was never the constraint.
- **If a job is projecting hours, suspect the query plan before the instance size.** Doubling the
  CPU does not fix a quadratic scan. `EXPLAIN` the statement that `pg_stat_activity` shows.
- Expect buffer-cache pollution regardless: the live database's index working set gets evicted,
  so search runs cold for a while afterwards. Prefer a low-traffic window.

## Monitoring

Committed row counts **lag**, because each shard's derive is a single transaction: `count(*)`
sits flat and then jumps. A flat count is not a stall. Check `pg_stat_activity` (state,
`wait_event_type`, oldest `xact_start`) and the `tqdm` bars in the log, which reflect in-flight
progress. Watch `CPUUtilization` and `CPUCreditBalance` alongside.

## Verification

Run `scripts/verify_coverage.sql` against the candidate database. Under
`psql -v ON_ERROR_STOP=1` a violated invariant exits non-zero, so a cutover wrapper can
gate on process status. What it checks:

- **Coverage** (enforced per dataset × scope): every compound carrying a structural identifier
  (SMILES/InChI/MolBlock) should have a derived SMILES row. Checked **per (dataset, attachment
  scope) cell** — `ord.compound` split by reaction input, **workup** input, and **product
  measurement** (three attachments easy to lose to a join that only follows
  `reaction_input.reaction_id`), plus `ord.product_compound` (reaction products), each grouped
  by the reaction's dataset. Two gates, both robust to the RDKit-unparseable residue every real
  load carries (organometallics, charged-N rings — a structural identifier the load-time RDKit
  cannot parse or reconstruct; ~3.2k across ~18.8M derived rows on a full load):
  - **Omission** (existence): `_update_compound_smiles` derives a dataset's compounds across
    all attachment paths in one unsegmented pass, so a cell it visited keeps at least one
    derived row (its parseable compounds) while a cell it skipped has **exactly zero**. Any
    cell holding ≥ `min_scope_size` structural compounds must have a derived row. No
    percentage tolerance — a healthy cell keeps its parseable derivations however many
    unparseable rows it also carries, and `min_scope_size` (default 50) absorbs a hypothetical
    tiny all-unparseable cell. Catches the single-dataset × single-path omission a corpus-wide
    or reaction-keyed count would dilute.
  - **Failed shard** (completeness, per **shard-hash bucket**): SMILES derivation is **sharded**
    — one `hashtextextended(id)` partition of a dataset's compound ids per worker, with
    `num_shards = min(32, ceil(reactions / 50000))`, so only datasets over ~50k reactions shard
    at all. The loader raises on a failed shard, but this gate does not trust that exit status.
    A failed shard leaves **exactly one hash bucket wholly underived** — parseable compounds
    included — while the unparseable residue spreads uniformly across every bucket (parseability
    is uncorrelated with the id hash) and so can never empty one. The same shard writes that
    dataset's **reaction** SMILES in the same transaction under the same hash, so its reaction
    bucket is empty too. The check recomputes each dataset's `num_shards`, buckets both tables by
    that hash, and judges an empty bucket against its **sibling buckets**: they measure how often
    this dataset's rows legitimately fail to derive, so the bucket is flagged when the siblings
    make "every one of mine missed by chance" implausible (combined probability <
    `max_false_positive`, default `1e-6`). Using **both** tables is what keeps a *sparse* hole
    visible: sharding is chosen by **reaction** count, so every bucket of a sharded dataset holds
    thousands of reactions even when a name-heavy one leaves only a handful of structural
    compounds there — the reaction evidence is strongest exactly where the compound evidence is
    weakest, and no absolute floor over either table alone spans both. `shard_size` /
    `shard_cap` mirror the loader's constants — keep them in step if those change.

  Name-only compounds are excluded, so the printed `compounds > derived` gap does not count.
  Override any threshold with `-v NAME=N`.
- **Per-dataset reaction derivation** (enforced): an independent guard on
  `derived.reaction_smiles`, the separate table the compound gate never touches. `reaction_smiles`
  carries the reaction's `dataset_id`, so it is a cheap hash aggregate over `ord.reaction` (no
  compound→reaction join) — flag any dataset with fewer than half its reactions derived. A dataset
  omitted from the derive pass has ~none; the reaction residual stays well under half.
- **No stray unlinked rows** (enforced): the only unlinked derived rows are `[Ti+5]`
  structures, which `_update_rdkit_mols` deliberately keeps out of `rdkit.mols`. Any other
  unlinked row fails the script.
- **Payload parity** (enforced): every `ord.reaction` has its `public.reactions` proto and the
  dataset counts agree.

Two further checks worth doing before a cutover:

- **Dataset freshness**, using the loader's own hash: `public.datasets.md5` must equal
  `ord_schema.parquet.DatasetView(path).md5()` for every Parquet file in `ord-data`.
- **Reaction-set equality** against the database you are replacing — an order-independent
  fingerprint beats comparing counts alone:

  ```sql
  SELECT count(*), count(DISTINCT reaction_id), sum(hashtext(reaction_id)::bigint)
  FROM ord.reaction;
  ```

## Cleanup

Stop the dev VM, close any SSM tunnels, and do not leave the master password on disk.

## Gotchas

- `derived.compound_smiles` legitimately has fewer rows than `ord.compound`. Compounds whose only
  identifier is a name have no derivable SMILES and get no row. That is not data loss.
- `rdkit.mols` is deduplicated by SMILES **string**, not by structure, so two spellings of one
  molecule become two rows. This is why the writer's RDKit version matters.
- Reaction classification is a separate opt-in pass; a database can be complete and correct with
  `derived.reaction_classes` nearly empty.
