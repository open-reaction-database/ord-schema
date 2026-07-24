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

`--stages derived` re-derives over already-ingested datasets. It runs the same per-dataset,
per-shard RDKit passes as an end-to-end load, so nothing special is required:

```sh
python -u -m ord_schema.orm.scripts.add_datasets \
  --pattern "$ORD_DATA/data/*/*.parquet" --stages derived \
  --dsn "postgresql+psycopg://" --n_jobs 16
```

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

- **Coverage** (enforced above a tolerance): every compound carrying a structural identifier
  (SMILES/InChI/MolBlock) should have a derived SMILES row — across `ord.compound` (reaction
  input, **workup** input, or **product measurement**, three attachments easy to lose to a join
  that only follows `reaction_input.reaction_id`) and `ord.product_compound` (reaction products).
  The count of structural-identifier compounds without a derived row is printed, and the script
  fails only when it exceeds a small tolerance (default 0.1% of derived rows) — enough to catch a
  skipped compound, dataset, or attachment path, while tolerating the RDKit-unparseable residual
  that every real load carries (~0.02%: organometallics, charged-N rings — a structural
  identifier the load-time RDKit cannot parse or reconstruct). Name-only compounds are excluded,
  so the printed `compounds > derived` gap does not count. Lower the tolerance to tighten the gate.
- **No stray unlinked rows** (enforced): the only unlinked derived rows are `[Ti+5]`
  structures, which `_update_rdkit_mols` deliberately keeps out of `rdkit.mols`. Any other
  unlinked row fails the script.
- **Payload parity** (enforced): every `ord.reaction` has its `public.reactions` proto and the
  dataset counts agree.

Two further checks worth doing before a cutover:

- **Dataset freshness**, using the loader's own hash: `public.datasets.md5` must equal
  `ord_schema.parquet.streaming_md5(path)` for every Parquet file in `ord-data`.
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
