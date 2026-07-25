-- Copyright 2026 Open Reaction Database Project Authors
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

-- Verifies derived/RDKit coverage for a freshly loaded or backfilled ORM database.
-- Read-only. The invariant checks at the end RAISE on violation, so pass ON_ERROR_STOP
-- for a failure to exit non-zero -- without it psql reports the error but still exits 0:
--   psql -v ON_ERROR_STOP=1 "$DSN" -f verify_coverage.sql

\echo === Table sizes ===
SELECT 'ord.compound' AS table, count(*) FROM ord.compound
UNION ALL SELECT 'derived.compound_smiles', count(*) FROM derived.compound_smiles
UNION ALL SELECT 'derived.product_compound_smiles', count(*) FROM derived.product_compound_smiles
UNION ALL SELECT 'derived.reaction_smiles', count(*) FROM derived.reaction_smiles
UNION ALL SELECT 'derived.reaction_classes', count(*) FROM derived.reaction_classes
UNION ALL SELECT 'rdkit.mols', count(*) FROM rdkit.mols
UNION ALL SELECT 'rdkit.reactions', count(*) FROM rdkit.reactions;

\echo
\echo === Coverage by how each compound reaches its reaction ===
-- A compound hangs off a reaction input, a workup's input (whose ord.reaction_input row has
-- reaction_id NULL), or a product measurement (authentic standards; reaction_input_id NULL).
-- A dataset-scoping join that follows only reaction_input.reaction_id silently drops the latter
-- two, so all three must appear here with `derived` covering every derivable compound.
--
-- `compounds` > `derived` is expected and healthy: compounds whose only identifier is a name have
-- no derivable SMILES and get no derived row. `unlinked` must be 0 -- [Ti+5] structures are
-- excluded because _update_rdkit_mols deliberately keeps them out of rdkit.mols (issue #672).
SELECT CASE
         WHEN ri.reaction_id IS NOT NULL THEN 'reaction_input'
         WHEN c.reaction_input_id IS NOT NULL THEN 'workup_input'
         ELSE 'product_measurement'
       END AS attachment,
       count(*) AS compounds,
       count(d.compound_id) AS derived,
       count(*) FILTER (
           WHERE d.compound_id IS NOT NULL
             AND d.rdkit_mol_id IS NULL
             AND d.smiles NOT LIKE '%[Ti+5]%'
       ) AS unlinked
FROM ord.compound c
LEFT JOIN ord.reaction_input ri ON c.reaction_input_id = ri.id
LEFT JOIN derived.compound_smiles d ON d.compound_id = c.id
GROUP BY 1
ORDER BY 2 DESC;

\echo
\echo === Unlinked derived rows (expect only [Ti+5] structures) ===
SELECT left(smiles, 60) AS smiles, count(*)
FROM derived.compound_smiles WHERE rdkit_mol_id IS NULL GROUP BY 1
UNION ALL
SELECT left(smiles, 60), count(*)
FROM derived.product_compound_smiles WHERE rdkit_mol_id IS NULL GROUP BY 1;

\echo
\echo === Reaction-set fingerprint (compare against the database being replaced) ===
-- Order-independent, so it catches a missing or duplicated reaction that equal counts would hide.
SELECT count(*) AS reactions,
       count(DISTINCT reaction_id) AS distinct_ids,
       sum(hashtext(reaction_id)::bigint) AS fingerprint
FROM ord.reaction;

\echo
\echo === Payload parity: every reaction has its served proto ===
SELECT (SELECT count(*) FROM ord.reaction) AS index_rows,
       (SELECT count(*) FROM public.reactions) AS payload_rows,
       (SELECT count(*) FROM ord.dataset) AS datasets,
       (SELECT count(*) FROM public.datasets) AS dataset_metadata,
       (SELECT sum(num_reactions) FROM public.datasets) AS sum_num_reactions;

\echo
\echo === Invariant checks (raise, so exit status gates a cutover) ===
-- The sections above only print; these raise on a violated invariant so a wrapper running under
-- `psql -v ON_ERROR_STOP=1` fails non-zero instead of accepting a bad candidate database. Still
-- read-only. The fingerprint stays print-only: it is a manual compare against the database being
-- replaced, not a self-contained invariant.
-- Compound linkage only. derived.reaction_smiles is intentionally excluded: a reaction
-- SMILES the cartridge cannot parse is legitimately unlinked, with no clean [Ti+5]-style
-- carve-out to distinguish it from a real gap.
DO $$
DECLARE
    unlinked bigint;
BEGIN
    SELECT count(*) INTO unlinked
    FROM (
        SELECT smiles FROM derived.compound_smiles WHERE rdkit_mol_id IS NULL
        UNION ALL
        SELECT smiles FROM derived.product_compound_smiles WHERE rdkit_mol_id IS NULL
    ) t
    WHERE smiles NOT LIKE '%[Ti+5]%';
    IF unlinked > 0 THEN
        RAISE EXCEPTION 'coverage: % unlinked derived row(s) that are not [Ti+5]', unlinked;
    END IF;
END $$;

-- Every compound carrying a structural identifier should derive a SMILES; a load that
-- skipped a compound leaves structural-identifier compounds with no derived row. NAME-only
-- compounds are excluded, so the documented `compounds > derived` gap does not trip this.
-- The type list mirrors message_helpers.STRUCTURAL_IDENTIFIER_TYPES.
--
-- Checked PER (dataset, attachment scope) cell: ord.compound split by its three attachment
-- paths (reaction input, workup input, product measurement) plus ord.product_compound, each
-- grouped by the reaction's dataset. Two gates over the same cells, both robust to the
-- RDKit-unparseable residue every real load carries (organometallics such as Ir photocatalysts
-- or C[Al](C)(C)([Li])..., charged-N rings -- identifiers the load-time RDKit cannot parse or
-- reconstruct; observed ~3.2k across ~18.8M derived rows on a full load):
--
--   1. OMISSION (derived>0 existence). _update_compound_smiles derives a dataset's compounds
--      across all attachment paths together in one unsegmented pass, so a cell it visited keeps
--      at least one derived row (its parseable compounds) while a cell it skipped -- a whole
--      dataset, a whole path, or their intersection -- has exactly zero. Any cell holding
--      >= min_scope_size structural compounds must have a derived row. No percentage tolerance:
--      a healthy cell keeps its parseable derivations regardless of how many unparseable rows it
--      also carries, and min_scope_size absorbs a hypothetical tiny all-unparseable cell (real
--      reagent/product cells always contain parseable organics). Catches the single dataset x
--      single path omission a corpus-wide or reaction-keyed count would dilute.
--
--   2. PARTIAL derivation (completeness), detected by SHARD SHAPE. SMILES derivation is sharded
--      (loading.py / ord_schema.orm.sharding: one hashtextextended(id) partition of a dataset's
--      compound ids per worker; num_shards = min(32, ceil(reactions / 50000)), so only datasets
--      over ~50k reactions shard at all). The loader raises on a failed shard, but a
--      defence-in-depth gate should not trust that exit status. A failed shard leaves exactly one
--      hash bucket entirely underived -- parseable compounds included -- while the unparseable
--      residue spreads uniformly across all buckets (unparseability is uncorrelated with the id
--      hash) and so never empties one. So recompute each dataset's num_shards from its reaction
--      count, bucket its structural compounds by the same ((h % n) + n) % n hash, and flag any
--      bucket of a sharded dataset with >= min_bucket_size structural compounds but no derived
--      row. Detecting the hole by its shape (a wholly-empty hash class) rather than its size needs
--      no residue-calibrated tolerance and catches the failure however small or name-heavy the
--      dataset is. A healthy database populates every bucket whatever num_shards is, so an
--      inaccurate recomputation can only miss a hole, never invent one.
--
-- The total underived count is printed. Defaults: min_scope_size 50, min_bucket_size 50, and
-- shard_size/shard_cap mirroring loading._DERIVE_SHARD_SIZE/_DERIVE_SHARD_CAP -- keep those two
-- in step with the loader if its constants change (override any with -v NAME=N).
\if :{?min_scope_size}
\else
  \set min_scope_size 50
\endif
\if :{?min_bucket_size}
\else
  \set min_bucket_size 50
\endif
\if :{?shard_size}
\else
  \set shard_size 50000
\endif
\if :{?shard_cap}
\else
  \set shard_cap 32
\endif
SELECT set_config('ord.verify_min_scope_size', :'min_scope_size', false) AS min_scope_size,
       set_config('ord.verify_min_bucket_size', :'min_bucket_size', false) AS min_bucket_size,
       set_config('ord.verify_shard_size', :'shard_size', false) AS shard_size,
       set_config('ord.verify_shard_cap', :'shard_cap', false) AS shard_cap;
DO $$
DECLARE
    total_missing bigint;
    bad_omitted bigint;
    example_omitted text;
    bad_shard bigint;
    example_shard text;
    min_scope_size int := current_setting('ord.verify_min_scope_size')::int;
    min_bucket_size int := current_setting('ord.verify_min_bucket_size')::int;
    shard_size int := current_setting('ord.verify_shard_size')::int;
    shard_cap int := current_setting('ord.verify_shard_cap')::int;
BEGIN
    WITH scoped AS (
        SELECT r.dataset_id AS dataset_id, 'reaction_input' AS scope,
               hashtextextended(c.id::text, 0) AS id_hash,
               EXISTS (
                   SELECT 1 FROM ord.compound_identifier ci
                   WHERE ci.compound_id = c.id
                     AND ci.type IN ('SMILES', 'INCHI', 'MOLBLOCK') AND ci.value <> ''
               ) AS structural,
               EXISTS (
                   SELECT 1 FROM derived.compound_smiles d WHERE d.compound_id = c.id
               ) AS derived
        FROM ord.compound c
        JOIN ord.reaction_input ri ON c.reaction_input_id = ri.id
        JOIN ord.reaction r ON ri.reaction_id = r.id
        UNION ALL
        SELECT r.dataset_id, 'workup_input',
               hashtextextended(c.id::text, 0),
               EXISTS (
                   SELECT 1 FROM ord.compound_identifier ci
                   WHERE ci.compound_id = c.id
                     AND ci.type IN ('SMILES', 'INCHI', 'MOLBLOCK') AND ci.value <> ''
               ),
               EXISTS (
                   SELECT 1 FROM derived.compound_smiles d WHERE d.compound_id = c.id
               )
        FROM ord.compound c
        JOIN ord.reaction_input ri ON c.reaction_input_id = ri.id
        JOIN ord.reaction_workup rw ON ri.reaction_workup_id = rw.id
        JOIN ord.reaction r ON rw.reaction_id = r.id
        UNION ALL
        SELECT r.dataset_id, 'product_measurement',
               hashtextextended(c.id::text, 0),
               EXISTS (
                   SELECT 1 FROM ord.compound_identifier ci
                   WHERE ci.compound_id = c.id
                     AND ci.type IN ('SMILES', 'INCHI', 'MOLBLOCK') AND ci.value <> ''
               ),
               EXISTS (
                   SELECT 1 FROM derived.compound_smiles d WHERE d.compound_id = c.id
               )
        FROM ord.compound c
        JOIN ord.product_measurement pm ON c.product_measurement_id = pm.id
        JOIN ord.product_compound pc ON pm.product_compound_id = pc.id
        JOIN ord.reaction_outcome ro ON pc.reaction_outcome_id = ro.id
        JOIN ord.reaction r ON ro.reaction_id = r.id
        UNION ALL
        SELECT r.dataset_id, 'product_compound',
               hashtextextended(pc.id::text, 0),
               EXISTS (
                   SELECT 1 FROM ord.compound_identifier ci
                   WHERE ci.product_compound_id = pc.id
                     AND ci.type IN ('SMILES', 'INCHI', 'MOLBLOCK') AND ci.value <> ''
               ),
               EXISTS (
                   SELECT 1 FROM derived.product_compound_smiles d
                   WHERE d.product_compound_id = pc.id
               )
        FROM ord.product_compound pc
        JOIN ord.reaction_outcome ro ON pc.reaction_outcome_id = ro.id
        JOIN ord.reaction r ON ro.reaction_id = r.id
    ),
    per_cell AS (
        SELECT dataset_id, scope,
               count(*) FILTER (WHERE structural) AS structural,
               count(*) FILTER (WHERE structural AND derived) AS derived
        FROM scoped
        GROUP BY dataset_id, scope
    ),
    -- num_shards the loader used for each dataset's SMILES derivation:
    -- min(shard_cap, ceil(size / shard_size)), at least 1. Mirrors loading._derive_shard_items,
    -- reading size from the same place it does (get_dataset_size -> public.datasets.num_reactions,
    -- 0 when the metadata row is missing) rather than counting ord.reaction -- if the two ever
    -- disagreed across a shard-size boundary, a recomputed num_shards would bucket the compounds
    -- differently than the loader did and smear a shard hole across buckets instead of exposing it.
    dataset_shards AS (
        SELECT d.id AS dataset_id,
               greatest(1, least(shard_cap,
                                 ceil(coalesce(pd.num_reactions, 0)::numeric / shard_size)::int))
                   AS num_shards
        FROM ord.dataset d
        LEFT JOIN public.datasets pd ON pd.dataset_id = d.dataset_id
    ),
    -- Bucket each dataset's structural compounds by the loader's shard hash, so a failed shard
    -- shows up as a bucket with structural rows but no derived row (the residue spreads across
    -- every bucket and cannot empty one). num_shards = 1 datasets never sharded, so their lone
    -- bucket only restates the omission gate and is excluded below.
    per_bucket AS (
        SELECT dataset_id, num_shards, bucket,
               count(*) FILTER (WHERE structural) AS structural,
               count(*) FILTER (WHERE structural AND derived) AS derived
        FROM (
            SELECT s.dataset_id, ds.num_shards,
                   ((s.id_hash % ds.num_shards) + ds.num_shards) % ds.num_shards AS bucket,
                   s.structural, s.derived
            FROM scoped s
            JOIN dataset_shards ds ON ds.dataset_id = s.dataset_id
        ) b
        GROUP BY dataset_id, num_shards, bucket
    )
    -- Scalar subqueries over the shared CTEs; scoped (hence the UNION) materializes once.
    SELECT
        (SELECT coalesce(sum(structural - derived), 0) FROM per_cell),
        (SELECT count(*) FROM per_cell
             WHERE structural >= min_scope_size AND derived = 0),
        (SELECT min(scope || ' in dataset ' || dataset_id) FROM per_cell
             WHERE structural >= min_scope_size AND derived = 0),
        (SELECT count(*) FROM per_bucket
             WHERE num_shards > 1 AND structural >= min_bucket_size AND derived = 0),
        (SELECT min('shard ' || bucket || '/' || num_shards || ' of dataset ' || dataset_id
                    || ' (' || structural || ' structural compounds, 0 derived)')
             FROM per_bucket
             WHERE num_shards > 1 AND structural >= min_bucket_size AND derived = 0)
      INTO total_missing, bad_omitted, example_omitted, bad_shard, example_shard;
    RAISE INFO 'coverage: % structural-identifier compound(s) lack a derived row (existence gate per (dataset, scope), min scope size %; shard-hole gate per hash bucket, min bucket size %)', total_missing, min_scope_size, min_bucket_size;
    IF bad_omitted > 0 THEN
        RAISE EXCEPTION 'coverage: % (dataset, attachment scope) cell(s) hold >= % structural-identifier compounds but no derived row (e.g. %) -- a dataset or attachment path omitted from the derive pass', bad_omitted, min_scope_size, example_omitted;
    END IF;
    IF bad_shard > 0 THEN
        RAISE EXCEPTION 'coverage: % shard-hash bucket(s) hold >= % structural-identifier compounds but no derived row (e.g. %) -- a failed derivation shard leaves exactly one bucket wholly underived', bad_shard, min_bucket_size, example_shard;
    END IF;
END $$;

-- Per-dataset reaction derivation completeness. Independent of the compound gate above: it
-- validates derived.reaction_smiles, a separate table the per-cell compound check never
-- touches. reaction_smiles carries the reaction's dataset_id, so this is cheap (a hash
-- aggregate over ord.reaction, no compound->reaction join) -- flag any dataset with fewer than
-- half its reactions derived. A dataset omitted from the derive pass has ~none derived; the
-- ~0.4% unparseable reaction residual and normal variation stay well under half.
DO $$
DECLARE
    omitted bigint;
BEGIN
    SELECT count(*) INTO omitted
    FROM (
        SELECT r.dataset_id, count(*) AS reactions, count(rs.reaction_id) AS derived
        FROM ord.reaction r
        LEFT JOIN derived.reaction_smiles rs ON rs.reaction_id = r.id
        GROUP BY r.dataset_id
    ) t
    WHERE derived * 2 < reactions;
    IF omitted > 0 THEN
        RAISE EXCEPTION 'coverage: % dataset(s) have fewer than half their reactions derived -- omitted from the derive pass', omitted;
    END IF;
END $$;

DO $$
DECLARE
    index_rows bigint := (SELECT count(*) FROM ord.reaction);
    payload_rows bigint := (SELECT count(*) FROM public.reactions);
    datasets bigint := (SELECT count(*) FROM ord.dataset);
    dataset_metadata bigint := (SELECT count(*) FROM public.datasets);
    declared bigint := (SELECT coalesce(sum(num_reactions), 0) FROM public.datasets);
BEGIN
    -- An empty candidate satisfies every equality below (0 = 0); reject it outright so a
    -- cutover gate never green-lights an unloaded database.
    IF index_rows = 0 OR datasets = 0 THEN
        RAISE EXCEPTION 'parity: candidate is empty (% reactions, % datasets)',
            index_rows, datasets;
    END IF;
    IF index_rows <> payload_rows THEN
        RAISE EXCEPTION 'parity: % index rows vs % payload rows', index_rows, payload_rows;
    END IF;
    IF datasets <> dataset_metadata THEN
        RAISE EXCEPTION 'parity: % datasets vs % metadata rows', datasets, dataset_metadata;
    END IF;
    IF declared <> index_rows THEN
        RAISE EXCEPTION 'parity: datasets declare % reactions but % are indexed', declared, index_rows;
    END IF;
END $$;
\echo All invariant checks passed.
