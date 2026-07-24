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
--   2. PARTIAL derivation (completeness). SMILES derivation is sharded (loading.py: one hash
--      partition of a >50k-reaction dataset's compound ids per worker, capped at 32), and the
--      loader raises on a failed shard -- but a defence-in-depth gate should not trust that exit
--      status. A failed shard leaves a cell with derived>0 yet ~1/num_shards of it underived:
--      tens of thousands of compounds at >= ~3% of the cell (num_shards <= 32). That dwarfs the
--      residue (a few thousand corpus-wide), so a cell is flagged when its underived count clears
--      BOTH an absolute floor (max_cell_missing, above any per-cell residue) AND a fraction
--      (min_gap_pct, below the smallest 1/32 shard hole). Requiring both keeps a small
--      organometallic cell (high fraction, tiny absolute) and a huge cell holding a little
--      residue (large absolute, tiny fraction) from tripping, while a shard hole clears both.
--
-- The total underived count is printed. Defaults: min_scope_size 50, max_cell_missing 10000,
-- min_gap_pct 1 (override with -v NAME=N).
\if :{?min_scope_size}
\else
  \set min_scope_size 50
\endif
\if :{?max_cell_missing}
\else
  \set max_cell_missing 10000
\endif
\if :{?min_gap_pct}
\else
  \set min_gap_pct 1
\endif
SELECT set_config('ord.verify_min_scope_size', :'min_scope_size', false) AS min_scope_size,
       set_config('ord.verify_max_cell_missing', :'max_cell_missing', false) AS max_cell_missing,
       set_config('ord.verify_min_gap_pct', :'min_gap_pct', false) AS min_gap_pct;
DO $$
DECLARE
    total_missing bigint;
    bad_omitted bigint;
    example_omitted text;
    bad_partial bigint;
    example_partial text;
    min_scope_size int := current_setting('ord.verify_min_scope_size')::int;
    max_cell_missing int := current_setting('ord.verify_max_cell_missing')::int;
    min_gap_pct numeric := current_setting('ord.verify_min_gap_pct')::numeric;
BEGIN
    WITH scoped AS (
        SELECT r.dataset_id AS dataset_id, 'reaction_input' AS scope,
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
    )
    SELECT
        coalesce(sum(structural - derived), 0),
        count(*) FILTER (WHERE structural >= min_scope_size AND derived = 0),
        min(scope || ' in dataset ' || dataset_id)
            FILTER (WHERE structural >= min_scope_size AND derived = 0),
        count(*) FILTER (
            WHERE derived > 0
              AND structural - derived >= max_cell_missing
              AND (structural - derived)::numeric >= min_gap_pct / 100 * structural
        ),
        min(scope || ' in dataset ' || dataset_id
            || ' (' || (structural - derived) || ' of ' || structural || ' missing)')
            FILTER (
                WHERE derived > 0
                  AND structural - derived >= max_cell_missing
                  AND (structural - derived)::numeric >= min_gap_pct / 100 * structural
            )
      INTO total_missing, bad_omitted, example_omitted, bad_partial, example_partial
    FROM per_cell;
    RAISE INFO 'coverage: % structural-identifier compound(s) lack a derived row (per (dataset, scope) gates: existence min scope size %, completeness max % missing and >= % pct per cell)', total_missing, min_scope_size, max_cell_missing, min_gap_pct;
    IF bad_omitted > 0 THEN
        RAISE EXCEPTION 'coverage: % (dataset, attachment scope) cell(s) hold >= % structural-identifier compounds but no derived row (e.g. %) -- a dataset or attachment path omitted from the derive pass', bad_omitted, min_scope_size, example_omitted;
    END IF;
    IF bad_partial > 0 THEN
        RAISE EXCEPTION 'coverage: % (dataset, attachment scope) cell(s) are partially derived beyond the unparseable residue (e.g. %) -- likely a failed derivation shard, which underives ~1/num_shards of a cell', bad_partial, example_partial;
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
