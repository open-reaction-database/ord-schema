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
-- skipped a compound -- or an entire attachment path (workup input, product measurement,
-- reaction product) -- leaves it with a structural identifier but no derived row. Covers
-- both ord.compound and ord.product_compound. NAME-only compounds are excluded, so the
-- documented `compounds > derived` gap does not trip this. The type list mirrors
-- message_helpers.STRUCTURAL_IDENTIFIER_TYPES.
--
-- This does not raise on the first missing row. Even on a load whose RDKit matches the
-- writer, a small residue of structural identifiers is unparseable/unreconstructable by
-- RDKit (organometallics such as Ir photocatalysts or C[Al](C)(C)([Li])..., charged-N
-- ring systems) and legitimately derives no row -- on the order of 0.02% of the full ORD
-- corpus. Raising on any miss would fail every real load. Instead it raises only when the
-- miss rate exceeds `tolerance_pct` of derived rows, which still catches a systematic gap
-- (a skipped dataset or attachment path, orders of magnitude larger) while tolerating that
-- residue. The count is always printed; lower the tolerance to tighten the gate.
DO $$
DECLARE
    missing bigint;
    derived_total bigint;
    tolerance_pct numeric := 0.1;
BEGIN
    SELECT count(*) INTO missing
    FROM (
        SELECT c.id
        FROM ord.compound c
        WHERE EXISTS (
                SELECT 1 FROM ord.compound_identifier ci
                WHERE ci.compound_id = c.id
                  AND ci.type IN ('SMILES', 'INCHI', 'MOLBLOCK')
                  AND ci.value <> ''
            )
          AND NOT EXISTS (
                SELECT 1 FROM derived.compound_smiles d WHERE d.compound_id = c.id
            )
        UNION ALL
        SELECT pc.id
        FROM ord.product_compound pc
        WHERE EXISTS (
                SELECT 1 FROM ord.compound_identifier ci
                WHERE ci.product_compound_id = pc.id
                  AND ci.type IN ('SMILES', 'INCHI', 'MOLBLOCK')
                  AND ci.value <> ''
            )
          AND NOT EXISTS (
                SELECT 1 FROM derived.product_compound_smiles d
                WHERE d.product_compound_id = pc.id
            )
    ) t;
    SELECT (SELECT count(*) FROM derived.compound_smiles)
         + (SELECT count(*) FROM derived.product_compound_smiles)
      INTO derived_total;
    RAISE INFO 'coverage: % structural-identifier compound(s) of % derived rows lack a derived row (tolerance % percent)', missing, derived_total, tolerance_pct;
    IF missing::numeric > tolerance_pct / 100 * derived_total THEN
        RAISE EXCEPTION 'coverage: % structural-identifier compounds lack a derived row, above the tolerance of % percent of % derived rows -- likely a systematic miss, not the RDKit-unparseable residual', missing, tolerance_pct, derived_total;
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
