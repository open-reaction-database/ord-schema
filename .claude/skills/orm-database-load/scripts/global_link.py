# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Links rdkit.mols ids into the derived SMILES tables in one global pass.

Equivalent to ``update_rdkit_tables`` + ``update_rdkit_ids`` over every dataset,
minus the per-dataset/per-shard scoping. That scoping makes the planner rescan the
whole unlinked partial index once per ``(dataset, shard)`` pair, because the shard
hash is a filter rather than an index condition. Backfilling every dataset at once
leaves millions of unlinked rows standing for the whole link phase, so a single scan
is dramatically cheaper: 372 s versus a >9 h projection, measured on 4.7M compounds
across 53 datasets.

Use this only for backfills. A normal end-to-end load drains the unlinked set
dataset by dataset and never hits the pathology. Idempotent: every statement is
guarded, so re-running is safe.

Connects via standard ``PG*`` environment variables. Run from a host whose Python
RDKit matches the one that wrote the database's SMILES.
"""

import argparse
import time

import psycopg

# Mirrors _update_rdkit_mols: dedupe candidate SMILES, skip structures already
# present, and keep the [Ti+5] carve-out (issue #672). mol_from_smiles returns NULL
# for unparseable input; the guard keeps NULL-mol rows out.
_INSERT_MOLS = """
WITH candidates AS (
    SELECT smiles FROM derived.compound_smiles
    WHERE rdkit_mol_id IS NULL AND smiles IS NOT NULL
    UNION
    SELECT smiles FROM derived.product_compound_smiles
    WHERE rdkit_mol_id IS NULL AND smiles IS NOT NULL
), new_smiles AS MATERIALIZED (
    SELECT smiles FROM candidates c
    WHERE c.smiles NOT LIKE '%[Ti+5]%'
      AND NOT EXISTS (SELECT 1 FROM rdkit.mols m WHERE m.smiles = c.smiles)
)
INSERT INTO rdkit.mols (smiles, mol, morgan_bfp, morgan_sfp)
SELECT smiles, mol, morganbv_fp(mol) AS morgan_bfp, morgan_fp(mol) AS morgan_sfp
FROM (SELECT smiles, mol_from_smiles(smiles::cstring) AS mol FROM new_smiles) computed
WHERE mol IS NOT NULL
ON CONFLICT (smiles) DO NOTHING
"""

# Chunked so each statement bounds its WAL and lock footprint. Skip the permanently
# unlinkable [Ti+5] rows (issue #672): they never match rdkit.mols, so leaving them in
# the window could yield a zero-row UPDATE and end the loop while linkable rows remain.
_LINK_CHUNK = """
UPDATE {table} d
SET rdkit_mol_id = m.id
FROM rdkit.mols m
WHERE m.smiles = d.smiles
  AND d.{id_column} IN (
      SELECT {id_column} FROM {table}
      WHERE rdkit_mol_id IS NULL AND smiles NOT LIKE '%[Ti+5]%'
      LIMIT {chunk}
  )
"""

_DERIVED_TABLES = (
    ("derived.compound_smiles", "compound_id"),
    ("derived.product_compound_smiles", "product_compound_id"),
)


def _positive_int(value: str) -> int:
    """Parses a strictly positive integer for argparse.

    A non-positive chunk renders ``LIMIT 0``, so every statement links nothing yet
    the run reports each table done -- a silent no-op that looks like a completed
    recovery.

    Args:
        value: The raw command-line string.

    Returns:
        The parsed positive integer.

    Raises:
        argparse.ArgumentTypeError: If ``value`` is not a positive integer.
    """
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments.

    Args:
        argv: Command-line arguments; defaults to ``sys.argv`` when None.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Links rdkit.mols ids into the derived SMILES tables in one global pass."
        )
    )
    parser.add_argument(
        "--chunk",
        type=_positive_int,
        default=500_000,
        help="Rows to link per statement; bounds WAL and lock footprint",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Inserts any missing structures, links every unlinked derived row, then ANALYZEs.

    Args:
        args: Parsed command-line arguments.
    """
    with psycopg.connect(autocommit=True) as connection:
        # _INSERT_MOLS and the ANALYZEs are single statements over millions of rows; a
        # statement timeout would only fire partway through work that is safe to finish.
        connection.execute("SET statement_timeout = 0")
        start = time.time()

        cursor = connection.execute(_INSERT_MOLS)
        print(
            f"[{time.time() - start:7.1f}s] rdkit.mols inserted: {cursor.rowcount}",
            flush=True,
        )

        for table, id_column in _DERIVED_TABLES:
            linked = 0
            while True:
                cursor = connection.execute(
                    _LINK_CHUNK.format(
                        table=table, id_column=id_column, chunk=args.chunk
                    )
                )
                if not cursor.rowcount:
                    break
                linked += cursor.rowcount
                print(
                    f"[{time.time() - start:7.1f}s] {table}: linked {linked}",
                    flush=True,
                )
            print(
                f"[{time.time() - start:7.1f}s] {table}: done, {linked} linked",
                flush=True,
            )

        # The link flips millions of rows from NULL to non-NULL; the planner needs to
        # know before anything else queries the unlinked partial index.
        for table in (*(table for table, _ in _DERIVED_TABLES), "rdkit.mols"):
            connection.execute(f"ANALYZE {table}")
            print(f"[{time.time() - start:7.1f}s] analyzed {table}", flush=True)


if __name__ == "__main__":
    main(parse_args())
