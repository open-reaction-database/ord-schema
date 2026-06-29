# Copyright 2022 Open Reaction Database Project Authors
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

"""CLI for loading ORD datasets into the ORM database.

Thin wrapper over ``ord_schema.orm.loading.load_datasets``; see that module for the
ingest/derivation staging the ``--stages`` flag selects.
"""

import argparse
import logging
import os

from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.orm import database, loading

logger = get_logger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Add datasets to the ORM database")
    parser.add_argument(
        "--pattern", required=True, help="Pattern for dataset filenames"
    )
    parser.add_argument(
        "--stages",
        default=",".join(loading.STAGES),
        help="Comma-separated stages to run: 'ingest' (ord.*/public.*) and/or 'derived' "
        "(derived.* SMILES, RDKit links, reaction classes). Derived-only runs over "
        "already-ingested datasets matching --pattern.",
    )
    parser.add_argument(
        "--classify_reactions",
        action="store_true",
        help="Assign reaction class/name labels in the derived stage "
        "(requires the 'reaction-class' extra)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Update changed datasets"
    )
    parser.add_argument("--dsn", default=None, help="Postgres connection string")
    parser.add_argument("--database", default="orm", help="Database")
    parser.add_argument("--username", default="postgres", help="Database username")
    parser.add_argument("--password", default=None, help="Database password")
    parser.add_argument("--host", default="localhost", help="Database host")
    parser.add_argument("--port", type=int, default=5432, help="Database port")
    parser.add_argument(
        "--n_jobs", type=int, default=1, help="Number of parallel workers"
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Runs the selected ingest and/or derived stages over the matching datasets."""
    silence_rdkit_logs()
    if args.debug:
        get_logger(database.__name__, level=logging.DEBUG)
    stages = [stage.strip() for stage in args.stages.split(",") if stage.strip()]
    unknown = set(stages) - set(loading.STAGES)
    if unknown:
        raise SystemExit(
            f"Unknown --stages values {sorted(unknown)}; choose from {list(loading.STAGES)}"
        )
    if not stages:
        raise SystemExit(f"--stages must select at least one of {list(loading.STAGES)}")
    if args.classify_reactions and database.update_reaction_classes is None:
        raise SystemExit(
            "--classify_reactions requires the 'reaction-class' extra: "
            "pip install ord-schema[reaction-class]"
        )
    if args.dsn:
        dsn = args.dsn
    else:
        dsn = database.get_connection_string(
            database=args.database,
            username=args.username,
            password=args.password or os.environ["PGPASSWORD"],
            host=args.host,
            port=args.port,
        )
    loading.load_datasets(
        args.pattern,
        dsn,
        stages=stages,
        overwrite=args.overwrite,
        classify_reactions=args.classify_reactions,
        n_jobs=args.n_jobs,
    )


if __name__ == "__main__":
    main(parse_args())
