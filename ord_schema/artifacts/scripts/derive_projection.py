# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Derives queryable projections from Parquet datasets.

Each input dataset yields one projection carrying every field of every message reachable
from Reaction; see ``ord_schema.artifacts.projection`` for what it normalizes and why.
Outputs mirror the inputs' directory layout beneath ``--output_dir``, rooted at the
leading components of ``--input_pattern`` that hold no wildcard::

    derive_projection.py --input_pattern='data/*/*.parquet' --output_dir=projections

    data/aa/ord_dataset-<id>.parquet  ->  projections/aa/ord_dataset-<id>.parquet

Projections whose footer already records the current source content, library version,
and artifact version are skipped, so re-running is cheap; ``--force`` rewrites them
anyway. A match that is itself a derived artifact is ignored rather than derived from,
so ``--output_dir`` may sit inside a recursive pattern's reach; an ``--output_dir`` that
would write over any source dataset is an error, as is a run that derives nothing.
"""

import argparse

from ord_schema.artifacts import base, projection
from ord_schema.logging import silence_rdkit_logs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_pattern", required=True, help="Input pattern for Parquet datasets"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to write projections beneath"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite projections that are already current",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Derives a projection for every dataset matching the input pattern.

    Args:
        args: Parsed arguments, read for ``input_pattern``, ``output_dir``, and
            ``force``.

    Raises:
        ValueError: If the pattern matched only derived artifacts, which means it was
            aimed at an output tree rather than a source tree. Silence there would let a
            pipeline step downstream proceed as though the projections had been built.
    """
    silence_rdkit_logs()
    written, skipped, ignored = base.derive_tree(
        args.input_pattern,
        args.output_dir,
        artifact=projection.ARTIFACT,
        write=projection.write_projection,
        schema=projection.SCHEMA,
        force=args.force,
    )
    if not written and not skipped:
        raise ValueError(
            f"no sources derived: all {ignored} matches for {args.input_pattern!r} are "
            "already derived artifacts"
        )


if __name__ == "__main__":
    main(parse_args())
