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

r"""Derives structure-search artifacts from projections.

Each projection yields one structures artifact: one row per distinct structure, with the
fingerprints a query engine screens against and the serialized molecule verification
reads; see ``ord_schema.artifacts.structures`` for the columns and why. Outputs mirror
the inputs' directory layout beneath ``--output_dir``::

    derive_structures.py --input_pattern='projections/*/*.parquet' \
        --output_dir=structures

    projections/aa/ord_dataset-<id>.parquet -> structures/aa/ord_dataset-<id>.parquet

Artifacts whose footer already records the current source content, library version, and
artifact version are skipped, so re-running is cheap; --force rewrites them anyway.

A match that is not a projection -- a source dataset, or a structures artifact from an
earlier run -- is ignored rather than derived from, so --output_dir may sit inside a
recursive pattern's reach. These are errors: an --output_dir that would write over any
input, a match that cannot be read as Parquet at all, and a run that finds no
projections at all.
"""

import argparse

from ord_schema.artifacts import base, projection, structures
from ord_schema.logging import silence_rdkit_logs


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input_pattern", required=True, help="Input pattern for Parquet projections"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to write artifacts beneath"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite artifacts that are already current",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Derives a structures artifact for every projection matching the input pattern.

    Args:
        args: Parsed arguments, read for ``input_pattern``, ``output_dir``, and
            ``force``.

    Raises:
        ValueError: If the pattern matched no projections, which usually means it was
            aimed at the source tree, or at an output tree, rather than at the
            projections. Silence there would let a pipeline step downstream proceed as
            though the artifacts had been built.
    """
    silence_rdkit_logs()
    written, skipped, ignored = base.derive_tree(
        args.input_pattern,
        args.output_dir,
        artifact=structures.ARTIFACT,
        write=structures.write_structures,
        force=args.force,
        parent_artifact=projection.ARTIFACT,
    )
    if not written and not skipped:
        raise ValueError(
            f"no structures derived: none of the {ignored} matches for "
            f"{args.input_pattern!r} are projections"
        )


if __name__ == "__main__":
    main(parse_args())
