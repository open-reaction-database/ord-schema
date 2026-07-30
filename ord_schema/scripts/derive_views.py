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

"""Derives tabular views from Parquet datasets.

Each input dataset yields one view (see ``ord_schema.views`` for what a view contains,
and what it deliberately does not). Outputs mirror the inputs' directory layout beneath
``--output_dir``, taking the part of ``--input_pattern`` before its first wildcard as
the root::

    derive_views.py --input_pattern='data/*/*.parquet' --output_dir=views

    data/aa/ord_dataset-<id>.parquet  ->  views/aa/ord_dataset-<id>.parquet

Views whose footer already records the current source content, library version, and view
definition are skipped, so re-running is cheap; ``--force`` rewrites them anyway.
"""

import argparse
import glob
import pathlib

from ord_schema import parquet, views
from ord_schema.logging import get_logger, silence_rdkit_logs

logger = get_logger(__name__)


def glob_root(pattern: str) -> pathlib.PurePath:
    """Returns the leading directories of a glob pattern that contain no wildcards.

    Output paths are built relative to this, so ``data/*/x.parquet`` under
    ``--output_dir=views`` lands at ``views/<subdir>/x.parquet``. A pattern naming a
    single file has no wildcard to stop at, so its own last component is the file
    rather than a directory and the root is the directory holding it.
    """
    parts = pathlib.PurePath(pattern).parts
    fixed = []
    for part in parts:
        if any(char in part for char in "*?["):
            break
        fixed.append(part)
    else:
        fixed = fixed[:-1]
    return pathlib.PurePath(*fixed)


def output_path(source: str, pattern: str, output_dir: str) -> pathlib.Path:
    """Maps a source path to its view path under ``output_dir``."""
    relative = pathlib.PurePath(source).relative_to(glob_root(pattern))
    return pathlib.Path(output_dir) / relative


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_pattern", required=True, help="Input pattern for Parquet datasets"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to write views beneath"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite views that are already current",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Derives a view for every dataset matching the input pattern.

    Raises:
        ValueError: If the input pattern matches nothing.
    """
    silence_rdkit_logs()
    sources = sorted(glob.glob(args.input_pattern, recursive=True))
    if not sources:
        raise ValueError(f"no datasets matched: {args.input_pattern}")
    logger.info("Found %d datasets", len(sources))
    written = skipped = 0
    for source in sources:
        destination = output_path(source, args.input_pattern, args.output_dir)
        source_md5, _ = parquet.streaming_md5(source)
        if not args.force and views.is_current(destination, source_md5):
            logger.info("%s is current; skipping", destination)
            skipped += 1
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        rows = views.write_view(source, destination, source_md5=source_md5)
        logger.info("Wrote %d rows to %s", rows, destination)
        written += 1
    logger.info("Derived %d views (%d already current)", written, skipped)


if __name__ == "__main__":
    main(parse_args())
