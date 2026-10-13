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

r"""Publishes a derived artifact tree to the ord-artifacts Hugging Face repository.

The tree holds one directory per artifact, as the derive scripts write them::

    upload_artifacts.py --root=artifacts \
        --message="Derive from ord-data <commit>" --dry_run

    artifacts/projections/aa/ord_dataset-<id>.parquet
        -> projections/projection-1/rdkit-<version>/aa/ord_dataset-<id>.parquet
    artifacts/pivots/inputs.components/aa/ord_dataset-<id>.parquet
        -> pivots/pivot-1_projection-1/rdkit-<version>/inputs.components/aa/...

Each file goes under the build its own footer names, everything lands in one commit,
and files a build no longer holds are deleted from it. A tree that spans more than one
build, or whose artifacts do not all restate the same source datasets, is refused.
Authenticates the way ``huggingface_hub`` does: ``hf auth login`` or ``$HF_TOKEN``.
"""

import argparse

from ord_schema import huggingface
from ord_schema.logging import get_logger

logger = get_logger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--root", required=True, help="Directory holding one directory per artifact"
    )
    parser.add_argument("--message", required=True, help="Commit message")
    parser.add_argument("--revision", default="main", help="Branch to commit to")
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Check the tree and report the commit without making it",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Checks the tree and publishes it."""
    args = parse_args(argv)
    upload = huggingface.upload_artifacts(
        args.root,
        commit_message=args.message,
        revision=args.revision,
        dry_run=args.dry_run,
    )
    for artifact, prefix in upload.prefixes.items():
        logger.info("%s -> %s", artifact, prefix)
    verb = "would be" if args.dry_run else "were"
    logger.info(
        "%d files %s uploaded (unchanged ones are dropped from the commit), %d deleted",
        len(upload.uploaded),
        verb,
        len(upload.deleted),
    )
    for path in upload.deleted:
        logger.info("deleting %s", path)
    if upload.commit_url is not None:
        logger.info("committed %s", upload.commit_url)


if __name__ == "__main__":
    main()
