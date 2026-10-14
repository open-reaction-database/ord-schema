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

"""Checks that a derived artifact tree covers exactly the datasets ord-data holds.

``upload_artifacts`` checks that a tree is one consistent build, but a build of only
some datasets is consistent too. This compares the tree's projections with the
ord-data checkout they were derived from, and that checkout with the ord-data Hugging
Face mirror:

* every Parquet dataset in the checkout has exactly one projection, stamped with that
  dataset's current content hash, and no projection restates anything else;
* the checkout holds the same dataset IDs as the mirror at ``--revision``, so a stale
  checkout cannot pass for the whole corpus.

``upload_artifacts`` carries the projections' dataset set to every other artifact.
Exits non-zero on any mismatch. Hashing the checkout takes under a minute over ORD.
"""

import argparse
import collections
import pathlib
import re
import sys

from huggingface_hub import list_repo_files

from ord_schema import huggingface, parquet
from ord_schema.artifacts import base, projection

_MIRROR_DATASET = re.compile(r"^data/\w{2}/(ord_dataset-\w+)\.(parquet|pb\.gz)$")


def _checkout(ord_data: pathlib.Path) -> dict[str, str]:
    """Returns each Parquet dataset's content hash in the checkout, keyed by ID."""
    hashes = {}
    for path in sorted(ord_data.glob("data/*/*.parquet")):
        view = parquet.DatasetView(path)
        if not view.dataset_id or view.dataset_id in hashes:
            raise ValueError(f"{path} has a missing or repeated dataset ID")
        hashes[view.dataset_id] = view.md5()
    return hashes


def _mirror(revision: str) -> set[str]:
    """Returns the dataset IDs on the ord-data mirror at ``revision``."""
    files = list_repo_files(
        huggingface.ORD_DATA_HF_REPO, repo_type="dataset", revision=revision
    )
    return {match[1] for name in files if (match := _MIRROR_DATASET.match(name))}


def _projections(root: pathlib.Path) -> dict[str, list[str]]:
    """Returns the projections under ``root``, keyed by the source hash they restate."""
    held = collections.defaultdict(list)
    for path in huggingface.parquet_files(root / "projections"):
        stamps = base.load_stamps(path)
        if stamps.artifact != projection.ARTIFACT:
            raise ValueError(f"{path} holds {stamps.artifact}, not a projection")
        held[stamps.source_md5].append(f"{stamps.source_dataset_id} ({path})")
    return held


def check(root: pathlib.Path, ord_data: pathlib.Path, revision: str) -> list[str]:
    """Returns every way the tree and the checkout fail to cover the same datasets."""
    checkout = _checkout(ord_data)
    mirror = _mirror(revision)
    held = _projections(root)
    problems = [
        f"on the mirror but not in the checkout: {dataset_id}"
        for dataset_id in sorted(mirror - set(checkout))
    ]
    problems.extend(
        f"in the checkout but not on the mirror: {dataset_id}"
        for dataset_id in sorted(set(checkout) - mirror)
    )
    expected = {md5: dataset_id for dataset_id, md5 in checkout.items()}
    problems.extend(
        f"no current projection of {dataset_id}"
        for dataset_id, md5 in sorted(checkout.items())
        if md5 not in held
    )
    for md5, projections in sorted(held.items()):
        if md5 not in expected:
            problems.extend(
                f"projection restates no dataset in the checkout: {name}"
                for name in projections
            )
        elif len(projections) > 1:
            problems.append(
                f"{len(projections)} projections restate one dataset: "
                + ", ".join(projections)
            )
    print(
        f"{len(checkout)} datasets in the checkout, {len(mirror)} on the mirror at "
        f"{revision}, {sum(map(len, held.values()))} projections"
    )
    return problems


def main() -> None:
    """Checks the tree against the checkout and the mirror."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, help="The local artifact tree")
    parser.add_argument("--ord_data", required=True, help="The ord-data checkout")
    parser.add_argument("--revision", default="main", help="Mirror revision to match")
    args = parser.parse_args()
    problems = check(
        pathlib.Path(args.root).expanduser(),
        pathlib.Path(args.ord_data).expanduser(),
        args.revision,
    )
    for problem in problems:
        print(f"MISMATCH {problem}")
    if problems:
        sys.exit(1)
    print("the tree covers exactly the checkout's datasets")


if __name__ == "__main__":
    main()
