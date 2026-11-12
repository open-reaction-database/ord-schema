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

"""Writes the review report bundle for a dataset.

Usage:
    python -m ord_schema.scripts.review_dataset --input=dataset.pb.gz
        --output_dir=review/ [--profile=uniform|varied]

Writes ``<dataset>_review.html``, ``tables/`` and the AI review ``packet/`` into the
output directory; see
:mod:`ord_schema.review.report`. Findings never make this fail: it exits 0 once the
report is written, and non-zero only when the dataset cannot be read. Needs the
``review`` extra.
"""

import argparse
import pathlib
import sys

from ord_schema import datasets
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.review import packet, report, structure
from ord_schema.validations import Severity

logger = get_logger(__name__)

_LFS_POINTER_HEADER = b"version https://git-lfs.github.com/spec/"


def describe_read_failure(filename: str, error: BaseException) -> str:
    """Says why a dataset file could not be read.

    A Git LFS pointer, the usual cause in an ord-data checkout, is named as one, with
    what to do about it. Kept here rather than imported from validate_dataset, so the
    review tooling does not depend on another script's internals.
    """
    try:
        with pathlib.Path(filename).open("rb") as f:
            is_pointer = f.read(len(_LFS_POINTER_HEADER)) == _LFS_POINTER_HEADER
    except OSError:
        is_pointer = False
    if is_pointer:
        return (
            f"{filename}: file is a Git LFS pointer, not dataset content; "
            "run `git lfs pull` to fetch it"
        )
    return f"{filename}: {type(error).__name__}: {error}"


def main(argv: list[str] | None = None) -> int:
    """Writes the bundle; returns the process exit status."""
    parser = argparse.ArgumentParser(
        description="Writes the review report bundle for a dataset."
    )
    parser.add_argument("--input", required=True, help="Path to a Dataset file.")
    parser.add_argument(
        "--output_dir", required=True, help="Where to write the bundle."
    )
    parser.add_argument(
        "--profile",
        choices=[profile.value for profile in structure.ValidationProfile],
        default=structure.ValidationProfile.UNIFORM.value,
        help="The dataset's profile; uniform unless its reactions share no template.",
    )
    parser.add_argument(
        "--packet_budget_tokens",
        type=int,
        default=packet.DEFAULT_BUDGET_TOKENS,
        help="The most estimated tokens each file of the AI review packet may hold.",
    )
    args = parser.parse_args(argv)
    silence_rdkit_logs()
    try:
        dataset = datasets.load_dataset(args.input)
    except Exception as error:  # noqa: BLE001 -- any read failure is reported the same way
        # Without the traceback: the description says what to do.
        logger.error(describe_read_failure(args.input, error))  # noqa: TRY400
        return 1
    review = report.build_report(dataset, structure.ValidationProfile(args.profile))
    page = report.write_bundle(
        review, args.output_dir, packet_budget_tokens=args.packet_budget_tokens
    )
    counts = review.severity_counts()
    logger.info(
        f"{review.reaction_count:,} reactions: {counts[Severity.ERROR]} errors, "
        f"{counts[Severity.WARNING]} warnings, {counts[Severity.INFO]} "
        f"notes; wrote {page}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
