# Copyright 2020 Open Reaction Database Project Authors
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
"""Validates a set of Dataset protocol buffers.

pb inputs are validated as a single per-file task. Parquet inputs are fanned out one
task per row group so ``--n_jobs`` actually saturates on a single large dataset; per-
file dataset-level checks (cross-references, scalar fields) are then run once after the
per-row-group results merge.
"""

import argparse
import dataclasses
import glob
import re
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from tqdm import tqdm

from ord_schema import message_helpers, parquet, validations
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)

# Header of an unsmudged Git LFS pointer. Datasets in ord-data are LFS objects,
# so a checkout that skipped `git lfs pull` leaves these small text files in
# place of content, and each format's parser then fails on them in its own
# confusing way (gzip reports a bad magic number, parquet a missing footer).
_LFS_POINTER_HEADER = b"version https://git-lfs.github.com/spec/"


def filter_filenames(filenames: Iterable[str], pattern: str) -> list[str]:
    """Filters filenames according to a regex pattern."""
    return [filename for filename in filenames if re.search(pattern, filename)]


def _is_lfs_pointer(filename: str) -> bool:
    """Tests whether a file holds a Git LFS pointer rather than dataset content."""
    try:
        with Path(filename).open("rb") as f:
            return f.read(len(_LFS_POINTER_HEADER)) == _LFS_POINTER_HEADER
    except OSError:
        return False


def _describe_failure(filename: str, error: BaseException) -> str:
    """Formats a file that could not be read or parsed at all."""
    if _is_lfs_pointer(filename):
        return (
            f"{filename}: file is a Git LFS pointer, not dataset content; "
            "run `git lfs pull` to fetch it"
        )
    return f"{filename}: {type(error).__name__}: {error}"


def _validate_pb(filename: str, options: validations.ValidationOptions) -> list[str]:
    """Validates a single pb dataset; returns formatted error lines."""
    silence_rdkit_logs()
    dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
    try:
        validations.validate_datasets({filename: dataset}, options=options)
    except validations.ValidationError as error:
        return [str(error)]
    return []


def _validate_row_group(
    filename: str, row_group: int, options: validations.ValidationOptions
) -> tuple[list[str], validations.DatasetCrossRefState]:
    """Validates one row group: per-reaction checks + cross-ref observations."""
    silence_rdkit_logs()
    errors: list[str] = []
    state = validations.DatasetCrossRefState()
    for i, (_, reaction) in enumerate(
        parquet.iter_reactions(filename, row_group=row_group)
    ):
        output = validations.validate_message(
            reaction, raise_on_error=False, options=options
        )
        errors.extend(
            f"row_group {row_group}, reaction {i}: {error}" for error in output.errors
        )
        state.observe(reaction)
    return errors, state


def _finalize_parquet(
    footer: parquet.ParquetFooter,
    state: validations.DatasetCrossRefState,
    options: validations.ValidationOptions,
) -> list[str]:
    """Runs dataset-level checks on aggregated parquet state; returns errors."""
    context = validations.ValidationContext(options=options)
    validations.validate_dataset_streaming(
        context=context,
        name=footer.dataset.name,
        description=footer.dataset.description,
        dataset_id=footer.dataset.dataset_id,
        reaction_ids=[],
        has_reactions=footer.num_rows > 0,
        state=state,
    )
    return [
        f"Dataset: {text}"
        for text, severity in context.findings
        if severity >= validations.Severity.ERROR
    ]


@dataclasses.dataclass
class _ParquetEntry:
    remaining: int
    footer: parquet.ParquetFooter
    state: validations.DatasetCrossRefState


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Validate Dataset protocol buffers")
    parser.add_argument(
        "--input_pattern",
        required=True,
        help="Input pattern for Dataset protos",
    )
    parser.add_argument("--filter", default=None, help="Regex filename filter")
    parser.add_argument(
        "--n_jobs", type=int, default=1, help="Number of parallel workers"
    )
    parser.add_argument(
        "--validate_ids",
        action="store_true",
        help="Require well-formed dataset and reaction ids. Off by default: ids are "
        "assigned when a dataset is submitted, so a draft does not have them yet. "
        "Corpus sweeps should pass this, matching what submission checks.",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Validates matching Dataset protos in parallel and reports any failures."""
    filenames = sorted(glob.glob(args.input_pattern, recursive=True))
    logger.info("Found %d datasets", len(filenames))
    if args.filter:
        filenames = filter_filenames(filenames, args.filter)
        logger.info("Filtered to %d datasets", len(filenames))

    options = validations.ValidationOptions(validate_ids=args.validate_ids)
    parquet_entries: dict[str, _ParquetEntry] = {}
    failures: list[str] = []

    with ProcessPoolExecutor(args.n_jobs) as executor:
        futures: dict = {}
        for filename in filenames:
            if filename.endswith(".parquet"):
                try:
                    footer = parquet.load_footer(filename)
                # A file this job cannot read at all is that file's failure, not
                # the run's: record it and keep validating the rest.
                except Exception as error:  # noqa: BLE001
                    failures.append(_describe_failure(filename, error))
                    continue
                entry = _ParquetEntry(
                    remaining=footer.num_row_groups,
                    footer=footer,
                    state=validations.DatasetCrossRefState(),
                )
                parquet_entries[filename] = entry
                if footer.num_row_groups == 0:
                    failures.extend(
                        f"{filename}: {e}"
                        for e in _finalize_parquet(footer, entry.state, options)
                    )
                    continue
                for row_group in range(footer.num_row_groups):
                    future = executor.submit(
                        _validate_row_group, filename, row_group, options
                    )
                    futures[future] = ("parquet", filename, row_group)
            else:
                future = executor.submit(_validate_pb, filename, options)
                futures[future] = ("pb", filename, None)

        total_tasks = len(futures)
        for future in tqdm(as_completed(futures), total=total_tasks):
            kind, filename, _ = futures[future]
            try:
                result = future.result()
            # As above: attribute an unreadable or unparseable file to itself.
            # A parquet row group that fails this way leaves `remaining` above
            # zero, which skips the dataset-level check for that file -- correct,
            # since the aggregated state is incomplete and the file already
            # counts as failed.
            except Exception as error:  # noqa: BLE001
                failures.append(_describe_failure(filename, error))
                continue
            if kind == "pb":
                failures.extend(f"{filename}: {error}" for error in result)
            else:
                errors, state = result
                failures.extend(f"{filename}: {error}" for error in errors)
                entry = parquet_entries[filename]
                entry.state.merge(state)
                entry.remaining -= 1
                if entry.remaining == 0:
                    failures.extend(
                        f"{filename}: {error}"
                        for error in _finalize_parquet(
                            entry.footer, entry.state, options
                        )
                    )

    if failures:
        text = "\n".join(failures)
        raise validations.ValidationError(f"Dataset(s) failed validation:\n{text}")


if __name__ == "__main__":
    main(parse_args())
