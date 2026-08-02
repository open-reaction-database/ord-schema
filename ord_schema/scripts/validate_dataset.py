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
import warnings
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm

from ord_schema import message_helpers, parquet, validations
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


def parse_shard(value: str | None) -> tuple[int, int] | None:
    """Parses a ``--shard`` value of the form ``I/N``.

    Args:
        value: Shard specification, e.g. ``0/4``, or None for no sharding.

    Returns:
        ``(index, count)``, or None when ``value`` is None.

    Raises:
        ValueError: If ``value`` is malformed or the index is out of range.
    """
    if value is None:
        return None
    try:
        index_text, count_text = value.split("/", 1)
        index, count = int(index_text), int(count_text)
    except ValueError:
        msg = f"--shard must look like I/N; got {value!r}"
        raise ValueError(msg) from None
    if count < 1 or not 0 <= index < count:
        msg = f"--shard needs 0 <= I < N; got {value!r}"
        raise ValueError(msg)
    return index, count


def filter_filenames(filenames: Iterable[str], pattern: str) -> list[str]:
    """Filters filenames according to a regex pattern."""
    return [filename for filename in filenames if re.search(pattern, filename)]


def _validate_pb(filename: str) -> list[str]:
    """Validates a single pb dataset; returns formatted error lines."""
    silence_rdkit_logs()
    dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
    try:
        validations.validate_datasets({filename: dataset})
    except validations.ValidationError as error:
        return [str(error)]
    return []


def _validate_row_group(
    filename: str, row_group: int
) -> tuple[list[str], validations.DatasetCrossRefState]:
    """Validates one row group: per-reaction checks + cross-ref observations."""
    silence_rdkit_logs()
    errors: list[str] = []
    state = validations.DatasetCrossRefState()
    for i, (_, reaction) in enumerate(
        parquet.iter_reactions(filename, row_group=row_group)
    ):
        output = validations.validate_message(reaction, raise_on_error=False)
        errors.extend(
            f"row_group {row_group}, reaction {i}: {error}" for error in output.errors
        )
        state.observe(reaction)
    return errors, state


def _cross_ref_state(filename: str) -> validations.DatasetCrossRefState:
    """Builds a file's complete cross-reference state, decoding but not validating.

    Cross-reference checks need every reaction in the file at once, which a shard
    holding a subset of row groups cannot supply. Reading them is cheap next to
    validating them -- on ord-data's 1.77M-reaction uspto parquet this pass takes
    about 20 seconds against the ~52 CPU-minutes of validate_message -- so one
    shard runs it over the whole file rather than the checks being skipped.

    Args:
        filename: Parquet dataset to scan.

    Returns:
        Aggregated cross-reference observations for every reaction in the file.
    """
    silence_rdkit_logs()
    state = validations.DatasetCrossRefState()
    for _, reaction in parquet.iter_reactions(filename):
        state.observe(reaction)
    return state


def _finalize_parquet(
    footer: parquet.ParquetFooter, state: validations.DatasetCrossRefState
) -> list[str]:
    """Runs dataset-level checks on aggregated parquet state; returns errors."""
    with warnings.catch_warnings(record=True) as tape:
        warnings.simplefilter("always", validations.ValidationError)
        validations.validate_dataset_streaming(
            name=footer.dataset.name,
            description=footer.dataset.description,
            dataset_id=footer.dataset.dataset_id,
            reaction_ids=[],
            has_reactions=footer.num_rows > 0,
            state=state,
        )
    return [
        f"Dataset: {w.message}"
        for w in tape
        if issubclass(w.category, validations.ValidationError)
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
        "--shard",
        default=None,
        help=(
            "Validate part of the work, as I/N (e.g. 0/4). Tasks -- one per "
            "parquet row group, one per pb file -- are dealt out round-robin, so "
            "N shards cover everything between them. Dataset-level "
            "cross-reference checks are skipped, since they need every row group "
            "of a file at once."
        ),
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Validates matching Dataset protos in parallel and reports any failures."""
    filenames = sorted(glob.glob(args.input_pattern, recursive=True))
    logger.info("Found %d datasets", len(filenames))
    if args.filter:
        filenames = filter_filenames(filenames, args.filter)
        logger.info("Filtered to %d datasets", len(filenames))

    shard = parse_shard(args.shard)
    if shard is not None:
        logger.info(
            "Shard %d/%d of the per-reaction checks%s",
            shard[0],
            shard[1],
            "; also running the dataset-level cross-reference pass"
            if shard[0] == 0
            else "",
        )
    # Dealt out round-robin over a counter that spans files, so shards stay even
    # whether the input is one big dataset or many small ones.
    task_index = 0

    def mine() -> bool:
        """Returns whether the next task is this shard's, advancing the counter."""
        nonlocal task_index
        index = task_index
        task_index += 1
        return shard is None or index % shard[1] == shard[0]

    parquet_entries: dict[str, _ParquetEntry] = {}
    failures: list[str] = []

    with ProcessPoolExecutor(args.n_jobs) as executor:
        futures: dict = {}
        for filename in filenames:
            if filename.endswith(".parquet"):
                footer = parquet.load_footer(filename)
                entry = _ParquetEntry(
                    remaining=footer.num_row_groups,
                    footer=footer,
                    state=validations.DatasetCrossRefState(),
                )
                parquet_entries[filename] = entry
                # Dataset-level checks need every reaction in the file at once.
                # Unsharded they fall out of merging the row-group states below;
                # sharded, no one shard holds them all, so shard 0 makes a single
                # cheap pass over the whole file instead.
                if shard is not None and shard[0] == 0:
                    future = executor.submit(_cross_ref_state, filename)
                    futures[future] = ("crossref", filename, None)
                if footer.num_row_groups == 0:
                    if shard is None:
                        failures.extend(
                            f"{filename}: {e}"
                            for e in _finalize_parquet(footer, entry.state)
                        )
                    continue
                submitted = 0
                for row_group in range(footer.num_row_groups):
                    if not mine():
                        continue
                    future = executor.submit(_validate_row_group, filename, row_group)
                    futures[future] = ("parquet", filename, row_group)
                    submitted += 1
                # Only the row groups this shard runs can count down towards
                # finalization; the rest belong to other shards.
                entry.remaining = submitted
            elif mine():
                future = executor.submit(_validate_pb, filename)
                futures[future] = ("pb", filename, None)

        total_tasks = len(futures)
        for future in tqdm(as_completed(futures), total=total_tasks):
            kind, filename, _ = futures[future]
            if kind == "pb":
                failures.extend(f"{filename}: {error}" for error in future.result())
            elif kind == "crossref":
                failures.extend(
                    f"{filename}: {error}"
                    for error in _finalize_parquet(
                        parquet_entries[filename].footer, future.result()
                    )
                )
            else:
                errors, state = future.result()
                failures.extend(f"{filename}: {error}" for error in errors)
                entry = parquet_entries[filename]
                entry.state.merge(state)
                entry.remaining -= 1
                if entry.remaining == 0 and shard is None:
                    failures.extend(
                        f"{filename}: {error}"
                        for error in _finalize_parquet(entry.footer, entry.state)
                    )

    if failures:
        text = "\n".join(failures)
        raise validations.ValidationError(f"Dataset(s) failed validation:\n{text}")


if __name__ == "__main__":
    main(parse_args())
