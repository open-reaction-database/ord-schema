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
"""Dataset processing script for database submissions.

This script is meant to be a one-stop shop for preparing submissions to the
Open Reaction Database.

By default, the script only validates the input Dataset messages. Validation
may introduce changes to the Reaction messages, such as the addition of SMILES
for compounds identified only by NAME. Users should frequently run these checks
as they are preparing a dataset for submission.

With the optional --update flag, the script also performs database-specific
updates (such as adding record IDs). These operations are meant to be run as
part of the submission process and not as part of the pre-submission validation
cycle.
"""

import argparse
import dataclasses
import glob
import gzip
import os
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Mapping

import github

from ord_schema import message_helpers, parquet_dataset, updates, validations
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


@dataclasses.dataclass(eq=True, frozen=True, order=True)
class FileStatus:
    """A filename and its status in Git."""

    filename: str
    status: str
    original_filename: str

    def __post_init__(self):
        if self.status[0] not in ["A", "D", "M", "R"]:
            raise ValueError(f"unsupported file status: {self.status}")


def _get_inputs(args) -> list[FileStatus]:
    """Gets a list of Dataset proto filenames to process.

    Returns:
        List of FileStatus objects.

    Raises:
        ValueError: If a git-diff status is not one of {'A', 'D', 'M', 'R'}.
    """
    if args.input_pattern:
        # Setting recursive=True allows recursive matching with '**'.
        filenames = glob.glob(args.input_pattern, recursive=True)
        return [FileStatus(filename, "A", "") for filename in filenames]
    if args.input_file:
        inputs = []
        with open(args.input_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) == 3:
                    status, original_filename, filename = fields
                    if not status.startswith("R"):
                        raise ValueError(f"malformed status line: {line.strip()}")
                else:
                    status, filename = fields
                    if not status.startswith(("A", "D", "M")):
                        raise ValueError(f"unsupported git-diff statue: {status}")
                    original_filename = ""
                inputs.append(FileStatus(filename, status, original_filename))
        return inputs
    raise ValueError("one of --input_pattern or --input_file is required")


def cleanup(filename: str, output_filename: str):
    """Removes and/or renames the input Dataset files.

    Args:
        filename: Original dataset filename.
        output_filename: Updated dataset filename.
    """
    if filename == output_filename:
        logger.info("editing an existing dataset; no cleanup needed")
        return  # Reuse the existing dataset ID.
    args = ["git", "mv", filename, output_filename]
    logger.info("Running command: %s", " ".join(args))
    subprocess.run(args, check=True)


def _get_reaction_ids(dataset: dataset_pb2.Dataset | parquet_dataset.DatasetView) -> set[str]:
    """Returns a set containing the reaction IDs in a Dataset.

    For ``DatasetView``, the ``reaction_id`` column is read directly from the
    Parquet file so we never decode Reaction blobs just to collect IDs.
    """
    if isinstance(dataset, parquet_dataset.DatasetView):
        return {rid for rid in parquet_dataset.iter_reaction_ids(dataset.path) if rid}
    return {reaction.reaction_id for reaction in dataset.reactions if reaction.reaction_id}


def _load_base_dataset(file_status: FileStatus, base: str) -> dataset_pb2.Dataset | parquet_dataset.DatasetView | None:
    """Loads a Dataset from another git branch.

    Parquet inputs are spilled to a temp file and wrapped in a ``DatasetView``
    so the diff path can scan the ``reaction_id`` column without decoding any
    Reaction blobs. The temp file outlives the view (process-lifetime leak),
    which is fine for this CLI script — the OS reclaims it on exit.
    """
    if file_status.status.startswith("A"):
        return None  # Dataset only exists in the submission.
    # NOTE(kearnes): Use --no-pager to avoid a non-zero exit code.
    git_args = ["git", "--no-pager", "show"]
    if file_status.status.startswith("R"):
        git_args.append(f"{base}:{file_status.original_filename}")
    else:
        git_args.append(f"{base}:{file_status.filename}")
    logger.info("Running command: %s", " ".join(git_args))
    serialized = subprocess.run(git_args, capture_output=True, check=True, text=False)
    if serialized.stdout.startswith(b"version"):
        # Convert Git LFS pointers to real data.
        serialized = subprocess.run(
            ["git", "lfs", "smudge"],
            input=serialized.stdout,
            capture_output=True,
            check=True,
            text=False,
        )
    if git_args[-1].endswith(".parquet"):
        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as temp:
            temp.write(serialized.stdout)
            temp_path = temp.name
        return parquet_dataset.DatasetView(temp_path)
    if git_args[-1].endswith(".gz"):
        value = gzip.decompress(serialized.stdout)
    else:
        value = serialized.stdout
    return dataset_pb2.Dataset.FromString(value)


def get_change_stats(
    datasets: Mapping[str, dataset_pb2.Dataset | parquet_dataset.DatasetView | None],
    inputs: Iterable[FileStatus],
    base: str,
) -> tuple[set[str], set[str], set[str]]:
    """Computes diff statistics for the submission.

    Args:
        datasets: Dict mapping filenames to Dataset messages. Values may be None only
            for deleted files; non-deleted inputs must have a Dataset for `filename`.
        inputs: List of FileStatus objects.
        base: Git branch to diff against.

    Returns:
        added: Set of added reaction IDs.
        removed: Set of deleted reaction IDs.
        changed: Set of changed reaction IDs.
    """
    old, new = set(), set()
    for file_status in inputs:
        if not file_status.status.startswith("D"):
            current = datasets[file_status.filename]
            if current is None:
                raise ValueError(f"missing dataset for non-deleted file: {file_status.filename}")
            new.update(_get_reaction_ids(current))
        dataset = _load_base_dataset(file_status, base)
        if dataset is not None:
            old.update(_get_reaction_ids(dataset))
    return new - old, old - new, new & old


def _run_updates(
    datasets: Mapping[str, dataset_pb2.Dataset | parquet_dataset.DatasetView],
    *,
    root: str,
    output_format: str,
    write_errors: bool,
    cleanup_files: bool,
) -> None:
    """Updates the submission files.

    When the input is a ``DatasetView`` and the output is also Parquet, the
    update runs as a streaming two-pass over the input file with an atomic
    temp-then-rename publish (validation runs against the temp before the
    rename). Otherwise the in-memory path mutates the Dataset in place,
    validates, and writes through ``message_helpers.write_dataset``.
    """
    options = validations.ValidationOptions(validate_ids=True, require_provenance=True)
    for input_filename, dataset in datasets.items():
        is_parquet_stream = isinstance(dataset, parquet_dataset.DatasetView) and output_format == ".parquet"
        if is_parquet_stream:
            # Resolve dataset_id up-front so the output filename is known
            # before we open the streaming writer.
            updates.assign_dataset_id(dataset)
            output_filename = os.path.join(
                root,
                message_helpers.id_filename(f"{dataset.dataset_id}{output_format}"),
            )
            os.makedirs(os.path.dirname(output_filename), exist_ok=True)
            if cleanup_files:
                cleanup(input_filename, output_filename)
            temp_filename = output_filename + ".tmp"
            try:
                updates.update_dataset_parquet(input_filename, temp_filename, dataset_id=dataset.dataset_id)
                validations.validate_datasets(
                    {input_filename: parquet_dataset.DatasetView(temp_filename)},
                    write_errors,
                    options=options,
                )
            except BaseException:
                if os.path.exists(temp_filename):
                    os.unlink(temp_filename)
                raise
            logger.info("writing Dataset to %s", output_filename)
            os.replace(temp_filename, output_filename)
            continue
        # In-memory path: materialize a Parquet input if the requested output
        # format is not Parquet (so we can mutate via update_dataset).
        if isinstance(dataset, parquet_dataset.DatasetView):
            dataset = parquet_dataset.read_dataset(input_filename)
        updates.update_dataset(dataset)
        validations.validate_datasets({input_filename: dataset}, write_errors, options=options)
        output_filename = os.path.join(
            root,
            message_helpers.id_filename(f"{dataset.dataset_id}{output_format}"),
        )
        os.makedirs(os.path.dirname(output_filename), exist_ok=True)
        if cleanup_files:
            cleanup(input_filename, output_filename)
        logger.info("writing Dataset to %s", output_filename)
        message_helpers.write_dataset(dataset, output_filename)


def run(args) -> tuple[set[str] | None, set[str] | None, set[str] | None]:
    """Main function that returns added/removed reaction ID sets.

    This function should be called directly by tests to get access to the
    return values. If main() returns something other than None it will break
    shell error code logic downstream.

    Returns:
        added: Set of added reaction IDs.
        removed: Set of deleted reaction IDs.
        changed: Set of changed reaction IDs.
    """
    inputs = sorted(_get_inputs(args))
    if not inputs:
        logger.info("nothing to do")
        return set(), set(), set()  # Nothing to do.
    # NOTE(kearnes): Process one dataset at a time to avoid OOM errors.
    change_stats = {}
    for file_status in inputs:
        dataset: dataset_pb2.Dataset | parquet_dataset.DatasetView | None
        if file_status.status == "D":
            dataset = None
        elif file_status.filename.endswith(".parquet"):
            dataset = parquet_dataset.DatasetView(file_status.filename)
            logger.info("%s: %d reactions", file_status.filename, len(dataset.reactions))
        else:
            dataset = message_helpers.load_message(file_status.filename, dataset_pb2.Dataset)
            logger.info("%s: %d reactions", file_status.filename, len(dataset.reactions))
        datasets: dict[str, dataset_pb2.Dataset | parquet_dataset.DatasetView | None] = {file_status.filename: dataset}
        datasets_checked: dict[str, dataset_pb2.Dataset | parquet_dataset.DatasetView] = (
            {file_status.filename: dataset} if dataset is not None else {}
        )
        if not args.no_validate and dataset is not None:
            # Note: this does not check if IDs are malformed.
            validations.validate_datasets(datasets_checked, args.write_errors)
            # Check reaction sizes.
            for reaction in dataset.reactions:
                reaction_size = sys.getsizeof(reaction.SerializeToString()) / 1e6
                if reaction_size > args.max_size:
                    raise ValueError(f"Reaction is larger than --max_size ({reaction_size} vs {args.max_size})")
        if args.base:
            added, removed, changed = get_change_stats(datasets, [file_status], base=args.base)
            change_stats[file_status.filename] = (added, removed, changed)
            logger.info(
                "Summary: +%d -%d Δ%d reaction IDs",
                len(added),
                len(removed),
                len(changed),
            )
        if args.update and dataset is not None:
            _run_updates(
                datasets_checked,
                root=args.root,
                output_format=args.output_format,
                write_errors=args.write_errors,
                cleanup_files=args.cleanup,
            )
    if change_stats:
        total_added, total_removed, total_changed = set(), set(), set()
        comment = [
            "Change summary:",
            "| Filename | Added | Removed | Changed |",
            "| -------- | ----- | ------- | ------- |",
        ]
        for filename, (added, removed, changed) in change_stats.items():
            comment.append(f"| {filename} | {len(added)} | {len(removed)} | {len(changed)} |")
            total_added |= added
            total_removed |= removed
            total_changed |= changed
        comment.append(f"| | **{len(total_added)}** | **{len(total_removed)}** | **{len(total_changed)}** |")
        if args.issue and args.token:
            client = github.Github(args.token)
            repo = client.get_repo(os.environ["GITHUB_REPOSITORY"])
            issue = repo.get_issue(int(args.issue))
            issue.create_comment("\n".join(comment))
    else:
        total_added, total_removed, total_changed = None, None, None
    return total_added, total_removed, total_changed


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Process datasets for database submissions")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input_pattern", default=None, help="Pattern matching input Dataset protos")
    group.add_argument("--input_file", default=None, help="File containing Dataset proto filenames")
    parser.add_argument("--root", default="", help="Root of the repository")
    parser.add_argument("--output_format", default=".pb.gz", help="Dataset output format")
    parser.add_argument("--write_errors", action="store_true", help="If True, errors will be written to *.error")
    parser.add_argument("--no-validate", action="store_true", help="If set, reactions will not be validated")
    parser.add_argument("--update", action="store_true", help="If True, update Reaction protos")
    parser.add_argument("--cleanup", action="store_true", help="If True, use git to clean up")
    parser.add_argument("--max_size", type=float, default=10.0, help="Maximum size (in MB) for any Reaction message")
    parser.add_argument("--base", default=None, help="Git branch to diff against")
    parser.add_argument("--issue", default=None, help="GitHub pull request number")
    parser.add_argument("--token", default=None, help="GitHub authentication token")
    return parser.parse_args(argv)


def main(args):
    silence_rdkit_logs()
    run(args)


if __name__ == "__main__":
    main(parse_args())
