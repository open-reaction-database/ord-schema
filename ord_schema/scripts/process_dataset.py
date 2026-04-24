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
import re
import subprocess
import sys
import tempfile
import uuid
from collections.abc import Iterable

import github

from ord_schema import message_helpers, parquet_dataset, updates, validations
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.proto import dataset_pb2, reaction_pb2

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


def _is_parquet(filename: str) -> bool:
    return filename.endswith(".parquet")


def _get_reaction_ids(dataset: dataset_pb2.Dataset) -> set[str]:
    """Returns a set containing the reaction IDs in a Dataset."""
    reaction_ids = set()
    for reaction in dataset.reactions:
        if reaction.reaction_id:
            reaction_ids.add(reaction.reaction_id)
    return reaction_ids


def _reaction_ids_from_path(filename: str) -> set[str]:
    """Streams reaction IDs from a Parquet file without deserializing Reactions."""
    return {rid for rid in parquet_dataset.iter_reaction_ids(filename) if rid}


def _load_base_reaction_ids(file_status: FileStatus, base: str) -> set[str] | None:
    """Returns the reaction IDs from the base-branch version of a file, or None if new."""
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
    git_path = git_args[-1]
    if _is_parquet(git_path) or git_path.endswith(".parquet.gz"):
        # Second clause is defensive; Parquet files are not gzip-wrapped.
        return {rid for rid in parquet_dataset.iter_reaction_ids(serialized.stdout) if rid}
    if git_path.endswith(".gz"):
        value = gzip.decompress(serialized.stdout)
    else:
        value = serialized.stdout
    return _get_reaction_ids(dataset_pb2.Dataset.FromString(value))


def get_change_stats(
    reaction_ids: Iterable[tuple[FileStatus, set[str] | None]],
    base: str,
) -> tuple[set[str], set[str], set[str]]:
    """Computes diff statistics for the submission.

    Args:
        reaction_ids: Iterable of (FileStatus, current_reaction_ids). The
            current set is None only for deleted files.
        base: Git branch to diff against.

    Returns:
        added: Set of added reaction IDs.
        removed: Set of deleted reaction IDs.
        changed: Set of changed reaction IDs.
    """
    old, new = set(), set()
    for file_status, current in reaction_ids:
        if not file_status.status.startswith("D"):
            if current is None:
                raise ValueError(f"missing reaction IDs for non-deleted file: {file_status.filename}")
            new.update(current)
        base_ids = _load_base_reaction_ids(file_status, base)
        if base_ids is not None:
            old.update(base_ids)
    return new - old, old - new, new & old


def _apply_id_substitutions(reaction: reaction_pb2.Reaction, id_substitutions: dict[str, str]) -> None:
    """Rewrites cross-referenced reaction IDs in-place using ``updates.update_dataset`` semantics."""
    for reaction_input in reaction.inputs.values():
        for component in reaction_input.components:
            for preparation in component.preparations:
                if preparation.reaction_id and preparation.reaction_id in id_substitutions:
                    preparation.reaction_id = id_substitutions[preparation.reaction_id]
        for crude_component in reaction_input.crude_components:
            if crude_component.reaction_id in id_substitutions:
                crude_component.reaction_id = id_substitutions[crude_component.reaction_id]


def _run_updates_eager(
    filename: str,
    dataset: dataset_pb2.Dataset,
    *,
    root: str,
    output_format: str,
    write_errors: bool,
    cleanup_files: bool,
) -> None:
    """Updates the submission file for a single non-Parquet Dataset."""
    updates.update_dataset(dataset)
    options = validations.ValidationOptions(validate_ids=True, require_provenance=True)
    validations.validate_datasets({filename: dataset}, write_errors, options=options)
    output_filename = os.path.join(
        root,
        message_helpers.id_filename(f"{dataset.dataset_id}{output_format}"),
    )
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    if cleanup_files:
        cleanup(filename, output_filename)
    logger.info("writing Dataset to %s", output_filename)
    message_helpers.write_message(dataset, output_filename)


def _run_updates_streaming(
    filename: str,
    *,
    root: str,
    output_format: str,
    write_errors: bool,
    cleanup_files: bool,
) -> None:
    """Updates a Parquet submission file via a streaming pipeline."""
    if output_format != ".parquet":
        raise ValueError(
            f"Parquet inputs require --output_format .parquet to preserve streaming; got {output_format!r}"
        )
    # Resolve the final output path up-front so the streaming pipeline can
    # write directly there. We peek at metadata to decide the dataset_id.
    metadata = parquet_dataset.read_metadata(filename)
    if re.fullmatch("^ord_dataset-[0-9a-f]{32}$", metadata.dataset_id):
        dataset_id = metadata.dataset_id
    else:
        dataset_id = f"ord_dataset-{uuid.uuid4().hex}"
    output_filename = os.path.join(
        root,
        message_helpers.id_filename(f"{dataset_id}{output_format}"),
    )
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    if filename == output_filename:
        # Streaming requires distinct input and output; write to a sibling and rename.
        tmp_output = f"{output_filename}.new"
        metadata.dataset_id = dataset_id  # Ensure consistency if we bypass the ID reassignment.
        _update_dataset_parquet_streaming_with_id(filename, tmp_output, dataset_id=dataset_id)
        os.replace(tmp_output, output_filename)
    else:
        _update_dataset_parquet_streaming_with_id(filename, output_filename, dataset_id=dataset_id)
        if cleanup_files:
            cleanup(filename, output_filename)
    options = validations.ValidationOptions(validate_ids=True, require_provenance=True)
    validations.validate_parquet_datasets([output_filename], write_errors=write_errors, options=options)
    logger.info("writing Dataset to %s", output_filename)


def _update_dataset_parquet_streaming_with_id(
    input_path: str,
    output_path: str,
    *,
    dataset_id: str,
) -> None:
    """Variant of the two-pass update pipeline with a pre-resolved dataset_id."""
    metadata = parquet_dataset.read_metadata(input_path)
    metadata.dataset_id = dataset_id
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    tmp_fd, tmp_path = tempfile.mkstemp(
        prefix=".update-",
        suffix=".parquet",
        dir=os.path.dirname(output_path) or None,
    )
    os.close(tmp_fd)
    try:
        id_substitutions: dict[str, str] = {}
        with parquet_dataset.DatasetWriter(
            tmp_path,
            name=metadata.name,
            description=metadata.description,
            dataset_id=metadata.dataset_id,
        ) as writer:
            for _, reaction in parquet_dataset.iter_reactions(input_path):
                id_substitutions.update(updates.update_reaction(reaction))
                writer.write(reaction)
        with parquet_dataset.DatasetWriter(
            output_path,
            name=metadata.name,
            description=metadata.description,
            dataset_id=metadata.dataset_id,
        ) as writer:
            for _, reaction in parquet_dataset.iter_reactions(tmp_path):
                _apply_id_substitutions(reaction, id_substitutions)
                writer.write(reaction)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


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
        is_parquet = _is_parquet(file_status.filename)
        deleted = file_status.status == "D"
        # Validation + size checks (stream for Parquet, load fully otherwise).
        current_reaction_ids: set[str] | None
        if deleted:
            current_reaction_ids = None
        elif is_parquet:
            if not args.no_validate:
                validations.validate_parquet_datasets([file_status.filename], write_errors=args.write_errors)
            current_reaction_ids = set()
            max_bytes = args.max_size * 1e6
            for reaction_id, reaction in parquet_dataset.iter_reactions(file_status.filename):
                blob = reaction.SerializeToString()
                if sys.getsizeof(blob) > max_bytes:
                    reaction_size = sys.getsizeof(blob) / 1e6
                    raise ValueError(f"Reaction is larger than --max_size ({reaction_size} vs {args.max_size})")
                if reaction_id:
                    current_reaction_ids.add(reaction_id)
            logger.info("%s: %d reactions", file_status.filename, len(current_reaction_ids))
        else:
            dataset = message_helpers.load_message(file_status.filename, dataset_pb2.Dataset)
            logger.info("%s: %d reactions", file_status.filename, len(dataset.reactions))
            if not args.no_validate:
                validations.validate_datasets({file_status.filename: dataset}, args.write_errors)
                for reaction in dataset.reactions:
                    reaction_size = sys.getsizeof(reaction.SerializeToString()) / 1e6
                    if reaction_size > args.max_size:
                        raise ValueError(f"Reaction is larger than --max_size ({reaction_size} vs {args.max_size})")
            current_reaction_ids = _get_reaction_ids(dataset)
        if args.base:
            added, removed, changed = get_change_stats(
                [(file_status, current_reaction_ids)],
                base=args.base,
            )
            change_stats[file_status.filename] = (added, removed, changed)
            logger.info(
                "Summary: +%d -%d Δ%d reaction IDs",
                len(added),
                len(removed),
                len(changed),
            )
        if args.update and not deleted:
            if is_parquet:
                _run_updates_streaming(
                    file_status.filename,
                    root=args.root,
                    output_format=args.output_format,
                    write_errors=args.write_errors,
                    cleanup_files=args.cleanup,
                )
            else:
                _run_updates_eager(
                    file_status.filename,
                    dataset,
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
