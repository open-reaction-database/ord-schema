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

Usage:
    process_dataset.py (--input_pattern=<str> | --input_file=<str>) [options]

Options:
    --input_pattern=<str>       Pattern matching input Dataset protos
    --input_file=<str>          File containing Dataset proto filenames
    --root=<str>                Root of the repository [default: ]
    --output_format=<str>       Dataset output format [default: .pb.gz]
    --write_errors              If True, errors will be written to *.error
    --no-validate               If set, reactions will not be validated
    --update                    If True, update Reaction protos
    --cleanup                   If True, use git to clean up
    --max_size=<float>          Maximum size (in MB) for any Reaction message [default: 10.0]
    --base=<str>                Git branch to diff against
    --issue=<str>               GitHub pull request number; if provided, a comment will be added
    --token=<str>               GitHub authentication token
"""
import dataclasses
import glob
import gzip
import os
import subprocess
import sys
from collections.abc import Iterable, Mapping
from typing import Optional

import docopt
import github
from rdkit import RDLogger

from ord_schema.logging import get_logger
from ord_schema import message_helpers
from ord_schema import updates
from ord_schema import validations
from ord_schema.proto import dataset_pb2

logger = get_logger(__name__)


# pylint: disable=too-many-branches,too-many-locals


@dataclasses.dataclass(eq=True, frozen=True, order=True)
class FileStatus:
    """A filename and its status in Git."""

    filename: str
    status: str
    original_filename: str

    def __post_init__(self):
        if self.status[0] not in ["A", "D", "M", "R"]:
            raise ValueError(f"unsupported file status: {self.status}")


def _get_inputs(kwargs) -> list[FileStatus]:
    """Gets a list of Dataset proto filenames to process.

    Returns:
        List of FileStatus objects.

    Raises:
        ValueError: If a git-diff status is not one of {'A', 'D', 'M', 'R'}.
    """
    if kwargs["--input_pattern"]:
        # Setting recursive=True allows recursive matching with '**'.
        filenames = glob.glob(kwargs["--input_pattern"], recursive=True)
        return [FileStatus(filename, "A", "") for filename in filenames]
    if kwargs["--input_file"]:
        inputs = []
        with open(kwargs["--input_file"]) as f:
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


def _get_reaction_ids(dataset: dataset_pb2.Dataset) -> set[str]:
    """Returns a set containing the reaction IDs in a Dataset."""
    reaction_ids = set()
    for reaction in dataset.reactions:
        if reaction.reaction_id:
            reaction_ids.add(reaction.reaction_id)
    return reaction_ids


def _load_base_dataset(file_status: FileStatus, base: str) -> dataset_pb2.Dataset:
    """Loads a Dataset message from another branch."""
    if file_status.status.startswith("A"):
        return None  # Dataset only exists in the submission.
    # NOTE(kearnes): Use --no-pager to avoid a non-zero exit code.
    args = ["git", "--no-pager", "show"]
    if file_status.status.startswith("R"):
        args.append(f"{base}:{file_status.original_filename}")
    else:
        args.append(f"{base}:{file_status.filename}")
    logger.info("Running command: %s", " ".join(args))
    serialized = subprocess.run(args, capture_output=True, check=True, text=False)
    if serialized.stdout.startswith(b"version"):
        # Convert Git LFS pointers to real data.
        serialized = subprocess.run(
            ["git", "lfs", "smudge"],
            input=serialized.stdout,
            capture_output=True,
            check=True,
            text=False,
        )
    if args[-1].endswith(".gz"):
        value = gzip.decompress(serialized.stdout)
    else:
        value = serialized.stdout
    return dataset_pb2.Dataset.FromString(value)


def get_change_stats(
    datasets: Mapping[str, dataset_pb2.Dataset], inputs: Iterable[FileStatus], base: str
) -> tuple[set[str], set[str], set[str]]:
    """Computes diff statistics for the submission.

    Args:
        datasets: Dict mapping filenames to Dataset messages.
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
            new.update(_get_reaction_ids(datasets[file_status.filename]))
        dataset = _load_base_dataset(file_status, base)
        if dataset is not None:
            old.update(_get_reaction_ids(dataset))
    return new - old, old - new, new & old


def _run_updates(datasets: Mapping[str, dataset_pb2.Dataset], kwargs) -> None:
    """Updates the submission files.

    Args:
        datasets: Dict mapping filenames to Dataset messages.

    Raises:
        ValueError: if any Reaction is larger than FLAGS.max_size.
    """
    for dataset in datasets.values():
        # Set reaction_ids, resolve names, fix cross-references, etc.
        updates.update_dataset(dataset)
    # Final validation to make sure we didn't break anything.
    options = validations.ValidationOptions(validate_ids=True, require_provenance=True)
    validations.validate_datasets(datasets, kwargs["--write_errors"], options=options)
    for filename, dataset in datasets.items():
        output_filename = os.path.join(
            kwargs["--root"],
            message_helpers.id_filename(f'{dataset.dataset_id}{kwargs["--output_format"]}'),
        )
        os.makedirs(os.path.dirname(output_filename), exist_ok=True)
        if kwargs["--cleanup"]:
            cleanup(filename, output_filename)
        logger.info("writing Dataset to %s", output_filename)
        message_helpers.write_message(dataset, output_filename)


def run(kwargs) -> tuple[Optional[set[str]], Optional[set[str]], Optional[set[str]]]:
    """Main function that returns added/removed reaction ID sets.

    This function should be called directly by tests to get access to the
    return values. If main() returns something other than None it will break
    shell error code logic downstream.

    Returns:
        added: Set of added reaction IDs.
        removed: Set of deleted reaction IDs.
        changed: Set of changed reaction IDs.
    """
    inputs = sorted(_get_inputs(kwargs))
    if not inputs:
        logger.info("nothing to do")
        return set(), set(), set()  # Nothing to do.
    # NOTE(kearnes): Process one dataset at a time to avoid OOM errors.
    change_stats = {}
    for file_status in inputs:
        if file_status.status == "D":
            dataset = None
        else:
            dataset = message_helpers.load_message(file_status.filename, dataset_pb2.Dataset)
            logger.info("%s: %d reactions", file_status.filename, len(dataset.reactions))
        datasets = {file_status.filename: dataset}
        if not kwargs["--no-validate"] and dataset is not None:
            # Note: this does not check if IDs are malformed.
            validations.validate_datasets(datasets, kwargs["--write_errors"])
            # Check reaction sizes.
            for reaction in dataset.reactions:
                reaction_size = sys.getsizeof(reaction.SerializeToString()) / 1e6
                if reaction_size > float(kwargs["--max_size"]):
                    raise ValueError(
                        "Reaction is larger than --max_size " f'({reaction_size} vs {kwargs["--max_size"]}'
                    )
        if kwargs["--base"]:
            added, removed, changed = get_change_stats(datasets, [file_status], base=kwargs["--base"])
            change_stats[file_status.filename] = (added, removed, changed)
            logger.info(
                "Summary: +%d -%d Δ%d reaction IDs",
                len(added),
                len(removed),
                len(changed),
            )
        if kwargs["--update"] and dataset is not None:
            _run_updates(datasets, kwargs)
    if change_stats:
        total_added, total_removed, total_changed = set(), set(), set()
        comment = [
            "Change summary:",
            "| Filename | Added | Removed | Changed |",
            "| -------- | ----- | ------- | ------- |",
        ]
        for filename, (added, removed, changed) in change_stats.items():
            comment.append(f"| {filename} | " f"{len(added)} | {len(removed)} | {len(changed)} |")
            total_added |= added
            total_removed |= removed
            total_changed |= changed
        comment.append(f"| | **{len(total_added)}** | " f"**{len(total_removed)}** | " f"**{len(total_changed)}** |")
        if kwargs["--issue"] and kwargs["--token"]:
            client = github.Github(kwargs["--token"])
            repo = client.get_repo(os.environ["GITHUB_REPOSITORY"])
            issue = repo.get_issue(int(kwargs["--issue"]))
            issue.create_comment("\n".join(comment))
    else:
        total_added, total_removed, total_changed = None, None, None
    return total_added, total_removed, total_changed


def main(kwargs):
    RDLogger.DisableLog("rdApp.*")  # Disable RDKit logging.
    run(kwargs)


if __name__ == "__main__":
    main(docopt.docopt(__doc__))
