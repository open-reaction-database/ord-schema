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

Example usage:

* For normal validation-only operation:
  $ python process_dataset.py --input_pattern=my_dataset.pb
* To write the validated protos to disk:
  $ python process_dataset.py --input_pattern=my_dataset.pb --update
* To write errors to disk:
  $ python process_dataset.py --input_pattern=my_dataset.pb --write_errors
* To process multiple Dataset protos:
  $ python process_dataset.py --input_pattern="my_dataset-*.pb"
"""

import dataclasses
import glob
import gzip
import os
import subprocess
import sys
from typing import Iterable, List, Mapping, Set, Tuple

from absl import app
from absl import flags
from absl import logging
import github
from rdkit import RDLogger

from ord_schema import message_helpers
from ord_schema import updates
from ord_schema import validations
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string("root", "", "Root of the repository.")
flags.DEFINE_string("input_pattern", None, "Pattern (glob) matching input Dataset protos.")
flags.DEFINE_string("input_file", None, "Filename containing Dataset proto filenames.")
flags.DEFINE_enum(
    "output_format",
    ".pb.gz",
    [".pb.gz", ".pb", ".pbtxt.gz", ".pbtxt"],
    "Dataset output format.",
)
flags.DEFINE_boolean("write_errors", False, "If True, errors will be written to <filename>.error.")
flags.DEFINE_boolean("validate", True, "If True, validate input Reaction protos.")
flags.DEFINE_boolean("update", False, "If True, update Reaction protos.")
flags.DEFINE_boolean("cleanup", False, "If True, use git to clean up.")
flags.DEFINE_float("max_size", 10.0, "Maximum size (in MB) for any Reaction message.")
flags.DEFINE_string("base", None, "Git branch to diff against.")
flags.DEFINE_integer("issue", None, "GitHub pull request number. If provided, a comment will be added.")
flags.DEFINE_string("token", None, "GitHub authentication token.")

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


def _get_inputs() -> List[FileStatus]:
    """Gets a list of Dataset proto filenames to process.

    Returns:
        List of FileStatus objects.

    Raises:
        ValueError: If a git-diff status is not one of {'A', 'D', 'M', 'R'}.
    """
    if FLAGS.input_pattern and FLAGS.input_file:
        raise ValueError("one of --input_pattern or --input_file is required, not both")
    if FLAGS.input_pattern:
        # Setting recursive=True allows recursive matching with '**'.
        filenames = glob.glob(FLAGS.input_pattern, recursive=True)
        return [FileStatus(filename, "A", "") for filename in filenames]
    if FLAGS.input_file:
        inputs = []
        with open(FLAGS.input_file) as f:
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
        logging.info("editing an existing dataset; no cleanup needed")
        return  # Reuse the existing dataset ID.
    args = ["git", "mv", filename, output_filename]
    logging.info("Running command: %s", " ".join(args))
    subprocess.run(args, check=True)


def _get_reaction_ids(dataset: dataset_pb2.Dataset) -> Set[str]:
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
    logging.info("Running command: %s", " ".join(args))
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
) -> Tuple[Set[str], Set[str], Set[str]]:
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


def _run_updates(datasets: Mapping[str, dataset_pb2.Dataset]):
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
    validations.validate_datasets(datasets, FLAGS.write_errors, options=options)
    for filename, dataset in datasets.items():
        output_filename = os.path.join(
            FLAGS.root,
            message_helpers.id_filename(f"{dataset.dataset_id}{FLAGS.output_format}"),
        )
        os.makedirs(os.path.dirname(output_filename), exist_ok=True)
        if FLAGS.cleanup:
            cleanup(filename, output_filename)
        logging.info("writing Dataset to %s", output_filename)
        message_helpers.write_message(dataset, output_filename)


def run() -> Tuple[Set[str], Set[str], Set[str]]:
    """Main function that returns added/removed reaction ID sets.

    This function should be called directly by tests to get access to the
    return values. If main() returns something other than None it will break
    shell error code logic downstream.

    Returns:
        added: Set of added reaction IDs.
        removed: Set of deleted reaction IDs.
        changed: Set of changed reaction IDs.
    """
    inputs = sorted(_get_inputs())
    if not inputs:
        logging.info("nothing to do")
        return set(), set(), set()  # Nothing to do.
    # NOTE(kearnes): Process one dataset at a time to avoid OOM errors.
    change_stats = {}
    for file_status in inputs:
        if file_status.status == "D":
            dataset = None
        else:
            dataset = message_helpers.load_message(file_status.filename, dataset_pb2.Dataset)
            logging.info("%s: %d reactions", file_status.filename, len(dataset.reactions))
        datasets = {file_status.filename: dataset}
        if FLAGS.validate and dataset is not None:
            # Note: this does not check if IDs are malformed.
            validations.validate_datasets(datasets, FLAGS.write_errors)
            # Check reaction sizes.
            for reaction in dataset.reactions:
                reaction_size = sys.getsizeof(reaction.SerializeToString()) / 1e6
                if reaction_size > FLAGS.max_size:
                    raise ValueError("Reaction is larger than --max_size " f"({reaction_size} vs {FLAGS.max_size}")
        if FLAGS.base:
            added, removed, changed = get_change_stats(datasets, [file_status], base=FLAGS.base)
            change_stats[file_status.filename] = (added, removed, changed)
            logging.info(
                "Summary: +%d -%d Δ%d reaction IDs",
                len(added),
                len(removed),
                len(changed),
            )
        if FLAGS.update and dataset is not None:
            _run_updates(datasets)
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
        if FLAGS.issue and FLAGS.token:
            client = github.Github(FLAGS.token)
            repo = client.get_repo(os.environ["GITHUB_REPOSITORY"])
            issue = repo.get_issue(FLAGS.issue)
            issue.create_comment("\n".join(comment))
    else:
        total_added, total_removed, total_changed = None, None, None
    return total_added, total_removed, total_changed


def main(argv):
    del argv  # Only used by app.run().
    RDLogger.DisableLog("rdApp.*")  # Disable RDKit logging.
    run()


if __name__ == "__main__":
    app.run(main)
