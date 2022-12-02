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
"""Tests for ord_schema.scripts.process_dataset."""
import glob
import os
import subprocess
from typing import Optional

import docopt
import pytest
from rdkit import RDLogger

from ord_schema.logging import get_logger
from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import process_dataset

logger = get_logger(__name__)


class TestProcessDataset:
    @pytest.fixture
    def setup(self, tmp_path) -> tuple[str, str]:
        # Suppress RDKit warnings to clean up the test output.
        RDLogger.logger().setLevel(RDLogger.CRITICAL)
        reaction1 = reaction_pb2.Reaction()
        dummy_input = reaction1.inputs["dummy_input"]
        dummy_component = dummy_input.components.add()
        dummy_component.identifiers.add(type="CUSTOM")
        dummy_component.identifiers[0].details = "custom_identifier"
        dummy_component.identifiers[0].value = "custom_value"
        dummy_component.is_limiting = True
        dummy_component.amount.mass.value = 1
        dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
        reaction1.outcomes.add().conversion.value = 75
        reaction1.provenance.record_created.time.value = "2020-01-01"
        reaction1.provenance.record_created.person.username = "test"
        reaction1.provenance.record_created.person.email = "test@example.com"
        dataset1 = dataset_pb2.Dataset(
            dataset_id="ord_dataset-00000000000000000000000000000000",
            reactions=[reaction1],
        )
        dataset1_filename = (tmp_path / "dataset1.pbtxt").as_posix()
        message_helpers.write_message(dataset1, dataset1_filename)
        # reaction2 is empty.
        reaction2 = reaction_pb2.Reaction()
        dataset2 = dataset_pb2.Dataset(reactions=[reaction1, reaction2])
        dataset2_filename = (tmp_path / "dataset2.pb").as_posix()
        message_helpers.write_message(dataset2, dataset2_filename)
        yield dataset1_filename, dataset2_filename

    def test_main_with_input_pattern(self, setup):
        dataset1_filename, _ = setup
        argv = ["--input_pattern", dataset1_filename, "--base", "main"]
        process_dataset.main(docopt.docopt(process_dataset.__doc__, argv))

    def test_main_with_input_file(self, setup, tmp_path):
        dataset1_filename, _ = setup
        input_file = (tmp_path / "input_file.txt").as_posix()
        with open(input_file, "w") as f:
            f.write(f"A\t{dataset1_filename}\n")
        argv = ["--input_file", input_file, "--base", "main"]
        process_dataset.main(docopt.docopt(process_dataset.__doc__, argv))

    def test_main_with_validation_errors(self, setup):
        _, dataset2_filename = setup
        argv = ["--input_pattern", dataset2_filename, "--write_errors"]
        with pytest.raises(validations.ValidationError, match="validation encountered errors"):
            process_dataset.main(docopt.docopt(process_dataset.__doc__, argv))
        error_filename = f"{dataset2_filename}.error"
        assert os.path.exists(error_filename)
        expected_output = [
            "Reaction: Reactions should have at least 1 reaction input\n",
            "Reaction: Reactions should have at least 1 reaction outcome\n",
        ]
        with open(error_filename) as f:
            assert f.readlines() == expected_output

    def test_main_with_updates(self, setup):
        dataset1_filename, _ = setup
        dirname = os.path.dirname(dataset1_filename)
        argv = [
            "--input_pattern",
            dataset1_filename,
            "--root",
            dirname,
            "--base",
            "main",
            "--update",
        ]
        process_dataset.main(docopt.docopt(process_dataset.__doc__, argv))
        expected_output = os.path.join(dirname, "data", "00", "ord_dataset-00000000000000000000000000000000.pb.gz")
        assert os.path.exists(expected_output)
        dataset = message_helpers.load_message(expected_output, dataset_pb2.Dataset)
        assert len(dataset.reactions) == 1
        assert dataset.reactions[0].reaction_id.startswith("ord-")


class TestSubmissionWorkflow:
    """Test suite for the ORD submission workflow.

    setUp() starts each test with a clean git environment containing some
    data. To create a new test, you should made a modification to the git
    repo (e.g. adding a new dataset or editing an existing one) and call
    self._run() to commit the changes and run process_datasets.py.
    """

    _DEFAULT_BRANCH = "main"

    @pytest.fixture
    def setup(self, tmp_path) -> tuple[str, str]:
        test_subdirectory = tmp_path.as_posix()
        os.chdir(test_subdirectory)
        subprocess.run(["git", "init", "-b", self._DEFAULT_BRANCH], check=True)
        subprocess.run(["git", "config", "--local", "user.email", "test@ord-schema"], check=True)
        subprocess.run(["git", "config", "--local", "user.name", "Test Runner"], check=True)
        # Add some initial data.
        reaction = reaction_pb2.Reaction()
        methylamine = reaction.inputs["methylamine"]
        component = methylamine.components.add()
        component.identifiers.add(type="SMILES", value="CN")
        component.is_limiting = True
        component.amount.moles.value = 1
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 75
        reaction.provenance.record_created.time.value = "2020-01-01"
        reaction.provenance.record_created.person.username = "test"
        reaction.provenance.record_created.person.email = "test@example.com"
        reaction.reaction_id = "ord-10aed8b5dffe41fab09f5b2cc9c58ad9"
        dataset_id = "ord_dataset-64b14868c5cd46dd8e75560fd3589a6b"
        dataset = dataset_pb2.Dataset(reactions=[reaction], dataset_id=dataset_id)
        # Make sure the initial dataset is valid.
        validations.validate_message(dataset)
        os.makedirs(os.path.join("data", "64"))
        dataset_filename = os.path.join(test_subdirectory, "data", "64", f"{dataset_id}.pb.gz")
        message_helpers.write_message(dataset, dataset_filename)
        subprocess.run(["git", "add", "data"], check=True)
        subprocess.run(["git", "commit", "-m", "Initial commit"], check=True)
        # Use a new branch for tests.
        subprocess.run(["git", "checkout", "-b", "test"], check=True)
        yield test_subdirectory, dataset_filename

    def _run(self, test_subdirectory: str, extra_argv: Optional[list[str]] = None):
        """Runs process_dataset.main().

        Args:
            test_subdirectory: Directory containing test inputs/outputs.
            extra_argv: Extra arguments.

        Returns:
            added: Set of added reaction IDs.
            removed: Set of deleted reaction IDs.
            changed: Set of changed reaction IDs.
            filenames: List of .pb filenames in the updated database.
        """
        # These commands will fail if there are no files to match for a given
        # pattern, so run them separately to make sure we pick up changes.
        try:
            subprocess.run(["git", "add", "*.pb*"], check=True)
        except subprocess.CalledProcessError as error:
            logger.info(error)
        try:
            subprocess.run(["git", "add", "data/*/*.pb*"], check=True)
        except subprocess.CalledProcessError as error:
            logger.info(error)
        changed = subprocess.run(
            ["git", "diff", "--name-status", self._DEFAULT_BRANCH],
            check=True,
            capture_output=True,
            text=True,
        )
        with open("changed.txt", "w") as f:
            f.write(changed.stdout)
        logger.info(f"Changed files:\n{changed.stdout}")
        subprocess.run(["git", "commit", "-m", "Submission"], check=True)
        argv = [
            "--input_file",
            "changed.txt",
            "--update",
            "--cleanup",
            "--base",
            "main",
        ]
        if extra_argv:
            argv.extend(extra_argv)
        added, removed, changed = process_dataset.run(docopt.docopt(process_dataset.__doc__, argv))
        filenames = glob.glob(os.path.join(test_subdirectory, "**/*.pb*"), recursive=True)
        return added, removed, changed, filenames

    def test_add_dataset(self, setup):
        test_subdirectory, dataset_filename = setup
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs["ethylamine"]
        component = ethylamine.components.add()
        component.identifiers.add(type="SMILES", value="CCN")
        component.is_limiting = True
        component.amount.moles.value = 2
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = "2020-01-01"
        reaction.provenance.record_created.person.username = "test"
        reaction.provenance.record_created.person.email = "test@example.com"
        reaction.reaction_id = "test"
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        this_dataset_filename = os.path.join(test_subdirectory, "test.pbtxt")
        message_helpers.write_message(dataset, this_dataset_filename)
        added, removed, changed, filenames = self._run(test_subdirectory)
        assert added == {"test"}
        assert not removed
        assert not changed
        assert len(filenames) == 2
        assert not os.path.exists(this_dataset_filename)
        # Check for assignment of dataset and reaction IDs.
        filenames.pop(filenames.index(dataset_filename))
        assert len(filenames) == 1
        dataset = message_helpers.load_message(filenames[0], dataset_pb2.Dataset)
        assert dataset.dataset_id
        assert len(dataset.reactions) == 1
        assert dataset.reactions[0].reaction_id
        # Check for binary output.
        assert filenames[0].endswith(".pb.gz")

    def test_add_sharded_dataset(self, setup):
        test_subdirectory, dataset_filename = setup
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs["ethylamine"]
        component = ethylamine.components.add()
        component.identifiers.add(type="SMILES", value="CCN")
        component.is_limiting = True
        component.amount.moles.value = 2
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = "2020-01-02"
        reaction.provenance.record_created.person.username = "test2"
        reaction.provenance.record_created.person.email = "test2@example.com"
        reaction.reaction_id = "test1"
        dataset1 = dataset_pb2.Dataset(reactions=[reaction])
        dataset1_filename = os.path.join(test_subdirectory, "test1.pbtxt")
        message_helpers.write_message(dataset1, dataset1_filename)
        reaction.provenance.record_created.time.value = "2020-01-03"
        reaction.provenance.record_created.person.username = "test3"
        reaction.provenance.record_created.person.email = "test3@example.com"
        reaction.reaction_id = "test2"
        dataset2 = dataset_pb2.Dataset(reactions=[reaction])
        dataset2_filename = os.path.join(test_subdirectory, "test2.pbtxt")
        message_helpers.write_message(dataset2, dataset2_filename)
        added, removed, changed, filenames = self._run(test_subdirectory)
        assert added == {"test1", "test2"}
        assert not removed
        assert not changed
        assert len(filenames) == 3
        filenames.pop(filenames.index(dataset_filename))
        assert len(filenames) == 2
        for filename in filenames:
            dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
            assert len(dataset.reactions) == 1
        assert not os.path.exists(dataset1_filename)
        assert not os.path.exists(dataset2_filename)

    def test_add_dataset_with_existing_reaction_ids(self, setup):
        test_subdirectory, dataset_filename = setup
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs["ethylamine"]
        component = ethylamine.components.add()
        component.identifiers.add(type="SMILES", value="CCN")
        component.is_limiting = True
        component.amount.moles.value = 2
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction_id = "ord-10aed8b5dffe41fab09f5b2cc9c58ad9"
        reaction.reaction_id = reaction_id
        reaction.provenance.record_created.time.value = "2020-01-01"
        reaction.provenance.record_created.person.username = "test"
        reaction.provenance.record_created.person.email = "test@example.com"
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        this_dataset_filename = os.path.join(test_subdirectory, "test.pbtxt")
        message_helpers.write_message(dataset, this_dataset_filename)
        added, removed, changed, filenames = self._run(test_subdirectory)
        assert added == {"ord-10aed8b5dffe41fab09f5b2cc9c58ad9"}
        assert not removed
        assert not changed
        assert len(filenames) == 2
        assert not os.path.exists(this_dataset_filename)
        filenames.pop(filenames.index(dataset_filename))
        assert len(filenames) == 1
        dataset = message_helpers.load_message(filenames[0], dataset_pb2.Dataset)
        # Check that existing record IDs for added datasets are not overridden.
        assert dataset.reactions[0].reaction_id == reaction_id
        assert len(dataset.reactions[0].provenance.record_modified) == 0

    def test_modify_dataset(self, setup):
        test_subdirectory, dataset_filename = setup
        dataset = message_helpers.load_message(dataset_filename, dataset_pb2.Dataset)
        # Modify the existing reaction...
        reaction1 = dataset.reactions[0]
        reaction1.inputs["methylamine"].components[0].amount.moles.value = 2
        # ...and add a new reaction.
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs["ethylamine"]
        component = ethylamine.components.add()
        component.identifiers.add(type="SMILES", value="CCN")
        component.is_limiting = True
        component.amount.moles.value = 2
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = "2020-01-01"
        reaction.provenance.record_created.person.username = "test"
        reaction.provenance.record_created.person.email = "test@example.com"
        reaction.reaction_id = "test"
        dataset.reactions.add().CopyFrom(reaction)
        message_helpers.write_message(dataset, dataset_filename)
        added, removed, changed, filenames = self._run(test_subdirectory)
        assert added == {"test"}
        assert not removed
        assert changed == {"ord-10aed8b5dffe41fab09f5b2cc9c58ad9"}
        assert filenames == [dataset_filename]
        # Check for preservation of dataset and record IDs.
        updated_dataset = message_helpers.load_message(dataset_filename, dataset_pb2.Dataset)
        assert len(updated_dataset.reactions) == 2
        assert dataset.dataset_id == updated_dataset.dataset_id
        assert dataset.reactions[0].reaction_id == updated_dataset.reactions[0].reaction_id
        assert updated_dataset.reactions[1].reaction_id

    def test_modify_reaction_id(self, setup):
        test_subdirectory, dataset_filename = setup
        dataset = message_helpers.load_message(dataset_filename, dataset_pb2.Dataset)
        dataset.reactions[0].reaction_id = "test_rename"
        message_helpers.write_message(dataset, dataset_filename)
        added, removed, changed, filenames = self._run(test_subdirectory)
        assert added == {"test_rename"}
        assert removed == {"ord-10aed8b5dffe41fab09f5b2cc9c58ad9"}
        assert not changed
        assert filenames == [dataset_filename]

    def test_add_dataset_with_validation_errors(self, setup):
        test_subdirectory, _ = setup
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs["ethylamine"]
        component = ethylamine.components.add()
        component.identifiers.add(type="SMILES", value="C#O")
        component.is_limiting = True
        component.amount.moles.value = 2
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(test_subdirectory, "test.pbtxt")
        message_helpers.write_message(dataset, dataset_filename)
        with pytest.raises(validations.ValidationError, match="could not validate SMILES"):
            self._run(test_subdirectory)

    def test_add_sharded_dataset_with_validation_errors(self, setup):
        test_subdirectory, _ = setup
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs["ethylamine"]
        component = ethylamine.components.add()
        component.identifiers.add(type="SMILES", value="CCN")
        component.is_limiting = True
        component.amount.moles.value = 2
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        reaction.provenance.record_created.time.value = "2021-02-09"
        reaction.provenance.record_created.person.username = "bob"
        reaction.provenance.record_created.person.email = "bob@bob.com"
        dataset1 = dataset_pb2.Dataset(reactions=[reaction])
        dataset1_filename = os.path.join(test_subdirectory, "test1.pbtxt")
        message_helpers.write_message(dataset1, dataset1_filename)
        reaction.inputs["ethylamine"].components[0].identifiers[0].value = "C#O"
        dataset2 = dataset_pb2.Dataset(reactions=[reaction])
        dataset2_filename = os.path.join(test_subdirectory, "test2.pbtxt")
        message_helpers.write_message(dataset2, dataset2_filename)
        with pytest.raises(validations.ValidationError, match="could not validate SMILES"):
            self._run(test_subdirectory)

    def test_modify_dataset_with_validation_errors(self, setup):
        test_subdirectory, dataset_filename = setup
        dataset = message_helpers.load_message(dataset_filename, dataset_pb2.Dataset)
        reaction = dataset.reactions[0]
        reaction.inputs["methylamine"].components[0].amount.moles.value = -2
        message_helpers.write_message(dataset, dataset_filename)
        with pytest.raises(validations.ValidationError, match="must be non-negative"):
            self._run(test_subdirectory)

    def test_add_dataset_with_too_large_reaction(self, setup):
        test_subdirectory, _ = setup
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs["ethylamine"]
        component = ethylamine.components.add()
        component.identifiers.add(type="SMILES", value="CCN")
        component.is_limiting = True
        component.amount.moles.value = 2
        component.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.outcomes.add().conversion.value = 25
        image = reaction.observations.add().image
        image.bytes_value = b"test data value"
        image.format = "png"
        dataset = dataset_pb2.Dataset(reactions=[reaction])
        dataset_filename = os.path.join(test_subdirectory, "test.pbtxt")
        message_helpers.write_message(dataset, dataset_filename)
        with pytest.raises(ValueError, match="larger than --max_size"):
            self._run(test_subdirectory, ["--max_size", "0.0"])

    def test_delete_dataset(self, setup):
        test_subdirectory, dataset_filename = setup
        subprocess.run(["git", "rm", dataset_filename], check=True)
        added, removed, changed, filenames = self._run(test_subdirectory)
        assert not added
        assert len(removed) == 1
        assert not changed
        assert not filenames

    def test_replace_dataset(self, setup):
        test_subdirectory, dataset_filename = setup
        dataset = message_helpers.load_message(dataset_filename, dataset_pb2.Dataset)
        this_dataset_filename = os.path.join(test_subdirectory, "test.pbtxt")
        message_helpers.write_message(dataset, this_dataset_filename)
        subprocess.run(["git", "rm", dataset_filename], check=True)
        added, removed, changed, filenames = self._run(test_subdirectory)
        assert len(added) == 1
        assert len(removed) == 1
        assert not changed
        assert len(filenames) == 1
