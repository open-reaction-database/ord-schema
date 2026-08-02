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
"""Tests for ord_schema.scripts.validate_dataset."""

import pathlib

import pytest

from ord_schema import datasets, parquet, validations
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.scripts import validate_dataset


@pytest.fixture(params=[".pbtxt", ".parquet"])
def setup(request, tmp_path) -> tuple[str, str]:
    """Writes dataset1 and dataset2 to ``tmp_path`` in the parametrized format.

    Returns ``(tmp_path, suffix)``.
    """
    suffix = request.param
    test_subdirectory = tmp_path.as_posix()
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
    reaction1.provenance.record_created.time.value = "2023-07-01"
    reaction1.provenance.record_created.person.name = "test"
    reaction1.provenance.record_created.person.email = "test@example.com"
    dataset1 = dataset_pb2.Dataset(
        name="test1", description="test1", reactions=[reaction1]
    )
    datasets.save_dataset(
        dataset1, pathlib.Path(test_subdirectory) / f"dataset1{suffix}"
    )
    # reaction2 is empty.
    reaction2 = reaction_pb2.Reaction()
    dataset2 = dataset_pb2.Dataset(
        name="test2", description="test2", reactions=[reaction1, reaction2]
    )
    datasets.save_dataset(
        dataset2, pathlib.Path(test_subdirectory) / f"dataset2{suffix}"
    )
    return test_subdirectory, suffix


def test_simple(setup):
    test_subdirectory, suffix = setup
    argv = [
        "--input_pattern",
        str(pathlib.Path(test_subdirectory) / f"dataset1{suffix}"),
    ]
    validate_dataset.main(validate_dataset.parse_args(argv))


def test_filter(setup):
    test_subdirectory, suffix = setup
    argv = [
        "--input_pattern",
        str(pathlib.Path(test_subdirectory) / f"dataset1{suffix}"),
        "--filter",
        "dataset",
    ]
    validate_dataset.main(validate_dataset.parse_args(argv))


def test_validation_errors(setup):
    test_subdirectory, suffix = setup
    argv = [
        "--input_pattern",
        str(pathlib.Path(test_subdirectory) / f"dataset*{suffix}"),
    ]
    with pytest.raises(
        validations.ValidationError,
        match="Reactions should have at least 1 reaction input",
    ):
        validate_dataset.main(validate_dataset.parse_args(argv))


def _valid_reaction(reaction_id: str = "") -> reaction_pb2.Reaction:
    """Builds a Reaction that passes per-reaction validation."""
    reaction = reaction_pb2.Reaction()
    if reaction_id:
        reaction.reaction_id = reaction_id
    dummy_input = reaction.inputs["dummy_input"]
    dummy_component = dummy_input.components.add()
    dummy_component.identifiers.add(type="CUSTOM")
    dummy_component.identifiers[0].details = "custom_identifier"
    dummy_component.identifiers[0].value = "custom_value"
    dummy_component.is_limiting = True
    dummy_component.amount.mass.value = 1
    dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
    reaction.outcomes.add().conversion.value = 75
    reaction.provenance.record_created.time.value = "2023-07-01"
    reaction.provenance.record_created.person.name = "test"
    reaction.provenance.record_created.person.email = "test@example.com"
    return reaction


def test_parquet_cross_row_group_duplicate_id(tmp_path):
    """A duplicate reaction_id split across row groups must still be caught.

    Exercises the row-group parallelism path: row_group_size=2 places the
    repeated id in row groups 0 and 2, so cross-reference detection only
    succeeds if per-row-group state is correctly merged.
    """
    repeated_id = "ord-deadbeef00000000000000000000a1"
    reactions = [
        _valid_reaction(reaction_id=repeated_id),
        _valid_reaction(reaction_id="ord-deadbeef00000000000000000000a2"),
        _valid_reaction(reaction_id="ord-deadbeef00000000000000000000a3"),
        _valid_reaction(reaction_id="ord-deadbeef00000000000000000000a4"),
        _valid_reaction(reaction_id=repeated_id),
    ]
    path = tmp_path / "cross_group_dup.parquet"
    with parquet.DatasetWriter(
        str(path), name="test", description="test", row_group_size=2
    ) as writer:
        writer.write_all(reactions)
    # Sanity-check the layout the test is asserting against.
    assert parquet.num_row_groups(str(path)) == 3
    with pytest.raises(
        validations.ValidationError,
        match="Multiple Reactions should never have the same IDs",
    ):
        validate_dataset.main(
            validate_dataset.parse_args(["--input_pattern", str(path)])
        )


def test_parquet_empty_file(tmp_path):
    """A zero-row-group parquet still hits the dataset-level finalize path.

    Exercises the synchronous fallthrough in ``main`` (no row-group tasks submitted) and
    asserts the standard 'requires reactions or reaction_ids' error surfaces.
    """
    path = tmp_path / "empty.parquet"
    with parquet.DatasetWriter(str(path), name="test", description="test") as writer:
        writer.write_all([])
    assert parquet.num_row_groups(str(path)) == 0
    with pytest.raises(
        validations.ValidationError,
        match="Dataset requires reactions or reaction_ids",
    ):
        validate_dataset.main(
            validate_dataset.parse_args(["--input_pattern", str(path)])
        )


def test_parquet_per_reaction_error_in_late_row_group(tmp_path):
    """Per-reaction errors must surface from any row group, not just the first."""
    reactions = [
        _valid_reaction(),
        _valid_reaction(),
        reaction_pb2.Reaction(),  # invalid: no inputs / no outcomes
    ]
    path = tmp_path / "late_invalid.parquet"
    with parquet.DatasetWriter(
        str(path), name="test", description="test", row_group_size=2
    ) as writer:
        writer.write_all(reactions)
    assert parquet.num_row_groups(str(path)) == 2
    with pytest.raises(
        validations.ValidationError,
        match="Reactions should have at least 1 reaction input",
    ):
        validate_dataset.main(
            validate_dataset.parse_args(["--input_pattern", str(path)])
        )


@pytest.mark.parametrize(
    ("pattern", "expected"),
    [
        (r"data/\d{2}", ["data/11/foo.pb"]),
        (r"data/\d[a-z]", ["data/1a/foo.pb"]),
        (r"data/[a-z]\d", ["data/a1/foo.pb"]),
        (r"data/[a-z]{2}", ["data/aa/foo.pb"]),
    ],
)
def test_filter_filenames(pattern, expected):
    filenames = [
        "dataset.pb",
        "data/a1/foo.pb",
        "data/aa/foo.pb",
        "data/11/foo.pb",
        "data/1a/foo.pb",
    ]
    assert expected == validate_dataset.filter_filenames(filenames, pattern)


@pytest.mark.parametrize(
    ("value", "expected"),
    [(None, None), ("0/1", (0, 1)), ("0/4", (0, 4)), ("3/4", (3, 4))],
)
def test_parse_shard(value, expected):
    assert validate_dataset.parse_shard(value) == expected


@pytest.mark.parametrize("value", ["", "1", "4/", "/4", "a/4", "4/4", "-1/4", "0/0"])
def test_parse_shard_rejects_malformed(value):
    with pytest.raises(ValueError, match="--shard"):
        validate_dataset.parse_shard(value)


def _sharded_dataset(tmp_path, name):
    """Writes a 4-row-group parquet whose third row group holds an invalid reaction."""
    reactions = [_valid_reaction() for _ in range(7)]
    reactions.append(reaction_pb2.Reaction())  # invalid: no inputs / no outcomes
    path = tmp_path / name
    with parquet.DatasetWriter(
        str(path), name="test", description="test", row_group_size=2
    ) as writer:
        writer.write_all(reactions)
    assert parquet.num_row_groups(str(path)) == 4
    return path


def test_shard_covers_every_row_group_between_them(tmp_path):
    """Some shard must catch an error, and exactly the shard holding it."""
    path = _sharded_dataset(tmp_path, "sharded.parquet")
    failing = []
    for index in range(4):
        try:
            validate_dataset.main(
                validate_dataset.parse_args(
                    ["--input_pattern", str(path), "--shard", f"{index}/4"]
                )
            )
        except validations.ValidationError:
            failing.append(index)
    # The invalid reaction is the 8th of 8, so it lands in the last row group,
    # which round-robin deals to shard 3.
    assert failing == [3]


def test_shard_zero_runs_dataset_level_checks(tmp_path):
    """A duplicate id spanning row groups must still be caught when sharded."""
    reaction = _valid_reaction()
    reaction.reaction_id = "ord-0000000000000000000000000000000a"
    path = tmp_path / "duplicate.parquet"
    with parquet.DatasetWriter(
        str(path), name="test", description="test", row_group_size=1
    ) as writer:
        writer.write_all([reaction, reaction])
    assert parquet.num_row_groups(str(path)) == 2
    # Unsharded, the duplicate is reported.
    with pytest.raises(validations.ValidationError, match=r"[Dd]uplicate"):
        validate_dataset.main(
            validate_dataset.parse_args(["--input_pattern", str(path)])
        )
    # No shard holds both copies, so the check cannot come from merging their
    # partial state. Shard 0 makes its own pass over the whole file and catches
    # it; the other shard reports only its own per-reaction findings.
    with pytest.raises(validations.ValidationError, match=r"[Dd]uplicate"):
        validate_dataset.main(
            validate_dataset.parse_args(
                ["--input_pattern", str(path), "--shard", "0/2"]
            )
        )
    validate_dataset.main(
        validate_dataset.parse_args(["--input_pattern", str(path), "--shard", "1/2"])
    )


def test_shard_zero_catches_cross_file_reference(tmp_path):
    """An undefined cross-reference is found even from another shard's slice."""
    defined = _valid_reaction()
    defined.reaction_id = "ord-0000000000000000000000000000000a"
    referring = _valid_reaction()
    referring.reaction_id = "ord-0000000000000000000000000000000b"
    component = referring.inputs["dummy_input"].components[0]
    component.preparations.add(
        type="SYNTHESIZED", reaction_id="ord-000000000000000000000000000000ff"
    )
    path = tmp_path / "crossref.parquet"
    with parquet.DatasetWriter(
        str(path), name="test", description="test", row_group_size=1
    ) as writer:
        writer.write_all([defined, referring])
    assert parquet.num_row_groups(str(path)) == 2
    # Shard 0 owns row group 0, so the referring reaction is not even in its
    # slice -- the finding can only come from its whole-file cross-ref pass.
    with pytest.raises(validations.ValidationError, match="undefined reaction_ids"):
        validate_dataset.main(
            validate_dataset.parse_args(
                ["--input_pattern", str(path), "--shard", "0/2"]
            )
        )


def test_sharded_and_unsharded_report_the_same_findings(tmp_path):
    """The merge path and shard 0's whole-file pass must not drift apart.

    Dataset-level state is assembled two ways -- merged from per-row-group
    observations when unsharded, and from one pass over the file on shard 0 --
    and both feed the same finalizer. Nothing else forces them to agree.
    """
    reaction = _valid_reaction()
    reaction.reaction_id = "ord-0000000000000000000000000000000a"
    referring = _valid_reaction()
    referring.reaction_id = "ord-0000000000000000000000000000000b"
    referring.inputs["dummy_input"].components[0].preparations.add(
        type="SYNTHESIZED", reaction_id="ord-000000000000000000000000000000ff"
    )
    path = tmp_path / "both.parquet"
    with parquet.DatasetWriter(
        str(path), name="test", description="test", row_group_size=1
    ) as writer:
        writer.write_all([reaction, reaction, referring])

    def findings(argv):
        try:
            validate_dataset.main(validate_dataset.parse_args(argv))
        except validations.ValidationError as error:
            return {
                line.split("Dataset: ", 1)[1]
                for line in str(error).splitlines()
                if "Dataset: " in line
            }
        return set()

    unsharded = findings(["--input_pattern", str(path)])
    sharded = findings(["--input_pattern", str(path), "--shard", "0/3"])
    assert unsharded, "expected the duplicate and the undefined reference"
    assert sharded == unsharded
