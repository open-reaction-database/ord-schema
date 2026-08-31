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

"""Tests for ord_schema.artifacts.scripts.derive_occurrences."""

import pathlib
import shutil

import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import occurrences, pivot, projection
from ord_schema.artifacts.scripts import derive_occurrences
from ord_schema.proto import dataset_pb2, reaction_pb2

_ROLE = reaction_pb2.ReactionRole.ReactionRoleType


@pytest.fixture
def pivots(tmp_path) -> pathlib.Path:
    """Pivots over every level the indexed paths are read from, in two shards."""
    root = tmp_path / "pivots"
    for shard in ("aa", "bb"):
        reaction = reaction_pb2.Reaction(reaction_id=f"ord-{shard}01")
        component = reaction.inputs["in"].components.add(reaction_role=_ROLE.SOLVENT)
        component.identifiers.add(type="SMILES", value="CCO")
        source = tmp_path / "data" / shard / f"ord_dataset-{shard}.parquet"
        source.parent.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            dataset_pb2.Dataset(
                dataset_id=f"ord_dataset-{shard}",
                name="test",
                description="test",
                reactions=[reaction],
            ),
            str(source),
        )
        projected = tmp_path / "projections" / shard / source.name
        projected.parent.mkdir(parents=True, exist_ok=True)
        projection.write_projection(source, projected)
        for level in {level.path for level, _ in occurrences.PATHS.values()}:
            target = root / level / shard / source.name
            target.parent.mkdir(parents=True, exist_ok=True)
            pivot.write_pivot(projected, target, level_path=level)
    return root


def _run(pivots: pathlib.Path, output: pathlib.Path, *extra: str) -> None:
    derive_occurrences.main(
        derive_occurrences.parse_args(
            [f"--pivots_dir={pivots}", f"--output_dir={output}", *extra]
        )
    )


def test_a_run_covers_every_path_and_mirrors_the_shards(pivots, tmp_path):
    # The layout is a contract with whatever reads it: a corpus globs a path's
    # directory recursively, so a file written anywhere else is one it cannot find.
    output = tmp_path / "occurrences"
    _run(pivots, output)
    for path in occurrences.PATHS:
        found = occurrences.artifact_paths(output, path)
        assert [artifact.parent.name for artifact in found] == ["aa", "bb"], path
    held = occurrences.artifact_paths(output, "inputs.components")
    assert sum(pq.read_metadata(artifact).num_rows for artifact in held) == 2


def test_a_second_run_writes_nothing(pivots, tmp_path, caplog):
    # Stamps say an artifact is current for its source, so a re-run over an unchanged
    # tree is footer reads. An interrupted run is finished by the next rather than
    # started again.
    output = tmp_path / "occurrences"
    _run(pivots, output)
    with caplog.at_level("INFO"):
        _run(pivots, output)
    messages = [record.getMessage() for record in caplog.records]
    assert any("0 written, 2 already current" in message for message in messages)


def test_forcing_rewrites_what_is_current(pivots, tmp_path, caplog):
    output = tmp_path / "occurrences"
    _run(pivots, output)
    with caplog.at_level("INFO"):
        _run(pivots, output, "--force")
    messages = [record.getMessage() for record in caplog.records]
    assert any("2 written, 0 already current" in message for message in messages)


def test_a_path_the_schema_carries_no_structure_at_is_refused(pivots, tmp_path):
    with pytest.raises(ValueError, match="no structure sits at"):
        _run(pivots, tmp_path / "occurrences", "--paths", "workups")


def test_a_path_named_twice_is_refused(pivots, tmp_path):
    # Deriving it twice would report the second pass as already current, which reads
    # like a tree that was half up to date.
    with pytest.raises(ValueError, match="named more than once"):
        _run(
            pivots,
            tmp_path / "occurrences",
            "--paths",
            "inputs.components",
            "inputs.components",
        )


def test_a_level_with_no_pivots_is_refused_rather_than_left_empty(tmp_path):
    # The silent case: a path with no artifacts leaves the structures sitting only
    # there reachable by no query, and a corpus reading the tree refuses it for having
    # lost them -- minutes into a load rather than here.
    with pytest.raises(ValueError, match="Derive them first with derive_pivots"):
        _run(tmp_path / "empty", tmp_path / "occurrences")


def test_the_levels_to_derive_first_are_printed(capsys):
    # The pivots have to exist before this script can read them, and which levels those
    # are is not obvious from the paths: three name a level of their own and the
    # authentic standard is read from the measurements'.
    derive_occurrences.main(derive_occurrences.parse_args(["--print_levels"]))
    printed = capsys.readouterr().out.split()
    assert printed == sorted({level.path for level, _ in occurrences.PATHS.values()})


def test_an_artifact_filed_under_another_path_is_refused(pivots, tmp_path):
    # Every path writes the same three columns, so a file copied from one to another
    # passes the currency check on every measure it has -- artifact name, source hash,
    # versions, columns -- and is skipped by every later run, answering that path's
    # quantifiers with another level's elements. The stamped path is the only thing
    # that tells them apart.
    output = tmp_path / "occurrences"
    _run(pivots, output)
    stranger = occurrences.artifact_paths(output, "inputs.components")[0]
    shutil.copy(stranger, output / "outcomes.products" / "aa" / stranger.name)
    with pytest.raises(ValueError, match="stamped with another path"):
        _run(pivots, output)


def test_a_path_whose_artifacts_all_hold_nothing_warns(pivots, tmp_path, caplog):
    # Ordinary for a path nothing in the corpus records, and identical on disk to a
    # derivation that filtered wrongly: an empty artifact is stamped current like any
    # other, so no later run revisits it and no reader tells the two apart.
    output = tmp_path / "occurrences"
    with caplog.at_level("WARNING"):
        _run(pivots, output, "--paths", "workups.input.components")
    messages = [record.getMessage() for record in caplog.records]
    assert any("every artifact is empty at this path" in m for m in messages)


def test_a_path_holding_rows_does_not_warn(pivots, tmp_path, caplog):
    output = tmp_path / "occurrences"
    with caplog.at_level("WARNING"):
        _run(pivots, output, "--paths", "inputs.components")
    assert not [record.getMessage() for record in caplog.records]
