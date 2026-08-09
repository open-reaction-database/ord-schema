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

"""Tests for ord_schema.scripts.derive_structures.

Artifact behavior is covered in ord_schema.structures_test; these cover the CLI: path
mapping, which matches count as sources, the skip-if-current shortcut, and --force.
"""

import pathlib

import pyarrow.parquet as pq
import pytest

from ord_schema import artifacts, parquet, structures
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.scripts import derive_projection, derive_structures


def _dataset(dataset_id: str):
    reaction = reaction_pb2.Reaction(reaction_id="ord-0001")
    component = reaction.inputs["a"].components.add()
    component.identifiers.add(type="SMILES", value="c1ccccc1")
    return dataset_pb2.Dataset(
        dataset_id=dataset_id, name="test", description="desc", reactions=[reaction]
    )


def _corpus(root: pathlib.Path, shards=("aa", "bb")) -> None:
    """Writes data/<shard>/ord_dataset-<shard>.parquet under root."""
    for shard in shards:
        directory = root / "data" / shard
        directory.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            _dataset(f"ord_dataset-{shard}"),
            str(directory / f"ord_dataset-{shard}.parquet"),
        )


def _project(root: pathlib.Path) -> None:
    """Projects the corpus, which is what structures are derived from."""
    derive_projection.main(
        derive_projection.parse_args(
            [
                f"--input_pattern={root / 'data' / '*' / '*.parquet'}",
                f"--output_dir={root / 'projections'}",
            ]
        )
    )


def _run(root: pathlib.Path, *extra: str) -> None:
    derive_structures.main(
        derive_structures.parse_args(
            [
                f"--input_pattern={root / 'projections' / '*' / '*.parquet'}",
                f"--output_dir={root / 'structures'}",
                *extra,
            ]
        )
    )


def test_main_writes_one_artifact_per_dataset(tmp_path):
    _corpus(tmp_path)
    _project(tmp_path)
    _run(tmp_path)
    for shard in ("aa", "bb"):
        artifact = tmp_path / "structures" / shard / f"ord_dataset-{shard}.parquet"
        assert artifact.exists()
        table = pq.read_table(artifact)
        assert table.schema.names == structures.SCHEMA.names
        assert table.column("smiles").to_pylist() == ["c1ccccc1"]


def test_main_stamps_artifacts_with_their_source(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    source = tmp_path / "data" / "aa" / "ord_dataset-aa.parquet"
    stamps = artifacts.load_stamps(
        tmp_path / "structures" / "aa" / "ord_dataset-aa.parquet"
    )
    assert stamps.artifact == structures.ARTIFACT
    assert stamps.source_dataset_id == "ord_dataset-aa"
    assert stamps.source_md5 == parquet.DatasetView(source).md5()


# write_structures publishes by renaming a fresh temp file into place, so the inode
# changes if and only if the artifact was rewritten. That is exact, unlike mtimes.
def test_main_skips_artifacts_that_are_current(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    artifact = tmp_path / "structures" / "aa" / "ord_dataset-aa.parquet"
    before = artifact.stat().st_ino
    _run(tmp_path)
    assert artifact.stat().st_ino == before


def test_force_rewrites_a_current_artifact(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    artifact = tmp_path / "structures" / "aa" / "ord_dataset-aa.parquet"
    before = artifact.stat().st_ino
    _run(tmp_path, "--force")
    assert artifact.stat().st_ino != before


def test_main_ignores_source_datasets(tmp_path):
    # The structures artifact featurizes a projection, so a pattern aimed at the
    # sources has nothing to read -- rather than silently deriving one the slow way.
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    with pytest.raises(ValueError, match="are projections"):
        derive_structures.main(
            derive_structures.parse_args(
                [
                    f"--input_pattern={tmp_path}/data/*/*.parquet",
                    f"--output_dir={tmp_path}/structures",
                ]
            )
        )


def test_main_raises_when_nothing_matches(tmp_path):
    with pytest.raises(ValueError, match="no datasets matched"):
        derive_structures.main(
            derive_structures.parse_args(
                [
                    f"--input_pattern={tmp_path / 'absent' / '*.parquet'}",
                    f"--output_dir={tmp_path / 'structures'}",
                ]
            )
        )
