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

"""Tests for ord_schema.scripts.derive_projection."""

import pyarrow.parquet as pq
import pytest

from ord_schema import artifacts, parquet, projection
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.scripts import derive_projection


def _dataset(tmp_path, shard: str, dataset_id: str):
    reaction = reaction_pb2.Reaction(reaction_id=f"ord-{shard}")
    reaction.identifiers.add(type="REACTION_SMILES", value="C>>CO")
    directory = tmp_path / "data" / shard
    directory.mkdir(parents=True)
    path = directory / f"{dataset_id}.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id=dataset_id,
            name="test",
            description="test dataset",
            reactions=[reaction],
        ),
        path,
    )
    return path


def test_main_mirrors_the_input_layout(tmp_path):
    _dataset(tmp_path, "aa", "ord_dataset-a")
    _dataset(tmp_path, "bb", "ord_dataset-b")
    output_dir = tmp_path / "projections"
    derive_projection.main(
        derive_projection.parse_args(
            [
                f"--input_pattern={tmp_path}/data/*/*.parquet",
                f"--output_dir={output_dir}",
            ]
        )
    )
    assert (output_dir / "aa" / "ord_dataset-a.parquet").exists()
    assert (output_dir / "bb" / "ord_dataset-b.parquet").exists()


def test_main_writes_a_readable_projection(tmp_path):
    _dataset(tmp_path, "aa", "ord_dataset-a")
    output_dir = tmp_path / "projections"
    derive_projection.main(
        derive_projection.parse_args(
            [
                f"--input_pattern={tmp_path}/data/*/*.parquet",
                f"--output_dir={output_dir}",
            ]
        )
    )
    output = output_dir / "aa" / "ord_dataset-a.parquet"
    table = pq.read_table(output)
    assert table.num_rows == 1
    assert table.column("reaction_id").to_pylist() == ["ord-aa"]
    assert artifacts.load_stamps(output).artifact == projection.ARTIFACT


def test_rerunning_skips_current_projections(tmp_path):
    _dataset(tmp_path, "aa", "ord_dataset-a")
    output_dir = tmp_path / "projections"
    args = derive_projection.parse_args(
        [f"--input_pattern={tmp_path}/data/*/*.parquet", f"--output_dir={output_dir}"]
    )
    derive_projection.main(args)
    output = output_dir / "aa" / "ord_dataset-a.parquet"
    # Inode rather than mtime: the atomic publish renames, so identity is exact.
    first = output.stat().st_ino
    derive_projection.main(args)
    assert output.stat().st_ino == first


def test_force_rewrites_a_current_projection(tmp_path):
    _dataset(tmp_path, "aa", "ord_dataset-a")
    output_dir = tmp_path / "projections"
    base = [
        f"--input_pattern={tmp_path}/data/*/*.parquet",
        f"--output_dir={output_dir}",
    ]
    derive_projection.main(derive_projection.parse_args(base))
    output = output_dir / "aa" / "ord_dataset-a.parquet"
    first = output.stat().st_ino
    derive_projection.main(derive_projection.parse_args([*base, "--force"]))
    assert output.stat().st_ino != first


def test_an_exact_filename_writes_a_file_not_a_directory(tmp_path):
    source = _dataset(tmp_path, "aa", "ord_dataset-a")
    output_dir = tmp_path / "projections"
    derive_projection.main(
        derive_projection.parse_args(
            [f"--input_pattern={source}", f"--output_dir={output_dir}"]
        )
    )
    output = output_dir / "ord_dataset-a.parquet"
    assert output.is_file()


def test_writing_into_the_source_tree_is_an_error(tmp_path):
    _dataset(tmp_path, "aa", "ord_dataset-a")
    with pytest.raises(ValueError, match="cannot be written over its own source"):
        derive_projection.main(
            derive_projection.parse_args(
                [
                    f"--input_pattern={tmp_path}/data/*/*.parquet",
                    f"--output_dir={tmp_path}/data",
                ]
            )
        )


def test_a_pattern_reaching_the_output_tree_stays_rerunnable(tmp_path):
    _dataset(tmp_path, "aa", "ord_dataset-a")
    output_dir = tmp_path / "data" / "projections"
    args = derive_projection.parse_args(
        [
            f"--input_pattern={tmp_path}/data/**/*.parquet",
            f"--output_dir={output_dir}",
        ]
    )
    derive_projection.main(args)
    output = output_dir / "aa" / "ord_dataset-a.parquet"
    first = output.stat().st_ino
    # The projection is now itself a match for the pattern; it is not a source.
    derive_projection.main(args)
    assert output.stat().st_ino == first
    assert artifacts.load_stamps(output).artifact == projection.ARTIFACT


def test_an_unmatched_pattern_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="no datasets matched"):
        derive_projection.main(
            derive_projection.parse_args(
                [
                    f"--input_pattern={tmp_path}/data/*/*.parquet",
                    f"--output_dir={tmp_path}/out",
                ]
            )
        )
