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

from datetime import datetime

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import base, first_commits, projection
from ord_schema.artifacts.scripts import derive_projection
from ord_schema.proto import dataset_pb2, reaction_pb2


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
    assert base.load_stamps(output).artifact == projection.ARTIFACT


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
    with pytest.raises(ValueError, match="would write over its own inputs"):
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
    assert base.load_stamps(output).artifact == projection.ARTIFACT


def test_a_pattern_aimed_at_the_output_tree_is_an_error(tmp_path):
    # Exits 0 with an empty output directory otherwise, and a downstream step proceeds
    # as though the projections had been built.
    _dataset(tmp_path, "aa", "ord_dataset-a")
    derive_projection.main(
        derive_projection.parse_args(
            [
                f"--input_pattern={tmp_path}/data/*/*.parquet",
                f"--output_dir={tmp_path}/projections",
            ]
        )
    )
    with pytest.raises(ValueError, match="no sources derived"):
        derive_projection.main(
            derive_projection.parse_args(
                [
                    f"--input_pattern={tmp_path}/projections/*/*.parquet",
                    f"--output_dir={tmp_path}/again",
                ]
            )
        )


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


def _write_index(path, chain, *reaction_ids):
    """Writes an index placing every reaction at the first commit of ``chain``."""
    rows = [
        {
            "reaction_id": reaction_id,
            "sha": chain[0],
            "timestamp": datetime(2021, 7, 28, 14, 7, 53),
            "utc_offset_minutes": -420,
            "dataset_id": "ord_dataset-a",
        }
        for reaction_id in reaction_ids
    ]
    table = pa.Table.from_pylist(rows, schema=first_commits.INDEX_SCHEMA)
    first_commits.write_index(first_commits.Index(table, chain), path)
    return path


class TestFirstCommits:
    def _main(self, tmp_path, *extra):
        derive_projection.main(
            derive_projection.parse_args(
                [
                    f"--input_pattern={tmp_path}/data/*/*.parquet",
                    f"--output_dir={tmp_path}/projections",
                    *extra,
                ]
            )
        )
        return tmp_path / "projections" / "aa" / "ord_dataset-a.parquet"

    def test_fills_first_commit_from_the_index(self, tmp_path):
        _dataset(tmp_path, "aa", "ord_dataset-a")
        index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa")
        output = self._main(tmp_path, f"--first_commits={index}")
        [arrival] = pq.read_table(output).column("first_commit").to_pylist()
        assert arrival["sha"] == "a" * 40
        assert first_commits.read_stamp(output) == (first_commits.VERSION, "a" * 40)

    def test_derives_again_a_projection_derived_without_an_index(self, tmp_path):
        _dataset(tmp_path, "aa", "ord_dataset-a")
        output = self._main(tmp_path)
        first = output.stat().st_ino
        index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa")
        self._main(tmp_path, f"--first_commits={index}")
        assert output.stat().st_ino != first
        assert first_commits.read_stamp(output) is not None

    def test_skips_a_projection_an_extended_index_covers(self, tmp_path):
        _dataset(tmp_path, "aa", "ord_dataset-a")
        index = tmp_path / "index.parquet"
        output = self._main(
            tmp_path, f"--first_commits={_write_index(index, ['a' * 40], 'ord-aa')}"
        )
        first = output.stat().st_ino
        _write_index(index, ["a" * 40, "b" * 40], "ord-aa")
        self._main(tmp_path, f"--first_commits={index}")
        assert output.stat().st_ino == first

    def test_derives_again_a_projection_from_another_history(self, tmp_path):
        _dataset(tmp_path, "aa", "ord_dataset-a")
        index = tmp_path / "index.parquet"
        output = self._main(
            tmp_path, f"--first_commits={_write_index(index, ['a' * 40], 'ord-aa')}"
        )
        _write_index(index, ["c" * 40], "ord-aa")
        self._main(tmp_path, f"--first_commits={index}")
        assert first_commits.read_stamp(output) == (first_commits.VERSION, "c" * 40)

    def test_without_an_index_skips_a_projection_derived_with_one(self, tmp_path):
        _dataset(tmp_path, "aa", "ord_dataset-a")
        index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa")
        output = self._main(tmp_path, f"--first_commits={index}")
        first = output.stat().st_ino
        self._main(tmp_path)
        assert output.stat().st_ino == first

    def test_refuses_a_dataset_the_index_does_not_reach(self, tmp_path):
        _dataset(tmp_path, "aa", "ord_dataset-a")
        index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-zz")
        with pytest.raises(ValueError, match=r"\['ord-aa'\]; walk ord-data"):
            self._main(tmp_path, f"--first_commits={index}")
