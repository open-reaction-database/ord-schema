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

"""Tests for ord_schema.artifacts.scripts.derive_arrivals."""

from datetime import datetime

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import arrivals, history, projection
from ord_schema.artifacts.scripts import derive_arrivals
from ord_schema.proto import dataset_pb2, reaction_pb2


def _projection(tmp_path, shard, reaction_id):
    source = tmp_path / "data" / shard / f"ord_dataset-{shard}.parquet"
    source.parent.mkdir(parents=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id=f"ord_dataset-{shard}",
            name="test",
            description="test",
            reactions=[reaction_pb2.Reaction(reaction_id=reaction_id)],
        ),
        source,
    )
    output = tmp_path / "projections" / shard / source.name
    output.parent.mkdir(parents=True)
    projection.write_projection(source, output)
    return output


def _write_index(path, chain, *reaction_ids):
    """Writes an index placing every reaction at the first commit of ``chain``."""
    rows = [
        {
            "reaction_id": reaction_id,
            "sha": chain[0],
            "timestamp": datetime(2021, 7, 28, 14, 7, 53),
            "utc_offset_minutes": -420,
            "dataset_id": "ord_dataset-aa",
        }
        for reaction_id in reaction_ids
    ]
    table = pa.Table.from_pylist(rows, schema=history.INDEX_SCHEMA)
    history.write_index(history.Index(table, chain), path)
    return path


def _main(tmp_path, index, *extra):
    derive_arrivals.main(
        derive_arrivals.parse_args(
            [
                f"--input_pattern={tmp_path}/projections/*/*.parquet",
                f"--output_dir={tmp_path}/arrivals",
                f"--index={index}",
                *extra,
            ]
        )
    )


def test_mirrors_the_projection_tree(tmp_path):
    _projection(tmp_path, "aa", "ord-aa")
    _projection(tmp_path, "bb", "ord-bb")
    index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa", "ord-bb")
    _main(tmp_path, index)
    for shard in ("aa", "bb"):
        output = tmp_path / "arrivals" / shard / f"ord_dataset-{shard}.parquet"
        [row] = pq.read_table(output).to_pylist()
        assert row["reaction_id"] == f"ord-{shard}"
        assert row[arrivals.COLUMN]["sha"] == "a" * 40


def test_skips_an_artifact_an_extended_index_covers(tmp_path):
    _projection(tmp_path, "aa", "ord-aa")
    index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa")
    _main(tmp_path, index)
    output = tmp_path / "arrivals" / "aa" / "ord_dataset-aa.parquet"
    # Inode rather than mtime: the atomic publish renames, so identity is exact.
    first = output.stat().st_ino
    _write_index(index, ["a" * 40, "b" * 40], "ord-aa")
    _main(tmp_path, index)
    assert output.stat().st_ino == first


def test_derives_again_an_artifact_from_another_history(tmp_path):
    _projection(tmp_path, "aa", "ord-aa")
    index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa")
    _main(tmp_path, index)
    _write_index(index, ["c" * 40], "ord-aa")
    _main(tmp_path, index)
    output = tmp_path / "arrivals" / "aa" / "ord_dataset-aa.parquet"
    assert history.read_stamp(output) == (history.VERSION, "c" * 40)


def test_force_rewrites_a_current_artifact(tmp_path):
    _projection(tmp_path, "aa", "ord-aa")
    index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa")
    _main(tmp_path, index)
    output = tmp_path / "arrivals" / "aa" / "ord_dataset-aa.parquet"
    first = output.stat().st_ino
    _main(tmp_path, index, "--force")
    assert output.stat().st_ino != first


def test_refuses_a_projection_the_index_does_not_reach(tmp_path):
    _projection(tmp_path, "aa", "ord-aa")
    index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-zz")
    with pytest.raises(ValueError, match=r"\['ord-aa'\]; walk ord-data"):
        _main(tmp_path, index)


def test_refuses_a_run_that_finds_no_projections(tmp_path):
    projected = _projection(tmp_path, "aa", "ord-aa")
    index = _write_index(tmp_path / "index.parquet", ["a" * 40], "ord-aa")
    with pytest.raises(ValueError, match=r"none of the 1 matches .* are projections"):
        derive_arrivals.main(
            derive_arrivals.parse_args(
                [
                    f"--input_pattern={tmp_path}/data/*/*.parquet",
                    f"--output_dir={tmp_path}/arrivals",
                    f"--index={index}",
                ]
            )
        )
    assert projected.exists()
