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

"""Tests for ord_schema.artifacts.arrivals."""

import dataclasses
from datetime import datetime

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import arrivals, base, history, projection
from ord_schema.proto import dataset_pb2, reaction_pb2

_EARLY = {
    "sha": "a" * 40,
    "timestamp": datetime(2022, 12, 6, 13, 51, 13),
    "utc_offset_minutes": -300,
    "dataset_id": "ord_dataset-0",
}
_LATE = {
    "sha": "b" * 40,
    "timestamp": datetime(2025, 1, 10, 0, 53, 10),
    "utc_offset_minutes": 780,
    "dataset_id": "ord_dataset-1",
}


def _index(rows, chain=("a" * 40, "b" * 40)):
    """Returns an index holding ``rows``, each a reaction ID and its first commit."""
    table = pa.Table.from_pylist(
        [{"reaction_id": reaction_id, **row} for reaction_id, row in rows.items()],
        schema=history.INDEX_SCHEMA,
    )
    return history.Index(table, list(chain))


@pytest.fixture(name="projected")
def projected_fixture(tmp_path):
    """A projection of three reactions, in the order its source holds them."""
    source = tmp_path / "source.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-1",
            name="test",
            description="test",
            reactions=[
                reaction_pb2.Reaction(reaction_id=reaction_id)
                for reaction_id in ("ord-3", "ord-1", "ord-2")
            ],
        ),
        source,
    )
    output = tmp_path / "projection.parquet"
    projection.write_projection(source, output)
    return output


def test_holds_each_reaction_s_first_commit_in_projection_order(projected, tmp_path):
    output = tmp_path / "arrivals.parquet"
    index = _index({"ord-1": _EARLY, "ord-2": _LATE, "ord-3": _LATE, "ord-9": _EARLY})
    assert arrivals.write_arrivals(projected, output, index=index) == 3
    table = pq.read_table(output)
    assert table.schema.remove_metadata() == arrivals.SCHEMA
    assert table.to_pylist() == [
        {"reaction_id": "ord-3", "first_commit": _LATE},
        {"reaction_id": "ord-1", "first_commit": _EARLY},
        {"reaction_id": "ord-2", "first_commit": _LATE},
    ]


def test_passes_the_projection_s_provenance_through(projected, tmp_path):
    output = tmp_path / "arrivals.parquet"
    index = _index({"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY})
    arrivals.write_arrivals(projected, output, index=index)
    parent = base.load_stamps(projected)
    stamps = base.load_stamps(output)
    assert stamps == dataclasses.replace(parent, artifact=arrivals.ARTIFACT)
    assert history.read_stamp(output) == (history.VERSION, "b" * 40)


def test_refuses_a_reaction_the_index_lacks(projected, tmp_path):
    output = tmp_path / "arrivals.parquet"
    index = _index({"ord-1": _EARLY})
    with pytest.raises(
        ValueError, match=r"2 reactions are not in .* among them \['ord-2', 'ord-3'\]"
    ):
        arrivals.write_arrivals(projected, output, index=index)
    assert not output.exists()


def test_refuses_a_parent_that_is_not_a_projection(projected, tmp_path):
    index = _index({"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY})
    first = tmp_path / "arrivals.parquet"
    arrivals.write_arrivals(projected, first, index=index)
    with pytest.raises(ValueError, match="holds arrivals, not a projection"):
        arrivals.write_arrivals(first, tmp_path / "again.parquet", index=index)


def test_refuses_a_stale_projection(projected, tmp_path, monkeypatch):
    index = _index({"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY})
    monkeypatch.setattr(base, "ARTIFACT_SCHEMA_VERSION", "2")
    with pytest.raises(ValueError, match="is a stale projection"):
        arrivals.write_arrivals(projected, tmp_path / "arrivals.parquet", index=index)


def test_writes_row_groups_of_the_size_asked(projected, tmp_path):
    output = tmp_path / "arrivals.parquet"
    index = _index({"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY})
    arrivals.write_arrivals(projected, output, index=index, row_group_size=2)
    with pq.ParquetFile(output) as parquet_file:
        assert [
            parquet_file.metadata.row_group(group).num_rows
            for group in range(parquet_file.num_row_groups)
        ] == [2, 1]


class TestIsCurrent:
    def test_is_current_under_an_index_extending_the_one_it_read(
        self, projected, tmp_path
    ):
        output = tmp_path / "arrivals.parquet"
        rows = {"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY}
        arrivals.write_arrivals(
            projected, output, index=_index(rows, chain=("a" * 40,))
        )
        md5 = base.load_stamps(projected).source_md5
        assert arrivals.is_current(output, md5, _index(rows))

    def test_is_stale_under_another_history(self, projected, tmp_path):
        output = tmp_path / "arrivals.parquet"
        rows = {"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY}
        arrivals.write_arrivals(projected, output, index=_index(rows))
        md5 = base.load_stamps(projected).source_md5
        assert not arrivals.is_current(output, md5, _index(rows, chain=("c" * 40,)))

    def test_is_stale_for_another_source(self, projected, tmp_path):
        output = tmp_path / "arrivals.parquet"
        rows = {"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY}
        arrivals.write_arrivals(projected, output, index=_index(rows))
        assert not arrivals.is_current(output, "0" * 32, _index(rows))


@pytest.mark.parametrize("size", [0, -1])
def test_refuses_a_row_group_size_no_row_can_fill(projected, tmp_path, size):
    output = tmp_path / "arrivals.parquet"
    index = _index({"ord-1": _EARLY, "ord-2": _EARLY, "ord-3": _EARLY})
    with pytest.raises(ValueError, match=f"row_group_size is {size}"):
        arrivals.write_arrivals(projected, output, index=index, row_group_size=size)
    assert not output.exists()
