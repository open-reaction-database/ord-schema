# Copyright 2026 Open Reaction Database Project Authors
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
"""Tests for ord_schema.parquet_dataset."""

import os

import pyarrow.parquet as pq
import pytest

from ord_schema import parquet_dataset as dataset
from ord_schema.proto import dataset_pb2, reaction_pb2


def _make_reaction(reaction_id: str, conversion: float = 50.0) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.outcomes.add().conversion.value = conversion
    return reaction


def _make_dataset(n: int = 3, *, dataset_id="ord_dataset-test123", name="test", description="desc"):
    return dataset_pb2.Dataset(
        dataset_id=dataset_id,
        name=name,
        description=description,
        reactions=[_make_reaction(f"ord-{i:04d}", conversion=float(i)) for i in range(n)],
    )


def test_round_trip(tmp_path):
    original = _make_dataset(n=5)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    loaded = dataset.read_dataset(path)
    assert loaded.dataset_id == original.dataset_id
    assert loaded.name == original.name
    assert loaded.description == original.description
    assert list(loaded.reactions) == list(original.reactions)


def test_write_rejects_empty_reactions(tmp_path):
    path = os.path.join(tmp_path, "empty.parquet")
    with pytest.raises(ValueError, match="no reactions"):
        dataset.write_dataset(dataset_pb2.Dataset(name="x"), path)


@pytest.mark.parametrize("missing", ["name", "description"])
def test_write_rejects_empty_name_or_description(tmp_path, missing):
    fields = {"name": "n", "description": "d"}
    fields[missing] = ""
    ds = dataset_pb2.Dataset(**fields, reactions=[_make_reaction("ord-0000")])
    path = os.path.join(tmp_path, "ds.parquet")
    with pytest.raises(ValueError, match=missing):
        dataset.write_dataset(ds, path)


def test_read_metadata_skips_reactions(tmp_path):
    original = _make_dataset(n=4)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    metadata = dataset.read_metadata(path)
    assert metadata.dataset_id == original.dataset_id
    assert metadata.name == original.name
    assert metadata.description == original.description
    assert list(metadata.reaction_ids) == [r.reaction_id for r in original.reactions]
    assert len(metadata.reactions) == 0


def test_iter_reactions_all(tmp_path):
    original = _make_dataset(n=3)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    pairs = list(dataset.iter_reactions(path))
    assert [rid for rid, _ in pairs] == [r.reaction_id for r in original.reactions]
    assert [r.outcomes[0].conversion.value for _, r in pairs] == [0.0, 1.0, 2.0]


def test_iter_reactions_filtered(tmp_path):
    original = _make_dataset(n=5)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    # Output order follows file order, not the order of ``reaction_ids``;
    # ord-9999 is absent and silently skipped.
    wanted = {"ord-0003", "ord-0001", "ord-9999"}
    pairs = list(dataset.iter_reactions(path, reaction_ids=wanted))
    assert [rid for rid, _ in pairs] == ["ord-0001", "ord-0003"]


def test_iter_reactions_empty_filter_yields_nothing(tmp_path):
    original = _make_dataset(n=3)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    assert list(dataset.iter_reactions(path, reaction_ids=[])) == []


def test_read_reaction_hit(tmp_path):
    original = _make_dataset(n=5)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    reaction = dataset.read_reaction(path, "ord-0002")
    assert reaction.reaction_id == "ord-0002"
    assert reaction.outcomes[0].conversion.value == 2.0


def test_read_reaction_miss(tmp_path):
    original = _make_dataset(n=2)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    with pytest.raises(KeyError):
        dataset.read_reaction(path, "ord-nope")


def test_row_group_boundaries(tmp_path):
    original = _make_dataset(n=7)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path, row_group_size=3)
    parquet_file = pq.ParquetFile(path)
    # 7 rows with row_group_size=3 -> groups of [3, 3, 1].
    assert parquet_file.num_row_groups == 3


def test_streaming_writer(tmp_path):
    path = os.path.join(tmp_path, "streamed.parquet")
    with dataset.DatasetWriter(
        path, dataset_id="ord_dataset-stream", name="stream", description="streamed", row_group_size=2
    ) as writer:
        for i in range(5):
            writer.write(_make_reaction(f"ord-s{i}", conversion=float(i)))
    loaded = dataset.read_dataset(path)
    assert loaded.dataset_id == "ord_dataset-stream"
    assert loaded.name == "stream"
    assert [r.reaction_id for r in loaded.reactions] == [f"ord-s{i}" for i in range(5)]


def test_streaming_writer_close_is_idempotent(tmp_path):
    path = os.path.join(tmp_path, "streamed.parquet")
    writer = dataset.DatasetWriter(path, name="x", description="y")
    writer.write(_make_reaction("ord-a"))
    writer.close()
    writer.close()  # Second close is a no-op.


def test_writer_requires_name_and_description(tmp_path):
    path = os.path.join(tmp_path, "ds.parquet")
    with pytest.raises(ValueError, match="name"):
        dataset.DatasetWriter(path, name="", description="d")
    with pytest.raises(ValueError, match="description"):
        dataset.DatasetWriter(path, name="n", description="")


def test_footer_omits_empty_dataset_id(tmp_path):
    ds = dataset_pb2.Dataset(name="n", description="d", reactions=[_make_reaction("ord-0000")])
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(ds, path)
    schema = pq.ParquetFile(path).schema_arrow
    keys = {key.decode() for key in (schema.metadata or {})}
    assert keys == {"ord.schema_version", "ord.name", "ord.description"}


def test_read_rejects_unknown_schema_version(tmp_path):
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(_make_dataset(n=1), path)
    # Rewrite the file with a bogus schema_version in the footer.
    table = pq.read_table(path)
    bad_metadata = dict(table.schema.metadata or {})
    bad_metadata[b"ord.schema_version"] = b"999"
    table = table.replace_schema_metadata(bad_metadata)
    bad_path = os.path.join(tmp_path, "bad.parquet")
    pq.write_table(table, bad_path)
    with pytest.raises(ValueError, match="schema version"):
        dataset.read_metadata(bad_path)


@pytest.mark.parametrize("missing_key", ["ord.schema_version", "ord.name", "ord.description"])
def test_read_rejects_missing_required_footer_keys(tmp_path, missing_key):
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(_make_dataset(n=1), path)
    table = pq.read_table(path)
    bad_metadata = {key: value for key, value in (table.schema.metadata or {}).items() if key != missing_key.encode()}
    table = table.replace_schema_metadata(bad_metadata)
    bad_path = os.path.join(tmp_path, "bad.parquet")
    pq.write_table(table, bad_path)
    with pytest.raises(ValueError, match=missing_key):
        dataset.read_metadata(bad_path)


def test_read_metadata_respects_include_reaction_ids_flag(tmp_path):
    original = _make_dataset(n=3)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    metadata = dataset.read_metadata(path, include_reaction_ids=False)
    assert metadata.name == original.name
    assert metadata.description == original.description
    assert list(metadata.reaction_ids) == []


def test_unicode_metadata_round_trips(tmp_path):
    original = dataset_pb2.Dataset(
        name="名前 🧪",
        description="déscription — über",
        dataset_id="ord_dataset-中文",
        reactions=[_make_reaction("ord-0000")],
    )
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path)
    loaded = dataset.read_metadata(path, include_reaction_ids=False)
    assert loaded.name == original.name
    assert loaded.description == original.description
    assert loaded.dataset_id == original.dataset_id


def test_row_group_size_larger_than_dataset(tmp_path):
    original = _make_dataset(n=3)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path, row_group_size=1000)
    parquet_file = pq.ParquetFile(path)
    assert parquet_file.num_row_groups == 1
    assert dataset.read_dataset(path).reactions == original.reactions


def test_multi_row_group_streaming_preserves_order(tmp_path):
    original = _make_dataset(n=2500)
    path = os.path.join(tmp_path, "ds.parquet")
    dataset.write_dataset(original, path, row_group_size=500)
    parquet_file = pq.ParquetFile(path)
    assert parquet_file.num_row_groups == 5
    streamed = [rid for rid, _ in dataset.iter_reactions(path)]
    assert streamed == [r.reaction_id for r in original.reactions]


def test_writer_close_propagates_flush_error(tmp_path):
    path = os.path.join(tmp_path, "ds.parquet")
    writer = dataset.DatasetWriter(path, name="n", description="d")
    writer.write(_make_reaction("ord-0000"))

    class FlushError(RuntimeError):
        pass

    class CloseError(RuntimeError):
        pass

    def fail_flush():
        raise FlushError("flush failed")

    def fail_close():
        raise CloseError("close failed")

    # Force both flush and the underlying writer's close to raise; the flush
    # error must surface and the secondary close error must be swallowed.
    original_parquet_close = writer._writer.close
    writer._flush = fail_flush  # ty: ignore[invalid-assignment]
    writer._writer.close = fail_close  # ty: ignore[invalid-assignment]
    try:
        with pytest.raises(FlushError, match="flush failed"):
            writer.close()
        # Second close is still a no-op even after the failure.
        writer.close()
    finally:
        # Restore so pyarrow's GC-time close doesn't resurrect CloseError.
        writer._writer.close = original_parquet_close  # ty: ignore[invalid-assignment]
        writer._writer.close()
