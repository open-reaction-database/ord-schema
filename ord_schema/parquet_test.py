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
"""Tests for ord_schema.parquet."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ord_schema import parquet as dataset
from ord_schema import updates
from ord_schema.proto import dataset_pb2, reaction_pb2


def _make_reaction(reaction_id: str, conversion: float = 50.0) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.outcomes.add().conversion.value = conversion
    return reaction


def _make_dataset(
    n: int = 3, *, dataset_id="ord_dataset-test123", name="test", description="desc"
):
    return dataset_pb2.Dataset(
        dataset_id=dataset_id,
        name=name,
        description=description,
        reactions=[
            _make_reaction(f"ord-{i:04d}", conversion=float(i)) for i in range(n)
        ],
    )


def test_round_trip(tmp_path):
    original = _make_dataset(n=5)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    loaded = dataset.load_dataset(path)
    assert loaded.dataset_id == original.dataset_id
    assert loaded.name == original.name
    assert loaded.description == original.description
    assert list(loaded.reactions) == list(original.reactions)


def test_save_rejects_empty_reactions(tmp_path):
    path = tmp_path / "empty.parquet"
    with pytest.raises(ValueError, match="no reactions"):
        dataset.save_dataset(dataset_pb2.Dataset(name="x"), path)


@pytest.mark.parametrize("missing", ["name", "description"])
def test_save_rejects_empty_name_or_description(tmp_path, missing):
    fields = {"name": "n", "description": "d"}
    fields[missing] = ""
    ds = dataset_pb2.Dataset(**fields, reactions=[_make_reaction("ord-0000")])
    path = tmp_path / "ds.parquet"
    with pytest.raises(ValueError, match=missing):
        dataset.save_dataset(ds, path)


def test_view_reads_scalars_without_column_data(tmp_path, monkeypatch):
    """Opening a view reads the footer only; the row count comes with it."""
    original = _make_dataset(n=4)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    reads = []
    unpatched = pq.ParquetFile.read_row_group

    def _record(self, row_group, *args, **kwargs):
        reads.append(row_group)
        return unpatched(self, row_group, *args, **kwargs)

    monkeypatch.setattr(pq.ParquetFile, "read_row_group", _record)
    view = dataset.DatasetView(path)
    assert view.dataset_id == original.dataset_id
    assert view.name == original.name
    assert view.description == original.description
    assert len(view.reactions) == 4
    assert list(view.reaction_ids) == []
    assert reads == []


def test_iter_reactions_all(tmp_path):
    original = _make_dataset(n=3)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    pairs = list(dataset.DatasetView(path).iter_reactions())
    assert [rid for rid, _ in pairs] == [r.reaction_id for r in original.reactions]
    assert [r.outcomes[0].conversion.value for _, r in pairs] == [0.0, 1.0, 2.0]


def test_iter_reactions_matches_reactions_stream(tmp_path):
    """``reactions`` is ``iter_reactions()`` with the IDs dropped."""
    original = _make_dataset(n=4)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    view = dataset.DatasetView(path)
    assert [reaction for _, reaction in view.iter_reactions()] == list(view.reactions)


def test_iter_reactions_filtered(tmp_path):
    original = _make_dataset(n=5)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    # Output order follows file order, not the order of ``reaction_ids``;
    # ord-9999 is absent and silently skipped.
    wanted = {"ord-0003", "ord-0001", "ord-9999"}
    pairs = list(dataset.DatasetView(path).iter_reactions(reaction_ids=wanted))
    assert [rid for rid, _ in pairs] == ["ord-0001", "ord-0003"]


def test_iter_reactions_row_group(tmp_path):
    original = _make_dataset(n=5)
    path = tmp_path / "ds.parquet"
    # row_group_size=2 -> row groups of [0,1], [2,3], [4]
    with dataset.DatasetWriter(
        path, name=original.name, description=original.description, row_group_size=2
    ) as writer:
        writer.write_all(original.reactions)
    view = dataset.DatasetView(path)
    assert view.num_row_groups == 3
    assert [rid for rid, _ in view.iter_reactions(row_group=0)] == [
        "ord-0000",
        "ord-0001",
    ]
    assert [rid for rid, _ in view.iter_reactions(row_group=1)] == [
        "ord-0002",
        "ord-0003",
    ]
    assert [rid for rid, _ in view.iter_reactions(row_group=2)] == ["ord-0004"]
    with pytest.raises(IndexError, match="out of range"):
        list(view.iter_reactions(row_group=3))


def test_iter_reactions_empty_filter_raises(tmp_path):
    original = _make_dataset(n=3)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    with pytest.raises(ValueError, match="non-empty"):
        list(dataset.DatasetView(path).iter_reactions(reaction_ids=[]))


def test_get_reaction_hit(tmp_path):
    original = _make_dataset(n=5)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    reaction = dataset.DatasetView(path).get_reaction("ord-0002")
    assert reaction.reaction_id == "ord-0002"
    assert reaction.outcomes[0].conversion.value == 2.0


def test_get_reaction_miss(tmp_path):
    original = _make_dataset(n=2)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    with pytest.raises(KeyError, match="ord-nope"):
        dataset.DatasetView(path).get_reaction("ord-nope")


def test_get_reaction_across_row_groups(tmp_path):
    """Every ID resolves to its own Reaction, including across row group bounds."""
    original = _make_dataset(n=7)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path, row_group_size=3)
    view = dataset.DatasetView(path)
    for reaction in original.reactions:
        assert view.get_reaction(reaction.reaction_id) == reaction


def test_get_reaction_builds_id_index_once(tmp_path, monkeypatch):
    """The reaction_id column is read once, not per lookup."""
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=6), path, row_group_size=2)
    view = dataset.DatasetView(path)
    calls = []
    unpatched = dataset.DatasetView.iter_reaction_ids

    def _record(self):
        calls.append(self.path)
        return unpatched(self)

    monkeypatch.setattr(dataset.DatasetView, "iter_reaction_ids", _record)
    for i in range(6):
        assert view.get_reaction(f"ord-{i:04d}").reaction_id == f"ord-{i:04d}"
    assert len(calls) == 1


def test_row_group_read_shared_by_iteration_and_lookup(tmp_path, monkeypatch):
    """Streaming a row group, indexing into it, and ID lookup share one cached read.

    All three go through the same _RowGroupReader, so once a row group is resolved the
    other two paths must not re-read it.
    """
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=6), path, row_group_size=3)
    view = dataset.DatasetView(path)
    reads = []
    unpatched = pq.ParquetFile.read_row_group

    def _record(self, row_group, *args, **kwargs):
        reads.append(row_group)
        return unpatched(self, row_group, *args, **kwargs)

    monkeypatch.setattr(pq.ParquetFile, "read_row_group", _record)
    assert [rid for rid, _ in view.iter_reactions(row_group=1)] == [
        "ord-0003",
        "ord-0004",
        "ord-0005",
    ]
    assert view.reactions[4].reaction_id == "ord-0004"
    assert view.get_reaction("ord-0005").reaction_id == "ord-0005"
    assert reads == [1]


def test_row_group_boundaries(tmp_path):
    original = _make_dataset(n=7)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path, row_group_size=3)
    parquet_file = pq.ParquetFile(path)
    # 7 rows with row_group_size=3 -> groups of [3, 3, 1].
    assert parquet_file.num_row_groups == 3


def test_streaming_writer(tmp_path):
    path = tmp_path / "streamed.parquet"
    with dataset.DatasetWriter(
        path,
        dataset_id="ord_dataset-stream",
        name="stream",
        description="streamed",
        row_group_size=2,
    ) as writer:
        for i in range(5):
            writer.write(_make_reaction(f"ord-s{i}", conversion=float(i)))
    loaded = dataset.load_dataset(path)
    assert loaded.dataset_id == "ord_dataset-stream"
    assert loaded.name == "stream"
    assert [r.reaction_id for r in loaded.reactions] == [f"ord-s{i}" for i in range(5)]


def test_streaming_writer_close_is_idempotent(tmp_path):
    path = tmp_path / "streamed.parquet"
    writer = dataset.DatasetWriter(path, name="x", description="y")
    writer.write(_make_reaction("ord-a"))
    writer.close()
    writer.close()  # Second close is a no-op.


def test_writer_requires_name_and_description(tmp_path):
    path = tmp_path / "ds.parquet"
    with pytest.raises(ValueError, match="name"):
        dataset.DatasetWriter(path, name="", description="d")
    with pytest.raises(ValueError, match="description"):
        dataset.DatasetWriter(path, name="n", description="")


def test_footer_omits_empty_dataset_id(tmp_path):
    ds = dataset_pb2.Dataset(
        name="n", description="d", reactions=[_make_reaction("ord-0000")]
    )
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(ds, path)
    schema = pq.ParquetFile(path).schema_arrow
    keys = {key.decode() for key in (schema.metadata or {})}
    assert keys == {"ord.schema_version", "ord.name", "ord.description"}


def test_load_rejects_unknown_schema_version(tmp_path):
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=1), path)
    # Rewrite the file with a bogus schema_version in the footer.
    table = pq.read_table(path)
    bad_metadata = dict(table.schema.metadata or {})
    bad_metadata[b"ord.schema_version"] = b"999"
    table = table.replace_schema_metadata(bad_metadata)
    bad_path = tmp_path / "bad.parquet"
    pq.write_table(table, bad_path)
    with pytest.raises(ValueError, match="schema version"):
        dataset.DatasetView(bad_path)


@pytest.mark.parametrize(
    "missing_key", ["ord.schema_version", "ord.name", "ord.description"]
)
def test_load_rejects_missing_required_footer_keys(tmp_path, missing_key):
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=1), path)
    table = pq.read_table(path)
    bad_metadata = {
        key: value
        for key, value in (table.schema.metadata or {}).items()
        if key != missing_key.encode()
    }
    table = table.replace_schema_metadata(bad_metadata)
    bad_path = tmp_path / "bad.parquet"
    pq.write_table(table, bad_path)
    with pytest.raises(ValueError, match=missing_key):
        dataset.DatasetView(bad_path)


def test_unicode_metadata_round_trips(tmp_path):
    original = dataset_pb2.Dataset(
        name="名前 🧪",
        description="déscription — über",
        dataset_id="ord_dataset-中文",
        reactions=[_make_reaction("ord-0000")],
    )
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    loaded = dataset.DatasetView(path)
    assert loaded.name == original.name
    assert loaded.description == original.description
    assert loaded.dataset_id == original.dataset_id


def test_row_group_size_larger_than_dataset(tmp_path):
    original = _make_dataset(n=3)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path, row_group_size=1000)
    parquet_file = pq.ParquetFile(path)
    assert parquet_file.num_row_groups == 1
    assert dataset.load_dataset(path).reactions == original.reactions


def test_multi_row_group_streaming_preserves_order(tmp_path):
    original = _make_dataset(n=2500)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path, row_group_size=500)
    parquet_file = pq.ParquetFile(path)
    assert parquet_file.num_row_groups == 5
    streamed = [rid for rid, _ in dataset.DatasetView(path).iter_reactions()]
    assert streamed == [r.reaction_id for r in original.reactions]


def test_writer_publishes_atomically(tmp_path):
    """Destination must not exist until close(); a sibling temp does."""
    path = tmp_path / "ds.parquet"
    writer = dataset.DatasetWriter(path, name="n", description="d")
    writer.write(_make_reaction("ord-0000"))
    # Mid-write: the destination doesn't exist yet, only the temp does.
    assert not path.exists()
    entries = sorted(q.name for q in tmp_path.iterdir())
    assert "ds.parquet" not in entries
    assert len(entries) == 1
    writer.close()
    # After close: the destination exists and is the only file.
    assert sorted(q.name for q in tmp_path.iterdir()) == ["ds.parquet"]


def test_writer_aborts_on_exception_in_context(tmp_path):
    """An exception in the with-body must leave neither destination nor temp behind."""
    path = tmp_path / "ds.parquet"
    # (abort-on-exception test: the multi-statement body in the context is under test)
    with (  # noqa: PT012
        pytest.raises(RuntimeError, match="boom"),
        dataset.DatasetWriter(path, name="n", description="d") as writer,
    ):
        writer.write(_make_reaction("ord-0000"))
        raise RuntimeError("boom")
    assert sorted(q.name for q in tmp_path.iterdir()) == []


def test_writer_aborts_on_keyboard_interrupt_in_context(tmp_path):
    """KeyboardInterrupt in the with-body triggers _abort and leaves no temp behind."""
    path = tmp_path / "ds.parquet"
    # (abort-on-exception test: the multi-statement body in the context is under test)
    with (  # noqa: PT012
        pytest.raises(KeyboardInterrupt),
        dataset.DatasetWriter(path, name="n", description="d") as writer,
    ):
        writer.write(_make_reaction("ord-0000"))
        raise KeyboardInterrupt
    assert sorted(q.name for q in tmp_path.iterdir()) == []


def test_writer_aborts_preserves_existing_destination(tmp_path):
    """A failed atomic publish must not clobber an existing file at the destination."""
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=2, name="old", description="old desc"), path)
    with path.open("rb") as f:
        original_bytes = f.read()
    # (abort-on-exception test: the multi-statement body in the context is under test)
    with (  # noqa: PT012
        pytest.raises(RuntimeError, match="boom"),
        dataset.DatasetWriter(path, name="new", description="new desc") as writer,
    ):
        writer.write(_make_reaction("ord-new"))
        raise RuntimeError("boom")
    with path.open("rb") as f:
        assert f.read() == original_bytes
    assert sorted(q.name for q in tmp_path.iterdir()) == ["ds.parquet"]


def test_md5_is_deterministic(tmp_path):
    original = _make_dataset(n=4)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    view = dataset.DatasetView(path)
    assert len(view.reactions) == 4
    # Same file hashes the same on a second pass, and across view instances.
    assert view.md5() == view.md5() == dataset.DatasetView(path).md5()


def test_md5_is_decoupled_from_row_group_size(tmp_path):
    """The same content rewritten with different row group sizes hashes the same."""
    original = _make_dataset(n=12)
    path_small = tmp_path / "small.parquet"
    path_large = tmp_path / "large.parquet"
    dataset.save_dataset(original, path_small, row_group_size=3)
    dataset.save_dataset(original, path_large, row_group_size=1000)
    assert (
        dataset.DatasetView(path_small).md5() == dataset.DatasetView(path_large).md5()
    )


def test_md5_changes_when_scalars_change(tmp_path):
    original = _make_dataset(n=2)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    baseline = dataset.DatasetView(path).md5()
    renamed = dataset_pb2.Dataset()
    renamed.CopyFrom(original)
    renamed.name = "different"
    other_path = tmp_path / "other.parquet"
    dataset.save_dataset(renamed, other_path)
    assert dataset.DatasetView(other_path).md5() != baseline


def test_md5_changes_when_a_reaction_changes(tmp_path):
    original = _make_dataset(n=3)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(original, path)
    baseline = dataset.DatasetView(path).md5()
    mutated = dataset_pb2.Dataset()
    mutated.CopyFrom(original)
    mutated.reactions[1].outcomes[0].conversion.value = 999.0
    other_path = tmp_path / "other.parquet"
    dataset.save_dataset(mutated, other_path)
    assert dataset.DatasetView(other_path).md5() != baseline


def test_writer_close_propagates_flush_error(tmp_path):
    path = tmp_path / "ds.parquet"
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


def test_dataset_view_exposes_all_dataset_fields(tmp_path):
    """DatasetView must expose every field declared on the Dataset proto.

    Fails if a new field is added to the proto without a matching attribute on
    DatasetView, so the view stays a drop-in stand-in for validation.
    """
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=1), path)
    view = dataset.DatasetView(path)
    descriptor = dataset_pb2.Dataset.DESCRIPTOR
    assert descriptor is not None  # Type hint.
    missing = {field.name for field in descriptor.fields} - set(dir(view))
    assert not missing, f"DatasetView is missing Dataset fields: {sorted(missing)}"


def test_dataset_view_reactions_is_read_only(tmp_path):
    """Rebinding ``view.reactions`` must raise so the stream can't be swapped."""
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=1), path)
    view = dataset.DatasetView(path)
    with pytest.raises(AttributeError):
        view.reactions = []  # ty: ignore[invalid-assignment]


def test_dataset_view_empty_parquet_is_falsy(tmp_path):
    """DatasetView.reactions reports length 0 / False for an empty Parquet.

    ``validate_dataset``'s "Dataset requires reactions or reaction_ids" warning uses
    ``if not message.reactions`` — this test guards against that branch going dead when
    a bare generator would be truthy.
    """
    path = tmp_path / "empty.parquet"
    # DatasetWriter (unlike save_dataset) permits zero reactions.
    with dataset.DatasetWriter(path, name="n", description="d"):
        pass
    view = dataset.DatasetView(path)
    assert not view.reactions
    assert len(view.reactions) == 0
    assert list(view.reactions) == []


def test_dataset_view_values_round_trip(tmp_path):
    """DatasetView scalars and re-iterated Reactions match the source Dataset.

    Scalars come from the footer; ``.reactions`` opens a fresh stream on each access, so
    two validation passes iterate it twice.
    """
    source = _make_dataset(n=3, dataset_id="ord_dataset-abc", name="n", description="d")
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(source, path)
    view = dataset.DatasetView(path)
    assert view.name == source.name
    assert view.description == source.description
    assert view.dataset_id == source.dataset_id
    assert list(view.reaction_ids) == []  # Not persisted in Parquet.
    assert len(view.reactions) == len(source.reactions)
    assert list(view.reactions) == list(source.reactions)
    # Re-iterating must yield the same sequence (each access re-streams).
    assert list(view.reactions) == list(source.reactions)


def _write_row_groups(path, group_sizes: list[int]) -> list[str]:
    """Writes a Parquet dataset with explicitly sized row groups.

    DatasetWriter flushes on a fixed threshold, so it can only ever produce a short
    group at the end; this drops to pyarrow for the interior short groups and zero-row
    groups that random access also has to resolve correctly.

    Args:
        path: Destination path for the Parquet file.
        group_sizes: Number of reactions to write into each row group, in order.

    Returns:
        The reaction IDs written, in file order.
    """
    schema = dataset._build_schema(name="n", description="d", dataset_id=None)
    reaction_ids = []
    with pq.ParquetWriter(path, schema) as writer:
        for group, size in enumerate(group_sizes):
            ids = [f"ord-{group:02d}{i:02d}" for i in range(size)]
            blobs = [
                _make_reaction(reaction_id).SerializeToString(deterministic=True)
                for reaction_id in ids
            ]
            writer.write_table(
                pa.table({"reaction_id": ids, "reaction": blobs}, schema=schema)
            )
            reaction_ids.extend(ids)
    return reaction_ids


def test_dataset_view_reactions_indexing_matches_source(tmp_path):
    """``view.reactions[i]`` must match the source Dataset across row group bounds."""
    source = _make_dataset(n=7)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(source, path, row_group_size=3)
    view = dataset.DatasetView(path)
    assert [view.reactions[i] for i in range(7)] == list(source.reactions)
    # Descending order exercises cache misses in the other direction.
    assert [view.reactions[i] for i in reversed(range(7))] == list(
        reversed(source.reactions)
    )


def test_dataset_view_reactions_negative_indices(tmp_path):
    source = _make_dataset(n=5)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(source, path, row_group_size=2)
    view = dataset.DatasetView(path)
    assert view.reactions[-1] == source.reactions[-1]
    assert view.reactions[-5] == source.reactions[-5]


@pytest.mark.parametrize("index", [5, 6, -6, -100])
def test_dataset_view_reactions_index_out_of_range(tmp_path, index):
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=5), path, row_group_size=2)
    view = dataset.DatasetView(path)
    with pytest.raises(IndexError, match=f"reaction index {index} out of range"):
        view.reactions[index]


def test_dataset_view_reactions_slices(tmp_path):
    """Slices materialize a list, matching protobuf repeated-field slicing."""
    source = _make_dataset(n=6)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(source, path, row_group_size=2)
    view = dataset.DatasetView(path)
    assert view.reactions[:] == list(source.reactions)
    assert view.reactions[1:4] == list(source.reactions[1:4])
    assert view.reactions[::2] == list(source.reactions[::2])
    assert view.reactions[::-1] == list(source.reactions[::-1])
    assert view.reactions[10:20] == []


def test_dataset_view_reactions_indexing_ragged_row_groups(tmp_path):
    """Index resolution must use footer row counts, not an assumed uniform size.

    Includes a zero-row group, which pyarrow does write and which shares a start
    position with the group that follows it.
    """
    path = tmp_path / "ragged.parquet"
    reaction_ids = _write_row_groups(path, [3, 0, 1, 4])
    # Assert the layout rather than trusting it: a pyarrow that dropped empty tables on
    # write would leave this test passing with its coverage silently gone.
    with pq.ParquetFile(path) as parquet_file:
        metadata = parquet_file.metadata
        counts = [
            metadata.row_group(i).num_rows for i in range(metadata.num_row_groups)
        ]
    assert counts == [3, 0, 1, 4]
    view = dataset.DatasetView(path)
    assert len(view.reactions) == len(reaction_ids)
    assert [view.reactions[i].reaction_id for i in range(len(reaction_ids))] == (
        reaction_ids
    )


def test_dataset_view_indexing_reads_each_row_group_once(tmp_path, monkeypatch):
    """Walking indices in order must cost one row-group read per group, not per row.

    Without the cache this degrades to a full decompress per element, which turns a
    ``for i in range(len(...))`` loop over a large dataset into an apparent hang.
    """
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=10), path, row_group_size=5)
    view = dataset.DatasetView(path)
    reads = []
    opens = []
    unpatched_read = pq.ParquetFile.read_row_group
    unpatched_init = pq.ParquetFile.__init__

    def _record_read(self, row_group, *args, **kwargs):
        reads.append(row_group)
        return unpatched_read(self, row_group, *args, **kwargs)

    def _record_init(self, *args, **kwargs):
        opens.append(1)
        return unpatched_init(self, *args, **kwargs)

    monkeypatch.setattr(pq.ParquetFile, "read_row_group", _record_read)
    monkeypatch.setattr(pq.ParquetFile, "__init__", _record_init)
    for i in range(10):
        view.reactions[i]
    assert reads == [0, 1]
    # One open to build the row-group index, then one per row-group read. A regression
    # that rebuilt the index per access would leave `reads` correct and inflate this.
    assert len(opens) == 3


def test_dataset_view_reactions_empty_dataset(tmp_path):
    """Indexing an empty dataset raises rather than reaching the row-group index.

    A zero-row-group file gives an empty prefix sum, where a bisect would resolve to
    row group -1 and silently read from the end.
    """
    path = tmp_path / "empty.parquet"
    with dataset.DatasetWriter(path, name="n", description="d"):
        pass
    view = dataset.DatasetView(path)
    assert view.reactions[:] == []
    with pytest.raises(IndexError, match="reaction index 0 out of range"):
        view.reactions[0]
    with pytest.raises(IndexError, match="reaction index -1 out of range"):
        view.reactions[-1]


@pytest.mark.parametrize("key", ["0", 1.0, 2.7, None])
def test_dataset_view_reactions_rejects_non_integer_key(tmp_path, key):
    """Non-integer keys raise TypeError instead of being coerced.

    Guards ``operator.index``: ``int(key)`` would truncate ``reactions[n / 2]`` to a
    silently wrong Reaction rather than failing.
    """
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=5), path)
    view = dataset.DatasetView(path)
    with pytest.raises(TypeError):
        view.reactions[key]


def test_dataset_view_detects_replaced_file(tmp_path):
    """A view must refuse to resolve indices against a file replaced underneath it.

    DatasetWriter publishes by renaming onto its destination, so a view outliving a
    republish would otherwise resolve positions against a stale row count.
    """
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=10), path, row_group_size=5)
    view = dataset.DatasetView(path)
    dataset.save_dataset(_make_dataset(n=3), path)
    with pytest.raises(RuntimeError, match="changed since this view was opened"):
        view.reactions[0]


def test_dataset_view_to_proto_round_trips(tmp_path):
    """``to_proto`` returns a real message supporting the protobuf surface."""
    source = _make_dataset(n=3)
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(source, path)
    materialized = dataset.DatasetView(path).to_proto()
    assert isinstance(materialized, dataset_pb2.Dataset)
    assert materialized.SerializeToString(
        deterministic=True
    ) == source.SerializeToString(deterministic=True)


def test_dataset_view_to_proto_carries_scalar_mutations(tmp_path):
    """Scalars assigned on the view must survive materialization.

    ``updates.assign_dataset_id`` accepts a DatasetView and writes ``dataset_id`` in
    place, so materializing from the footer alone would silently drop the assignment.
    """
    source = _make_dataset(n=2, dataset_id="")
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(source, path)
    view = dataset.DatasetView(path)
    assigned = updates.assign_dataset_id(view)
    view.name = "renamed"
    materialized = view.to_proto()
    assert materialized.dataset_id == assigned
    assert materialized.name == "renamed"
    assert list(materialized.reactions) == list(source.reactions)


def test_dataset_view_to_proto_carries_every_scalar_field(tmp_path):
    """Every scalar Dataset field set on the view must survive materialization.

    ``to_proto`` names the scalars one by one, so a field added to the proto would be
    dropped silently. Driving off the descriptor fails here until it is carried.
    """
    path = tmp_path / "ds.parquet"
    dataset.save_dataset(_make_dataset(n=1), path)
    view = dataset.DatasetView(path)
    descriptor = dataset_pb2.Dataset.DESCRIPTOR
    assert descriptor is not None  # Type hint.
    expected = {}
    for field in descriptor.fields:
        if field.label == field.LABEL_REPEATED:
            continue
        assert field.type == field.TYPE_STRING, (
            f"non-string scalar {field.name!r} added to Dataset; extend this test"
        )
        expected[field.name] = f"mutated-{field.name}"
        setattr(view, field.name, expected[field.name])
    assert expected, "Dataset declares no scalar fields; this test is vacuous"
    materialized = view.to_proto()
    assert {name: getattr(materialized, name) for name in expected} == expected
