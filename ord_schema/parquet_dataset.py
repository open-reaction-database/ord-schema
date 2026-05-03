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
"""Parquet-based serialization for Dataset messages.

Each Parquet file is a serialized ``Dataset``:

* Rows are ``Reaction`` messages with columns ``reaction_id`` (string) and
  ``reaction`` (serialized wire-format bytes).
* The ``Dataset`` scalar fields (``dataset_id``, ``name``, ``description``)
  are stored in the Parquet footer key-value metadata under ``ord.*`` keys.
* ``Dataset.reaction_ids`` is not persisted — the ``reaction_id`` column is
  the source of truth.

This layout supports random access (by row group) and streaming iteration
without loading the full dataset into memory.
"""

import hashlib
from collections.abc import Iterable, Iterator

import pyarrow as pa
import pyarrow.parquet as pq

from ord_schema.proto import dataset_pb2, reaction_pb2

SCHEMA_VERSION = "1"

_META_PREFIX = "ord."
_META_SCHEMA_VERSION = _META_PREFIX + "schema_version"
_META_DATASET_ID = _META_PREFIX + "dataset_id"
_META_NAME = _META_PREFIX + "name"
_META_DESCRIPTION = _META_PREFIX + "description"

_SCHEMA = pa.schema(
    [
        pa.field("reaction_id", pa.string(), nullable=False),
        pa.field("reaction", pa.binary(), nullable=False),
    ]
)


class DatasetWriter:
    """Streaming writer for a Parquet dataset.

    Reactions are buffered and flushed to a new row group every
    ``row_group_size`` rows, so peak memory stays bounded regardless of the
    total number of reactions written.

    ``name`` and ``description`` are required; ``dataset_id`` is optional
    because it is assigned during the submission process.

    Example:
        with DatasetWriter(path, name="big dataset", description="...") as writer:
            for reaction in source:
                writer.write(reaction)
    """

    def __init__(
        self,
        path: str,
        *,
        name: str,
        description: str,
        dataset_id: str | None = None,
        compression: str = "zstd",
        row_group_size: int = 1000,
    ):
        if not name:
            raise ValueError("DatasetWriter requires a non-empty name")
        if not description:
            raise ValueError("DatasetWriter requires a non-empty description")
        self._row_group_size = row_group_size
        self._schema = _build_schema(name=name, description=description, dataset_id=dataset_id)
        self._writer = pq.ParquetWriter(path, self._schema, compression=compression)
        self._buffer_ids: list[str] = []
        self._buffer_blobs: list[bytes] = []
        self._closed = False

    def write(self, reaction: reaction_pb2.Reaction) -> None:
        """Buffers a Reaction; flushes a row group when the buffer is full."""
        self._buffer_ids.append(reaction.reaction_id)
        self._buffer_blobs.append(reaction.SerializeToString(deterministic=True))
        if len(self._buffer_ids) >= self._row_group_size:
            self._flush()

    def write_all(self, reactions: Iterable[reaction_pb2.Reaction]) -> None:
        """Writes an iterable of Reactions."""
        for reaction in reactions:
            self.write(reaction)

    def close(self) -> None:
        """Flushes the final row group and closes the underlying Parquet writer.

        If ``_flush`` raises, the Parquet writer is still closed best-effort
        (any secondary close error is suppressed) so the original error is the
        one propagated.
        """
        if self._closed:
            return
        self._closed = True
        try:
            self._flush()
        except Exception:
            try:
                self._writer.close()
            except Exception:
                pass
            raise
        self._writer.close()

    def _flush(self) -> None:
        if not self._buffer_ids:
            return
        table = pa.table(
            {"reaction_id": self._buffer_ids, "reaction": self._buffer_blobs},
            schema=self._schema,
        )
        self._writer.write_table(table, row_group_size=self._row_group_size)
        self._buffer_ids.clear()
        self._buffer_blobs.clear()

    def __enter__(self) -> "DatasetWriter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()


def write_dataset(
    dataset: dataset_pb2.Dataset,
    path: str,
    *,
    compression: str = "zstd",
    row_group_size: int = 1000,
) -> None:
    """Writes a Dataset as a Parquet file.

    ``Dataset.reaction_ids`` is not persisted; the ``reaction_id`` column
    derived from ``Dataset.reactions`` is the source of truth.

    Args:
        dataset: Dataset message. Must have at least one Reaction.
        path: Output filename.
        compression: Parquet column compression codec.
        row_group_size: Maximum number of rows per row group. Smaller groups
            give finer-grained random access; larger groups compress better.

    Raises:
        ValueError: If ``dataset.reactions`` is empty or required scalar
            fields (``name``, ``description``) are unset.
    """
    if not dataset.reactions:
        raise ValueError("Dataset has no reactions; Parquet datasets require at least one Reaction.")
    with DatasetWriter(
        path,
        name=dataset.name,
        description=dataset.description,
        dataset_id=dataset.dataset_id or None,
        compression=compression,
        row_group_size=row_group_size,
    ) as writer:
        writer.write_all(dataset.reactions)


class _ReactionStream:
    """Re-iterable Reaction stream with a known length.

    Each iteration opens a fresh read over the backing Parquet file so
    callers that iterate more than once stay memory-bounded. ``__len__``
    and ``__bool__`` come from the footer row count, so emptiness checks
    (e.g., ``if not dataset.reactions``) behave like a list instead of
    always being truthy the way a bare generator would be.
    """

    def __init__(self, path: str, num_reactions: int):
        self._path = path
        self._num_reactions = num_reactions

    def __iter__(self) -> Iterator[reaction_pb2.Reaction]:
        for _, reaction in iter_reactions(self._path):
            yield reaction

    def __len__(self) -> int:
        return self._num_reactions

    def __bool__(self) -> bool:
        return self._num_reactions > 0


class DatasetView:
    """A read-only, streaming view of a Parquet-serialized Dataset.

    Quacks like a ``dataset_pb2.Dataset`` for the read-only attributes used
    during validation: ``name``, ``description``, and ``dataset_id`` come
    from the Parquet footer; ``reaction_ids`` is always empty (Parquet does
    not persist it); and ``reactions`` is a re-iterable ``_ReactionStream``
    that opens a fresh read on each iteration and reports its length from
    the footer, so emptiness/length checks behave like a list. ``reactions``
    is exposed as a read-only property so accidental rebinding raises.
    """

    def __init__(self, path: str):
        self._path = path
        with pq.ParquetFile(path) as parquet_file:
            scalars = _dataset_from_metadata(parquet_file.schema_arrow.metadata)
            num_rows = parquet_file.metadata.num_rows
        self.name = scalars.name
        self.description = scalars.description
        self.dataset_id = scalars.dataset_id
        self.reaction_ids: list[str] = []
        self._reactions = _ReactionStream(path, num_rows)

    @property
    def path(self) -> str:
        """Path to the backing Parquet file (read-only)."""
        return self._path

    @property
    def reactions(self) -> _ReactionStream:
        return self._reactions


def read_dataset(path: str) -> dataset_pb2.Dataset:
    """Reads a full Dataset from a Parquet file.

    All reactions are deserialized into memory; prefer ``iter_reactions`` or
    ``read_reaction`` for large datasets.
    """
    with pq.ParquetFile(path) as parquet_file:
        dataset = _dataset_from_metadata(parquet_file.schema_arrow.metadata)
        for _, reaction in _iter_reactions(parquet_file, filter_set=None):
            dataset.reactions.append(reaction)
    return dataset


def read_metadata(path: str) -> dataset_pb2.Dataset:
    """Reads Dataset scalar fields (``name``, ``description``, ``dataset_id``)
    from the Parquet footer. No column data is read. The returned ``Dataset``
    has no ``reactions`` or ``reaction_ids`` populated.
    """
    with pq.ParquetFile(path) as parquet_file:
        return _dataset_from_metadata(parquet_file.schema_arrow.metadata)


def iter_reactions(
    path: str,
    reaction_ids: Iterable[str] | None = None,
) -> Iterator[tuple[str, reaction_pb2.Reaction]]:
    """Yields ``(reaction_id, Reaction)`` pairs from a Parquet dataset.

    Args:
        path: Parquet file path.
        reaction_ids: Optional iterable of reaction IDs to restrict iteration
            to. IDs not present in the file are silently skipped. Filtering
            is applied row-wise; row groups are still read in full. Pass
            ``None`` (the default) to iterate all reactions; an explicitly
            empty iterable raises ``ValueError`` since it is almost always a
            bug.
    """
    if reaction_ids is None:
        filter_set = None
    else:
        filter_set = set(reaction_ids)
        if not filter_set:
            raise ValueError("reaction_ids must be non-empty when provided; pass None to iterate all")
    with pq.ParquetFile(path) as parquet_file:
        yield from _iter_reactions(parquet_file, filter_set=filter_set)


def _iter_reactions(
    parquet_file: pq.ParquetFile,
    *,
    filter_set: set[str] | None,
) -> Iterator[tuple[str, reaction_pb2.Reaction]]:
    for batch in parquet_file.iter_batches(columns=["reaction_id", "reaction"]):
        ids = batch.column("reaction_id").to_pylist()
        blobs = batch.column("reaction").to_pylist()
        for reaction_id, blob in zip(ids, blobs, strict=True):
            if filter_set is not None and reaction_id not in filter_set:
                continue
            yield reaction_id, reaction_pb2.Reaction.FromString(blob)


def streaming_md5(path: str) -> tuple[str, int]:
    """Returns ``(md5_hexdigest, num_reactions)`` for a Parquet dataset.

    Streams the file row-group at a time so peak memory stays bounded. The
    hash is fed in this order:

    * ``name=<value>\\n``, ``description=<value>\\n``, ``dataset_id=<value>\\n``
      (omitted when unset) from the footer scalars.
    * For each Reaction in iteration order: the 8-byte big-endian length of
      the serialized Reaction blob, then the blob bytes. The reaction_id is
      not hashed separately since it is already inside the blob (field 1 of
      ``Reaction``); the length prefix disambiguates blob boundaries.

    Different from ``md5(Dataset.SerializeToString(deterministic=True))`` —
    re-ingesting an existing ``.pb.gz`` dataset as ``.parquet`` will look
    like a content change once. That is correct: the on-disk bytes really
    did change, and the hash is only a cheap "did the content change since I
    last saw it?" indicator. The scheme here is decoupled from the Parquet
    file layout (row group sizes, compression) so the same logical content
    rewritten with different writer settings still hashes the same.
    """
    hasher = hashlib.md5()
    metadata = read_metadata(path)
    if metadata.name:
        hasher.update(f"name={metadata.name}\n".encode())
    if metadata.description:
        hasher.update(f"description={metadata.description}\n".encode())
    if metadata.dataset_id:
        hasher.update(f"dataset_id={metadata.dataset_id}\n".encode())
    num_reactions = 0
    with pq.ParquetFile(path) as parquet_file:
        for batch in parquet_file.iter_batches(columns=["reaction"]):
            for blob in batch.column("reaction").to_pylist():
                hasher.update(len(blob).to_bytes(8, "big"))
                hasher.update(blob)
                num_reactions += 1
    return hasher.hexdigest(), num_reactions


def iter_reaction_ids(path: str) -> Iterator[str]:
    """Yields ``reaction_id`` values from a Parquet dataset without decoding Reaction blobs.

    Reads only the ``reaction_id`` column row-group at a time, so this is
    cheap even for very large files. Iteration order matches ``iter_reactions``.
    """
    with pq.ParquetFile(path) as parquet_file:
        for batch in parquet_file.iter_batches(columns=["reaction_id"]):
            yield from batch.column("reaction_id").to_pylist()


def read_reaction(path: str, reaction_id: str) -> reaction_pb2.Reaction:
    """Reads a single Reaction by ID.

    Passes ``reaction_id`` as a Parquet predicate, which lets row groups be
    pruned via per-group min/max statistics **when** ``reaction_id`` values
    are clustered within groups (e.g., sorted). For unsorted writes every
    row group's range typically covers the query, so this falls back to a
    full scan of the ``reaction_id``/``reaction`` columns.

    Raises:
        KeyError: If ``reaction_id`` is not present in the file.
    """
    table = pq.read_table(
        path,
        columns=["reaction_id", "reaction"],
        filters=[("reaction_id", "=", reaction_id)],
    )
    ids = table.column("reaction_id").to_pylist()
    blobs = table.column("reaction").to_pylist()
    for candidate_id, blob in zip(ids, blobs, strict=True):
        if candidate_id == reaction_id:
            return reaction_pb2.Reaction.FromString(blob)
    raise KeyError(reaction_id)


def _build_schema(*, name: str, description: str, dataset_id: str | None) -> pa.Schema:
    metadata = {
        _META_SCHEMA_VERSION: SCHEMA_VERSION,
        _META_NAME: name,
        _META_DESCRIPTION: description,
    }
    if dataset_id:
        metadata[_META_DATASET_ID] = dataset_id
    return _SCHEMA.with_metadata(metadata)


def _dataset_from_metadata(raw_metadata: dict[bytes, bytes] | None) -> dataset_pb2.Dataset:
    metadata = raw_metadata or {}

    def _get(key: str) -> str | None:
        value = metadata.get(key.encode())
        return value.decode() if value is not None else None

    version = _get(_META_SCHEMA_VERSION)
    if version is None:
        raise ValueError(f"missing required Parquet footer key: {_META_SCHEMA_VERSION!r}")
    if version != SCHEMA_VERSION:
        raise ValueError(f"unsupported Parquet dataset schema version: {version!r}")
    name = _get(_META_NAME)
    description = _get(_META_DESCRIPTION)
    if not name:
        raise ValueError(f"missing or empty required Parquet footer key: {_META_NAME!r}")
    if not description:
        raise ValueError(f"missing or empty required Parquet footer key: {_META_DESCRIPTION!r}")
    dataset = dataset_pb2.Dataset(name=name, description=description)
    dataset_id = _get(_META_DATASET_ID)
    if dataset_id:
        dataset.dataset_id = dataset_id
    return dataset
