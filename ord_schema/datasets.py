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
"""Whole-dataset file I/O, dispatching on filename suffix.

The entry point for reading and writing a complete ``Dataset``, whatever it is
serialized as. ``.parquet`` reads back as a ``parquet.DatasetView``, which
streams from the file instead of materializing it; every other suffix reads back
as a ``dataset_pb2.Dataset`` through ``message_helpers`` single-message I/O.
"""

import os
import pathlib
from typing import Literal, overload

from ord_schema import message_helpers, parquet
from ord_schema.proto import dataset_pb2


@overload
def load_dataset(
    filename: str | os.PathLike[str], *, as_dataset: Literal[False] = False
) -> dataset_pb2.Dataset | parquet.DatasetView: ...


@overload
def load_dataset(
    filename: str | os.PathLike[str], *, as_dataset: Literal[True]
) -> dataset_pb2.Dataset: ...


def load_dataset(
    filename: str | os.PathLike[str], *, as_dataset: bool = False
) -> dataset_pb2.Dataset | parquet.DatasetView:
    """Loads a Dataset from disk, dispatching on filename suffix.

    A Parquet dataset reads back as a ``parquet.DatasetView``, which takes its
    scalars and row count from the file footer and reads Reactions on demand.
    The view is a read-only stand-in for the message: iteration, ``len``,
    emptiness, indexing, and slicing over ``reactions`` behave the same, so code
    that only reads needs no change. Materializing every Reaction is the
    exception rather than the default, which matters because published datasets
    run to millions of rows.

    Other suffixes have no streaming form and read back as a real message.

    Args:
        filename: Path to a serialized Dataset (``.parquet`` or any suffix
            understood by ``message_helpers.load_message``).
        as_dataset: If True, always return a ``dataset_pb2.Dataset``,
            materializing a Parquet dataset in full. Use this when the caller
            needs the protobuf surface -- serialization, JSON conversion,
            mutation -- and cannot know which format it was handed; a caller
            that knows it holds a view can call ``to_proto`` instead.

    Returns:
        A ``parquet.DatasetView`` for ``.parquet`` unless ``as_dataset`` is set,
        and a ``dataset_pb2.Dataset`` otherwise.
    """
    if pathlib.Path(filename).suffix == ".parquet":
        view = parquet.DatasetView(filename)
        return view.to_proto() if as_dataset else view
    return message_helpers.load_message(filename, dataset_pb2.Dataset)


def save_dataset(
    dataset: dataset_pb2.Dataset, filename: str | os.PathLike[str]
) -> None:
    """Writes a Dataset to disk, dispatching on filename suffix.

    ``.parquet`` routes to ``parquet.save_dataset``; other suffixes go through
    ``message_helpers.save_message``.

    Args:
        dataset: Dataset message to write.
        filename: Output path; its suffix selects the serialization.
    """
    if pathlib.Path(filename).suffix == ".parquet":
        parquet.save_dataset(dataset, filename)
        return
    message_helpers.save_message(dataset, filename)
