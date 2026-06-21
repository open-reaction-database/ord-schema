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

These are the convenience entry points for reading and writing a complete
``Dataset``: ``.parquet`` routes to the (lower-level) ``parquet``
serialization, and other suffixes go through ``message_helpers`` single-message
I/O. For large Parquet datasets prefer the streaming ``parquet``
interfaces (``DatasetView`` / ``iter_reactions``) over loading the whole thing
into memory.
"""

import os
import pathlib
import warnings

from ord_schema import message_helpers, parquet
from ord_schema.proto import dataset_pb2


def load_dataset(filename: str | os.PathLike[str]) -> dataset_pb2.Dataset:
    """Loads a Dataset from disk, dispatching on filename suffix.

    The dataset-level counterpart to ``message_helpers.load_message``:
    ``.parquet`` routes to ``parquet.load_dataset``, which deserializes
    every reaction into memory; other suffixes go through ``load_message``.

    Loading an entire Parquet dataset into memory defeats the format's
    streaming benefits, so this path warns and recommends the streaming
    ``parquet.DatasetView`` loader for large datasets.

    Args:
        filename: Path to a serialized Dataset (``.parquet`` or any suffix
            understood by ``load_message``).

    Returns:
        The in-memory Dataset proto.
    """
    if pathlib.Path(filename).suffix == ".parquet":
        warnings.warn(
            f"Loading the entire Parquet dataset {filename} into memory; for large datasets prefer the "
            "streaming loader ord_schema.parquet.DatasetView (or iter_reactions/load_reaction).",
            UserWarning,
            stacklevel=2,
        )
        return parquet.load_dataset(filename)
    return message_helpers.load_message(filename, dataset_pb2.Dataset)


def save_dataset(
    dataset: dataset_pb2.Dataset, filename: str | os.PathLike[str]
) -> None:
    """Writes a Dataset to disk, dispatching on filename suffix.

    ``.parquet`` routes to ``parquet.save_dataset``; other suffixes go
    through ``message_helpers.save_message``.
    """
    if pathlib.Path(filename).suffix == ".parquet":
        parquet.save_dataset(dataset, filename)
        return
    message_helpers.save_message(dataset, filename)
