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
"""Tests for ord_schema.datasets."""

import pytest

from ord_schema import datasets, message_helpers, parquet
from ord_schema.proto import dataset_pb2, reaction_pb2


class TestLoadAndSaveDataset:
    @pytest.mark.parametrize(
        "suffix",
        [
            ".pbtxt",
            ".pb",
            ".pb.gz",
            ".json",
            ".parquet",
            ".txtpb",
            ".binpb",
            ".binpb.gz",
            ".txtpb.gz",
        ],
    )
    def test_save_dataset(self, suffix, tmp_path):
        dataset = dataset_pb2.Dataset(
            name="n",
            description="d",
            reactions=[reaction_pb2.Reaction(reaction_id="ord-0")],
        )
        path = (tmp_path / f"ds{suffix}").as_posix()
        datasets.save_dataset(dataset, path)
        # For .parquet, exercise the DatasetView entry point callers will use;
        # for other formats, use the generic load_message.
        if suffix == ".parquet":
            loaded = parquet.DatasetView(path)
        else:
            loaded = message_helpers.load_message(path, dataset_pb2.Dataset)
        assert loaded.name == "n"
        assert loaded.description == "d"
        assert list(loaded.reactions) == list(dataset.reactions)

    @pytest.mark.parametrize(
        "suffix",
        [
            ".pbtxt",
            ".pb",
            ".pb.gz",
            ".json",
            ".txtpb",
            ".binpb",
            ".binpb.gz",
            ".txtpb.gz",
        ],
    )
    def test_load_dataset(self, suffix, tmp_path):
        dataset = dataset_pb2.Dataset(
            name="n",
            description="d",
            reactions=[reaction_pb2.Reaction(reaction_id="ord-0")],
        )
        path = (tmp_path / f"ds{suffix}").as_posix()
        datasets.save_dataset(dataset, path)
        loaded = datasets.load_dataset(path)
        assert loaded.name == "n"
        assert loaded.description == "d"
        assert list(loaded.reactions) == list(dataset.reactions)

    def test_load_dataset_parquet_returns_a_view(self, tmp_path):
        """Parquet reads back as a streaming view, not a materialized message."""
        dataset = dataset_pb2.Dataset(
            name="n",
            description="d",
            reactions=[reaction_pb2.Reaction(reaction_id="ord-0")],
        )
        path = (tmp_path / "ds.parquet").as_posix()
        datasets.save_dataset(dataset, path)
        loaded = datasets.load_dataset(path)
        assert isinstance(loaded, parquet.DatasetView)
        # The read-only surface a caller would otherwise have had to change.
        assert loaded.name == "n"
        assert loaded.description == "d"
        assert len(loaded.reactions) == 1
        assert loaded.reactions[0] == dataset.reactions[0]
        assert list(loaded.reactions) == list(dataset.reactions)

    def test_load_dataset_parquet_as_dataset(self, tmp_path):
        """``as_dataset`` materializes, for callers that need the message surface."""
        dataset = dataset_pb2.Dataset(
            name="n",
            description="d",
            reactions=[reaction_pb2.Reaction(reaction_id="ord-0")],
        )
        path = (tmp_path / "ds.parquet").as_posix()
        datasets.save_dataset(dataset, path)
        loaded = datasets.load_dataset(path, as_dataset=True)
        assert isinstance(loaded, dataset_pb2.Dataset)
        assert loaded.SerializeToString(
            deterministic=True
        ) == dataset.SerializeToString(deterministic=True)

    def test_load_dataset_as_dataset_is_a_no_op_for_proto_formats(self, tmp_path):
        """``as_dataset`` must not change what the non-streaming formats return."""
        dataset = dataset_pb2.Dataset(
            name="n",
            description="d",
            reactions=[reaction_pb2.Reaction(reaction_id="ord-0")],
        )
        path = (tmp_path / "ds.pbtxt").as_posix()
        datasets.save_dataset(dataset, path)
        assert datasets.load_dataset(path) == datasets.load_dataset(
            path, as_dataset=True
        )
