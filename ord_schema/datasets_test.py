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

from ord_schema import datasets, message_helpers, parquet_dataset
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
            loaded = parquet_dataset.DatasetView(path)
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

    def test_load_dataset_parquet_warns(self, tmp_path):
        dataset = dataset_pb2.Dataset(
            name="n",
            description="d",
            reactions=[reaction_pb2.Reaction(reaction_id="ord-0")],
        )
        path = (tmp_path / "ds.parquet").as_posix()
        datasets.save_dataset(dataset, path)
        with pytest.warns(UserWarning, match="DatasetView"):
            loaded = datasets.load_dataset(path)
        assert loaded.name == "n"
        assert loaded.description == "d"
        assert list(loaded.reactions) == list(dataset.reactions)
