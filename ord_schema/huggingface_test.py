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
"""Tests for ord_schema.huggingface."""

import pytest
from huggingface_hub.errors import EntryNotFoundError

from ord_schema import huggingface, message_helpers

_DATASET_ID = "ord_dataset-35a5a513f1dd44a3a97c88da99f81a00"


class TestFetchDataset:
    def test_fetch_dataset_prefers_parquet(self, tmp_path, monkeypatch):
        parquet_path = (tmp_path / "ds.parquet").as_posix()
        requested = []

        def fake_download(*, filename, **_):
            requested.append(filename)
            return parquet_path

        monkeypatch.setattr(huggingface, "hf_hub_download", fake_download)
        path = huggingface.fetch_dataset(_DATASET_ID)
        # Parquet is tried first, so a single download resolves it.
        assert requested == [message_helpers.id_filename(f"{_DATASET_ID}.parquet")]
        assert path == parquet_path

    def test_fetch_dataset_falls_back_to_pb_gz(self, tmp_path, monkeypatch):
        pb_path = (tmp_path / "ds.pb.gz").as_posix()
        requested = []

        def fake_download(*, filename, **_):
            requested.append(filename)
            if filename.endswith(".parquet"):
                raise EntryNotFoundError("missing parquet")
            return pb_path

        monkeypatch.setattr(huggingface, "hf_hub_download", fake_download)
        path = huggingface.fetch_dataset(_DATASET_ID)
        # Parquet is attempted before .pb.gz.
        assert requested == [
            message_helpers.id_filename(f"{_DATASET_ID}.parquet"),
            message_helpers.id_filename(f"{_DATASET_ID}.pb.gz"),
        ]
        assert path == pb_path

    def test_fetch_dataset_forwards_cache_options(self, monkeypatch):
        captured = {}

        def fake_download(**kwargs):
            captured.update(kwargs)
            return "/some/path.parquet"

        monkeypatch.setattr(huggingface, "hf_hub_download", fake_download)
        huggingface.fetch_dataset(
            _DATASET_ID, revision="abc123", cache_dir="/c", local_dir="/l"
        )
        assert captured["revision"] == "abc123"
        assert captured["cache_dir"] == "/c"
        assert captured["local_dir"] == "/l"

    def test_fetch_dataset_not_found(self, monkeypatch):
        def fake_download(**_):
            raise EntryNotFoundError("missing")

        monkeypatch.setattr(huggingface, "hf_hub_download", fake_download)
        with pytest.raises(RuntimeError, match="not found"):
            huggingface.fetch_dataset(_DATASET_ID)

    def test_fetch_dataset_invalid_id(self):
        with pytest.raises(ValueError, match="Invalid dataset ID"):
            huggingface.fetch_dataset("not-a-dataset-id")
