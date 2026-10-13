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

import glob
import re
import types

import pytest
from huggingface_hub.errors import EntryNotFoundError
from rdkit import rdBase

from ord_schema import huggingface, message_helpers
from ord_schema.artifacts import base

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


def _published(levels=("inputs", "inputs.components")):
    """Returns repository paths holding this library's build of every artifact."""
    shard = "aa/ord_dataset-aa.parquet"
    pivots = huggingface.current_prefix("pivot")
    return [
        ".gitattributes",
        "README.md",
        *(
            f"{huggingface.current_prefix(artifact)}/{shard}"
            for artifact in ("projection", "structures", "occurrences")
        ),
        *(f"{pivots}/{level}/{shard}" for level in levels),
    ]


def _touch(*paths):
    """Creates empty files, and their directories."""
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()


class _FakeHub:
    """Stands in for the Hub's read side: a head commit, its files, the downloads."""

    def __init__(self, monkeypatch, published, head="0123abcd"):
        self.published = list(published)
        self.head = head
        self.listed = []
        self.downloads = []
        monkeypatch.setattr(huggingface, "repo_info", self._repo_info)
        monkeypatch.setattr(huggingface, "list_repo_files", self._list_repo_files)
        monkeypatch.setattr(huggingface, "snapshot_download", self._snapshot_download)

    def _repo_info(self, repo_id, *, repo_type, revision):
        assert repo_id == huggingface.ORD_ARTIFACTS_HF_REPO
        assert repo_type == "dataset"
        return types.SimpleNamespace(sha=self.head)

    def _list_repo_files(self, repo_id, *, repo_type, revision):
        self.listed.append(revision)
        return list(self.published)

    def _snapshot_download(self, **kwargs):
        self.downloads.append(kwargs)
        return str(kwargs["local_dir"])


class TestArtifactPrefix:
    def test_keys_by_lineage_and_rdkit(self):
        prefix = huggingface.artifact_prefix(
            "structures", "structures=1,projection=2", "2026.03.5"
        )
        assert prefix == "structures/structures-1_projection-2/rdkit-2026.03.5"

    def test_spells_no_hive_partition(self):
        # DuckDB reads a key=value directory as a partition column.
        assert "=" not in huggingface.current_prefix("occurrences")

    def test_current_prefix_names_this_library(self):
        assert huggingface.current_prefix("pivot") == huggingface.artifact_prefix(
            "pivot", base.lineage("pivot"), rdBase.rdkitVersion
        )


class TestFetchArtifacts:
    def test_fetches_this_librarys_build_at_one_commit(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published())
        fetched = huggingface.fetch_artifacts(tmp_path, revision="main")
        [download] = hub.downloads
        assert download["repo_id"] == huggingface.ORD_ARTIFACTS_HF_REPO
        assert download["local_dir"] == tmp_path
        assert sorted(download["allow_patterns"]) == sorted(
            f"{huggingface.current_prefix(artifact)}/**"
            for artifact in huggingface.ARTIFACT_DIRECTORIES
        )
        # The branch is resolved once, so the build checked is the build downloaded.
        assert hub.listed == [hub.head]
        assert download["revision"] == hub.head
        assert fetched.revision == hub.head
        assert fetched.pivots == tmp_path / huggingface.current_prefix("pivot")
        assert fetched.occurrences == tmp_path / huggingface.current_prefix(
            "occurrences"
        )

    def test_patterns_find_what_was_fetched(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        fetched = huggingface.fetch_artifacts(tmp_path / "a[1]")
        shard = fetched.projections / "aa" / "ord_dataset-aa.parquet"
        _touch(shard)
        # Escaped, so a bracket in the directory is not a character class.
        assert glob.glob(fetched.projection_pattern, recursive=True) == [str(shard)]

    def test_deletes_what_the_build_no_longer_holds(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        projections = tmp_path / huggingface.current_prefix("projection")
        kept = projections / "aa" / "ord_dataset-aa.parquet"
        removed = projections / "cc" / "ord_dataset-cc.parquet"
        other_build = tmp_path / "projections/projection-0/rdkit-1/cc/x.parquet"
        _touch(kept, removed, other_build)
        huggingface.fetch_artifacts(tmp_path)
        assert kept.exists()
        assert not removed.exists()
        # Only what was fetched is pruned; another build's files are not its to judge.
        assert other_build.exists()

    def test_leaves_a_linked_directory_alone(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        shared = tmp_path / "shared" / "ord_dataset-cc.parquet"
        _touch(shared)
        projections = tmp_path / "tree" / huggingface.current_prefix("projection")
        projections.mkdir(parents=True)
        (projections / "cc").symlink_to(shared.parent)
        with pytest.raises(ValueError, match="linked from outside the tree"):
            huggingface.fetch_artifacts(tmp_path / "tree")
        assert shared.exists()

    def test_removes_a_linked_file_as_a_link(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        target = tmp_path / "elsewhere.parquet"
        _touch(target)
        link = (
            tmp_path / "tree" / huggingface.current_prefix("projection") / "cc.parquet"
        )
        link.parent.mkdir(parents=True)
        link.symlink_to(target)
        huggingface.fetch_artifacts(tmp_path / "tree")
        assert not link.is_symlink()
        assert target.exists()

    def test_fetches_named_pivot_levels(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published())
        huggingface.fetch_artifacts(tmp_path, pivot_levels=["inputs.components"])
        pivots = huggingface.current_prefix("pivot")
        [download] = hub.downloads
        assert f"{pivots}/inputs.components/**" in download["allow_patterns"]
        assert f"{pivots}/**" not in download["allow_patterns"]
        assert f"{pivots}/inputs/**" not in download["allow_patterns"]

    def test_refreshes_levels_already_held(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published())
        pivots = huggingface.current_prefix("pivot")
        held = tmp_path / pivots / "inputs" / "aa" / "ord_dataset-aa.parquet"
        dropped = tmp_path / pivots / "inputs" / "cc" / "ord_dataset-cc.parquet"
        _touch(held, dropped)
        huggingface.fetch_artifacts(tmp_path, pivot_levels=["inputs.components"])
        # A held level left at an older commit is one a Corpus refuses.
        [download] = hub.downloads
        assert f"{pivots}/inputs/**" in download["allow_patterns"]
        assert held.exists()
        assert not dropped.exists()

    def test_refuses_missing_pivot_level(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        with pytest.raises(LookupError, match=r"no pivots for \['workups'\]"):
            huggingface.fetch_artifacts(tmp_path, pivot_levels=["workups"])

    def test_refuses_a_revision_without_a_commit(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published(), head=None)
        with pytest.raises(LookupError, match="reports no commit for main"):
            huggingface.fetch_artifacts(tmp_path)
        assert hub.listed == []

    def test_refuses_when_no_build_matches(self, tmp_path, monkeypatch):
        other = huggingface.artifact_prefix(
            "projection", base.lineage("projection"), "1999.09.1"
        )
        published = [
            path for path in _published() if not path.startswith("projections/")
        ]
        published.append(f"{other}/ord_dataset-aa.parquet")
        hub = _FakeHub(monkeypatch, published)
        with pytest.raises(LookupError, match=re.escape(f"it holds ['{other}']")):
            huggingface.fetch_artifacts(tmp_path)
        assert hub.downloads == []
