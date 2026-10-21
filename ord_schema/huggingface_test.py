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
import io
import re
import types

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from huggingface_hub.errors import EntryNotFoundError
from rdkit import rdBase

from ord_schema import huggingface, message_helpers
from ord_schema.artifacts import (
    base,
    first_commits,
    occurrences,
    pivot,
    projection,
    structures,
)

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
    """Returns repository paths holding a build of every artifact."""
    shard = "aa/ord_dataset-aa.parquet"
    return [
        ".gitattributes",
        "README.md",
        f"projections/{shard}",
        f"structures/{shard}",
        f"first_commits/{shard}",
        f"occurrences/inputs.components/{shard}",
        *(f"pivots/{level}/{shard}" for level in levels),
    ]


def _footer(artifact_schema_version=None, rdkit_version=None):
    """Returns a Parquet file's bytes, stamped as an occurrences artifact."""
    stamps = base.Stamps(
        artifact="occurrences",
        source_dataset_id=None,
        source_md5="aa",
        artifact_schema_version=artifact_schema_version or base.ARTIFACT_SCHEMA_VERSION,
        rdkit_version=rdkit_version or rdBase.rdkitVersion,
    )
    table = pa.table({"x": [1]}).replace_schema_metadata(base.to_metadata(stamps))
    sink = io.BytesIO()
    pq.write_table(table, sink)
    return sink.getvalue()


def _touch(*paths):
    """Creates empty files, and their directories."""
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()


class _FakeHub:
    """Stands in for the Hub's read side: a head commit, its files, the downloads."""

    def __init__(self, monkeypatch, published, head="0123abcd", footer=None):
        self.published = list(published)
        self.head = head
        self.footer = footer or _footer()
        self.listed = []
        self.opened = []
        self.downloads = []
        monkeypatch.setattr(huggingface, "repo_info", self._repo_info)
        monkeypatch.setattr(huggingface, "list_repo_files", self._list_repo_files)
        monkeypatch.setattr(huggingface, "HfFileSystem", lambda: self)
        monkeypatch.setattr(huggingface, "snapshot_download", self._snapshot_download)

    def _repo_info(self, repo_id, *, repo_type, revision):
        assert repo_id == huggingface.ORD_ARTIFACTS_HF_REPO
        assert repo_type == "dataset"
        return types.SimpleNamespace(sha=self.head)

    def _list_repo_files(self, repo_id, *, repo_type, revision):
        self.listed.append(revision)
        return list(self.published)

    def open(self, location, *, block_size):
        self.opened.append(location)
        return io.BytesIO(self.footer)

    def _snapshot_download(self, **kwargs):
        self.downloads.append(kwargs)
        return str(kwargs["local_dir"])


def test_directories_cover_every_artifact():
    # An artifact missing here would be skipped by fetch and upload alike.
    assert set(huggingface.ARTIFACT_DIRECTORIES) == {
        projection.ARTIFACT,
        structures.ARTIFACT,
        first_commits.ARTIFACT,
        pivot.ARTIFACT,
        occurrences.ARTIFACT,
    }


class TestPublishedBuild:
    def test_reads_one_occurrences_footer_at_the_commit(self, monkeypatch):
        hub = _FakeHub(
            monkeypatch, _published(), footer=_footer(artifact_schema_version="7")
        )
        stamps = huggingface.published_build(hub.head, _published())
        assert stamps is not None
        assert stamps.artifact_schema_version == "7"
        assert hub.opened == [
            f"datasets/{huggingface.ORD_ARTIFACTS_HF_REPO}@{hub.head}/"
            "occurrences/inputs.components/aa/ord_dataset-aa.parquet"
        ]

    def test_reports_no_build_without_occurrences(self, monkeypatch):
        published = [p for p in _published() if not p.startswith("occurrences/")]
        hub = _FakeHub(monkeypatch, published)
        assert huggingface.published_build(hub.head, published) is None
        assert hub.opened == []


class TestFetchArtifacts:
    def test_fetches_the_build_at_one_commit(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published())
        fetched = huggingface.fetch_artifacts(tmp_path, revision="main")
        [download] = hub.downloads
        assert download["repo_id"] == huggingface.ORD_ARTIFACTS_HF_REPO
        assert download["local_dir"] == tmp_path
        assert sorted(download["allow_patterns"]) == sorted(
            f"{directory}/**" for directory in huggingface.ARTIFACT_DIRECTORIES.values()
        )
        # The branch is resolved once, so the build checked is the build downloaded.
        assert hub.listed == [hub.head]
        assert all(f"@{hub.head}/" in location for location in hub.opened)
        assert download["revision"] == hub.head
        assert fetched.revision == hub.head
        assert {
            "projection": fetched.projections,
            "structures": fetched.structures,
            "first_commits": fetched.first_commits,
            "pivot": fetched.pivots,
            "occurrences": fetched.occurrences,
        } == {
            artifact: tmp_path / directory
            for artifact, directory in huggingface.ARTIFACT_DIRECTORIES.items()
        }

    def test_patterns_find_what_was_fetched(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        fetched = huggingface.fetch_artifacts(tmp_path / "a[1]")
        projected = fetched.projections / "aa" / "ord_dataset-aa.parquet"
        structured = fetched.structures / "aa" / "ord_dataset-aa.parquet"
        _touch(projected, structured)
        # Escaped, so a bracket in the directory is not a character class.
        assert glob.glob(fetched.projection_pattern, recursive=True) == [str(projected)]
        assert glob.glob(fetched.structures_pattern, recursive=True) == [
            str(structured)
        ]

    def test_deletes_what_the_build_no_longer_holds(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        kept = tmp_path / "projections" / "aa" / "ord_dataset-aa.parquet"
        removed = tmp_path / "projections" / "cc" / "ord_dataset-cc.parquet"
        elsewhere = tmp_path / "notes" / "cc" / "x.parquet"
        _touch(kept, removed, elsewhere)
        huggingface.fetch_artifacts(tmp_path)
        assert kept.exists()
        assert not removed.exists()
        # Only what was fetched is pruned; files outside it are not its to judge.
        assert elsewhere.exists()

    @pytest.mark.parametrize(
        ("artifact_schema_version", "rdkit_version", "message"),
        [
            ("0", None, "holds artifact schema version 0 derived under RDKit"),
            (None, "1999.09.1", "derived under RDKit 1999.09.1"),
        ],
        ids=["schema", "rdkit"],
    )
    def test_refuses_a_build_this_library_cannot_read(
        self, tmp_path, monkeypatch, artifact_schema_version, rdkit_version, message
    ):
        footer = _footer(
            artifact_schema_version=artifact_schema_version, rdkit_version=rdkit_version
        )
        hub = _FakeHub(monkeypatch, _published(), footer=footer)
        held = tmp_path / "projections" / "cc" / "ord_dataset-cc.parquet"
        _touch(held)
        with pytest.raises(LookupError, match=re.escape(message)):
            huggingface.fetch_artifacts(tmp_path)
        # Refused before anything is downloaded or pruned, so a tree holding a build
        # this library reads keeps it.
        assert hub.downloads == []
        assert held.exists()

    def test_refuses_a_missing_artifact(self, tmp_path, monkeypatch):
        published = [p for p in _published() if not p.startswith("structures/")]
        hub = _FakeHub(monkeypatch, published)
        with pytest.raises(LookupError, match="holds no structures artifacts"):
            huggingface.fetch_artifacts(tmp_path)
        assert hub.opened == []
        assert hub.downloads == []

    def test_a_longer_directory_name_is_not_the_artifact(self, tmp_path, monkeypatch):
        published = [p for p in _published() if not p.startswith("projections/")]
        # The same name with a digit more: a prefix of its path, not the artifact.
        published.append("projections0/aa/ord_dataset-aa.parquet")
        _FakeHub(monkeypatch, published)
        with pytest.raises(LookupError, match="holds no projection artifacts"):
            huggingface.fetch_artifacts(tmp_path)

    def test_leaves_a_linked_directory_alone(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        shared = tmp_path / "shared" / "ord_dataset-cc.parquet"
        _touch(shared)
        projections = tmp_path / "tree" / "projections"
        projections.mkdir(parents=True)
        (projections / "cc").symlink_to(shared.parent)
        with pytest.raises(ValueError, match="linked from outside the tree"):
            huggingface.fetch_artifacts(tmp_path / "tree")
        assert shared.exists()

    def test_leaves_a_linked_scope_alone(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        tree = tmp_path / "tree"
        shared = tmp_path / "shared"
        stale = tree / "projections" / "cc" / "x.parquet"
        # Lay the shared storage out as the tree would be beneath the link.
        _touch(shared / "cc" / "x.parquet")
        tree.mkdir()
        (tree / "projections").symlink_to(shared)
        with pytest.raises(ValueError, match="linked from outside the tree"):
            huggingface.fetch_artifacts(tree)
        assert stale.exists()

    def test_prunes_beneath_a_linked_root(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        real = tmp_path / "real"
        removed = real / "projections" / "cc" / "x.parquet"
        _touch(removed)
        (tmp_path / "tree").symlink_to(real)
        # The root is the caller's choice, so a link there is not a foreign directory.
        huggingface.fetch_artifacts(tmp_path / "tree")
        assert not removed.exists()

    def test_removes_a_linked_file_as_a_link(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        target = tmp_path / "elsewhere.parquet"
        _touch(target)
        link = tmp_path / "tree" / "projections" / "cc.parquet"
        link.parent.mkdir(parents=True)
        link.symlink_to(target)
        huggingface.fetch_artifacts(tmp_path / "tree")
        assert not link.is_symlink()
        assert target.exists()

    def test_fetches_named_pivot_levels(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published())
        huggingface.fetch_artifacts(tmp_path, pivot_levels=["inputs.components"])
        [download] = hub.downloads
        assert "pivots/inputs.components/**" in download["allow_patterns"]
        assert "pivots/**" not in download["allow_patterns"]
        assert "pivots/inputs/**" not in download["allow_patterns"]

    def test_refreshes_levels_already_held(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published())
        held = tmp_path / "pivots" / "inputs" / "aa" / "ord_dataset-aa.parquet"
        dropped = tmp_path / "pivots" / "inputs" / "cc" / "ord_dataset-cc.parquet"
        _touch(held, dropped)
        huggingface.fetch_artifacts(tmp_path, pivot_levels=["inputs.components"])
        # A held level left at an older commit is one a Corpus refuses.
        [download] = hub.downloads
        assert "pivots/inputs/**" in download["allow_patterns"]
        assert held.exists()
        assert not dropped.exists()

    def test_leaves_a_dotted_directory_alone(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published())
        kept = tmp_path / "pivots" / ".checkpoints" / "x.parquet"
        _touch(kept)
        huggingface.fetch_artifacts(tmp_path, pivot_levels=["inputs.components"])
        [download] = hub.downloads
        assert "pivots/.checkpoints/**" not in download["allow_patterns"]
        assert kept.exists()

    def test_refuses_missing_pivot_level(self, tmp_path, monkeypatch):
        _FakeHub(monkeypatch, _published())
        with pytest.raises(LookupError, match=r"no pivots for \['workups'\]"):
            huggingface.fetch_artifacts(tmp_path, pivot_levels=["workups"])

    def test_refuses_a_revision_without_a_commit(self, tmp_path, monkeypatch):
        hub = _FakeHub(monkeypatch, _published(), head=None)
        with pytest.raises(LookupError, match="reports no commit for main"):
            huggingface.fetch_artifacts(tmp_path)
        assert hub.listed == []
