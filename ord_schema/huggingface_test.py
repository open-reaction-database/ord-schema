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
import shutil
import types

import pyarrow as pa
import pyarrow.parquet as pq
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


def _write_artifact(
    path,
    artifact,
    source_md5,
    *,
    lineage=None,
    rdkit_version=None,
    stamp_rdkit=True,
):
    """Writes a one-row Parquet file carrying an artifact's footer stamps."""
    stamps = base.Stamps(
        artifact=artifact,
        source_dataset_id=None,
        source_md5=source_md5,
        ord_schema_version="0",
        artifact_lineage=lineage or base.lineage(artifact),
        rdkit_version=(rdkit_version or rdBase.rdkitVersion) if stamp_rdkit else None,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.table({"x": [1]}).replace_schema_metadata(base.to_metadata(stamps))
    pq.write_table(table, path)


def _write_tree(root, sources=("aa", "bb"), levels=("inputs", "inputs.components")):
    """Writes a complete artifact tree restating ``sources``, one shard per source."""
    for source in sources:
        name = f"{source}/ord_dataset-{source}.parquet"
        _write_artifact(root / "projections" / name, "projection", source)
        _write_artifact(root / "structures" / name, "structures", source)
        for level in levels:
            _write_artifact(root / "pivots" / level / name, "pivot", source)
        _write_artifact(
            root / "occurrences" / "inputs.components" / name, "occurrences", source
        )


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
    def test_fetches_this_librarys_build(self, tmp_path, monkeypatch):
        captured = {}

        def fake_snapshot(**kwargs):
            captured.update(kwargs)
            return str(tmp_path)

        monkeypatch.setattr(
            huggingface, "list_repo_files", lambda *_, **__: _published()
        )
        monkeypatch.setattr(huggingface, "snapshot_download", fake_snapshot)
        fetched = huggingface.fetch_artifacts(tmp_path, revision="abc123")
        assert captured["repo_id"] == huggingface.ORD_ARTIFACTS_HF_REPO
        assert captured["revision"] == "abc123"
        assert captured["local_dir"] == tmp_path
        assert sorted(captured["allow_patterns"]) == sorted(
            f"{huggingface.current_prefix(artifact)}/**"
            for artifact in huggingface.ARTIFACT_DIRECTORIES
        )
        assert fetched.pivots == tmp_path / huggingface.current_prefix("pivot")
        assert fetched.occurrences == tmp_path / huggingface.current_prefix(
            "occurrences"
        )

    def test_patterns_find_what_was_fetched(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            huggingface, "list_repo_files", lambda *_, **__: _published()
        )
        monkeypatch.setattr(huggingface, "snapshot_download", lambda **_: str(tmp_path))
        fetched = huggingface.fetch_artifacts(tmp_path / "a[1]")
        shard = fetched.projections / "aa" / "ord_dataset-aa.parquet"
        shard.parent.mkdir(parents=True)
        shard.touch()
        # Escaped, so a bracket in the directory is not a character class.
        assert glob.glob(fetched.projection_pattern, recursive=True) == [str(shard)]

    def test_fetches_named_pivot_levels(self, tmp_path, monkeypatch):
        captured = {}
        monkeypatch.setattr(
            huggingface, "list_repo_files", lambda *_, **__: _published()
        )
        monkeypatch.setattr(
            huggingface, "snapshot_download", lambda **kwargs: captured.update(kwargs)
        )
        huggingface.fetch_artifacts(tmp_path, pivot_levels=["inputs.components"])
        pivots = huggingface.current_prefix("pivot")
        assert f"{pivots}/inputs.components/**" in captured["allow_patterns"]
        assert f"{pivots}/**" not in captured["allow_patterns"]

    def test_refuses_missing_pivot_level(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            huggingface, "list_repo_files", lambda *_, **__: _published()
        )
        with pytest.raises(LookupError, match=r"no pivots for \['workups'\]"):
            huggingface.fetch_artifacts(tmp_path, pivot_levels=["workups"])

    def test_refuses_when_no_build_matches(self, tmp_path, monkeypatch):
        other = huggingface.artifact_prefix(
            "projection", base.lineage("projection"), "1999.09.1"
        )
        published = [
            path for path in _published() if not path.startswith("projections/")
        ]
        published.append(f"{other}/aa/ord_dataset-aa.parquet")
        monkeypatch.setattr(huggingface, "list_repo_files", lambda *_, **__: published)
        with pytest.raises(LookupError, match=re.escape(f"it holds ['{other}']")):
            huggingface.fetch_artifacts(tmp_path)


class TestUploadArtifacts:
    def _commit(self, monkeypatch, published=()):
        """Fakes the repository; returns the commits made."""
        commits = []

        def fake_commit(repo_id, operations, **kwargs):
            commits.append((repo_id, operations, kwargs))
            return types.SimpleNamespace(commit_url="https://example/commit/1")

        monkeypatch.setattr(
            huggingface, "list_repo_files", lambda *_, **__: list(published)
        )
        monkeypatch.setattr(huggingface, "create_commit", fake_commit)
        return commits

    def test_commits_each_file_under_its_build(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        projections = huggingface.current_prefix("projection")
        stale = f"{projections}/cc/ord_dataset-cc.parquet"
        elsewhere = "projections/projection-0/rdkit-1999.09.1/cc/ord_dataset-cc.parquet"
        commits = self._commit(monkeypatch, published=[stale, elsewhere, "README.md"])
        upload = huggingface.upload_artifacts(tmp_path, commit_message="derive")
        assert upload.commit_url == "https://example/commit/1"
        assert upload.prefixes == {
            artifact: huggingface.current_prefix(artifact)
            for artifact in huggingface.ARTIFACT_DIRECTORIES
        }
        # 2 sources x (projection, structures, 2 pivot levels, 1 indexed path).
        assert len(upload.uploaded) == 10
        assert f"{projections}/aa/ord_dataset-aa.parquet" in upload.uploaded
        pivots = huggingface.current_prefix("pivot")
        assert f"{pivots}/inputs/bb/ord_dataset-bb.parquet" in upload.uploaded
        # Only what this build no longer holds goes; another build is left alone.
        assert upload.deleted == [stale]
        [(repo_id, operations, kwargs)] = commits
        assert repo_id == huggingface.ORD_ARTIFACTS_HF_REPO
        assert kwargs["commit_message"] == "derive"
        assert kwargs["repo_type"] == "dataset"
        deletes = [
            op.path_in_repo
            for op in operations
            if isinstance(op, huggingface.CommitOperationDelete)
        ]
        assert deletes == [stale]

    def test_dry_run_commits_nothing(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        commits = self._commit(monkeypatch)
        upload = huggingface.upload_artifacts(
            tmp_path, commit_message="derive", dry_run=True
        )
        assert commits == []
        assert upload.commit_url is None
        assert len(upload.uploaded) == 10

    def test_refuses_a_tree_spanning_builds(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _write_artifact(
            tmp_path / "structures" / "bb" / "ord_dataset-bb.parquet",
            "structures",
            "bb",
            rdkit_version="1999.09.1",
        )
        commits = self._commit(monkeypatch)
        with pytest.raises(ValueError, match="spans more than one build"):
            huggingface.upload_artifacts(tmp_path, commit_message="derive")
        assert commits == []

    def test_refuses_a_level_missing_a_source(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        (tmp_path / "pivots" / "inputs" / "bb" / "ord_dataset-bb.parquet").unlink()
        self._commit(monkeypatch)
        with pytest.raises(ValueError, match="pivots/inputs lacks 1 and adds 0"):
            huggingface.upload_artifacts(tmp_path, commit_message="derive")

    def test_refuses_the_wrong_artifact(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _write_artifact(
            tmp_path / "structures" / "bb" / "ord_dataset-bb.parquet",
            "projection",
            "bb",
        )
        self._commit(monkeypatch)
        with pytest.raises(ValueError, match="holds projection, not structures"):
            huggingface.upload_artifacts(tmp_path, commit_message="derive")

    def test_refuses_a_file_without_an_rdkit_stamp(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _write_artifact(
            tmp_path / "projections" / "bb" / "ord_dataset-bb.parquet",
            "projection",
            "bb",
            stamp_rdkit=False,
        )
        self._commit(monkeypatch)
        with pytest.raises(ValueError, match="carries no RDKit stamp"):
            huggingface.upload_artifacts(tmp_path, commit_message="derive")

    def test_refuses_a_missing_artifact(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        shutil.rmtree(tmp_path / "occurrences")
        self._commit(monkeypatch)
        with pytest.raises(ValueError, match="no occurrences artifacts under"):
            huggingface.upload_artifacts(tmp_path, commit_message="derive")
