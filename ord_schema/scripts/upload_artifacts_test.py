# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for ord_schema.scripts.upload_artifacts."""

import fnmatch
import pathlib
import shutil
import types
from datetime import datetime

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from huggingface_hub import CommitOperationDelete
from rdkit import rdBase

from ord_schema import huggingface, parquet
from ord_schema.artifacts import (
    base,
    first_commits,
    occurrences,
    pivot,
    projection,
    structures,
)
from ord_schema.artifacts.scripts import derive_occurrences, derive_pivots
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.scripts import upload_artifacts
from ord_schema.search import execute, query

_LEVELS = ("inputs", "inputs.components")
_PATHS = ("inputs.components",)
_GROUP_KEYS = {
    "pivot": pivot.META_PIVOT_PATH,
    "occurrences": occurrences.META_OCCURRENCE_PATH,
}


def _history(*reaction_ids):
    """Returns an index placing every reaction at one commit."""
    rows = [
        {
            "reaction_id": reaction_id,
            "sha": "a" * 40,
            "timestamp": datetime(2021, 3, 3, 13, 33, 15),
            "utc_offset_minutes": -420,
            "dataset_id": "ord_dataset-aa",
        }
        for reaction_id in reaction_ids
    ]
    return first_commits.Index(
        pa.Table.from_pylist(rows, schema=first_commits.INDEX_SCHEMA), ["a" * 40]
    )


def _write_artifact(
    path,
    artifact,
    source_md5,
    *,
    group=None,
    artifact_schema_version=None,
    rdkit_version=None,
    stamp_rdkit=True,
    stamp_first_commits=True,
):
    """Writes a one-row Parquet file carrying an artifact's footer stamps.

    A projection also carries the stamp of the index of first commits it was derived
    with, unless ``stamp_first_commits`` is False.
    """
    stamps = base.Stamps(
        artifact=artifact,
        source_dataset_id=None,
        source_md5=source_md5,
        artifact_schema_version=artifact_schema_version or base.ARTIFACT_SCHEMA_VERSION,
        rdkit_version=(rdkit_version or rdBase.rdkitVersion) if stamp_rdkit else None,
    )
    metadata = base.to_metadata(stamps)
    if artifact == projection.ARTIFACT and stamp_first_commits:
        metadata |= _history().stamps()
    if group is not None:
        metadata[_GROUP_KEYS[artifact]] = group
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.table({"x": [1]}).replace_schema_metadata(metadata), path)


def _write_tree(root, sources=("aa", "bb")):
    """Writes a complete artifact tree restating ``sources``, one shard per source."""
    for source in sources:
        name = f"{source}/ord_dataset-{source}.parquet"
        _write_artifact(root / "projections" / name, "projection", source)
        _write_artifact(root / "structures" / name, "structures", source)
        for level in _LEVELS:
            _write_artifact(
                root / "pivots" / level / name, "pivot", source, group=level
            )
        for path in _PATHS:
            _write_artifact(
                root / "occurrences" / path / name, "occurrences", source, group=path
            )


class _FakeHub:
    """Stands in for the Hub: a head commit, its files and build, the commits made."""

    def __init__(
        self,
        monkeypatch,
        published=(),
        head="0123abcd",
        made="4567ef01",
        replaced=None,
    ):
        self.published = list(published)
        self.head = head
        self.made = made
        self.replaced = replaced
        self.listed = []
        self.read = []
        self.commits = []
        monkeypatch.setattr(
            huggingface,
            "repo_info",
            lambda *_, **__: types.SimpleNamespace(sha=self.head),
        )
        monkeypatch.setattr(upload_artifacts, "list_repo_files", self._list_repo_files)
        monkeypatch.setattr(huggingface, "published_build", self._published_build)
        monkeypatch.setattr(upload_artifacts, "create_commit", self._create_commit)

    def _list_repo_files(self, repo_id, *, repo_type, revision):
        self.listed.append(revision)
        return list(self.published)

    def _published_build(self, commit, published):
        self.read.append((commit, list(published)))
        return self.replaced

    def _create_commit(self, repo_id, operations, **kwargs):
        self.commits.append((repo_id, operations, kwargs))
        return types.SimpleNamespace(
            commit_url=f"https://example/commit/{self.made}", oid=self.made
        )


def _upload(root, **kwargs):
    """Uploads ``root``, requiring the levels and paths ``_write_tree`` writes."""
    kwargs.setdefault("pivot_levels", _LEVELS)
    kwargs.setdefault("occurrence_paths", _PATHS)
    return upload_artifacts.upload_artifacts(root, commit_message="derive", **kwargs)


class TestUploadArtifacts:
    @pytest.fixture(autouse=True)
    def _stamps_only(self, monkeypatch):
        # These trees carry stamps and no columns, which no Corpus opens; the real
        # opening is TestCorpusCheck's.
        monkeypatch.setattr(upload_artifacts, "_check_opens", lambda *_, **__: None)

    def test_commits_the_tree_at_the_repository_root(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        stale = "projections/cc/ord_dataset-cc.parquet"
        # A directory whose name extends an artifact's is not that artifact's.
        longer = "projections0/cc/ord_dataset-cc.parquet"
        hub = _FakeHub(monkeypatch, [stale, longer, "README.md"])
        upload = _upload(tmp_path, allow_deletes=True)
        assert upload.commit_url == f"https://example/commit/{hub.made}"
        assert upload.artifact_schema_version == base.ARTIFACT_SCHEMA_VERSION
        assert upload.rdkit_version == rdBase.rdkitVersion
        # 2 sources x (projection, structures, 2 pivot levels, 1 indexed path).
        assert len(upload.uploaded) == 10
        assert "projections/aa/ord_dataset-aa.parquet" in upload.uploaded
        assert "pivots/inputs/bb/ord_dataset-bb.parquet" in upload.uploaded
        # Only what the tree no longer holds goes; the card and the rest are left alone.
        assert upload.deleted == [stale]
        [(repo_id, operations, kwargs)] = hub.commits
        assert repo_id == huggingface.ORD_ARTIFACTS_HF_REPO
        assert kwargs["commit_message"] == "derive"
        # What a reader searches the history by for an earlier build.
        assert kwargs["commit_description"] == (
            f"Artifact schema version {base.ARTIFACT_SCHEMA_VERSION}, "
            f"RDKit {rdBase.rdkitVersion}."
        )
        assert kwargs["repo_type"] == "dataset"
        # Planned against the head it lists, so a commit landing in between fails it.
        assert hub.listed == [hub.head]
        assert kwargs["parent_commit"] == hub.head
        deletes = [
            op.path_in_repo
            for op in operations
            if isinstance(op, CommitOperationDelete)
        ]
        assert deletes == [stale]

    def test_dry_run_commits_nothing(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        hub = _FakeHub(monkeypatch)
        upload = _upload(tmp_path, dry_run=True)
        assert hub.commits == []
        assert upload.commit_url is None
        assert len(upload.uploaded) == 10

    @pytest.mark.parametrize("dry_run", [False, True])
    def test_refuses_to_delete_unless_allowed(self, tmp_path, monkeypatch, dry_run):
        # A tree holding only a new dataset is a complete build of that one dataset,
        # and every check passes it; publishing it as is would delete the rest.
        _write_tree(tmp_path, sources=("cc",))
        published = [
            f"projections/{source}/ord_dataset-{source}.parquet"
            for source in ("aa", "bb")
        ]
        hub = _FakeHub(monkeypatch, published)
        with pytest.raises(
            ValueError,
            match=(
                r"would delete 2 published files the tree does not hold: "
                r"2 under projections; pass --allow_deletes"
            ),
        ):
            _upload(tmp_path, dry_run=dry_run)
        assert hub.commits == []

    def test_dry_run_lists_allowed_deletes(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        stale = "structures/cc/ord_dataset-cc.parquet"
        hub = _FakeHub(monkeypatch, [stale])
        upload = _upload(tmp_path, dry_run=True, allow_deletes=True)
        assert upload.deleted == [stale]
        assert hub.commits == []

    def test_reports_the_build_it_replaces(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        published = ["occurrences/inputs.components/aa/ord_dataset-aa.parquet"]
        replaced = base.Stamps(
            artifact="occurrences",
            source_dataset_id=None,
            source_md5="aa",
            artifact_schema_version="0",
            rdkit_version="1999.09.1",
        )
        hub = _FakeHub(monkeypatch, published, replaced=replaced)
        upload = _upload(tmp_path, dry_run=True)
        # Read at the head the upload is planned against.
        assert hub.read == [(hub.head, published)]
        assert upload.replaced == replaced

    def test_reads_the_replaced_build_before_committing(self, tmp_path, monkeypatch):
        # A failed read after the commit would report an upload that landed as failed.
        _write_tree(tmp_path)
        hub = _FakeHub(monkeypatch)

        def unreadable(commit, published):
            raise OSError("footer unreadable")

        monkeypatch.setattr(huggingface, "published_build", unreadable)
        with pytest.raises(OSError, match="footer unreadable"):
            _upload(tmp_path)
        assert hub.commits == []

    def test_reports_no_commit_when_nothing_changed(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        # The Hub skips a commit that changes nothing and answers with the head.
        hub = _FakeHub(monkeypatch, made="0123abcd")
        upload = _upload(tmp_path)
        assert len(hub.commits) == 1
        assert upload.commit_url is None

    def test_requires_every_level_and_path_by_default(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match=r"pivots lacks \['conditions"):
            upload_artifacts.upload_artifacts(tmp_path, commit_message="derive")

    def test_refuses_a_missing_pivot_level(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        hub = _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match=r"pivots lacks \['workups'\]"):
            _upload(tmp_path, pivot_levels=[*_LEVELS, "workups"])
        assert hub.commits == []

    def test_refuses_an_unrequired_indexed_path(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _FakeHub(monkeypatch)
        with pytest.raises(
            ValueError, match=r"occurrences lacks \[\] and adds \['inputs.components'\]"
        ):
            _upload(tmp_path, occurrence_paths=[])

    def test_refuses_a_pivot_filed_under_another_level(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        # Coverage by directory name still holds; only the footer says otherwise.
        _write_artifact(
            tmp_path / "pivots" / "inputs" / "bb" / "ord_dataset-bb.parquet",
            "pivot",
            "bb",
            group="inputs.components",
        )
        _FakeHub(monkeypatch)
        with pytest.raises(
            ValueError, match=r"filed under inputs but holds inputs\.components"
        ):
            _upload(tmp_path)

    def test_refuses_artifacts_from_different_rdkit(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        for source in ("aa", "bb"):
            _write_artifact(
                tmp_path / "projections" / source / f"ord_dataset-{source}.parquet",
                "projection",
                source,
                rdkit_version="1999.09.1",
            )
        _FakeHub(monkeypatch)
        with pytest.raises(
            ValueError, match="different artifact schema versions or RDKit"
        ):
            _upload(tmp_path)

    def test_refuses_artifacts_of_different_schemas(self, tmp_path, monkeypatch):
        # Structures derived before their projection was, each internally whole.
        _write_tree(tmp_path)
        for source in ("aa", "bb"):
            _write_artifact(
                tmp_path / "structures" / source / f"ord_dataset-{source}.parquet",
                "structures",
                source,
                artifact_schema_version="0",
            )
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match=r"'structures': \('0', '[^']+'\)"):
            _upload(tmp_path)

    def test_refuses_a_tree_spanning_builds(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _write_artifact(
            tmp_path / "structures" / "bb" / "ord_dataset-bb.parquet",
            "structures",
            "bb",
            rdkit_version="1999.09.1",
        )
        hub = _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match="spans more than one build"):
            _upload(tmp_path)
        assert hub.commits == []

    def test_refuses_a_level_missing_a_source(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        (tmp_path / "pivots" / "inputs" / "bb" / "ord_dataset-bb.parquet").unlink()
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match="pivots/inputs lacks 1 and adds 0"):
            _upload(tmp_path)

    def test_refuses_a_dataset_restated_twice(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        # The same source under a second shard: the set of datasets is unchanged.
        _write_artifact(
            tmp_path / "structures" / "cc" / "ord_dataset-cc.parquet",
            "structures",
            "bb",
        )
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match="restates a source dataset another file"):
            _upload(tmp_path)

    def test_refuses_the_wrong_artifact(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _write_artifact(
            tmp_path / "structures" / "bb" / "ord_dataset-bb.parquet",
            "projection",
            "bb",
        )
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match="holds projection, not structures"):
            _upload(tmp_path)

    def test_refuses_a_file_without_an_rdkit_stamp(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        _write_artifact(
            tmp_path / "projections" / "bb" / "ord_dataset-bb.parquet",
            "projection",
            "bb",
            stamp_rdkit=False,
        )
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match="carries no RDKit stamp"):
            _upload(tmp_path)

    def test_refuses_a_projection_derived_without_first_commits(
        self, tmp_path, monkeypatch
    ):
        _write_tree(tmp_path)
        _write_artifact(
            tmp_path / "projections" / "bb" / "ord_dataset-bb.parquet",
            "projection",
            "bb",
            stamp_first_commits=False,
        )
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match=r"bb\.parquet was not derived with a"):
            _upload(tmp_path)

    def test_refuses_a_projection_of_another_walk(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        monkeypatch.setattr(first_commits, "VERSION", "2")
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match="version 2 index of first commits"):
            _upload(tmp_path)

    def test_refuses_a_missing_artifact(self, tmp_path, monkeypatch):
        _write_tree(tmp_path)
        shutil.rmtree(tmp_path / "occurrences")
        _FakeHub(monkeypatch)
        with pytest.raises(ValueError, match="no occurrences artifacts under"):
            _upload(tmp_path)


@pytest.fixture(scope="module")
def derived_tree(tmp_path_factory):
    """Derives every artifact, at every level and path, for two small datasets."""
    root = tmp_path_factory.mktemp("derived")
    for shard, smiles in (("aa", "c1ccncc1"), ("bb", "CCO")):
        reaction = reaction_pb2.Reaction(reaction_id=f"ord-{shard}01")
        component = reaction.inputs["input"].components.add()
        component.identifiers.add(type="SMILES", value=smiles)
        component.reaction_role = reaction_pb2.ReactionRole.REACTANT
        source = root / "data" / shard / f"ord_dataset-{shard}.parquet"
        source.parent.mkdir(parents=True)
        parquet.save_dataset(
            dataset_pb2.Dataset(
                dataset_id=f"ord_dataset-{shard}",
                name=shard,
                description=shard,
                reactions=[reaction],
            ),
            str(source),
        )
        projected = root / "artifacts" / "projections" / shard / source.name
        projected.parent.mkdir(parents=True)
        projection.write_projection(
            source, projected, history=_history("ord-aa01", "ord-bb01")
        )
        structured = root / "artifacts" / "structures" / shard / source.name
        structured.parent.mkdir(parents=True)
        structures.write_structures(projected, structured)
    artifacts = root / "artifacts"
    derive_pivots.main(
        derive_pivots.parse_args(
            [
                f"--input_pattern={artifacts}/projections/*/*.parquet",
                f"--output_dir={artifacts}/pivots",
            ]
        )
    )
    derive_occurrences.main(
        derive_occurrences.parse_args(
            [
                f"--pivots_dir={artifacts}/pivots",
                f"--output_dir={artifacts}/occurrences",
            ]
        )
    )
    return artifacts


class TestCorpusCheck:
    def test_publishes_a_tree_a_corpus_opens(self, derived_tree, monkeypatch):
        _FakeHub(monkeypatch)
        upload = upload_artifacts.upload_artifacts(
            derived_tree, commit_message="derive", dry_run=True
        )
        assert upload.artifact_schema_version == base.ARTIFACT_SCHEMA_VERSION
        assert upload.rdkit_version == rdBase.rdkitVersion

    def test_opens_a_subset_without_requiring_the_rest(
        self, derived_tree, tmp_path, monkeypatch
    ):
        root = tmp_path / "artifacts"
        shutil.copytree(derived_tree, root)
        levels = sorted(pivot.LEVELS)
        paths = sorted(occurrences.PATHS)
        shutil.rmtree(root / "pivots" / levels[-1])
        shutil.rmtree(root / "occurrences" / paths[-1])
        _FakeHub(monkeypatch)
        # A Corpus requiring every level and path would refuse the tree just named.
        upload = upload_artifacts.upload_artifacts(
            root,
            commit_message="derive",
            pivot_levels=levels[:-1],
            occurrence_paths=paths[:-1],
            dry_run=True,
        )
        assert upload.uploaded

    def test_refuses_occurrences_that_miss_a_structure(
        self, derived_tree, tmp_path, monkeypatch
    ):
        # Stamps and columns intact, rows gone: only building the index can tell.
        root = tmp_path / "artifacts"
        shutil.copytree(derived_tree, root)
        target = (
            root / "occurrences" / "inputs.components" / "aa" / "ord_dataset-aa.parquet"
        )
        with pq.ParquetFile(target) as artifact:
            table = artifact.read()
            schema = artifact.schema_arrow
        pq.write_table(table.slice(0, 0).cast(schema), target)
        hub = _FakeHub(monkeypatch)
        with pytest.raises(execute.PairingError, match="the occurrence index reached"):
            upload_artifacts.upload_artifacts(
                root, commit_message="derive", dry_run=True
            )
        assert hub.listed == []

    def test_refuses_structures_the_library_cannot_load(
        self, derived_tree, tmp_path, monkeypatch
    ):
        # A molecule without its fingerprint: only building the library can tell.
        root = tmp_path / "artifacts"
        shutil.copytree(derived_tree, root)
        target = root / "structures" / "aa" / "ord_dataset-aa.parquet"
        with pq.ParquetFile(target) as artifact:
            table = artifact.read()
            schema = artifact.schema_arrow
        column = schema.get_field_index("pattern_fp")
        emptied = pa.array([None] * table.num_rows, type=schema.field(column).type)
        pq.write_table(table.set_column(column, schema.field(column), emptied), target)
        hub = _FakeHub(monkeypatch)
        with pytest.raises(execute.PairingError, match="molecule but no fingerprint"):
            upload_artifacts.upload_artifacts(
                root, commit_message="derive", dry_run=True
            )
        assert hub.listed == []

    def test_refuses_structures_short_of_their_projection(
        self, derived_tree, tmp_path, monkeypatch
    ):
        # Stamps cannot see this: the short artifact keeps its source hash and schema.
        root = tmp_path / "artifacts"
        shutil.copytree(derived_tree, root)
        target = root / "structures" / "aa" / "ord_dataset-aa.parquet"
        with pq.ParquetFile(target) as artifact:
            table = artifact.read()
            schema = artifact.schema_arrow
        pq.write_table(table.slice(0, table.num_rows - 1).cast(schema), target)
        hub = _FakeHub(monkeypatch)
        with pytest.raises(execute.PairingError, match="would join to another dataset"):
            upload_artifacts.upload_artifacts(
                root, commit_message="derive", dry_run=True
            )
        # Refused before the Hub is asked anything.
        assert hub.listed == []


class _DirectoryHub:
    """A Hub backed by a local directory, so an upload can be fetched back."""

    def __init__(self, monkeypatch, directory):
        self.directory = directory
        directory.mkdir()
        sha = types.SimpleNamespace(sha="0123abcd")
        monkeypatch.setattr(huggingface, "repo_info", lambda *_, **__: sha)
        for module in (huggingface, upload_artifacts):
            monkeypatch.setattr(module, "list_repo_files", self._list_repo_files)
        monkeypatch.setattr(upload_artifacts, "create_commit", self._create_commit)
        monkeypatch.setattr(huggingface, "HfFileSystem", lambda: self)
        monkeypatch.setattr(huggingface, "snapshot_download", self._snapshot_download)

    def _list_repo_files(self, *_, **__):
        return [
            path.relative_to(self.directory).as_posix()
            for path in sorted(self.directory.rglob("*"))
            if path.is_file()
        ]

    def _create_commit(self, repo_id, operations, **kwargs):
        for operation in operations:
            target = self.directory / operation.path_in_repo
            if isinstance(operation, CommitOperationDelete):
                target.unlink()
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(operation.path_or_fileobj, target)
        return types.SimpleNamespace(commit_url="https://example/commit/1", oid="1")

    def open(self, location, *, block_size):
        repository = f"datasets/{huggingface.ORD_ARTIFACTS_HF_REPO}@0123abcd/"
        assert location.startswith(repository)
        return (self.directory / location.removeprefix(repository)).open("rb")

    def _snapshot_download(self, *, local_dir, allow_patterns, **_):
        for path in self._list_repo_files():
            if any(fnmatch.fnmatch(path, pattern) for pattern in allow_patterns):
                target = pathlib.Path(local_dir) / path
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(self.directory / path, target)
        return str(local_dir)


def test_a_published_build_fetches_into_a_corpus(derived_tree, tmp_path, monkeypatch):
    # Upload lays the tree out at the repository root and fetch reads its stamps back
    # before downloading it; the two have to meet, and what lands has to be what a
    # Corpus opens.
    _DirectoryHub(monkeypatch, tmp_path / "hub")
    upload_artifacts.upload_artifacts(derived_tree, commit_message="derive")
    fetched = huggingface.fetch_artifacts(tmp_path / "local")
    with execute.Corpus(
        fetched.projection_pattern,
        fetched.structures_pattern,
        pivots_dir=str(fetched.pivots),
        require_pivots=True,
        occurrences_dir=str(fetched.occurrences),
        require_occurrences=True,
        resolver={}.__getitem__,
    ) as corpus:
        table = corpus.search(
            query.Query.model_validate(
                {
                    "where": {
                        "op": "exists",
                        "path": "inputs.components",
                        "where": {
                            "op": "substructure",
                            "path": "smiles",
                            "smarts": "c1ccncc1",
                        },
                    }
                }
            )
        )
    assert table.column("reaction_id").to_pylist() == ["ord-aa01"]


def _report(commit_url=None):
    """Returns what an upload of nothing would report."""
    return upload_artifacts.ArtifactUpload(
        artifact_schema_version=base.ARTIFACT_SCHEMA_VERSION,
        rdkit_version=rdBase.rdkitVersion,
        replaced=None,
        uploaded=[],
        deleted=[],
        commit_url=commit_url,
    )


def test_main_requires_every_level_and_path(monkeypatch):
    calls = []

    def fake_upload(root, **kwargs):
        calls.append((root, kwargs))
        return _report()

    monkeypatch.setattr(upload_artifacts, "upload_artifacts", fake_upload)
    upload_artifacts.main(
        [
            "--root",
            "artifacts",
            "--message",
            "derive",
            "--revision",
            "next",
            "--dry_run",
        ]
    )
    # Levels and paths are left to their defaults, which are everything derived.
    assert calls == [
        (
            "artifacts",
            {
                "commit_message": "derive",
                "revision": "next",
                "dry_run": True,
                "allow_deletes": False,
            },
        )
    ]


def test_main_passes_allow_deletes(monkeypatch):
    calls = []

    def fake_upload(root, **kwargs):
        calls.append(kwargs["allow_deletes"])
        return _report()

    monkeypatch.setattr(upload_artifacts, "upload_artifacts", fake_upload)
    upload_artifacts.main(
        ["--root", "artifacts", "--message", "derive", "--allow_deletes"]
    )
    assert calls == [True]


@pytest.mark.parametrize(
    ("commit_url", "dry_run", "reported"),
    [
        (None, True, []),
        (None, False, ["nothing changed, so nothing was committed"]),
        ("https://example/commit/1", False, ["committed https://example/commit/1"]),
    ],
    ids=["dry", "unchanged", "committed"],
)
def test_main_reports_a_commit_only_when_made(
    monkeypatch, caplog, commit_url, dry_run, reported
):
    monkeypatch.setattr(
        upload_artifacts,
        "upload_artifacts",
        lambda *_, **__: _report(commit_url=commit_url),
    )
    argv = ["--root", "artifacts", "--message", "derive"]
    with caplog.at_level("INFO"):
        upload_artifacts.main([*argv, "--dry_run"] if dry_run else argv)
    assert [
        record.getMessage()
        for record in caplog.records
        if "committed" in record.getMessage()
    ] == reported


def test_message_is_required():
    with pytest.raises(SystemExit, match="2"):
        upload_artifacts.parse_args(["--root", "artifacts"])
