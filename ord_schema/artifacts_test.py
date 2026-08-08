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

"""Tests for ord_schema.artifacts."""

import pathlib
from importlib import metadata

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ord_schema import artifacts, parquet
from ord_schema.proto import dataset_pb2, reaction_pb2


def _write(path, metadata):
    schema = pa.schema([pa.field("x", pa.int32())]).with_metadata(metadata)
    pq.write_table(pa.table({"x": [1]}, schema=schema), path)


def _current_metadata(**overrides):
    """Stamps a parent this library considers current, which a chain requires."""
    return _valid_metadata(
        **{"ord.ord_schema_version": metadata.version("ord-schema"), **overrides}
    )


def _valid_metadata(**overrides):
    metadata = {
        "ord.artifact": "view",
        "ord.artifact_version": artifacts.ARTIFACT_VERSION,
        "ord.source_md5": "0" * 32,
        "ord.ord_schema_version": "9.9.9",
        "ord.source_dataset_id": "ord_dataset-1",
    }
    metadata.update(overrides)
    return metadata


# Stamps


def test_stamps_carries_the_current_versions():
    value = artifacts.current_stamps("view", "ord_dataset-1", "abc")
    assert value.artifact == "view"
    assert value.artifact_version == artifacts.ARTIFACT_VERSION
    assert value.ord_schema_version


def test_to_metadata_omits_a_missing_dataset_id():
    metadata = artifacts.to_metadata(artifacts.current_stamps("view", None, "abc"))
    assert "ord.source_dataset_id" not in metadata
    assert metadata["ord.source_md5"] == "abc"


def test_load_stamps_round_trips(tmp_path):
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata())
    stamps = artifacts.load_stamps(path)
    assert stamps.artifact == "view"
    assert stamps.source_dataset_id == "ord_dataset-1"
    assert stamps.source_md5 == "0" * 32


def test_load_stamps_tolerates_a_missing_dataset_id(tmp_path):
    path = tmp_path / "artifact.parquet"
    metadata = _valid_metadata()
    del metadata["ord.source_dataset_id"]
    _write(path, metadata)
    assert artifacts.load_stamps(path).source_dataset_id is None


@pytest.mark.parametrize(
    "missing",
    [
        "ord.artifact",
        "ord.artifact_version",
        "ord.source_md5",
        "ord.ord_schema_version",
    ],
)
def test_load_stamps_names_every_missing_key(tmp_path, missing):
    path = tmp_path / "artifact.parquet"
    metadata = _valid_metadata()
    del metadata[missing]
    _write(path, metadata)
    with pytest.raises(ValueError, match=missing):
        artifacts.load_stamps(path)


def test_load_stamps_rejects_a_file_with_no_metadata(tmp_path):
    path = tmp_path / "plain.parquet"
    pq.write_table(pa.table({"x": [1]}), path)
    with pytest.raises(ValueError, match="not a derived artifact"):
        artifacts.load_stamps(path)


def test_is_artifact_recognizes_only_stamped_files(tmp_path):
    _write(tmp_path / "artifact.parquet", _valid_metadata())
    pq.write_table(pa.table({"x": [1]}), tmp_path / "plain.parquet")
    assert artifacts.is_artifact(tmp_path / "artifact.parquet")
    assert not artifacts.is_artifact(tmp_path / "plain.parquet")
    assert not artifacts.is_artifact(tmp_path / "absent.parquet")


def test_is_artifact_refuses_to_call_an_unreadable_file_a_source(tmp_path):
    # ArrowInvalid subclasses ValueError, so catching ValueError to mean "no stamps"
    # would hand a truncated file to a reader that fails later without naming it.
    path = tmp_path / "truncated.parquet"
    _write(path, _valid_metadata())
    path.write_bytes(path.read_bytes()[:200])
    with pytest.raises(ValueError, match="not readable as Parquet"):
        artifacts.is_artifact(path)


def test_is_current_requires_every_stamp_to_match(tmp_path, monkeypatch):
    # Each of these alone can change a projected value, so each alone must force a
    # rebuild; a check that quietly stopped working would serve obsolete artifacts.
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata(**{"ord.ord_schema_version": "9.9.9"}))
    monkeypatch.setattr(artifacts.metadata, "version", lambda _: "9.9.9")
    assert artifacts.is_current(path, "view", "0" * 32)
    monkeypatch.setattr(artifacts.metadata, "version", lambda _: "9.9.10")
    assert not artifacts.is_current(path, "view", "0" * 32)
    monkeypatch.setattr(artifacts.metadata, "version", lambda _: "9.9.9")
    monkeypatch.setattr(artifacts, "ARTIFACT_VERSION", "99")
    assert not artifacts.is_current(path, "view", "0" * 32)


# Output paths


@pytest.mark.parametrize(
    ("pattern", "expected"),
    [
        ("data/*/*.parquet", "data"),
        ("data/**/*.parquet", "data"),
        ("a/b/c/*.parquet", "a/b/c"),
        # No wildcard: the last component is the file, so the root holds it.
        ("data/aa/one.parquet", "data/aa"),
        ("one.parquet", ""),
    ],
)
def test_glob_root(pattern, expected):
    assert artifacts.glob_root(pattern) == pathlib.PurePath(expected)


def test_output_path_mirrors_the_input_layout():
    assert artifacts.output_path(
        "data/aa/ord_dataset-x.parquet", "data/*/*.parquet", "views"
    ) == pathlib.Path("views/aa/ord_dataset-x.parquet")


def test_output_path_for_an_exact_filename_writes_into_the_directory():
    # A pattern naming one file must land inside output_dir, not become it.
    assert artifacts.output_path(
        "data/aa/one.parquet", "data/aa/one.parquet", "views"
    ) == pathlib.Path("views/one.parquet")


# Driver


def _fake_source(path):
    """Writes a minimal but genuine ORD dataset; derive_tree hashes it via parquet."""
    reaction = reaction_pb2.Reaction(reaction_id="ord-1")
    reaction.identifiers.add(type="REACTION_SMILES", value="C>>CO")
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-1",
            name="test",
            description="test dataset",
            reactions=[reaction],
        ),
        path,
    )


def test_derive_tree_raises_when_nothing_matches(tmp_path):
    with pytest.raises(ValueError, match="no datasets matched"):
        artifacts.derive_tree(
            str(tmp_path / "*.parquet"),
            str(tmp_path / "out"),
            artifact="view",
            write=lambda *args, **kwargs: 0,
        )


def test_derive_tree_writes_one_artifact_per_source(tmp_path):
    for name in ("aa", "bb"):
        (tmp_path / name).mkdir()
        _fake_source(tmp_path / name / "source.parquet")
    written_paths = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        del source, source_md5, source_dataset_id  # Unused.
        written_paths.append(output)
        pathlib.Path(output).write_bytes(b"")
        return 1

    written, skipped, ignored = artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=_write_one,
    )
    assert (written, skipped, ignored) == (2, 0, 0)
    assert {path.parent.name for path in written_paths} == {"aa", "bb"}


def test_derive_tree_hands_the_writer_the_source_and_its_provenance(tmp_path):
    (tmp_path / "aa").mkdir()
    source = tmp_path / "aa" / "source.parquet"
    _fake_source(source)
    seen = []

    def _write_one(path, output, *, source_md5, source_dataset_id):
        view = parquet.DatasetView(path)
        seen.append((source_dataset_id, view.md5() == source_md5))
        pathlib.Path(output).write_bytes(b"")
        return 1

    artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=_write_one,
    )
    assert seen == [("ord_dataset-1", True)]


def test_derive_tree_refuses_to_write_over_its_own_sources(tmp_path):
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    calls = []
    with pytest.raises(ValueError, match="would write over its inputs"):
        artifacts.derive_tree(
            str(tmp_path / "*" / "*.parquet"),
            str(tmp_path),
            artifact="view",
            write=lambda *args, **kwargs: calls.append(args) or 1,
        )
    assert not calls


def test_derive_tree_refuses_to_write_over_a_different_source(tmp_path):
    # data/source.parquet maps onto data/aa/source.parquet, another source this run
    # has not reached. A per-source check passes it and destroys the neighbor.
    _fake_source(tmp_path / "source.parquet")
    (tmp_path / "aa").mkdir()
    victim = tmp_path / "aa" / "source.parquet"
    _fake_source(victim)
    original = victim.read_bytes()
    with pytest.raises(ValueError, match="would write over its inputs"):
        artifacts.derive_tree(
            str(tmp_path / "**" / "*.parquet"),
            str(tmp_path / "aa"),
            artifact="view",
            write=lambda *args, **kwargs: 1,
        )
    assert victim.read_bytes() == original


def test_derive_tree_ignores_matches_that_are_themselves_artifacts(tmp_path):
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    # What a previous run would have left behind inside a recursive pattern's reach.
    _write(tmp_path / "aa" / "derived.parquet", _valid_metadata())
    written_paths = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        del source, source_md5, source_dataset_id  # Unused.
        written_paths.append(output)
        pathlib.Path(output).write_bytes(b"")
        return 1

    assert artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=_write_one,
    ) == (1, 0, 1)
    assert [path.name for path in written_paths] == ["source.parquet"]


def test_derive_tree_skips_current_artifacts_unless_forced(tmp_path, monkeypatch):
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    calls = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        del source, source_md5, source_dataset_id  # Unused.
        calls.append(output)
        pathlib.Path(output).write_bytes(b"")
        return 1

    monkeypatch.setattr(artifacts, "is_current", lambda *args: True)
    assert artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=_write_one,
    ) == (0, 1, 0)
    assert not calls
    assert artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=_write_one,
        force=True,
    ) == (1, 0, 0)
    assert len(calls) == 1


def test_derive_tree_reads_the_named_parent_artifact(tmp_path):
    (tmp_path / "aa").mkdir()
    _write(
        tmp_path / "aa" / "projected.parquet",
        _current_metadata(**{"ord.artifact": "projection", "ord.source_md5": "a" * 32}),
    )
    seen = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        seen.append((pathlib.Path(source).name, source_md5, source_dataset_id))
        pathlib.Path(output).write_bytes(b"")
        return 1

    assert artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=_write_one,
        parent_artifact="projection",
    ) == (1, 0, 0)
    # The source dataset's hash and ID pass through, so a view names the dataset it
    # reflects rather than the artifact it read.
    assert seen == [("projected.parquet", "a" * 32, "ord_dataset-1")]


def test_derive_tree_ignores_a_source_dataset_when_a_parent_artifact_is_named(tmp_path):
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    _write(
        tmp_path / "aa" / "projected.parquet",
        _current_metadata(**{"ord.artifact": "projection"}),
    )
    written_paths = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        del source, source_md5, source_dataset_id  # Unused.
        written_paths.append(output)
        pathlib.Path(output).write_bytes(b"")
        return 1

    assert artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=_write_one,
        parent_artifact="projection",
    ) == (1, 0, 1)
    assert [path.name for path in written_paths] == ["projected.parquet"]


def test_derive_tree_ignores_an_artifact_of_the_wrong_kind(tmp_path):
    # What an earlier run of this same command left in a recursive pattern's reach.
    (tmp_path / "aa").mkdir()
    _write(
        tmp_path / "aa" / "already.parquet", _valid_metadata(**{"ord.artifact": "view"})
    )
    assert artifacts.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="view",
        write=lambda *args, **kwargs: 1,
        parent_artifact="projection",
    ) == (0, 0, 1)


def test_derive_tree_refuses_a_stale_parent(tmp_path, monkeypatch):
    # Passing the parent's hash through carries the source content across the hop but
    # not the version stamps. A view written from a stale projection would stamp itself
    # with the current versions, and the dataset hash it inherits does not change when
    # the projection is rebuilt -- so nothing would ever mark it stale again.
    (tmp_path / "aa").mkdir()
    parent = tmp_path / "aa" / "projected.parquet"
    _write(parent, _current_metadata(**{"ord.artifact": "projection"}))
    monkeypatch.setattr(artifacts, "ARTIFACT_VERSION", "next")
    with pytest.raises(ValueError, match="stale projection"):
        artifacts.derive_tree(
            str(tmp_path / "*" / "*.parquet"),
            str(tmp_path / "out"),
            artifact="view",
            write=lambda *args, **kwargs: 1,
            parent_artifact="projection",
        )


def test_derive_tree_refuses_to_write_over_a_file_it_did_not_derive(tmp_path):
    # Comparing destinations against the run's own parents cannot catch this once a
    # derivation reads one tree and writes another: the source datasets are in neither
    # set, and replacing a corpus that cannot be regenerated is the one unrecoverable
    # mistake this driver can make.
    (tmp_path / "projections").mkdir()
    _write(
        tmp_path / "projections" / "ds.parquet",
        _current_metadata(**{"ord.artifact": "projection"}),
    )
    (tmp_path / "data").mkdir()
    source = tmp_path / "data" / "ds.parquet"
    _fake_source(source)
    before = source.read_bytes()
    with pytest.raises(ValueError, match="would write over its inputs"):
        artifacts.derive_tree(
            str(tmp_path / "projections" / "*.parquet"),
            str(tmp_path / "data"),
            artifact="view",
            write=lambda *args, **kwargs: 1,
            parent_artifact="projection",
        )
    assert source.read_bytes() == before


def test_derive_tree_still_rewrites_its_own_artifacts(tmp_path):
    # The guard above must not stop a re-run: an artifact this driver wrote is exactly
    # what --force is for.
    (tmp_path / "projections").mkdir()
    _write(
        tmp_path / "projections" / "ds.parquet",
        _current_metadata(**{"ord.artifact": "projection"}),
    )
    (tmp_path / "views").mkdir()
    _write(tmp_path / "views" / "ds.parquet", _valid_metadata())
    written, skipped, ignored = artifacts.derive_tree(
        str(tmp_path / "projections" / "*.parquet"),
        str(tmp_path / "views"),
        artifact="view",
        write=lambda *args, **kwargs: 1,
        force=True,
        parent_artifact="projection",
    )
    assert (written, skipped, ignored) == (1, 0, 0)
