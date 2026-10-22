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

"""Tests for ord_schema.artifacts.base."""

import collections
import dataclasses
import glob
import hashlib
import importlib.metadata
import json
import pathlib
from datetime import datetime

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import rdBase

from ord_schema import message_helpers, parquet
from ord_schema.artifacts import (
    arrivals,
    base,
    history,
    occurrences,
    pivot,
    projection,
    structures,
)
from ord_schema.proto import dataset_pb2, reaction_pb2


def _write(path, metadata):
    schema = pa.schema([pa.field("x", pa.int32())]).with_metadata(metadata)
    pq.write_table(pa.table({"x": [1]}, schema=schema), path)


def _valid_metadata(**overrides):
    metadata = {
        "ord.artifact": "structures",
        "ord.source_md5": "0" * 32,
        "ord.rdkit_version": rdBase.rdkitVersion,
        "ord.source_dataset_id": "ord_dataset-1",
        "ord.artifact_schema_version": base.ARTIFACT_SCHEMA_VERSION,
    }
    metadata.update(overrides)
    return metadata


# Stamps


def test_stamps_carries_the_current_versions():
    value = base.current_stamps("structures", "ord_dataset-1", "abc")
    assert value.artifact == "structures"
    assert value.artifact_schema_version == base.ARTIFACT_SCHEMA_VERSION
    assert value.rdkit_version == rdBase.rdkitVersion


def test_to_metadata_names_no_library_release():
    # A release need not change any artifact, so a key naming one would make a rebuild
    # of an unchanged artifact write different bytes.
    stamps = base.current_stamps("structures", "ord_dataset-1", "abc")
    assert set(base.to_metadata(stamps)) == {
        "ord.artifact",
        "ord.artifact_schema_version",
        "ord.source_md5",
        "ord.source_dataset_id",
        "ord.rdkit_version",
    }


def test_to_metadata_omits_a_missing_dataset_id():
    metadata = base.to_metadata(base.current_stamps("structures", None, "abc"))
    assert "ord.source_dataset_id" not in metadata
    assert metadata["ord.source_md5"] == "abc"


def test_load_stamps_round_trips(tmp_path):
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata())
    stamps = base.load_stamps(path)
    assert stamps.artifact == "structures"
    assert stamps.source_dataset_id == "ord_dataset-1"
    assert stamps.source_md5 == "0" * 32


def test_load_stamps_tolerates_a_missing_dataset_id(tmp_path):
    path = tmp_path / "artifact.parquet"
    metadata = _valid_metadata()
    del metadata["ord.source_dataset_id"]
    _write(path, metadata)
    assert base.load_stamps(path).source_dataset_id is None


@pytest.mark.parametrize(
    "missing",
    [
        "ord.artifact",
        "ord.artifact_schema_version",
        "ord.source_md5",
    ],
)
def test_load_stamps_names_every_missing_key(tmp_path, missing):
    path = tmp_path / "artifact.parquet"
    metadata = _valid_metadata()
    del metadata[missing]
    _write(path, metadata)
    with pytest.raises(ValueError, match=missing):
        base.load_stamps(path)


def test_load_stamps_rejects_a_file_with_no_metadata(tmp_path):
    path = tmp_path / "plain.parquet"
    pq.write_table(pa.table({"x": [1]}), path)
    with pytest.raises(ValueError, match="not a derived artifact"):
        base.load_stamps(path)


def test_is_artifact_recognizes_only_stamped_files(tmp_path):
    _write(tmp_path / "artifact.parquet", _valid_metadata())
    pq.write_table(pa.table({"x": [1]}), tmp_path / "plain.parquet")
    assert base.is_artifact(tmp_path / "artifact.parquet")
    assert not base.is_artifact(tmp_path / "plain.parquet")
    assert not base.is_artifact(tmp_path / "absent.parquet")


def test_is_artifact_refuses_to_call_an_unreadable_file_a_source(tmp_path):
    # ArrowInvalid subclasses ValueError, so catching ValueError to mean "no stamps"
    # would hand a truncated file to a reader that fails later without naming it.
    path = tmp_path / "truncated.parquet"
    _write(path, _valid_metadata())
    path.write_bytes(path.read_bytes()[:200])
    with pytest.raises(ValueError, match="not readable as Parquet"):
        base.is_artifact(path)


def test_is_current_requires_every_compared_stamp_to_match(tmp_path, monkeypatch):
    # Each of these alone can change a projected value, so each alone must force a
    # rebuild; a check that quietly stopped working would serve obsolete artifacts.
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata())
    assert base.is_current(path, "structures", "0" * 32)
    assert not base.is_current(path, "structures", "1" * 32)
    monkeypatch.setattr(base, "ARTIFACT_SCHEMA_VERSION", "99")
    assert not base.is_current(path, "structures", "0" * 32)


def test_a_footer_key_the_stamps_do_not_define_leaves_an_artifact_current(tmp_path):
    # A release number in the footer is read past rather than compared: a release need
    # not change any artifact, and ARTIFACT_SCHEMA_VERSION is what says the columns
    # moved.
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata(**{"ord.ord_schema_version": "0.0.1"}))
    assert base.is_current(path, "structures", "0" * 32)
    assert base.stamps_are_current(base.load_stamps(path), "structures")


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
        # glob.escape writes a wildcard character as a bracket holding only it, which
        # matches that character alone, so the root reads it as the character.
        (glob.escape("data[1]") + "/*/*.parquet", "data[1]"),
        (glob.escape("a*b?c") + "/x/*.parquet", "a*b?c/x"),
        (glob.escape("data[1]/aa/one.parquet"), "data[1]/aa"),
        # A bracket matching more than one character is still a wildcard, including one
        # that holds an escape-shaped run.
        ("data/[ab]/*.parquet", "data"),
        ("data/[a[*]/x.parquet", "data"),
    ],
)
def test_glob_root(pattern, expected):
    assert base.glob_root(pattern) == pathlib.PurePath(expected)


def test_output_path_mirrors_the_input_layout():
    assert base.output_path(
        "data/aa/ord_dataset-x.parquet", "data/*/*.parquet", "structures"
    ) == pathlib.Path("structures/aa/ord_dataset-x.parquet")


def test_output_path_for_an_exact_filename_writes_into_the_directory():
    # A pattern naming one file must land inside output_dir, not become it.
    assert base.output_path(
        "data/aa/one.parquet", "data/aa/one.parquet", "structures"
    ) == pathlib.Path("structures/one.parquet")


def test_output_path_mirrors_from_an_escaped_directory_name():
    pattern = glob.escape("data[1]/aa") + "/*.parquet"
    assert base.output_path(
        "data[1]/aa/ord_dataset-x.parquet", pattern, "structures"
    ) == pathlib.Path("structures/ord_dataset-x.parquet")


def test_an_explicit_root_mirrors_from_where_the_caller_says():
    assert base.output_path(
        "data/aa/ord_dataset-x.parquet",
        "data/*/*.parquet",
        "structures",
        root="data/aa",
    ) == pathlib.Path("structures/ord_dataset-x.parquet")


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
        base.derive_tree(
            str(tmp_path / "*.parquet"),
            str(tmp_path / "out"),
            artifact="structures",
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

    written, skipped, ignored = base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
    )
    assert (written, skipped, ignored) == (2, 0, 0)
    assert {path.parent.name for path in written_paths} == {"aa", "bb"}


def test_derive_tree_mirrors_from_a_directory_whose_name_holds_a_wildcard(tmp_path):
    data = tmp_path / "data[1]"
    for name in ("aa", "bb"):
        (data / name).mkdir(parents=True)
        _fake_source(data / name / "source.parquet")
    written_paths = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        del source, source_md5, source_dataset_id  # Unused.
        written_paths.append(pathlib.Path(output))
        pathlib.Path(output).write_bytes(b"")
        return 1

    output_dir = tmp_path / "out"
    written, _, _ = base.derive_tree(
        glob.escape(str(data)) + "/*/*.parquet",
        str(output_dir),
        artifact="structures",
        write=_write_one,
    )
    assert written == 2
    assert {path.parent for path in written_paths} == {
        output_dir / "aa",
        output_dir / "bb",
    }


def test_derive_tree_derives_a_source_reached_twice_only_once(tmp_path):
    # A pattern with more than one ** reaches the same file by more than one route.
    # Derived twice, the second pass reads what the first wrote and the counts returned
    # say two sources where the tree holds one -- which a caller comparing them against
    # the artifacts on disk reads as artifacts having gone missing.
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    seen = []

    def _write_one(path, output, *, source_md5, source_dataset_id):
        seen.append(str(path))
        pathlib.Path(output).write_bytes(b"")
        return 1

    written, skipped, _ = base.derive_tree(
        str(tmp_path / "**" / "**" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
    )
    assert len(seen) == 1
    assert (written, skipped) == (1, 0)


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

    base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
    )
    assert seen == [("ord_dataset-1", True)]


def test_derive_tree_refuses_to_write_over_its_own_sources(tmp_path):
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    calls = []
    with pytest.raises(ValueError, match="would write over its own inputs"):
        base.derive_tree(
            str(tmp_path / "*" / "*.parquet"),
            str(tmp_path),
            artifact="structures",
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
    with pytest.raises(ValueError, match="would write over its own inputs"):
        base.derive_tree(
            str(tmp_path / "**" / "*.parquet"),
            str(tmp_path / "aa"),
            artifact="structures",
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

    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
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

    monkeypatch.setattr(base, "is_current", lambda *args: True)
    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
    ) == (0, 1, 0)
    assert not calls
    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
        force=True,
    ) == (1, 0, 0)
    assert len(calls) == 1


def test_derive_tree_derives_again_what_the_caller_finds_out_of_date(
    tmp_path, monkeypatch
):
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    calls = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        del source, source_md5, source_dataset_id  # Unused.
        calls.append(output)
        pathlib.Path(output).write_bytes(b"")
        return 1

    monkeypatch.setattr(base, "is_current", lambda *args: True)
    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
        current=lambda path: True,
    ) == (0, 1, 0)
    asked = []
    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
        current=lambda path: asked.append(path) or False,
    ) == (1, 0, 0)
    assert asked == calls == [tmp_path / "out" / "aa" / "source.parquet"]


def test_derive_tree_reads_the_named_parent_artifact(tmp_path):
    (tmp_path / "aa").mkdir()
    _write(
        tmp_path / "aa" / "projected.parquet",
        _valid_metadata(**{"ord.artifact": "projection", "ord.source_md5": "a" * 32}),
    )
    seen = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        seen.append((pathlib.Path(source).name, source_md5, source_dataset_id))
        pathlib.Path(output).write_bytes(b"")
        return 1

    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
        parent_artifact="projection",
    ) == (1, 0, 0)
    # The source dataset's hash and ID pass through, so an artifact names the dataset it
    # reflects rather than the artifact it read.
    assert seen == [("projected.parquet", "a" * 32, "ord_dataset-1")]


def test_derive_tree_ignores_a_source_dataset_when_a_parent_artifact_is_named(tmp_path):
    (tmp_path / "aa").mkdir()
    _fake_source(tmp_path / "aa" / "source.parquet")
    _write(
        tmp_path / "aa" / "projected.parquet",
        _valid_metadata(**{"ord.artifact": "projection"}),
    )
    written_paths = []

    def _write_one(source, output, *, source_md5, source_dataset_id):
        del source, source_md5, source_dataset_id  # Unused.
        written_paths.append(output)
        pathlib.Path(output).write_bytes(b"")
        return 1

    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=_write_one,
        parent_artifact="projection",
    ) == (1, 0, 1)
    assert [path.name for path in written_paths] == ["projected.parquet"]


def test_derive_tree_ignores_an_artifact_of_the_wrong_kind(tmp_path):
    # What an earlier run of this same command left in a recursive pattern's reach.
    (tmp_path / "aa").mkdir()
    _write(
        tmp_path / "aa" / "already.parquet",
        _valid_metadata(**{"ord.artifact": "structures"}),
    )
    assert base.derive_tree(
        str(tmp_path / "*" / "*.parquet"),
        str(tmp_path / "out"),
        artifact="structures",
        write=lambda *args, **kwargs: 1,
        parent_artifact="projection",
    ) == (0, 0, 1)


def test_derive_tree_refuses_a_stale_parent(tmp_path, monkeypatch):
    # Passing the parent's hash through carries the source content across the hop but
    # not the version stamps. An artifact written from a stale projection would stamp
    # itself with the current versions, and the dataset hash it inherits does not change
    # when the projection is rebuilt -- so nothing would ever mark it stale again.
    (tmp_path / "aa").mkdir()
    parent = tmp_path / "aa" / "projected.parquet"
    _write(parent, _valid_metadata(**{"ord.artifact": "projection"}))
    monkeypatch.setattr(base, "ARTIFACT_SCHEMA_VERSION", "next")
    with pytest.raises(ValueError, match="stale projection"):
        base.derive_tree(
            str(tmp_path / "*" / "*.parquet"),
            str(tmp_path / "out"),
            artifact="structures",
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
        _valid_metadata(**{"ord.artifact": "projection"}),
    )
    (tmp_path / "data").mkdir()
    source = tmp_path / "data" / "ds.parquet"
    _fake_source(source)
    before = source.read_bytes()
    with pytest.raises(ValueError, match="files this library did not write"):
        base.derive_tree(
            str(tmp_path / "projections" / "*.parquet"),
            str(tmp_path / "data"),
            artifact="structures",
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
        _valid_metadata(**{"ord.artifact": "projection"}),
    )
    (tmp_path / "structures").mkdir()
    _write(tmp_path / "structures" / "ds.parquet", _valid_metadata())
    written, skipped, ignored = base.derive_tree(
        str(tmp_path / "projections" / "*.parquet"),
        str(tmp_path / "structures"),
        artifact="structures",
        write=lambda *args, **kwargs: 1,
        force=True,
        parent_artifact="projection",
    )
    assert (written, skipped, ignored) == (1, 0, 0)


def test_derive_tree_rewrites_an_artifact_stamped_by_an_older_library(tmp_path):
    # The one case the guard must not refuse and cannot see from the stamp set: a tree
    # written before a required key existed. Adding or renaming one leaves every file on
    # disk unreadable as stamps, and a guard that asks for the whole set then reads a
    # corpus this driver wrote as a stranger's -- refusing the re-derive that is the
    # only way to bring it up to date, with no --force to override it.
    (tmp_path / "projections").mkdir()
    _write(
        tmp_path / "projections" / "ds.parquet",
        _valid_metadata(**{"ord.artifact": "projection"}),
    )
    (tmp_path / "structures").mkdir()
    superseded = _valid_metadata()
    del superseded["ord.artifact_schema_version"]
    superseded["ord.artifact_lineage"] = "structures=1,projection=1"
    _write(tmp_path / "structures" / "ds.parquet", superseded)
    written, skipped, ignored = base.derive_tree(
        str(tmp_path / "projections" / "*.parquet"),
        str(tmp_path / "structures"),
        artifact="structures",
        write=lambda *args, **kwargs: 1,
        parent_artifact="projection",
    )
    assert (written, skipped, ignored) == (1, 0, 0)


def test_stamps_record_the_rdkit_version():
    value = base.current_stamps("structures", "ord_dataset-1", "abc")
    assert value.rdkit_version == rdBase.rdkitVersion
    assert base.stamps_are_current(value, "structures")


def test_an_artifact_from_a_different_rdkit_reads_stale():
    # Fingerprints and canonical SMILES are both functions of RDKit, so an upgrade
    # must mark artifacts stale or the screen's completeness rests on cross-version
    # fingerprint compatibility nothing enforces.
    value = dataclasses.replace(
        base.current_stamps("structures", "ord_dataset-1", "abc"),
        rdkit_version="0000.00.0",
    )
    assert not base.stamps_are_current(value, "structures")


def test_an_artifact_with_no_rdkit_stamp_reads_stale():
    # Files stamped before the key existed still load, and read stale: rebuilding is
    # the safe answer for an artifact whose RDKit nobody recorded.
    value = dataclasses.replace(
        base.current_stamps("structures", "ord_dataset-1", "abc"),
        rdkit_version=None,
    )
    assert not base.stamps_are_current(value, "structures")


# Every artifact this library derives.
_ARTIFACTS = frozenset(
    {
        projection.ARTIFACT,
        structures.ARTIFACT,
        arrivals.ARTIFACT,
        pivot.ARTIFACT,
        occurrences.ARTIFACT,
    }
)


# The artifact schema version


def test_a_schema_bump_makes_every_artifact_stale(tmp_path, monkeypatch):
    # One schema for every artifact: a build is published and fetched as a unit, so a
    # bump reaches all four rather than only the artifact whose columns moved.
    paths = {artifact: tmp_path / f"{artifact}.parquet" for artifact in _ARTIFACTS}
    for artifact, path in paths.items():
        _write(path, _valid_metadata(**{"ord.artifact": artifact}))
        assert base.is_current(path, artifact, "0" * 32), artifact
    monkeypatch.setattr(base, "ARTIFACT_SCHEMA_VERSION", "next")
    for artifact, path in paths.items():
        assert not base.is_current(path, artifact, "0" * 32), artifact


# What each artifact holds when derived from the fixtures below, pinned beside the
# schema version it was derived under. The stamps compare ARTIFACT_SCHEMA_VERSION and
# nothing else about this library, so a change to what an artifact holds reads as
# current everywhere; this is where that change surfaces for the decision. A change a
# reader of the pinned schema would misread bumps ARTIFACT_SCHEMA_VERSION. One it
# reads correctly -- a normalization fixed, a column added -- is derived again (with
# --force where the stamps still read current) and published in place.
_GOLDEN_SCHEMA_VERSION = "2"
_GOLDEN_OUTPUT = {
    "projection": "08d030fbeecbdc505fce84696fc69bffb6df744df045c83c6447c3bb517fc00a",
    "structures": "3c397cfc3ba2b8e0debd9b725ada6d119551e42aa7c20bf9791983e679ce7158",
    "arrivals": "d30f6710e8fbf5f2b1d11793c4de66b8d7b3bc3195e9fd83771627ad78f0d432",
    "pivot": "4951020b810d68ac68a6a909b2fb37bee2bf145fec279adc87ad5632db04b465",
    "occurrences": "d0617c0972b1c054c5f4507760cc53e055eb5cfd1336b0279e700a3c0a8c1118",
}
# The RDKit those hashes were derived under. An RDKit upgrade can change canonical
# SMILES and fingerprints, but ord.rdkit_version already makes every artifact stale
# then, so its hashes move without a schema bump.
_GOLDEN_RDKIT = "2026.03.5"

# One reaction setting nearly every field, and 80 real ones with slash-format dates.
_GOLDEN_FIXTURES = ("full.pbtxt", "ord-nielsen-example.pbtxt")

# Footer keys recording who derived a file rather than what it holds: the RDKit and the
# schema, pinned on their own above, the hash of how the fixture happened to serialize,
# and the index of first commits an arrivals artifact read, whose values its rows
# carry. The rest -- the artifact name, the pivot or occurrence path a reader checks a
# file against -- are part of what the artifact says.
_PROVENANCE_KEYS = frozenset(
    {
        b"ord.rdkit_version",
        b"ord.artifact_schema_version",
        b"ord.source_md5",
        b"ord.arrivals_version",
        b"ord.arrivals_head",
    }
)


def _content_hash(paths: list[pathlib.Path]) -> str:
    """Returns a hash of the footer, columns, and rows the artifacts at ``paths`` hold.

    Provenance stamps are left out (see ``_PROVENANCE_KEYS``); Parquet encoding is left
    out by hashing decoded rows; and row order is left out by sorting them, so the
    digest pins what the artifacts hold. The order they hold it in is pinned by the
    pivot and occurrence tests, and the bytes by a rebuild below.

    Args:
        paths: Artifact files, in a fixed order.

    Returns:
        A SHA-256 hex digest.
    """
    digest = hashlib.sha256()
    for path in paths:
        table = pq.read_table(path)
        footer = sorted(
            (key, value)
            for key, value in (table.schema.metadata or {}).items()
            if key not in _PROVENANCE_KEYS
        )
        digest.update(repr(footer).encode())
        digest.update(str(table.schema.remove_metadata()).encode())
        for row in sorted(
            json.dumps(row, sort_keys=True, default=repr) for row in table.to_pylist()
        ):
            digest.update(row.encode())
    return digest.hexdigest()


def _derive_fixtures(root: pathlib.Path) -> dict[str, list[pathlib.Path]]:
    """Derives every artifact from the golden fixtures beneath ``root``.

    Args:
        root: An empty directory to write the sources and artifacts beneath.

    Returns:
        Each artifact's files, in a fixed order.
    """
    testdata = pathlib.Path(__file__).parents[1] / "orm" / "testdata"
    derived: dict[str, list[pathlib.Path]] = collections.defaultdict(list)
    for position, name in enumerate(_GOLDEN_FIXTURES):
        source = root / "data" / f"{name}.parquet"
        source.parent.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            message_helpers.load_message(str(testdata / name), dataset_pb2.Dataset),
            source,
        )
        projected = root / "projections" / source.name
        projected.parent.mkdir(parents=True, exist_ok=True)
        projection.write_projection(source, projected)
        derived["projection"].append(projected)
        structured = root / "structures" / source.name
        structured.parent.mkdir(parents=True, exist_ok=True)
        structures.write_structures(projected, structured)
        derived["structures"].append(structured)
        # Each fixture first reached ord-data at a commit of its own, with an offset
        # each side of UTC, so the pin covers both signs.
        sha = f"{position + 1:040x}"
        index = history.Index(
            pa.Table.from_pylist(
                [
                    {
                        "reaction_id": reaction_id,
                        "sha": sha,
                        "timestamp": datetime(2021, 7, 28 + position, 14, 7, 53),
                        "utc_offset_minutes": -420 if position else 60,
                        "dataset_id": f"ord_dataset-{position}",
                    }
                    for reaction_id in parquet.DatasetView(source).iter_reaction_ids()
                ],
                schema=history.INDEX_SCHEMA,
            ),
            [sha],
        )
        arrived = root / "arrivals" / source.name
        arrived.parent.mkdir(parents=True, exist_ok=True)
        arrivals.write_arrivals(projected, arrived, index=index)
        derived["arrivals"].append(arrived)
        for level in sorted(pivot.LEVELS):
            pivoted = root / "pivots" / level / source.name
            pivoted.parent.mkdir(parents=True, exist_ok=True)
            pivot.write_pivot(projected, pivoted, level_path=level)
            derived["pivot"].append(pivoted)
        for path, (level, *_) in sorted(occurrences.PATHS.items()):
            indexed = root / "occurrences" / path / source.name
            indexed.parent.mkdir(parents=True, exist_ok=True)
            occurrences.write_occurrences(
                root / "pivots" / level.path / source.name, indexed, path=path
            )
            derived["occurrences"].append(indexed)
    return derived


@pytest.fixture(scope="module")
def golden_output(tmp_path_factory) -> dict[str, str]:
    """Derives every artifact from the fixtures and hashes each artifact's files."""
    derived = _derive_fixtures(tmp_path_factory.mktemp("golden"))
    return {artifact: _content_hash(paths) for artifact, paths in derived.items()}


def test_a_rebuild_under_another_release_writes_the_same_bytes(tmp_path, monkeypatch):
    # What lets an upload drop the files a build left unchanged: a release need not
    # change any artifact, so nothing it reports about itself may reach the bytes.
    first = _derive_fixtures(tmp_path / "first")
    monkeypatch.setattr(importlib.metadata, "version", lambda name: "0.0.1")
    second = _derive_fixtures(tmp_path / "second")
    assert first.keys() == second.keys() == _ARTIFACTS
    for artifact, paths in first.items():
        assert len(paths) == len(second[artifact])
        for path, rebuilt in zip(paths, second[artifact], strict=True):
            assert path.read_bytes() == rebuilt.read_bytes(), path


def test_every_artifact_has_pinned_output():
    assert set(_GOLDEN_OUTPUT) == _ARTIFACTS


@pytest.mark.parametrize("artifact", sorted(_GOLDEN_OUTPUT))
def test_an_artifact_holds_what_its_schema_pins(artifact, golden_output):
    expected = _GOLDEN_OUTPUT[artifact]
    found = golden_output[artifact]
    pin = f"_GOLDEN_OUTPUT[{artifact!r}] = {found!r}"
    assert _GOLDEN_SCHEMA_VERSION == base.ARTIFACT_SCHEMA_VERSION, (
        "ARTIFACT_SCHEMA_VERSION was bumped; set _GOLDEN_SCHEMA_VERSION = "
        f"{base.ARTIFACT_SCHEMA_VERSION!r} and pin each artifact's output again"
    )
    if found != expected and rdBase.rdkitVersion != _GOLDEN_RDKIT:
        pytest.fail(
            f"RDKit is {rdBase.rdkitVersion}, not the {_GOLDEN_RDKIT} these hashes "
            "were derived under, and every artifact already reads stale through "
            f"ord.rdkit_version: set _GOLDEN_RDKIT = {rdBase.rdkitVersion!r} and "
            f"{pin}, with no schema bump."
        )
    assert found == expected, (
        f"the {artifact} artifact derived from the fixtures changed, and artifacts "
        "already derived still read as current. Bump ARTIFACT_SCHEMA_VERSION if a "
        "reader of the current schema would misread the new files; otherwise derive "
        f"them again with --force. Then {pin}."
    )


def test_missing_columns_sees_a_field_added_inside_a_struct(tmp_path):
    # Nearly every field the projection has is nested, so a top-level comparison would
    # call a file current that is missing everything added since it was written -- and a
    # reader binding the new column would fail on a file nothing marked stale.
    path = tmp_path / "artifact.parquet"
    inner = pa.struct([pa.field("value", pa.string())])
    pq.write_table(pa.table({"provenance": pa.array([{"value": "x"}], inner)}), path)
    widened = pa.schema(
        [
            pa.field(
                "provenance",
                pa.struct(
                    [
                        pa.field("value", pa.string()),
                        pa.field("timestamp", pa.timestamp("us")),
                    ]
                ),
            )
        ]
    )
    assert base.missing_columns(path, widened) == ["provenance.timestamp"]


def test_missing_columns_sees_a_leaf_whose_type_changed(tmp_path):
    # A path comparison alone cannot see this: moving a timestamp from milliseconds to
    # microseconds leaves every name where it was, so the file would read as current
    # and keep serving values truncated against the current definition.
    path = tmp_path / "old.parquet"
    millis = pa.schema(
        [
            pa.field("reaction_id", pa.string()),
            pa.field(
                "provenance", pa.struct([pa.field("recorded", pa.timestamp("ms"))])
            ),
        ]
    )
    pq.write_table(millis.empty_table(), path)
    micros = pa.schema(
        [
            pa.field("reaction_id", pa.string()),
            pa.field(
                "provenance", pa.struct([pa.field("recorded", pa.timestamp("us"))])
            ),
        ]
    )
    assert base.missing_columns(path, micros) == ["provenance.recorded: timestamp[us]"]
    # The same file against its own schema is current, so the check is about the type
    # rather than about a nested field being reported twice.
    assert base.missing_columns(path, millis) == []


@pytest.mark.parametrize(
    ("name", "schema"),
    [("projection", projection.SCHEMA), ("structures", structures.SCHEMA)],
)
def test_every_declared_type_survives_a_parquet_round_trip(tmp_path, name, schema):
    # missing_columns compares leaf types, and Parquet does not store every Arrow type
    # it is given -- a second-resolution timestamp comes back as milliseconds. A type
    # the writer coerces would make every file of that artifact read as stale forever,
    # which is a re-derive that never converges rather than a wrong answer.
    path = tmp_path / f"{name}.parquet"
    pq.write_table(schema.empty_table(), path)
    assert base.missing_columns(path, schema) == []


def test_missing_columns_descends_through_lists_and_maps(tmp_path):
    # A list's own name is not a field a reader binds; what it holds is.
    path = tmp_path / "artifact.parquet"
    held = pa.list_(pa.struct([pa.field("a", pa.string())]))
    pq.write_table(pa.table({"items": pa.array([[{"a": "x"}]], held)}), path)
    widened = pa.schema(
        [
            pa.field(
                "items",
                pa.list_(
                    pa.struct([pa.field("a", pa.string()), pa.field("b", pa.int64())])
                ),
            )
        ]
    )
    assert base.missing_columns(path, widened) == ["items.b"]


def test_a_file_carrying_every_nested_field_lacks_nothing(tmp_path):
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata())
    assert base.missing_columns(path, pa.schema([pa.field("x", pa.int32())])) == []
