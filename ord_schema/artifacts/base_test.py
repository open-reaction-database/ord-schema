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

import dataclasses
import glob
import pathlib
from importlib import metadata

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import rdBase

from ord_schema import parquet
from ord_schema.artifacts import base
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
        "ord.artifact": "structures",
        "ord.source_md5": "0" * 32,
        "ord.ord_schema_version": "9.9.9",
        "ord.rdkit_version": rdBase.rdkitVersion,
        "ord.source_dataset_id": "ord_dataset-1",
    }
    metadata.update(overrides)
    # Derived from whichever artifact the caller asked for, so a test that overrides one
    # does not silently get another's chain and pass for the wrong reason.
    metadata.setdefault("ord.artifact_lineage", base.lineage(metadata["ord.artifact"]))
    return metadata


# Stamps


def test_stamps_carries_the_current_versions():
    value = base.current_stamps("structures", "ord_dataset-1", "abc")
    assert value.artifact == "structures"
    assert value.artifact_lineage == base.lineage("structures")
    assert value.ord_schema_version


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
        "ord.artifact_lineage",
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


def test_is_current_requires_every_stamp_to_match(tmp_path, monkeypatch):
    # Each of these alone can change a projected value, so each alone must force a
    # rebuild; a check that quietly stopped working would serve obsolete base.
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata(**{"ord.ord_schema_version": "9.9.9"}))
    monkeypatch.setattr(base.metadata, "version", lambda _: "9.9.9")
    assert base.is_current(path, "structures", "0" * 32)
    monkeypatch.setattr(base.metadata, "version", lambda _: "9.9.10")
    assert not base.is_current(path, "structures", "0" * 32)
    monkeypatch.setattr(base.metadata, "version", lambda _: "9.9.9")
    monkeypatch.setattr(
        base, "ARTIFACT_VERSIONS", base.ARTIFACT_VERSIONS | {"structures": "99"}
    )
    assert not base.is_current(path, "structures", "0" * 32)


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


def test_an_explicit_root_mirrors_from_where_the_caller_says():
    # A caller building a pattern around a directory name has to escape it, and the
    # escape is spelled with the wildcard it is escaping: glob_root stops at "data[[]1]"
    # and mirrors the source layout a second time beneath output_dir. The caller knows
    # the directory it meant, so it says so rather than spelling it as a pattern.
    pattern = glob.escape("data[1]/aa") + "/*.parquet"
    source = "data[1]/aa/ord_dataset-x.parquet"
    assert base.output_path(source, pattern, "structures") == pathlib.Path(
        "structures/data[1]/aa/ord_dataset-x.parquet"
    )
    assert base.output_path(
        source, pattern, "structures", root="data[1]/aa"
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
        _current_metadata(**{"ord.artifact": "projection"}),
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
    _write(parent, _current_metadata(**{"ord.artifact": "projection"}))
    monkeypatch.setattr(
        base, "ARTIFACT_VERSIONS", base.ARTIFACT_VERSIONS | {"projection": "next"}
    )
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
        _current_metadata(**{"ord.artifact": "projection"}),
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
        _current_metadata(**{"ord.artifact": "projection"}),
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
        _current_metadata(**{"ord.artifact": "projection"}),
    )
    (tmp_path / "structures").mkdir()
    superseded = _valid_metadata()
    del superseded["ord.artifact_lineage"]
    superseded["ord.artifact_version"] = "1"
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


# Per-artifact versions


def test_a_lineage_names_the_artifact_and_every_ancestor():
    assert base.lineage("projection") == "projection=1"
    assert base.lineage("structures") == "structures=1,projection=1"
    assert base.lineage("occurrences") == "occurrences=1,pivot=1,projection=1"


def test_every_artifact_reaches_a_source_dataset():
    # The two tables have to agree about which artifacts exist, and every chain has to
    # terminate: lineage runs for every artifact written and every currency check made.
    assert set(base.DERIVED_FROM) == set(base.ARTIFACT_VERSIONS)
    assert {
        parent for parent in base.DERIVED_FROM.values() if parent is not None
    } <= set(base.ARTIFACT_VERSIONS)
    roots = {name for name, parent in base.DERIVED_FROM.items() if parent is None}
    for artifact in base.ARTIFACT_VERSIONS:
        # Named by the root rather than by its version, which would make this fail on
        # the next bump of an artifact it is not about.
        last = base.lineage(artifact).split(",")[-1].split("=")[0]
        assert last in roots, artifact


@pytest.mark.timeout(15)
def test_a_cycle_is_refused_rather_than_walked_forever(monkeypatch):
    # An infinite loop where an error belongs: a typo in a four-line literal would
    # otherwise hang every write and every currency check rather than saying so.
    #
    # Bounded by the marker, because the failure this guards against is a hang: without
    # it a regression stops CI rather than failing it, and the run has to be killed by
    # hand to find out which test never came back.
    monkeypatch.setattr(
        base, "DERIVED_FROM", base.DERIVED_FROM | {"projection": "occurrences"}
    )
    with pytest.raises(ValueError, match="has a cycle"):
        base.lineage("occurrences")


def test_bumping_an_artifact_leaves_what_it_was_derived_from_alone(monkeypatch):
    # The whole point of splitting the version. Occurrences are 1.9 seconds to derive
    # and the projection under them 34.6 minutes, so a change to the cheap one must not
    # charge the expensive one.
    before = base.lineage("projection")
    monkeypatch.setattr(
        base, "ARTIFACT_VERSIONS", base.ARTIFACT_VERSIONS | {"occurrences": "2"}
    )
    assert base.lineage("projection") == before
    assert base.lineage("pivot") == "pivot=1,projection=1"
    assert base.lineage("occurrences") == "occurrences=2,pivot=1,projection=1"


def test_bumping_an_artifact_reaches_everything_derived_from_it(monkeypatch):
    # And the other half: a projection redefined has to invalidate the occurrences three
    # derivations down, which a stamp naming only the artifact's own version cannot say
    # and one naming only its parent's cannot either.
    monkeypatch.setattr(
        base, "ARTIFACT_VERSIONS", base.ARTIFACT_VERSIONS | {"projection": "2"}
    )
    assert base.lineage("occurrences") == "occurrences=1,pivot=1,projection=2"
    assert base.lineage("structures") == "structures=1,projection=2"


def test_a_grandparent_bump_makes_a_grandchild_stale(tmp_path, monkeypatch):
    # The case a direct-parent stamp misses: occurrences record the pivot they came
    # from, and a projection bump moves neither the occurrences' own version nor the
    # pivot's, so nothing but the lineage can see it.
    path = tmp_path / "artifact.parquet"
    _write(path, _valid_metadata(**{"ord.artifact": "occurrences"}))
    monkeypatch.setattr(base.metadata, "version", lambda _: "9.9.9")
    assert base.is_current(path, "occurrences", "0" * 32)
    monkeypatch.setattr(
        base, "ARTIFACT_VERSIONS", base.ARTIFACT_VERSIONS | {"projection": "2"}
    )
    assert not base.is_current(path, "occurrences", "0" * 32)


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
