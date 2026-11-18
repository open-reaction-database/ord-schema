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

"""Tests for ord_schema.artifacts.history."""

import gzip
import os
import pathlib
import shutil
import subprocess
from datetime import datetime

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from google.protobuf import text_format

from ord_schema import parquet
from ord_schema.artifacts import history
from ord_schema.artifacts.scripts import index_arrivals
from ord_schema.proto import dataset_pb2, reaction_pb2

# Isolated from the machine's own git configuration -- signing, hooks, a global LFS
# filter -- so the repository holds exactly what each test writes.
_ENV = {
    **os.environ,
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_AUTHOR_NAME": "ORD",
    "GIT_AUTHOR_EMAIL": "ord@example.com",
    "GIT_COMMITTER_NAME": "ORD",
    "GIT_COMMITTER_EMAIL": "ord@example.com",
}

# The author dates differ from the committer dates in day and offset, so a walk reading
# the wrong one is visible.
_AUTHOR_DATE = "2019-01-01 00:00:00 +0900"


def _git(repository, *args, committed=None):
    """Runs git in ``repository``, committing at ``committed`` where given."""
    env = dict(_ENV)
    if committed is not None:
        env["GIT_COMMITTER_DATE"] = committed
        env["GIT_AUTHOR_DATE"] = _AUTHOR_DATE
    result = subprocess.run(
        ["git", "-C", str(repository), *args],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def _dataset(dataset_id, *reaction_ids):
    reactions = []
    for reaction_id in reaction_ids:
        reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
        component = reaction.inputs["input"].components.add()
        component.identifiers.add(type="SMILES", value="CCO")
        reactions.append(reaction)
    return dataset_pb2.Dataset(
        dataset_id=dataset_id, name="test", description="test", reactions=reactions
    )


def _write(repository, path, dataset):
    """Writes ``dataset`` to ``path`` in the format its suffix names."""
    target = repository / path
    target.parent.mkdir(parents=True, exist_ok=True)
    if path.endswith(".parquet"):
        parquet.save_dataset(dataset, target)
    elif path.endswith(".pb.gz"):
        target.write_bytes(gzip.compress(dataset.SerializeToString()))
    elif path.endswith(".pb"):
        target.write_bytes(dataset.SerializeToString())
    else:
        target.write_text(text_format.MessageToString(dataset))


def _commit(repository, message, committed):
    _git(repository, "add", "-A")
    _git(repository, "commit", "-q", "-m", message, committed=committed)
    return _git(repository, "rev-parse", "HEAD")


def _init(repository):
    repository.mkdir()
    _git(repository, "init", "-q", "-b", "main")
    _git(repository, "lfs", "install", "--local")
    # As ord-data stores them: Parquet and binary datasets in LFS, text in git.
    (repository / ".gitattributes").write_text(
        "data/**/*.parquet filter=lfs diff=lfs merge=lfs -text\n"
        "data/**/*.pb.gz filter=lfs diff=lfs merge=lfs -text\n"
    )


@pytest.fixture(scope="module")
def ord_data(tmp_path_factory):
    """A repository whose main holds each shape ord-data's history does.

    Returns a namespace-like dict of the repository and the commits a test asks about.
    """
    repository = tmp_path_factory.mktemp("history") / "ord-data"
    _init(repository)
    (repository / "README.md").write_text("ord-data\n")
    commits = {"root": _commit(repository, "root", "2020-07-30 09:50:35 -0600")}
    # Text, with a reaction a submission left without an ID.
    _write(
        repository,
        "data/aa/ord_dataset-aa.pbtxt",
        _dataset("ord_dataset-aa", "ord-r1", "ord-r2", ""),
    )
    commits["text"] = _commit(repository, "add aa", "2020-09-29 12:52:22 -0600")
    # A submission merged in two steps: first into a branch, then into main.
    _git(repository, "checkout", "-q", "-b", "submission")
    _write(
        repository, "data/bb/ord_dataset-bb.pb.gz", _dataset("ord_dataset-bb", "ord-r3")
    )
    commits["side"] = _commit(repository, "add bb", "2021-02-25 06:35:30 -0700")
    _git(repository, "checkout", "-q", "main")
    _git(
        repository,
        "merge",
        "-q",
        "--no-ff",
        "-m",
        "Merge #70 in main",
        "submission",
        committed="2021-03-03 13:33:15 -0500",
    )
    commits["merge"] = _git(repository, "rev-parse", "HEAD")
    # A rename into another format, holding the same reactions.
    (repository / "data/aa/ord_dataset-aa.pbtxt").unlink()
    _write(
        repository,
        "data/aa/ord_dataset-aa.parquet",
        _dataset("ord_dataset-aa", "ord-r1", "ord-r2"),
    )
    commits["convert"] = _commit(repository, "to parquet", "2026-06-01 10:00:00 +0100")
    # An update that adds a reaction to a dataset already held.
    _write(
        repository,
        "data/aa/ord_dataset-aa.parquet",
        _dataset("ord_dataset-aa", "ord-r1", "ord-r2", "ord-r4"),
    )
    commits["update"] = _commit(repository, "add r4", "2026-07-24 23:30:00 +1300")
    # A reshard: bb's reaction moves to cc, beside a new one.
    (repository / "data/bb/ord_dataset-bb.pb.gz").unlink()
    _write(
        repository,
        "data/cc/ord_dataset-cc.parquet",
        _dataset("ord_dataset-cc", "ord-r3", "ord-r5"),
    )
    commits["split"] = _commit(repository, "split bb", "2026-10-15 08:00:00 +0000")
    return {"repository": repository, "commits": commits}


def _row(sha, timestamp, offset, dataset_id):
    return {
        "sha": sha,
        "timestamp": timestamp,
        "utc_offset_minutes": offset,
        "dataset_id": dataset_id,
    }


def _index(rows, chain):
    """Returns an index of ``rows``, each a reaction ID mapped to its first commit."""
    table = pa.Table.from_pylist(
        [{"reaction_id": reaction_id, **row} for reaction_id, row in rows.items()],
        schema=history.INDEX_SCHEMA,
    )
    return history.Index(table, chain)


class TestWalk:
    def test_each_reaction_reads_its_first_commit_on_main(self, ord_data):
        commits = ord_data["commits"]
        index = history.walk(ord_data["repository"])
        found = index.lookup(["ord-r1", "ord-r2", "ord-r3", "ord-r4", "ord-r5"])
        assert found == {
            # Unchanged by the rename into Parquet.
            "ord-r1": _row(
                commits["text"],
                datetime(2020, 9, 29, 12, 52, 22),
                -360,
                "ord_dataset-aa",
            ),
            "ord-r2": _row(
                commits["text"],
                datetime(2020, 9, 29, 12, 52, 22),
                -360,
                "ord_dataset-aa",
            ),
            # The merge into main rather than the earlier commit on the branch, and the
            # dataset it arrived in rather than the one holding it now.
            "ord-r3": _row(
                commits["merge"],
                datetime(2021, 3, 3, 13, 33, 15),
                -300,
                "ord_dataset-bb",
            ),
            "ord-r4": _row(
                commits["update"], datetime(2026, 7, 24, 23, 30), 780, "ord_dataset-aa"
            ),
            "ord-r5": _row(
                commits["split"], datetime(2026, 10, 15, 8), 0, "ord_dataset-cc"
            ),
        }

    def test_an_empty_id_names_no_reaction(self, ord_data):
        index = history.walk(ord_data["repository"])
        assert sorted(index.table["reaction_id"].to_pylist()) == [
            "ord-r1",
            "ord-r2",
            "ord-r3",
            "ord-r4",
            "ord-r5",
        ]

    def test_the_chain_follows_first_parents(self, ord_data):
        commits = ord_data["commits"]
        index = history.walk(ord_data["repository"])
        assert index.chain == tuple(
            commits[name]
            for name in ("root", "text", "merge", "convert", "update", "split")
        )
        assert index.head == commits["split"]

    def test_walks_to_a_named_revision(self, ord_data):
        commits = ord_data["commits"]
        index = history.walk(ord_data["repository"], revision=commits["merge"])
        assert index.head == commits["merge"]
        assert sorted(index.table["reaction_id"].to_pylist()) == [
            "ord-r1",
            "ord-r2",
            "ord-r3",
        ]

    def test_a_revision_is_never_read_as_an_option(self, ord_data, tmp_path):
        written = tmp_path / "written"
        with pytest.raises(ValueError, match="failed in"):
            history.walk(ord_data["repository"], revision=f"--output={written}")
        assert not written.exists()

    @pytest.mark.parametrize("walked_to", ["root", "text", "merge", "update", "split"])
    def test_extending_an_index_matches_walking_from_the_root(
        self, ord_data, walked_to
    ):
        repository = ord_data["repository"]
        previous = history.walk(repository, revision=ord_data["commits"][walked_to])
        extended = history.walk(repository, previous=previous)
        whole = history.walk(repository)
        assert extended.chain == whole.chain
        assert extended.table.equals(whole.table)

    def test_extending_reads_only_the_commits_after_the_head(
        self, ord_data, monkeypatch
    ):
        repository = ord_data["repository"]
        commits = ord_data["commits"]
        previous = history.walk(repository, revision=commits["update"])
        read = []
        original = history.read_ids
        monkeypatch.setattr(
            history,
            "read_ids",
            lambda path, content: read.append(path) or original(path, content),
        )
        history.walk(repository, previous=previous)
        assert read == ["data/cc/ord_dataset-cc.parquet"]

    def test_a_reaction_keeps_its_first_commit_across_a_later_walk(self, ord_data):
        # The index holds every ID it has met, so one an extension meets again is not
        # new, whatever the files between said.
        repository = ord_data["repository"]
        commits = ord_data["commits"]
        previous = _index(
            {"ord-r5": _row("f" * 40, datetime(2020, 1, 1), 0, "ord_dataset-zz")},
            [commits["root"], commits["text"], commits["merge"], commits["convert"]],
        )
        extended = history.walk(repository, previous=previous)
        assert extended.lookup(["ord-r5"])["ord-r5"]["dataset_id"] == "ord_dataset-zz"

    def test_refuses_a_shallow_clone(self, ord_data, tmp_path):
        # A shallow clone presents its oldest commit as a root, so every reaction it
        # holds would read as arriving there.
        shallow = tmp_path / "shallow"
        _git(
            tmp_path,
            "clone",
            "-q",
            "--depth=2",
            f"file://{ord_data['repository']}",
            str(shallow),
        )
        # Every LFS object fetched, as `git lfs pull` leaves a checkout, so nothing but
        # the depth stands between this clone and a walk.
        objects = pathlib.Path(".git", "lfs", "objects")
        shutil.copytree(ord_data["repository"] / objects, shallow / objects)
        with pytest.raises(
            ValueError, match=r"shallow clone.* git -C .* fetch --unshallow"
        ):
            history.walk(shallow)

    def _clone(self, ord_data, destination):
        """Clones the history with every LFS object, as `git lfs pull` leaves one."""
        _git(
            destination.parent,
            "clone",
            "-q",
            f"file://{ord_data['repository']}",
            str(destination),
        )
        objects = pathlib.Path(".git", "lfs", "objects")
        shutil.copytree(ord_data["repository"] / objects, destination / objects)
        return destination

    def test_refuses_a_grafted_history(self, ord_data, tmp_path):
        clone = self._clone(ord_data, tmp_path / "grafted")
        (clone / ".git" / "info").mkdir(exist_ok=True)
        (clone / ".git" / "info" / "grafts").write_text(
            ord_data["commits"]["merge"] + "\n"
        )
        with pytest.raises(ValueError, match="has a grafts file"):
            history.walk(clone)

    def test_ignores_replace_refs(self, ord_data, tmp_path):
        clone = self._clone(ord_data, tmp_path / "replaced")
        # Makes the merge read as a root, as a shallow boundary would.
        _git(clone, "replace", "--graft", ord_data["commits"]["merge"])
        assert history.walk(clone).table.equals(
            history.walk(ord_data["repository"]).table
        )

    def test_refuses_an_index_walked_along_another_chain(self, ord_data):
        commits = ord_data["commits"]
        previous = _index({}, [commits["root"], commits["side"]])
        with pytest.raises(ValueError, match="is not on the first-parent chain"):
            history.walk(ord_data["repository"], previous=previous)

    def test_refuses_a_missing_lfs_object(self, tmp_path):
        repository = tmp_path / "ord-data"
        _init(repository)
        _write(
            repository,
            "data/aa/ord_dataset-aa.parquet",
            _dataset("ord_dataset-aa", "ord-r1"),
        )
        commit = _commit(repository, "add aa", "2026-06-01 10:00:00 +0100")
        objects = repository / ".git" / "lfs" / "objects"
        for path in [path for path in objects.rglob("*") if path.is_file()]:
            path.unlink()
        with pytest.raises(
            ValueError,
            match=rf"1 LFS objects .* first {commit}:data/aa/ord_dataset-aa\.parquet",
        ):
            history.walk(repository)

    def test_refuses_a_file_recording_no_dataset_id(self, tmp_path):
        repository = tmp_path / "ord-data"
        _init(repository)
        _write(repository, "data/aa/ord_dataset-aa.pb", _dataset("", "ord-r1"))
        commit = _commit(repository, "add aa", "2026-06-01 10:00:00 +0100")
        with pytest.raises(
            ValueError, match=f"{commit}:data/aa/ord_dataset-aa.pb records no dataset"
        ):
            history.walk(repository)

    def test_refuses_a_reaction_first_appearing_in_two_datasets(self, tmp_path):
        repository = tmp_path / "ord-data"
        _init(repository)
        _write(repository, "data/aa/a.pb", _dataset("ord_dataset-aa", "ord-r1"))
        _write(repository, "data/bb/b.pb", _dataset("ord_dataset-bb", "ord-r1"))
        _commit(repository, "add both", "2026-06-01 10:00:00 +0100")
        with pytest.raises(
            ValueError,
            match="ord-r1 first appears in ord_dataset-aa and ord_dataset-bb",
        ):
            history.walk(repository)

    def test_refuses_an_lfs_pointer_naming_no_object(self, tmp_path):
        repository = tmp_path / "ord-data"
        _init(repository)
        target = repository / "data" / "aa" / "ord_dataset-aa.pb"
        target.parent.mkdir(parents=True)
        target.write_bytes(b"version https://git-lfs.github.com/spec/v1\nsize 12\n")
        _commit(repository, "add aa", "2026-06-01 10:00:00 +0100")
        with pytest.raises(ValueError, match="names no sha256 object"):
            history.walk(repository)

    def test_refuses_without_git(self, ord_data, monkeypatch):
        monkeypatch.setattr(history.shutil, "which", lambda name: None)
        with pytest.raises(ValueError, match="needs git on PATH"):
            history.walk(ord_data["repository"])

    def test_one_dataset_in_two_formats_at_once_is_one_arrival(self, tmp_path):
        repository = tmp_path / "ord-data"
        _init(repository)
        dataset = _dataset("ord_dataset-aa", "ord-r1")
        _write(repository, "data/aa/ord_dataset-aa.pb", dataset)
        _write(repository, "data/aa/ord_dataset-aa.pbtxt", dataset)
        _commit(repository, "add aa", "2020-12-01 16:33:43 -0700")
        index = history.walk(repository)
        assert index.table["reaction_id"].to_pylist() == ["ord-r1"]

    def test_files_outside_data_do_not_count(self, tmp_path):
        repository = tmp_path / "ord-data"
        _init(repository)
        _write(repository, "submission.pb", _dataset("ord_dataset-aa", "ord-r1"))
        _write(repository, "data/aa/ord_dataset-aa.pb", _dataset("ord_dataset-aa"))
        _commit(repository, "stage", "2021-04-20 10:00:00 -0400")
        assert history.walk(repository).table.num_rows == 0


class TestReadIds:
    @pytest.mark.parametrize(
        "name",
        ["d.pbtxt", "d.pb", "d.pb.gz", "d.parquet"],
    )
    def test_reads_every_format(self, tmp_path, name):
        _write(tmp_path, name, _dataset("ord_dataset-aa", "ord-r1", "ord-r2"))
        assert history.read_ids(name, tmp_path / name) == (
            "ord_dataset-aa",
            ["ord-r1", "ord-r2"],
        )

    @pytest.mark.parametrize("name", ["d.pbtxt", "d.pb", "d.pb.gz", "d.parquet"])
    def test_reads_bytes_as_well_as_files(self, tmp_path, name):
        _write(tmp_path, name, _dataset("ord_dataset-aa", "ord-r1"))
        content = (tmp_path / name).read_bytes()
        assert history.read_ids(name, content) == ("ord_dataset-aa", ["ord-r1"])

    def test_reads_text_a_retired_schema_wrote(self):
        # Fields today's Dataset no longer declares are skipped rather than refused.
        text = (
            'dataset_id: "ord_dataset-aa"\n'
            'retired_name: "x"\n'
            'reactions { reaction_id: "ord-r1" retired { value: 1 } }\n'
        )
        with pytest.raises(text_format.ParseError, match="retired_name"):
            text_format.Parse(text, dataset_pb2.Dataset())
        assert history.read_ids("d.pbtxt", text.encode()) == (
            "ord_dataset-aa",
            ["ord-r1"],
        )

    def test_refuses_another_suffix(self):
        with pytest.raises(ValueError, match=r"d\.json is not a dataset file"):
            history.read_ids("d.json", b"{}")


class TestCommitTime:
    @pytest.mark.parametrize(
        ("offset", "wall", "minutes"),
        [
            ("+0000", datetime(2021, 3, 3, 18, 33, 15), 0),
            ("-0500", datetime(2021, 3, 3, 13, 33, 15), -300),
            ("+0130", datetime(2021, 3, 3, 20, 3, 15), 90),
            ("+1300", datetime(2021, 3, 4, 7, 33, 15), 780),
        ],
    )
    def test_keeps_the_wall_clock_and_the_offset(self, offset, wall, minutes):
        # 2021-03-03 18:33:15 UTC.
        assert history.commit_time(1614796395, offset) == (wall, minutes)

    @pytest.mark.parametrize("offset", ["0500", "-05:00", "+05", "Z"])
    def test_refuses_another_spelling(self, offset):
        with pytest.raises(ValueError, match="not a git UTC offset"):
            history.commit_time(0, offset)


class TestIndex:
    def test_round_trips_through_a_file(self, ord_data, tmp_path):
        index = history.walk(ord_data["repository"])
        path = tmp_path / "history.parquet"
        history.write_index(index, path)
        read = history.read_index(path)
        assert read.chain == index.chain
        assert read.table.equals(index.table)

    def test_refuses_an_index_of_another_version(self, ord_data, tmp_path, monkeypatch):
        path = tmp_path / "history.parquet"
        history.write_index(history.walk(ord_data["repository"]), path)
        monkeypatch.setattr(history, "VERSION", "2")
        with pytest.raises(ValueError, match="walked under version 1, not 2"):
            history.read_index(path)

    def test_refuses_a_file_that_is_not_an_index(self, tmp_path):
        path = tmp_path / "other.parquet"
        pq.write_table(pa.table({"reaction_id": ["ord-r1"]}), path)
        with pytest.raises(ValueError, match="is not an arrivals index"):
            history.read_index(path)

    def test_refuses_an_empty_chain(self):
        with pytest.raises(ValueError, match="at least one commit"):
            _index({}, [])

    def test_lookup_leaves_out_what_it_does_not_hold(self):
        row = _row("a" * 40, datetime(2021, 1, 1), 0, "ord_dataset-aa")
        index = _index({"ord-r1": row}, ["a" * 40])
        assert index.lookup(["ord-r1", "ord-r9"]) == {"ord-r1": row}
        assert index.lookup([]) == {}


class TestCovers:
    def _stamped(self, path, stamps):
        pq.write_table(pa.table({"x": [1]}).replace_schema_metadata(stamps), path)
        return path

    def test_covers_a_projection_stamped_with_a_head_on_its_chain(self, tmp_path):
        earlier = _index({}, ["a" * 40])
        later = _index({}, ["a" * 40, "b" * 40])
        path = self._stamped(tmp_path / "p.parquet", earlier.stamps())
        assert earlier.covers(path)
        assert later.covers(path)

    def test_does_not_cover_a_head_off_its_chain(self, tmp_path):
        path = self._stamped(tmp_path / "p.parquet", _index({}, ["c" * 40]).stamps())
        assert not _index({}, ["a" * 40, "b" * 40]).covers(path)

    def test_does_not_cover_another_walk_version(self, tmp_path, monkeypatch):
        index = _index({}, ["a" * 40])
        path = self._stamped(tmp_path / "p.parquet", index.stamps())
        monkeypatch.setattr(history, "VERSION", "2")
        assert not index.covers(path)

    def test_does_not_cover_a_projection_without_a_stamp(self, tmp_path):
        index = _index({}, ["a" * 40])
        assert not index.covers(self._stamped(tmp_path / "p.parquet", {}))
        assert not index.covers(tmp_path / "absent.parquet")
        assert history.read_stamp(tmp_path / "absent.parquet") is None


class TestScript:
    def test_writes_then_extends_the_index(self, ord_data, tmp_path, monkeypatch):
        repository = ord_data["repository"]
        commits = ord_data["commits"]
        output = tmp_path / "history.parquet"
        index_arrivals.main(
            index_arrivals.parse_args(
                [
                    f"--ord_data={repository}",
                    f"--output={output}",
                    f"--revision={commits['merge']}",
                ]
            )
        )
        assert history.read_index(output).head == commits["merge"]
        walked = []
        original = history.walk
        monkeypatch.setattr(
            history,
            "walk",
            lambda *args, **kwargs: (
                walked.append(kwargs["previous"]) or original(*args, **kwargs)
            ),
        )
        index_arrivals.main(
            index_arrivals.parse_args(
                [f"--ord_data={repository}", f"--output={output}"]
            )
        )
        assert walked[0].head == commits["merge"]
        assert history.read_index(output).head == commits["split"]

    def test_rebuild_ignores_the_index_there(self, ord_data, tmp_path):
        repository = ord_data["repository"]
        output = tmp_path / "history.parquet"
        # Not an index at all, which only --rebuild may write over.
        pq.write_table(pa.table({"x": [1]}), output)
        args = [f"--ord_data={repository}", f"--output={output}"]
        with pytest.raises(ValueError, match="is not an arrivals index"):
            index_arrivals.main(index_arrivals.parse_args(args))
        index_arrivals.main(index_arrivals.parse_args([*args, "--rebuild"]))
        assert history.read_index(output).head == ord_data["commits"]["split"]
