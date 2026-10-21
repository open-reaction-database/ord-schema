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

"""Where and when each reaction first reached ord-data, read from its git history.

A record says when its depositor created it and when the submission pipeline modified
it, but not when it entered ORD's public corpus: the pipeline stamps its own
``record_modified`` event when it processes a submission, which can precede the merge
by days. ord-data's history holds the answer, so this module walks it into an index,
from which :mod:`ord_schema.artifacts.first_commits` derives one artifact per dataset.
The fact stays out of the records, which hold what the depositor and the pipeline
wrote, and out of the projection, which restates them.

A reaction's first commit is the first commit on ``main``, followed by first parent,
whose version of a dataset file under ``data/`` holds the reaction's ID:

* **First parent.** Early submissions merged in two steps, into a submission branch
  and then into ``main``, and the reaction entered the corpus at the second.
* **Every version of every file.** A dataset file's own first commit is not the
  answer: every dataset has been rewritten into another format at least once, and
  USPTO has been resharded more than once, so every version of every file the corpus
  has held is read, in whichever format it was written.
* **Per reaction.** A dataset updated later holds reactions that arrived at different
  commits.
* **The committer date.** When the commit was made on ``main``; for a merged pull
  request, the merge. ``timestamp`` holds its wall clock as git recorded it and
  ``utc_offset_minutes`` the offset git stored with it, the pair every projected time
  is read as, so nothing git records is dropped.
* **The dataset the file recorded,** which need not be the one holding the reaction
  now.

A reaction whose ID a migration reassigned first appears under its current ID at that
migration, since nothing in either file links the two IDs. ord-data holds two such
migrations: the Nanomole flow dataset merged in ord-data#41, given new IDs in #73; and
the USPTO grants merged in #66 through #69, given new IDs by the reshard in #151, #153,
and #155.

The index holds every reaction ID ``main`` has held, not only the current ones, and
records the chain of commits it walked. History only grows, so a later walk extends an
earlier index from its head rather than starting over, and a reaction that leaves the
corpus and returns keeps the commit it first arrived at.

Only the two IDs are read from each version, through a message declaring nothing else
rather than today's ``Dataset``: a version written under an older schema need not parse
under today's, and a full parse decodes every field of every reaction to reach one.
Neither ID has changed its field number or name since before ord-data held a dataset.

The walk needs the whole history as it was committed. A shallow clone presents its
oldest commit as a root, and a grafts file gives a commit whatever parents it names, so
either would place every reaction such a commit holds at that commit; both refuse the
walk. Replace refs, which do the same, are ignored.

Dataset files live in git LFS, and the walk reads only objects already fetched. One
that is not refuses the walk before any dataset file is read, naming the command that
fetches it, rather than fetching from whichever store the checkout points at.
"""

import dataclasses
import gzip
import os
import pathlib
import re
import shutil
import subprocess
from collections.abc import Collection, Iterator, Sequence
from datetime import UTC, datetime, timedelta
from typing import Any, cast

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from google.protobuf import (
    descriptor_pb2,
    descriptor_pool,
    message_factory,
    text_format,
)
from google.protobuf.message import Message

from ord_schema import atomic_io
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Bumped when the walk's definition changes -- which files count, which date is read --
# so an index walked under the old definition, and every artifact stamped with one,
# reads as stale.
VERSION = "1"

# What a first_commits artifact holds per reaction. The time is named as every projected
# time is:
# the wall clock as written, and the offset it was written under in minutes east of UTC.
TYPE = pa.struct(
    [
        pa.field("sha", pa.string()),
        pa.field("timestamp", pa.timestamp("us")),
        pa.field("utc_offset_minutes", pa.int16()),
        pa.field("dataset_id", pa.string()),
    ]
)

# Every reaction ID the walk has met, with its first commit. Nothing is null: a commit
# always has a time, and a file that records no dataset ID refuses the walk.
INDEX_SCHEMA = pa.schema(
    [
        pa.field("reaction_id", pa.string(), nullable=False),
        *(pa.field(field.name, field.type, nullable=False) for field in TYPE),
    ]
)

_META_VERSION = "ord.first_commits_version"
# In an index: every commit walked, oldest first, separated by spaces.
_META_CHAIN = "ord.first_commits_chain"
# In a first_commits artifact: the head of the index it was derived with.
_META_HEAD = "ord.first_commits_head"

# The footer key the Parquet datasets in ord-data's history record their ID under. Named
# here rather than taken from ord_schema.parquet: what matters is what those versions
# wrote, which a later change to the writer would not change.
_PARQUET_DATASET_ID = b"ord.dataset_id"

_DATASET_FILE = re.compile(r"\.(pbtxt|pb|pb\.gz|parquet)$")
_LFS_POINTER = b"version https://git-lfs.github.com/spec/v1\n"
_LFS_OID = re.compile(rb"^oid sha256:([0-9a-f]{64})$", re.MULTILINE)
# A pointer is a few short lines; anything longer is a blob holding its own content.
_LFS_POINTER_MAX_BYTES = 1024


def _ids_message() -> type[Message]:
    """Returns a message class that reads a Dataset's reaction IDs and its own ID alone.

    Every other field is unknown to it, so it parses any version of a Dataset: the
    binary form keeps what it does not declare as unknown bytes, and the text form
    skips it. The numbers are the ones every Dataset in ord-data's history was written
    with.
    """
    # descriptor_pb2 is generated, and ships without the stubs a checker would read.
    field = descriptor_pb2.FieldDescriptorProto  # ty: ignore[unresolved-attribute]
    proto = descriptor_pb2.FileDescriptorProto(  # ty: ignore[unresolved-attribute]
        name="ord_schema/artifacts/first_commits.proto",
        package="ord_schema.first_commits",
        syntax="proto3",
    )
    reaction = proto.message_type.add(name="Reaction")
    reaction.field.add(
        name="reaction_id",
        number=10,
        type=field.TYPE_STRING,
        label=field.LABEL_OPTIONAL,
    )
    dataset = proto.message_type.add(name="Dataset")
    dataset.field.add(
        name="reactions",
        number=3,
        type=field.TYPE_MESSAGE,
        label=field.LABEL_REPEATED,
        type_name=".ord_schema.first_commits.Reaction",
    )
    dataset.field.add(
        name="dataset_id",
        number=5,
        type=field.TYPE_STRING,
        label=field.LABEL_OPTIONAL,
    )
    pool = descriptor_pool.DescriptorPool()  # ty: ignore[possibly-missing-implicit-call]
    pool.Add(proto)
    return message_factory.GetMessageClass(
        pool.FindMessageTypeByName("ord_schema.first_commits.Dataset")
    )


_IDS_MESSAGE = _ids_message()


class Index:
    """Every reaction ID ord-data's ``main`` has held, and its first commit.

    Attributes:
        table: One row per reaction ID, in ``INDEX_SCHEMA``, in the order the walk met
            them.
        chain: Every commit walked, oldest first; the last is the head.
    """

    def __init__(self, table: pa.Table, chain: Sequence[str]) -> None:
        """Initializes an index.

        Args:
            table: Rows in ``INDEX_SCHEMA``.
            chain: Every commit walked, oldest first.

        Raises:
            ValueError: If ``chain`` is empty, since an index of no commits has no head
                to stamp.
        """
        if not chain:
            raise ValueError("an index walks at least one commit")
        self.table = table.cast(INDEX_SCHEMA)
        self.chain = tuple(chain)
        self._walked = frozenset(self.chain)

    @property
    def head(self) -> str:
        """The last commit walked."""
        return self.chain[-1]

    def lookup(self, reaction_ids: Collection[str]) -> dict[str, dict[str, Any]]:
        """Returns the first commit of each reaction ID the index holds.

        Args:
            reaction_ids: The IDs to look up.

        Returns:
            Rows in ``TYPE``, keyed by reaction ID. An ID the index does not hold is
            absent rather than mapped to None.
        """
        if not reaction_ids:
            return {}
        wanted = pa.array(list(reaction_ids), pa.string())
        # pyarrow.compute builds its functions at import; a checker sees none of them.
        held = pc.is_in(  # ty: ignore[unresolved-attribute]
            self.table["reaction_id"], value_set=wanted
        )
        found = self.table.filter(held)
        return {row.pop("reaction_id"): row for row in found.to_pylist()}

    def stamps(self) -> dict[str, str]:
        """Returns the footer keys an artifact derived with this index carries."""
        return {_META_VERSION: VERSION, _META_HEAD: self.head}

    def covers(self, path: str | os.PathLike[str]) -> bool:
        """Returns whether ``path`` was derived with this index or one it extends.

        Every first commit an earlier index records is one this index records too,
        since a walk only appends, so an artifact stamped with any head on this chain
        holds what this index would give it. One stamped with a head off the chain was
        derived from a history that has since been rewritten.

        Args:
            path: A first_commits artifact. Need not exist.

        Returns:
            Whether its stamp names this walk version and a commit on this chain. A
            file stamped by no index, or unreadable, is not covered.
        """
        stamp = read_stamp(path)
        return stamp is not None and stamp[0] == VERSION and stamp[1] in self._walked


def read_stamp(path: str | os.PathLike[str]) -> tuple[str, str] | None:
    """Returns the walk version and index head an artifact was derived with.

    Args:
        path: A first_commits artifact.

    Returns:
        ``(version, head)``, or None where ``path`` carries no index stamp: derived
        without an index, missing, or not Parquet.
    """
    try:
        with pq.ParquetFile(path) as parquet_file:
            metadata = parquet_file.schema_arrow.metadata or {}
    except (OSError, ValueError, pa.ArrowInvalid):
        return None
    version = metadata.get(_META_VERSION.encode())
    head = metadata.get(_META_HEAD.encode())
    if version is None or head is None:
        return None
    return version.decode(), head.decode()


def write_index(index: Index, path: str | os.PathLike[str]) -> None:
    """Writes ``index`` to ``path``, atomically.

    Args:
        index: The index to write.
        path: Where to write it.
    """
    table = index.table.replace_schema_metadata(
        {_META_VERSION: VERSION, _META_CHAIN: " ".join(index.chain)}
    )
    with atomic_io.atomic_path(path) as temp_path:
        pq.write_table(table, temp_path, compression="zstd")


def read_index(path: str | os.PathLike[str]) -> Index:
    """Reads an index written by ``write_index``.

    Args:
        path: The index file.

    Returns:
        The index.

    Raises:
        ValueError: If ``path`` is not an index, or was walked under another
            ``VERSION``, whose first commits this library does not stand behind.
    """
    table = pq.read_table(path)
    metadata = table.schema.metadata or {}
    version = metadata.get(_META_VERSION.encode())
    chain = metadata.get(_META_CHAIN.encode())
    if version is None or chain is None:
        raise ValueError(f"{path} is not a first-commits index")
    if version.decode() != VERSION:
        raise ValueError(
            f"{path} was walked under version {version.decode()}, not {VERSION}; walk "
            "ord-data again from the start"
        )
    return Index(table.replace_schema_metadata(None), chain.decode().split())


@dataclasses.dataclass(frozen=True)
class _Commit:
    """A commit on the walked chain.

    Attributes:
        sha: The commit.
        timestamp: The committer date's wall clock, as git recorded it.
        utc_offset_minutes: The offset git recorded with it, in minutes east of UTC.
    """

    sha: str
    timestamp: datetime
    utc_offset_minutes: int


@dataclasses.dataclass(frozen=True)
class _Version:
    """One version of a dataset file, as one commit changed it.

    Attributes:
        path: The file's path in the repository.
        blob: The blob git stores at that path, which for a file in LFS is its pointer.
        lfs_object: The local LFS object holding the content, or None for a file
            stored in git itself.
    """

    path: str
    blob: str
    lfs_object: pathlib.Path | None


def _git(repository: str | os.PathLike[str], *args: str, stdin: bytes = b"") -> bytes:
    """Runs git in ``repository`` and returns what it printed.

    Args:
        repository: The checkout to run in.
        *args: The git command and its arguments.
        stdin: Bytes to write to the command's standard input.

    Returns:
        The command's standard output.

    Raises:
        ValueError: If git is not installed, or the command fails.
    """
    executable = shutil.which("git")
    if executable is None:
        raise ValueError("walking ord-data's history needs git on PATH")
    # S603: the arguments are this module's own, plus a revision the caller names after
    # --end-of-options, passed as an argument list rather than through a shell.
    result = subprocess.run(  # noqa: S603
        # Replace refs rewrite the history git reports; the walk reads what was
        # committed.
        [executable, "--no-replace-objects", "-C", os.fspath(repository), *args],
        input=stdin,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise ValueError(
            f"git {' '.join(args)} failed in {repository}: "
            f"{result.stderr.decode(errors='replace').strip()}"
        )
    return result.stdout


def _require_whole_history(repository: str | os.PathLike[str]) -> None:
    """Refuses a checkout that reports less history, or other history, than it holds.

    Args:
        repository: The ord-data checkout.

    Raises:
        ValueError: If the checkout is a shallow clone, whose oldest commit reads as a
            root, or holds a grafts file, which gives commits the parents it names.
            Either would date every reaction such a commit holds to that commit.
    """
    if _git(repository, "rev-parse", "--is-shallow-repository").strip() == b"true":
        raise ValueError(
            f"{repository} is a shallow clone, so its oldest commit would read as "
            "where every reaction it holds first arrived; fetch the whole history with "
            f"git -C {repository} fetch --unshallow"
        )
    relative = _git(repository, "rev-parse", "--git-path", "info/grafts").decode()
    grafts = pathlib.Path(repository) / relative.strip()
    if grafts.exists():
        raise ValueError(
            f"{repository} has a grafts file, {grafts}, which rewrites the parents git "
            "reports; convert it with git replace --convert-graft-file, whose replace "
            "refs the walk ignores"
        )


def commit_time(seconds: int, offset: str) -> tuple[datetime, int]:
    """Returns the wall clock and UTC offset of a time git stores as seconds and offset.

    Args:
        seconds: Seconds since the epoch, as ``--date=raw`` prints them.
        offset: The offset git recorded, as ``-0400`` or ``+0130``.

    Returns:
        The wall clock at that offset, naive, and the offset in minutes east of UTC.
        Together they are everything git records, so the instant stays recoverable.

    Raises:
        ValueError: If ``offset`` is not a signed four-digit offset.
    """
    if not re.fullmatch(r"[+-]\d{4}", offset):
        raise ValueError(f"not a git UTC offset: {offset!r}")
    minutes = int(offset[1:3]) * 60 + int(offset[3:5])
    if offset[0] == "-":
        minutes = -minutes
    wall = datetime.fromtimestamp(seconds, UTC) + timedelta(minutes=minutes)
    return wall.replace(tzinfo=None), minutes


def _first_parent_chain(
    repository: str | os.PathLike[str], revision: str
) -> list[_Commit]:
    """Returns the first-parent chain ending at ``revision``, oldest first."""
    lines = _git(
        repository,
        "log",
        "--first-parent",
        "--reverse",
        "--format=%H %cd",
        "--date=raw",
        # A revision beginning with a dash is a revision, not an option.
        "--end-of-options",
        revision,
        "--",
    )
    chain = []
    for line in lines.decode().splitlines():
        sha, seconds, offset = line.split()
        timestamp, minutes = commit_time(int(seconds), offset)
        chain.append(_Commit(sha, timestamp, minutes))
    return chain


def _changed_files(
    repository: str | os.PathLike[str], parent: str | None, commit: str
) -> list[tuple[str, str]]:
    """Returns the dataset files ``commit`` adds or modifies against its first parent.

    Args:
        repository: The ord-data checkout.
        parent: The commit's first parent, or None for the root commit.
        commit: The commit.

    Returns:
        ``(path, blob)`` pairs in path order. A rename reads as the new path added, so
        its content is read like any other version.
    """
    trees = ["--root", commit] if parent is None else [parent, commit]
    output = _git(
        repository,
        "diff-tree",
        "-r",
        "-z",
        "--no-commit-id",
        "--no-renames",
        "--diff-filter=AM",
        *trees,
        "--",
        "data/",
    )
    fields = output.decode().split("\0")
    changed = []
    # Each change is a ":<modes> <old> <new> <status>" record followed by its path.
    for record, path in zip(fields[0:-1:2], fields[1::2], strict=True):
        if _DATASET_FILE.search(path):
            changed.append((path, record.split()[3]))
    return changed


def _read_blobs(
    repository: str | os.PathLike[str], blobs: Collection[str]
) -> Iterator[tuple[str, bytes]]:
    """Yields each blob's content, through one ``git cat-file`` for all of them."""
    if not blobs:
        return
    ordered = sorted(blobs)
    output = _git(
        repository,
        "cat-file",
        "--batch",
        stdin="".join(f"{b}\n" for b in ordered).encode(),
    )
    position = 0
    for blob in ordered:
        header_end = output.index(b"\n", position)
        size = int(output[position:header_end].split()[2])
        start = header_end + 1
        yield blob, output[start : start + size]
        # Each content is followed by a newline of its own.
        position = start + size + 1


def _blob_sizes(
    repository: str | os.PathLike[str], blobs: Collection[str]
) -> dict[str, int]:
    """Returns each blob's size, without reading its content."""
    if not blobs:
        return {}
    output = _git(
        repository,
        "cat-file",
        "--batch-check=%(objectname) %(objectsize)",
        stdin="".join(f"{b}\n" for b in blobs).encode(),
    )
    sizes = {}
    for line in output.decode().splitlines():
        blob, size = line.split()
        sizes[blob] = int(size)
    return sizes


def _lfs_objects_directory(repository: str | os.PathLike[str]) -> pathlib.Path:
    """Returns where the checkout keeps fetched LFS objects, as git-lfs reports it."""
    for line in _git(repository, "lfs", "env").decode().splitlines():
        if line.startswith("LocalMediaDir="):
            return pathlib.Path(line.removeprefix("LocalMediaDir="))
    raise ValueError(f"git lfs env names no LocalMediaDir in {repository}")


def _plan(
    repository: str | os.PathLike[str], chain: Sequence[_Commit], start: int
) -> list[list[_Version]]:
    """Lists the versions each commit from ``start`` on adds, and where each is held.

    Every LFS object is located before any is read, so a walk short of one refuses at
    once rather than after reading everything before it.

    Args:
        repository: The ord-data checkout.
        chain: The whole first-parent chain.
        start: The position of the first commit to walk.

    Returns:
        The versions each walked commit adds or modifies, one list per commit.

    Raises:
        ValueError: If an LFS object a version needs has not been fetched.
    """
    changed = [
        _changed_files(
            repository,
            chain[position - 1].sha if position else None,
            chain[position].sha,
        )
        for position in range(start, len(chain))
    ]
    blobs = {blob for files in changed for _, blob in files}
    sizes = _blob_sizes(repository, blobs)
    small = [blob for blob in blobs if sizes[blob] <= _LFS_POINTER_MAX_BYTES]
    oids: dict[str, str] = {}
    for blob, content in _read_blobs(repository, small):
        if content.startswith(_LFS_POINTER):
            match = _LFS_OID.search(content)
            if match is None:
                raise ValueError(f"LFS pointer {blob} names no sha256 object")
            oids[blob] = match.group(1).decode()
    directory = _lfs_objects_directory(repository) if oids else None
    missing: dict[str, str] = {}
    plan = []
    for offset, files in enumerate(changed):
        versions = []
        for path, blob in files:
            lfs_object = None
            if blob in oids:
                assert directory is not None  # Type hint.
                oid = oids[blob]
                lfs_object = directory / oid[0:2] / oid[2:4] / oid
                if not lfs_object.exists():
                    missing.setdefault(oid, f"{chain[start + offset].sha}:{path}")
            versions.append(_Version(path, blob, lfs_object))
        plan.append(versions)
    if missing:
        first = next(iter(missing.values()))
        raise ValueError(
            f"{len(missing)} LFS objects the walk reads are not in {repository}, the "
            f"first {first}; fetch every commit's objects with "
            f"git -C {repository} lfs fetch <remote> "
            f"$(git -C {repository} rev-list --first-parent {chain[-1].sha})"
        )
    return plan


def read_ids(path: str, content: bytes | pathlib.Path) -> tuple[str, list[str]]:
    """Returns a dataset file's dataset ID and the reaction IDs it holds.

    Args:
        path: The file's name, whose suffix says how it is written.
        content: Its bytes, or the local file holding them.

    Returns:
        The dataset ID the file records, empty where it records none, and every
        reaction ID in file order, empty ones included.

    Raises:
        ValueError: If the suffix is none a dataset file has been written with.
    """
    if path.endswith(".parquet"):
        source = pa.BufferReader(content) if isinstance(content, bytes) else content
        with pq.ParquetFile(source) as parquet_file:
            metadata = parquet_file.schema_arrow.metadata or {}
            ids = parquet_file.read(columns=["reaction_id"]).column(0).to_pylist()
        return metadata.get(_PARQUET_DATASET_ID, b"").decode(), ids
    data = content if isinstance(content, bytes) else content.read_bytes()
    message = _IDS_MESSAGE()
    if path.endswith(".pbtxt"):
        text_format.Parse(data.decode(), message, allow_unknown_field=True)
    elif path.endswith(".pb.gz"):
        message.ParseFromString(gzip.decompress(data))
    elif path.endswith(".pb"):
        message.ParseFromString(data)
    else:
        raise ValueError(f"{path} is not a dataset file")
    # The message class is built at runtime, so a checker sees none of its fields.
    dataset = cast(Any, message)
    return dataset.dataset_id, [reaction.reaction_id for reaction in dataset.reactions]


def walk(
    repository: str | os.PathLike[str],
    *,
    revision: str = "main",
    previous: Index | None = None,
) -> Index:
    """Walks ord-data's first-parent history into an index of first commits.

    Args:
        repository: An ord-data checkout holding every LFS object the walk reads.
        revision: The branch or commit to walk to, followed by first parent from the
            root.
        previous: An index to extend, walked to a commit on the same chain. Only the
            commits after its head are read.

    Returns:
        An index of every reaction ID the chain's dataset files have held, each at the
        first commit that held it, walked to ``revision``.

    Raises:
        ValueError: If the checkout's history is shallow or grafted, if ``previous``
            was walked along a chain that is not the start of this one, if an LFS object
            is missing, if a file records no dataset ID, or
            if a reaction first appears in two datasets in one commit, which would
            leave its first dataset unknowable.
    """
    _require_whole_history(repository)
    chain = _first_parent_chain(repository, revision)
    start = 0
    seen: set[str] = set()
    if previous is not None:
        walked = [commit.sha for commit in chain[: len(previous.chain)]]
        if tuple(walked) != previous.chain:
            raise ValueError(
                f"the index was walked to {previous.head}, which is not on the "
                f"first-parent chain of {revision} in {repository}; walk again from "
                "the start"
            )
        start = len(previous.chain)
        seen.update(previous.table["reaction_id"].to_pylist())
    plan = _plan(repository, chain, start)
    reaction_ids: list[str] = []
    positions: list[int] = []
    dataset_ids: list[str] = []
    read: set[str] = set()
    for offset, versions in enumerate(plan):
        commit = chain[start + offset]
        arrived: dict[str, str] = {}
        pending = [version for version in versions if version.blob not in read]
        stored = {version.blob for version in pending if version.lfs_object is None}
        contents = dict(_read_blobs(repository, stored))
        for version in pending:
            read.add(version.blob)
            content = (
                contents[version.blob]
                if version.lfs_object is None
                else version.lfs_object
            )
            dataset_id, ids = read_ids(version.path, content)
            if not dataset_id:
                raise ValueError(f"{commit.sha}:{version.path} records no dataset ID")
            for reaction_id in ids:
                # A submission staged before IDs were assigned holds empty ones, which
                # name no reaction.
                if not reaction_id or reaction_id in seen:
                    continue
                first = arrived.setdefault(reaction_id, dataset_id)
                if first != dataset_id:
                    raise ValueError(
                        f"{reaction_id} first appears in {first} and {dataset_id}, "
                        f"both at {commit.sha}"
                    )
        if arrived:
            logger.info("%s: %d reactions arrive", commit.sha[:12], len(arrived))
        for reaction_id, dataset_id in arrived.items():
            reaction_ids.append(reaction_id)
            positions.append(offset)
            dataset_ids.append(dataset_id)
        seen.update(arrived)
    walked_commits = chain[start:]
    take = pa.array(positions, pa.int64())
    added = pa.Table.from_arrays(
        [
            pa.array(reaction_ids, pa.string()),
            pa.array([c.sha for c in walked_commits], pa.string()).take(take),
            pa.array([c.timestamp for c in walked_commits], pa.timestamp("us")).take(
                take
            ),
            pa.array([c.utc_offset_minutes for c in walked_commits], pa.int16()).take(
                take
            ),
            pa.array(dataset_ids, pa.string()),
        ],
        schema=INDEX_SCHEMA,
    )
    table = added if previous is None else pa.concat_tables([previous.table, added])
    return Index(table, [commit.sha for commit in chain])
