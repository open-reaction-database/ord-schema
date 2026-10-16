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
"""Fetches datasets and derived artifacts from Hugging Face.

Source datasets come from the ord-data mirror. The derived artifacts of
:mod:`ord_schema.artifacts` live in a repository of their own, ord-artifacts, so that
rebuilding them leaves the mirror's history alone; they are published by
``ord_schema.scripts.upload_artifacts``.

Importing this module requires the optional ``huggingface`` extra (``pip install 'ord-
schema[huggingface]'``); ``huggingface_hub`` is not a default dependency because it
pulls in a heavy chain (typer, fsspec, the hf-xet native binary).
"""

import dataclasses
import glob
import os
import pathlib
from collections.abc import Collection, Iterable

import pyarrow.parquet as pq
from huggingface_hub import (
    HfFileSystem,
    hf_hub_download,
    list_repo_files,
    repo_info,
    snapshot_download,
)
from huggingface_hub.errors import EntryNotFoundError
from rdkit import rdBase

from ord_schema import validations
from ord_schema.artifacts import base
from ord_schema.message_helpers import id_filename

# Hugging Face dataset that mirrors the ord-data repository; see
# https://huggingface.co/datasets/open-reaction-database/ord-data.
ORD_DATA_HF_REPO = "open-reaction-database/ord-data"

# Hugging Face dataset holding the artifacts derived from ord-data; see
# https://huggingface.co/datasets/open-reaction-database/ord-artifacts.
ORD_ARTIFACTS_HF_REPO = "open-reaction-database/ord-artifacts"

# The directory each artifact is kept under, both in the repository and in a local tree.
# The repository holds one build at a time, laid out exactly as a local tree is.
ARTIFACT_DIRECTORIES = {
    "projection": "projections",
    "structures": "structures",
    "pivot": "pivots",
    "occurrences": "occurrences",
}

# Small, so reading a footer in place fetches little beyond it.
_FOOTER_BLOCK_SIZE = 64 * 1024


def fetch_dataset(
    dataset_id: str,
    *,
    revision: str = "main",
    cache_dir: str | None = None,
    local_dir: str | None = None,
) -> str:
    """Downloads a dataset file from the ord-data Hugging Face mirror.

    Prefers the Parquet serialization and falls back to the legacy ``.pb.gz``
    binary format for datasets that have not been converted yet. The file is
    not parsed; dispatch on the returned suffix to choose a reader, e.g.
    ``parquet.load_dataset``/``parquet.DatasetView`` for
    ``.parquet`` and ``load_message`` for ``.pb.gz``.

    Downloads are cached and content-verified by ``huggingface_hub``: repeat
    calls reuse the cached file (re-validating against the remote for moving
    revisions like ``"main"``) instead of re-downloading. By default the cache
    lives under ``~/.cache/huggingface/hub``; set ``$HF_HOME`` (or pass
    ``cache_dir``) to relocate it, or set ``$HF_HUB_OFFLINE=1`` to reuse an
    already-warm cache without network access.

    Args:
        dataset_id: Dataset ID.
        revision: Branch, tag, or commit SHA to download (default: ``"main"``).
        cache_dir: Override the Hugging Face cache directory. Files are stored
            content-addressed and returned via a symlink into the cache.
        local_dir: If set, materialize the file as a real (non-symlink) copy
            under this directory instead of returning a cache symlink; useful
            for staging files into a project tree.

    Returns:
        Local filesystem path to the downloaded dataset file.

    Raises:
        RuntimeError: If no file exists for the dataset.
        ValueError: If the dataset ID is invalid.
    """
    if not validations.is_valid_dataset_id(dataset_id):
        raise ValueError(f"Invalid dataset ID: {dataset_id}")
    for suffix in (".parquet", ".pb.gz"):
        try:
            return hf_hub_download(
                repo_id=ORD_DATA_HF_REPO,
                repo_type="dataset",
                revision=revision,
                filename=id_filename(f"{dataset_id}{suffix}"),
                cache_dir=cache_dir,
                local_dir=local_dir,
            )
        except EntryNotFoundError:
            continue
    raise RuntimeError(
        f"Dataset {dataset_id} not found in {ORD_DATA_HF_REPO}@{revision}"
    )


@dataclasses.dataclass(frozen=True)
class FetchedArtifacts:
    """Where ``fetch_artifacts`` put each artifact, in the shape ``Corpus`` reads.

    Attributes:
        revision: The commit the artifacts were fetched from.
        projections: Directory holding the projections.
        structures: Directory holding the structures artifacts.
        pivots: Directory holding one subdirectory per pivoted level.
        occurrences: Directory holding one subdirectory per indexed path.
    """

    revision: str
    projections: pathlib.Path
    structures: pathlib.Path
    pivots: pathlib.Path
    occurrences: pathlib.Path

    @property
    def projection_pattern(self) -> str:
        """A glob matching every projection, for ``Corpus``."""
        return parquet_pattern(self.projections)

    @property
    def structures_pattern(self) -> str:
        """A glob matching every structures artifact, for ``Corpus``."""
        return parquet_pattern(self.structures)


def parquet_pattern(directory: pathlib.Path) -> str:
    """Returns a recursive glob for the Parquet files under ``directory``, escaped."""
    return f"{glob.escape(str(directory))}/**/*.parquet"


def parquet_files(directory: pathlib.Path) -> list[pathlib.Path]:
    """Returns the Parquet files under ``directory`` the way a ``Corpus`` finds them.

    Through ``glob``, as the readers do, so a linked shard directory is descended and a
    dotted file is passed over.

    Args:
        directory: The directory to search beneath.

    Returns:
        The paths, sorted. Empty if the directory does not exist.
    """
    return sorted(
        pathlib.Path(found)
        for found in glob.glob(parquet_pattern(directory), recursive=True)
    )


def resolve_revision(revision: str) -> str:
    """Returns the commit ``revision`` names in ord-artifacts.

    Args:
        revision: Branch, tag, or commit SHA.

    Returns:
        The commit SHA.

    Raises:
        LookupError: If the Hub reports no commit for it.
    """
    commit = repo_info(
        ORD_ARTIFACTS_HF_REPO, repo_type="dataset", revision=revision
    ).sha
    if commit is None:
        raise LookupError(f"{ORD_ARTIFACTS_HF_REPO} reports no commit for {revision}")
    return commit


def published_build(commit: str, published: Iterable[str]) -> base.Stamps | None:
    """Returns the stamps of the build ord-artifacts holds at ``commit``, read in place.

    Every file of a published build carries the same artifact schema version and RDKit
    version, which ``upload_artifacts`` checks, so one footer speaks for the build.
    It is an occurrences file's, whose few columns make for a small footer, and
    nothing else of the file is downloaded.

    Args:
        commit: The commit SHA ``published`` was listed at.
        published: Every repository path the commit holds.

    Returns:
        The stamps, or None if the commit holds no occurrences artifacts.
    """
    paths = sorted(_parquet_beneath(published, ARTIFACT_DIRECTORIES["occurrences"]))
    if not paths:
        return None
    location = f"datasets/{ORD_ARTIFACTS_HF_REPO}@{commit}/{paths[0]}"
    with (
        HfFileSystem().open(location, block_size=_FOOTER_BLOCK_SIZE) as handle,
        pq.ParquetFile(handle) as parquet_file,
    ):
        return base.stamps_from_metadata(parquet_file.schema_arrow.metadata or {})


def fetch_artifacts(
    local_dir: str | os.PathLike[str],
    *,
    revision: str = "main",
    pivot_levels: Collection[str] | None = None,
) -> FetchedArtifacts:
    """Downloads the derived artifacts ord-artifacts holds, if this library reads them.

    The repository holds one build at a time, laid out as a local tree. A ``Corpus``
    refuses artifacts whose schema version or RDKit differs from its own, so the build's
    stamps are read in place first and a mismatch is refused before anything is
    downloaded: a tree holding a build this library reads is not overwritten with one it
    cannot. An earlier build is fetched by passing the commit that holds it. Files are
    written under ``local_dir`` rather than linked into the Hugging Face cache, because
    a ``Corpus`` confined to its own directories refuses a link leading out of them.

    Repeat calls download only what changed, and delete the Parquet files under what
    they fetched that the build no longer holds -- a dataset removed upstream would
    otherwise stay in every local tree that once fetched it. Pivot levels already held
    under ``local_dir`` are refreshed with the requested ones, since a ``Corpus``
    refuses a level left stale beside refreshed projections.

    Args:
        local_dir: Directory to download beneath; the repository's layout is kept.
        revision: Branch, tag, or commit SHA to download (default: ``"main"``). A
            branch is resolved to its commit first, so the build is checked and
            downloaded from the same one.
        pivot_levels: The pivoted levels to fetch. None fetches every level, which is
            what ``Corpus(require_pivots=True)`` needs; a level left out is built in
            process when a query first ranges over it.

    Returns:
        The commit fetched, and where each artifact landed.

    Raises:
        LookupError: If the revision holds no artifacts of some kind, holds a build of
            another artifact schema version or RDKit, or lacks a requested pivot level.
            The message names what it does hold.
        ValueError: If a file the build no longer holds sits in a directory reached
            through a link beneath ``local_dir``, which is not this function's to
            delete.
    """
    root = pathlib.Path(local_dir)
    commit = resolve_revision(revision)
    published = list_repo_files(
        ORD_ARTIFACTS_HF_REPO, repo_type="dataset", revision=commit
    )
    missing = [
        artifact
        for artifact, directory in ARTIFACT_DIRECTORIES.items()
        if not _parquet_beneath(published, directory)
    ]
    # With every artifact present, the occurrences have a footer to read.
    stamps = None if missing else published_build(commit, published)
    if stamps is None:
        raise LookupError(
            f"{ORD_ARTIFACTS_HF_REPO}@{revision} holds no "
            f"{', '.join(missing)} artifacts"
        )
    if (stamps.artifact_schema_version, stamps.rdkit_version) != (
        base.ARTIFACT_SCHEMA_VERSION,
        rdBase.rdkitVersion,
    ):
        raise LookupError(
            f"{ORD_ARTIFACTS_HF_REPO}@{revision} holds artifact schema version "
            f"{stamps.artifact_schema_version} derived under RDKit "
            f"{stamps.rdkit_version}, and this library reads schema version "
            f"{base.ARTIFACT_SCHEMA_VERSION} under RDKit "
            f"{rdBase.rdkitVersion}. Install the ord-schema and RDKit that build was "
            "derived with, or pass revision= an earlier commit holding a build this "
            "library reads, if one was published; each upload's commit description "
            "names the schema version and RDKit it holds."
        )
    pivots = ARTIFACT_DIRECTORIES["pivot"]
    scopes = [
        directory
        for artifact, directory in ARTIFACT_DIRECTORIES.items()
        if artifact != "pivot"
    ]
    if pivot_levels is None:
        scopes.append(pivots)
    else:
        levels = {
            path[len(pivots) + 1 :].split("/")[0]
            for path in _beneath(published, pivots)
        }
        missing = sorted(set(pivot_levels) - levels)
        if missing:
            raise LookupError(
                f"{ORD_ARTIFACTS_HF_REPO}@{revision} has no pivots for {missing}; it "
                f"has {sorted(levels)}"
            )
        # A dotted directory is some tool's, not a level; no reader looks inside one.
        held = (
            {
                entry.name
                for entry in (root / pivots).iterdir()
                if entry.is_dir() and not entry.name.startswith(".")
            }
            if (root / pivots).is_dir()
            else set()
        )
        scopes.extend(f"{pivots}/{level}" for level in sorted(set(pivot_levels) | held))
    snapshot_download(
        repo_id=ORD_ARTIFACTS_HF_REPO,
        repo_type="dataset",
        revision=commit,
        local_dir=root,
        allow_patterns=[f"{scope}/**" for scope in scopes],
    )
    _prune(root, scopes, set(published))
    return FetchedArtifacts(
        revision=commit,
        projections=root / ARTIFACT_DIRECTORIES["projection"],
        structures=root / ARTIFACT_DIRECTORIES["structures"],
        pivots=root / pivots,
        occurrences=root / ARTIFACT_DIRECTORIES["occurrences"],
    )


def _prune(root: pathlib.Path, scopes: Iterable[str], published: set[str]) -> None:
    """Deletes the Parquet files under ``scopes`` that ``published`` does not hold.

    A reader follows a linked directory, so a stale file inside one is still read; but
    deleting through the link would reach storage this tree does not own. Such a file
    is refused instead. A file that is itself a link is removed as a link.

    Args:
        root: The local tree, laid out as the repository is.
        scopes: Repository directories that were fetched.
        published: Every repository path the fetched commit holds.

    Raises:
        ValueError: If a file to delete sits in a directory reached through a link
            anywhere beneath ``root``.
    """
    foreign = []
    for scope in scopes:
        # Resolved only as far as root, which the caller chose: a link below it, the
        # scope's own directories included, leads to storage this tree does not own.
        directory = root.resolve() / scope
        for path in parquet_files(root / scope):
            if path.relative_to(root).as_posix() in published:
                continue
            if path.parent.resolve().is_relative_to(directory):
                path.unlink()
            else:
                foreign.append(str(path))
    if foreign:
        raise ValueError(
            "the fetched build no longer holds these files, which sit in directories "
            f"linked from outside the tree and were left in place: {foreign}"
        )


def _beneath(paths: Iterable[str], prefix: str) -> list[str]:
    """Returns the repository paths inside the directory ``prefix``."""
    return [path for path in paths if path.startswith(f"{prefix}/")]


def _parquet_beneath(paths: Iterable[str], prefix: str) -> list[str]:
    """Returns the repository's Parquet files inside the directory ``prefix``."""
    return [path for path in _beneath(paths, prefix) if path.endswith(".parquet")]
