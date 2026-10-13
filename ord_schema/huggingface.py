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
"""Fetches datasets and derived artifacts from Hugging Face, and publishes artifacts.

Source datasets come from the ord-data mirror. The derived artifacts of
:mod:`ord_schema.artifacts` live in a repository of their own, ord-artifacts, so that
rebuilding them leaves the mirror's history alone.

Importing this module requires the optional ``huggingface`` extra (``pip install 'ord-
schema[huggingface]'``); ``huggingface_hub`` is not a default dependency because it
pulls in a heavy chain (typer, fsspec, the hf-xet native binary).
"""

import collections
import dataclasses
import glob
import os
import pathlib
from collections.abc import Collection, Iterable

from huggingface_hub import (
    CommitOperationAdd,
    CommitOperationDelete,
    create_commit,
    hf_hub_download,
    list_repo_files,
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
ARTIFACT_DIRECTORIES = {
    "projection": "projections",
    "structures": "structures",
    "pivot": "pivots",
    "occurrences": "occurrences",
}


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


def artifact_prefix(artifact: str, lineage: str, rdkit_version: str) -> str:
    """Returns the repository directory holding one artifact as one library derives it.

    Keyed by the two stamps a reader compares besides the source: the artifact lineage
    and the RDKit version. Builds from different library versions therefore sit side by
    side, and a reader fetches only the one it would accept. The lineage's ``=`` and
    ``,`` are respelled because DuckDB reads a ``key=value`` directory as a Hive
    partition and adds a column for it to every file beneath.

    Args:
        artifact: Artifact name, a key of ``ARTIFACT_DIRECTORIES``.
        lineage: The artifact lineage, as ``base.lineage`` spells it.
        rdkit_version: The RDKit version, as ``rdBase.rdkitVersion`` spells it.

    Returns:
        A path relative to the repository root, e.g.
        ``structures/structures-1_projection-1/rdkit-2026.03.5``.
    """
    key = lineage.replace("=", "-").replace(",", "_")
    return f"{ARTIFACT_DIRECTORIES[artifact]}/{key}/rdkit-{rdkit_version}"


def current_prefix(artifact: str) -> str:
    """Returns the repository directory holding ``artifact`` as this library derives it.

    Args:
        artifact: Artifact name, a key of ``ARTIFACT_DIRECTORIES``.

    Returns:
        The prefix for this library's lineage and the installed RDKit.
    """
    return artifact_prefix(artifact, base.lineage(artifact), rdBase.rdkitVersion)


@dataclasses.dataclass(frozen=True)
class FetchedArtifacts:
    """Where ``fetch_artifacts`` put each artifact, in the shape ``Corpus`` reads.

    Attributes:
        projections: Directory holding the projections.
        structures: Directory holding the structures artifacts.
        pivots: Directory holding one subdirectory per pivoted level.
        occurrences: Directory holding one subdirectory per indexed path.
    """

    projections: pathlib.Path
    structures: pathlib.Path
    pivots: pathlib.Path
    occurrences: pathlib.Path

    @property
    def projection_pattern(self) -> str:
        """A glob matching every projection, for ``Corpus``."""
        return _parquet_pattern(self.projections)

    @property
    def structures_pattern(self) -> str:
        """A glob matching every structures artifact, for ``Corpus``."""
        return _parquet_pattern(self.structures)


def _parquet_pattern(directory: pathlib.Path) -> str:
    """Returns a recursive glob for the Parquet files under ``directory``."""
    return f"{glob.escape(str(directory))}/**/*.parquet"


def fetch_artifacts(
    local_dir: str | os.PathLike[str],
    *,
    revision: str = "main",
    pivot_levels: Collection[str] | None = None,
) -> FetchedArtifacts:
    """Downloads the derived artifacts this library and RDKit read as current.

    A ``Corpus`` refuses artifacts whose lineage or RDKit stamp differs from its own, so
    only the build matching both is fetched. Files are written under ``local_dir``
    rather than linked into the Hugging Face cache, because a ``Corpus`` confined to its
    own directories refuses a link leading out of them. Repeat calls download only what
    changed.

    Args:
        local_dir: Directory to download beneath; the repository's layout is kept.
        revision: Branch, tag, or commit SHA to download (default: ``"main"``).
        pivot_levels: The pivoted levels to fetch. None fetches every level published,
            which is what ``Corpus(require_pivots=True)`` needs; a level left out is
            built in process when a query first ranges over it.

    Returns:
        Where each artifact landed.

    Raises:
        LookupError: If the revision holds no build of an artifact for this library and
            RDKit, or lacks a requested pivot level. The message lists the builds and
            levels it does hold.
    """
    root = pathlib.Path(local_dir)
    prefixes = {artifact: current_prefix(artifact) for artifact in ARTIFACT_DIRECTORIES}
    published = list_repo_files(
        ORD_ARTIFACTS_HF_REPO, repo_type="dataset", revision=revision
    )
    builds = _published_builds(published)
    for artifact, prefix in prefixes.items():
        if prefix not in builds:
            available = sorted(
                build
                for build in builds
                if build.startswith(f"{ARTIFACT_DIRECTORIES[artifact]}/")
            )
            raise LookupError(
                f"{ORD_ARTIFACTS_HF_REPO}@{revision} holds no {artifact} artifacts for "
                f"this library and RDKit {rdBase.rdkitVersion} (wanted {prefix}); it "
                f"holds {available or 'none'}"
            )
    pivots = prefixes["pivot"]
    patterns = [
        f"{prefixes[artifact]}/**" for artifact in prefixes if artifact != "pivot"
    ]
    if pivot_levels is None:
        patterns.append(f"{pivots}/**")
    else:
        levels = _subdirectories(published, pivots)
        missing = sorted(set(pivot_levels) - levels)
        if missing:
            raise LookupError(
                f"{ORD_ARTIFACTS_HF_REPO}@{revision} has no pivots for {missing} under "
                f"{pivots}; it has {sorted(levels)}"
            )
        patterns.extend(f"{pivots}/{level}/**" for level in sorted(pivot_levels))
    snapshot_download(
        repo_id=ORD_ARTIFACTS_HF_REPO,
        repo_type="dataset",
        revision=revision,
        local_dir=root,
        allow_patterns=patterns,
    )
    return FetchedArtifacts(
        projections=root / prefixes["projection"],
        structures=root / prefixes["structures"],
        pivots=root / pivots,
        occurrences=root / prefixes["occurrences"],
    )


def _published_builds(paths: Iterable[str]) -> set[str]:
    """Returns the ``<artifact>/<lineage>/rdkit-<version>`` prefixes among ``paths``."""
    directories = set(ARTIFACT_DIRECTORIES.values())
    builds = set()
    for path in paths:
        parts = path.split("/")
        if len(parts) > 3 and parts[0] in directories:
            builds.add("/".join(parts[:3]))
    return builds


def _subdirectories(paths: Iterable[str], prefix: str) -> set[str]:
    """Returns the names of the directories immediately beneath ``prefix``."""
    found = set()
    for path in paths:
        if path.startswith(f"{prefix}/"):
            rest = path[len(prefix) + 1 :].split("/")
            if len(rest) > 1:
                found.add(rest[0])
    return found


@dataclasses.dataclass(frozen=True)
class ArtifactUpload:
    """What ``upload_artifacts`` sent, or would send.

    Attributes:
        prefixes: The repository directory each artifact went under.
        uploaded: Repository paths written. Unchanged files are among them;
            ``huggingface_hub`` drops those from the commit itself.
        deleted: Repository paths removed because the local tree no longer holds them.
        commit_url: The commit made, the head commit when nothing changed, or None
            for a dry run.
    """

    prefixes: dict[str, str]
    uploaded: list[str]
    deleted: list[str]
    commit_url: str | None


def upload_artifacts(
    root: str | os.PathLike[str],
    *,
    commit_message: str,
    revision: str = "main",
    dry_run: bool = False,
) -> ArtifactUpload:
    """Publishes a local artifact tree to ord-artifacts in one commit.

    ``root`` holds one directory per artifact, named as in ``ARTIFACT_DIRECTORIES``,
    each laid out as its derive script wrote it. Every file goes under the prefix its
    own footer names, so a file cannot be published as a build it does not belong to.

    Each artifact has to be one build -- one lineage and one RDKit version across its
    files -- and every artifact, pivot level, and indexed path has to restate the same
    source datasets as the projections. Anything else would publish a build that a
    ``Corpus`` refuses or answers from partially. Within each prefix uploaded, a file
    the local tree no longer holds is deleted; other builds are left alone.

    Args:
        root: The local tree.
        commit_message: Summary for the commit, e.g. naming the ord-data revision the
            artifacts were derived from.
        revision: Branch to commit to.
        dry_run: Check the tree and work out the commit without making it.

    Returns:
        What was, or would be, uploaded and deleted.

    Raises:
        ValueError: If an artifact directory is missing or empty, holds a file that is
            not that artifact or carries no RDKit stamp, spans more than one build, or
            restates different source datasets than the projections.
    """
    root = pathlib.Path(root)
    prefixes: dict[str, str] = {}
    files: dict[str, dict[str, pathlib.Path]] = {}
    sources: dict[str, set[str]] = {}
    for artifact, directory in ARTIFACT_DIRECTORIES.items():
        prefix, found, covered = _check_build(root / directory, artifact)
        prefixes[artifact] = prefix
        files[artifact] = found
        sources.update(
            {
                f"{directory}/{group}" if group else directory: md5s
                for group, md5s in covered.items()
            }
        )
    expected = sources[ARTIFACT_DIRECTORIES["projection"]]
    mismatched = sorted(name for name, md5s in sources.items() if md5s != expected)
    if mismatched:
        details = "; ".join(
            f"{name} lacks {len(expected - sources[name])} and adds "
            f"{len(sources[name] - expected)}"
            for name in mismatched
        )
        raise ValueError(
            "every artifact must restate the same source datasets as the "
            f"{len(expected)} projections: {details}"
        )
    published = set(
        list_repo_files(ORD_ARTIFACTS_HF_REPO, repo_type="dataset", revision=revision)
    )
    operations: list[CommitOperationAdd | CommitOperationDelete] = []
    uploaded: list[str] = []
    deleted: list[str] = []
    for artifact, prefix in prefixes.items():
        wanted = {
            f"{prefix}/{relative}": path for relative, path in files[artifact].items()
        }
        for path_in_repo, path in sorted(wanted.items()):
            operations.append(
                CommitOperationAdd(path_in_repo=path_in_repo, path_or_fileobj=str(path))
            )
            uploaded.append(path_in_repo)
        for path_in_repo in sorted(published):
            if path_in_repo.startswith(f"{prefix}/") and path_in_repo not in wanted:
                operations.append(CommitOperationDelete(path_in_repo=path_in_repo))
                deleted.append(path_in_repo)
    commit_url = None
    if not dry_run:
        info = create_commit(
            ORD_ARTIFACTS_HF_REPO,
            operations,
            commit_message=commit_message,
            repo_type="dataset",
            revision=revision,
        )
        commit_url = info.commit_url
    return ArtifactUpload(
        prefixes=prefixes, uploaded=uploaded, deleted=deleted, commit_url=commit_url
    )


def _check_build(
    directory: pathlib.Path, artifact: str
) -> tuple[str, dict[str, pathlib.Path], dict[str, set[str]]]:
    """Checks that one artifact's tree is a single build, and says what it covers.

    Args:
        directory: The artifact's local directory.
        artifact: The artifact it should hold.

    Returns:
        The build's repository prefix; its files keyed by path relative to
        ``directory``; and the source MD5s each group covers, where a group is the first
        directory beneath ``directory`` for the artifacts that have one per pivot level
        or indexed path, and the empty string for the others.

    Raises:
        ValueError: If the tree is missing or empty, holds a file that is not
            ``artifact`` or carries no RDKit stamp, or spans more than one build.
    """
    paths = sorted(
        pathlib.Path(found)
        for found in glob.glob(_parquet_pattern(directory), recursive=True)
    )
    if not paths:
        raise ValueError(f"no {artifact} artifacts under {directory}")
    grouped = artifact in ("pivot", "occurrences")
    builds: collections.Counter[str] = collections.Counter()
    files: dict[str, pathlib.Path] = {}
    covered: dict[str, set[str]] = collections.defaultdict(set)
    for path in paths:
        stamps = base.load_stamps(path)
        if stamps.artifact != artifact:
            raise ValueError(f"{path} holds {stamps.artifact}, not {artifact}")
        if stamps.rdkit_version is None:
            raise ValueError(f"{path} carries no RDKit stamp; derive it again")
        builds[
            artifact_prefix(artifact, stamps.artifact_lineage, stamps.rdkit_version)
        ] += 1
        relative = path.relative_to(directory)
        files[relative.as_posix()] = path
        covered[relative.parts[0] if grouped else ""].add(stamps.source_md5)
    if len(builds) > 1:
        raise ValueError(
            f"{directory} spans more than one build, so part of it is stale: "
            f"{dict(builds.most_common())}"
        )
    return next(iter(builds)), files, dict(covered)
