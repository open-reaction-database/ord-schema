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

r"""Publishes a derived artifact tree to the ord-artifacts Hugging Face repository.

The tree holds one directory per artifact, as the derive scripts write them::

    upload_artifacts.py --root=artifacts \
        --message="Derive from ord-data <commit>" --dry_run

    artifacts/projections/aa/ord_dataset-<id>.parquet
        -> projections/projection-1/rdkit-<version>/aa/ord_dataset-<id>.parquet
    artifacts/pivots/inputs.components/aa/ord_dataset-<id>.parquet
        -> pivots/pivot-1_projection-1/rdkit-<version>/inputs.components/aa/...

Each file goes under the build its own footer names, everything lands in one commit,
and files a build no longer holds are deleted from it. The tree has to be one build a
``Corpus(require_pivots=True, require_occurrences=True)`` opens, so it is refused if:

* any artifact spans more than one lineage or RDKit version, the artifacts disagree on
  the RDKit version, or one was derived from another version of its parent;
* a pivot or occurrences file sits under a level or path other than the one its footer
  names, or a level or path this library derives is missing;
* an artifact, pivot level, or indexed path restates a source dataset twice, or other
  datasets than the projections;
* a ``Corpus`` opened on it refuses it, which is what sees a structures artifact too
  short for its projection's IDs and a file this library and RDKit would not write.
  Run the upload under the library that derived the tree. Over ORD this takes about
  30 seconds and 4 GiB.

Needs the ``search`` and ``huggingface`` extras. Authenticates the way
``huggingface_hub`` does: ``hf auth login`` or ``$HF_TOKEN``.
"""

import argparse
import collections
import dataclasses
import os
import pathlib
from collections.abc import Collection

from huggingface_hub import (
    CommitOperationAdd,
    CommitOperationDelete,
    create_commit,
    list_repo_files,
)

from ord_schema import huggingface
from ord_schema.artifacts import base, occurrences, pivot
from ord_schema.logging import get_logger
from ord_schema.search import execute

logger = get_logger(__name__)


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


@dataclasses.dataclass(frozen=True)
class _Build:
    """One artifact's local tree, checked to be a single build.

    Attributes:
        lineage: The artifact lineage every file carries.
        rdkit_version: The RDKit version every file carries.
        files: The files, keyed by path relative to the artifact's directory.
        covered: The source MD5s each group restates. A group is the level or indexed
            path for pivots and occurrences, and the empty string otherwise.
    """

    lineage: str
    rdkit_version: str
    files: dict[str, pathlib.Path]
    covered: dict[str, set[str]]


# Where a grouped artifact's footer names the group its directory has to match.
_GROUP_STAMPS = {"pivot": pivot.pivot_path, "occurrences": occurrences.occurrence_path}


def upload_artifacts(
    root: str | os.PathLike[str],
    *,
    commit_message: str,
    revision: str = "main",
    pivot_levels: Collection[str] = tuple(sorted(pivot.LEVELS)),
    occurrence_paths: Collection[str] = tuple(sorted(occurrences.PATHS)),
    dry_run: bool = False,
) -> ArtifactUpload:
    """Publishes a local artifact tree to ord-artifacts in one commit.

    ``root`` holds one directory per artifact, named as in
    ``huggingface.ARTIFACT_DIRECTORIES``, each laid out as its derive script wrote it.
    Every file goes under the prefix its own footer names, so a file cannot be
    published as a build it does not belong to. Within each prefix uploaded, a file the
    local tree no longer holds is deleted; other builds are left alone.

    The commit names the head it was planned against as its parent, so one that lands
    in between fails this one rather than having its files deleted.

    Args:
        root: The local tree.
        commit_message: Summary for the commit, e.g. naming the ord-data revision the
            artifacts were derived from.
        revision: Branch to commit to.
        pivot_levels: The levels the pivots must cover, exactly; every level this
            library derives by default.
        occurrence_paths: The indexed paths the occurrences must cover, exactly; every
            path this library derives by default.
        dry_run: Check the tree and work out the commit without making it.

    Returns:
        What was, or would be, uploaded and deleted.

    Raises:
        ValueError: If the tree is not one complete build; the module docstring lists
            what that requires. A ``Corpus`` refusing it raises ``PairingError``, a
            ``ValueError``.
    """
    root = pathlib.Path(root)
    required = {"pivot": pivot_levels, "occurrences": occurrence_paths}
    builds = {
        artifact: _check_build(root / directory, artifact)
        for artifact, directory in huggingface.ARTIFACT_DIRECTORIES.items()
    }
    _check_consistent(builds)
    sources: dict[str, set[str]] = {}
    for artifact, build in builds.items():
        directory = huggingface.ARTIFACT_DIRECTORIES[artifact]
        groups = required.get(artifact)
        held = set(build.covered)
        if groups is not None and held != set(groups):
            raise ValueError(
                f"{directory} lacks {sorted(set(groups) - held)} and adds "
                f"{sorted(held - set(groups))}"
            )
        for group, md5s in build.covered.items():
            sources[f"{directory}/{group}" if group else directory] = md5s
    expected = sources[huggingface.ARTIFACT_DIRECTORIES["projection"]]
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
    _check_opens(
        root,
        require_pivots=set(pivot_levels) == set(pivot.LEVELS),
        require_occurrences=set(occurrence_paths) == set(occurrences.PATHS),
    )
    head = huggingface.resolve_revision(revision)
    published = list_repo_files(
        huggingface.ORD_ARTIFACTS_HF_REPO, repo_type="dataset", revision=head
    )
    prefixes = {
        artifact: huggingface.artifact_prefix(
            artifact, build.lineage, build.rdkit_version
        )
        for artifact, build in builds.items()
    }
    operations: list[CommitOperationAdd | CommitOperationDelete] = []
    uploaded: list[str] = []
    deleted: list[str] = []
    for artifact, prefix in prefixes.items():
        wanted = {
            f"{prefix}/{relative}": path
            for relative, path in builds[artifact].files.items()
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
            huggingface.ORD_ARTIFACTS_HF_REPO,
            operations,
            commit_message=commit_message,
            repo_type="dataset",
            revision=revision,
            parent_commit=head,
        )
        commit_url = info.commit_url
    return ArtifactUpload(
        prefixes=prefixes, uploaded=uploaded, deleted=deleted, commit_url=commit_url
    )


def _check_build(directory: pathlib.Path, artifact: str) -> _Build:
    """Checks that one artifact's tree is a single build, and says what it covers.

    Args:
        directory: The artifact's local directory.
        artifact: The artifact it should hold.

    Returns:
        The build.

    Raises:
        ValueError: If the tree is missing or empty, holds a file that is not
            ``artifact`` or carries no RDKit stamp, spans more than one build, files a
            pivot or occurrences artifact under a group its footer does not name, or
            restates a source dataset twice within a group.
    """
    paths = huggingface.parquet_files(directory)
    if not paths:
        raise ValueError(f"no {artifact} artifacts under {directory}")
    group_stamp = _GROUP_STAMPS.get(artifact)
    stamped: collections.Counter[tuple[str, str]] = collections.Counter()
    files: dict[str, pathlib.Path] = {}
    covered: dict[str, set[str]] = collections.defaultdict(set)
    for path in paths:
        stamps = base.load_stamps(path)
        if stamps.artifact != artifact:
            raise ValueError(f"{path} holds {stamps.artifact}, not {artifact}")
        if stamps.rdkit_version is None:
            raise ValueError(f"{path} carries no RDKit stamp; derive it again")
        stamped[stamps.artifact_lineage, stamps.rdkit_version] += 1
        relative = path.relative_to(directory)
        group = ""
        if group_stamp is not None:
            group = relative.parts[0]
            # A Corpus reads the group from the footer and refuses a file filed under
            # another, so the directory alone cannot be trusted to name it.
            held = group_stamp(path)
            if held != group:
                raise ValueError(f"{path} is filed under {group} but holds {held}")
        # A set cannot see a dataset restated twice, which a Corpus refuses.
        if stamps.source_md5 in covered[group]:
            raise ValueError(
                f"{path} restates a source dataset another file under "
                f"{directory / group} already does"
            )
        files[relative.as_posix()] = path
        covered[group].add(stamps.source_md5)
    if len(stamped) > 1:
        raise ValueError(
            f"{directory} spans more than one build, so part of it is stale: "
            f"{dict(stamped.most_common())}"
        )
    [(lineage, rdkit_version)] = stamped
    return _Build(
        lineage=lineage, rdkit_version=rdkit_version, files=files, covered=dict(covered)
    )


def _check_consistent(builds: dict[str, _Build]) -> None:
    """Checks that the artifacts were derived together, each from its parent's build.

    A build checked one artifact at a time can still pair a rederived projection with
    structures derived from the one before: each is internally whole, and each is
    published under a prefix no single library version fetches beside the other.

    Args:
        builds: Every artifact's build.

    Raises:
        ValueError: If the artifacts carry different RDKit versions, or an artifact's
            lineage does not continue its parent's.
    """
    versions = {artifact: build.rdkit_version for artifact, build in builds.items()}
    if len(set(versions.values())) > 1:
        raise ValueError(
            f"the artifacts were derived under different RDKit: {versions}"
        )
    for artifact, build in builds.items():
        parent = base.DERIVED_FROM[artifact]
        if parent is None:
            continue
        inherited = build.lineage.partition(",")[2]
        if inherited != builds[parent].lineage:
            raise ValueError(
                f"{artifact} was derived from {parent} {inherited!r}, but the tree's "
                f"{parent} is {builds[parent].lineage!r}"
            )


def _check_opens(
    root: pathlib.Path, *, require_pivots: bool, require_occurrences: bool
) -> None:
    """Opens a ``Corpus`` on the tree and publishes every pivot, as a reader will.

    The stamps cannot see a structures artifact short of the IDs its projection
    carries, or a file this library and RDKit would not write; opening can. Opened
    without warming, so the substructure library is not built.

    Args:
        root: The local tree.
        require_pivots: Require every level this library derives.
        require_occurrences: Require every indexed path this library derives.

    Raises:
        PairingError: If the corpus refuses the tree.
    """
    with execute.Corpus(
        huggingface.parquet_pattern(
            root / huggingface.ARTIFACT_DIRECTORIES["projection"]
        ),
        huggingface.parquet_pattern(
            root / huggingface.ARTIFACT_DIRECTORIES["structures"]
        ),
        pivots_dir=str(root / huggingface.ARTIFACT_DIRECTORIES["pivot"]),
        require_pivots=require_pivots,
        occurrences_dir=str(root / huggingface.ARTIFACT_DIRECTORIES["occurrences"]),
        require_occurrences=require_occurrences,
        warm=False,
    ) as corpus:
        corpus.check_pivots()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--root", required=True, help="Directory holding one directory per artifact"
    )
    parser.add_argument("--message", required=True, help="Commit message")
    parser.add_argument("--revision", default="main", help="Branch to commit to")
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Check the tree and report the commit without making it",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Checks the tree and publishes it."""
    args = parse_args(argv)
    upload = upload_artifacts(
        args.root,
        commit_message=args.message,
        revision=args.revision,
        dry_run=args.dry_run,
    )
    for artifact, prefix in upload.prefixes.items():
        logger.info("%s -> %s", artifact, prefix)
    verb = "would be" if args.dry_run else "were"
    logger.info(
        "%d files %s uploaded (unchanged ones are dropped from the commit), %d deleted",
        len(upload.uploaded),
        verb,
        len(upload.deleted),
    )
    for path in upload.deleted:
        logger.info("deleting %s", path)
    if upload.commit_url is not None:
        logger.info("committed %s", upload.commit_url)


if __name__ == "__main__":
    main()
