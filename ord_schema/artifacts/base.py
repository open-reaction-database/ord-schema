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

"""Shared mechanics for derived artifacts: footer stamps and the derive-a-tree driver.

A derived artifact restates a source dataset in a queryable shape. It adds no
information, so it is wrong if it disagrees with its source and stale if its source has
moved on. These stamps make both conditions detectable without reading column data.

More than one artifact is expected -- a projection carrying every field, structures
pulled out of it for search, and pivoted indexes over the paths that turn out to be hot
-- so the parts that do not vary between them live here: the stamps, the output-path
mapping, and the driver that walks a glob, refuses to write over what it reads, and
skips what is current. An artifact module supplies only its schema, its projection, and
its name.

An artifact may derive from another rather than from a source dataset, which is how the
structures artifact reaches the projection's columns without recomputing them. Only the
reading differs: the stamps describe the originating dataset either way, so a chain is
invisible to a consumer checking whether a file is current.

Every artifact carries these keys, in the ``ord.`` namespace that ``ord_schema.parquet``
already uses for source datasets:

* ``ord.artifact`` -- which artifact this is, so a reader can tell one from another
  without inspecting the schema.
* ``ord.artifact_lineage`` -- this artifact's version and every ancestor's, as one
  string. Per artifact, because the derivations cost wildly different amounts and a
  shared version charged the most expensive rebuild for a change to the cheapest. The
  whole chain rather than the artifact's own version, because what a version has to
  invalidate is everything built *from* it: occurrences derived from a pivot derived
  from a projection go stale when any of the three is redefined, and a stamp naming only
  one of them cannot say so.
* ``ord.source_md5`` -- the content of the source dataset it restates. An artifact
  derived from another artifact passes this through rather than hashing its parent, so
  every artifact names the dataset it reflects however many derivations away it sits,
  and one comparison answers "is this current for that dataset?"
* ``ord.ord_schema_version`` -- what derived it.
* ``ord.rdkit_version`` -- the RDKit that derived it. Every artifact leans on RDKit --
  canonical SMILES in the projection, fingerprints in the structures artifact -- and
  RDKit releases have changed both canonicalization and pattern fingerprints, either of
  which silently breaks an artifact whose consumers run the new RDKit against files
  built by the old. Optional to *read* so that a file stamped before the key existed
  still loads; absent reads as stale, which rebuilds it.
* ``ord.source_dataset_id`` -- the source's ID, written only when the source records
  one, and so the only key not required to read stamps back.

``is_current`` requires the artifact name, the source content, the library version, the
artifact lineage, and the RDKit version to all match, since any of the five changing can
change a value or mean the file answers a different question.
"""

import dataclasses
import glob
import os
import pathlib
from importlib import metadata
from typing import Protocol

import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import rdBase

from ord_schema import parquet
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Bumped when an artifact's definition changes -- a column's meaning, how it is
# populated, or a shared helper that feeds it. Per artifact, because a re-derive is not
# one price: measured over the source that is 96% of ORD, a projection is 34.6 minutes
# and its structures 29.6, while the occurrences at every indexed path are 1.9 seconds.
# Sharing one version across them would charge the 81-minute rebuild for a change to the
# 1.9-second artifact.
ARTIFACT_VERSIONS = {
    "projection": "1",
    "structures": "1",
    "pivot": "1",
    "occurrences": "1",
}

# What each artifact is derived from, or None where it is derived from a source dataset.
# Here rather than in the artifact modules because a version has to invalidate what was
# built *from* it, which is a fact about the chain rather than about either end: a
# projection bump has to reach the occurrences three derivations down, and neither
# module knows the other exists.
DERIVED_FROM: dict[str, str | None] = {
    "projection": None,
    "structures": "projection",
    "pivot": "projection",
    "occurrences": "pivot",
}


def lineage(artifact: str) -> str:
    """Returns an artifact's version and every ancestor's, as one comparable string.

    What is stamped, so that a bump anywhere up the chain makes a file stale. Stamping
    only the artifact's own version would leave the occurrences current when the pivot
    they came from was redefined; stamping only the parent's would leave them current
    when the *projection* under that pivot was, which is the case a chain of three makes
    reachable.

    Args:
        artifact: Artifact name, a key of ``ARTIFACT_VERSIONS``.

    Returns:
        The chain from this artifact to its root, e.g. ``occurrences=1,pivot=1,
        projection=1``.

    Raises:
        KeyError: If ``artifact`` is not one this library derives, which is a typo
            rather than a condition to handle.
        ValueError: If ``DERIVED_FROM`` describes a cycle. Bounded rather than walked to
            exhaustion because this runs on every artifact written and every currency
            check made, and a typo in a four-line literal should not hang all of them.
    """
    parts = []
    name: str | None = artifact
    for _ in range(len(ARTIFACT_VERSIONS)):
        if name is None:
            return ",".join(parts)
        parts.append(f"{name}={ARTIFACT_VERSIONS[name]}")
        name = DERIVED_FROM[name]
    raise ValueError(
        f"the chain from {artifact} visits {len(parts)} artifacts without reaching a "
        f"source dataset, so DERIVED_FROM has a cycle: {' -> '.join(parts)}"
    )


_META_ARTIFACT = "ord.artifact"
_META_ARTIFACT_LINEAGE = "ord.artifact_lineage"
_META_RDKIT_VERSION = "ord.rdkit_version"
_META_SOURCE_DATASET_ID = "ord.source_dataset_id"
_META_SOURCE_MD5 = "ord.source_md5"
_META_ORD_SCHEMA_VERSION = "ord.ord_schema_version"

_REQUIRED = (
    _META_ARTIFACT,
    _META_ARTIFACT_LINEAGE,
    _META_SOURCE_MD5,
    _META_ORD_SCHEMA_VERSION,
)


@dataclasses.dataclass(frozen=True)
class Stamps:
    """Footer stamps recording what an artifact was derived from, and by what.

    ``rdkit_version`` is None only when read from a file stamped before the key
    existed; such a file reads as stale.
    """

    artifact: str
    source_dataset_id: str | None
    source_md5: str
    ord_schema_version: str
    artifact_lineage: str
    rdkit_version: str | None = None


class ArtifactWriter(Protocol):
    """Writes one derived artifact from the Parquet file it derives from.

    The parent is passed as a path rather than an open reader, because parents are not
    all the same kind of file: a projection reads a source dataset, and the structures
    artifact reads a projection. Each writer opens what it knows how to read.

    The first two arguments are positional-only, so an implementation is free to name
    them for its own artifact.
    """

    def __call__(
        self,
        source: pathlib.Path,
        output: pathlib.Path,
        /,
        *,
        source_md5: str,
        source_dataset_id: str | None,
    ) -> int:
        """Returns the number of rows written to ``output``."""


def current_stamps(
    artifact: str, source_dataset_id: str | None, source_md5: str
) -> Stamps:
    """Returns the stamps to write for an artifact derived by the current library.

    Args:
        artifact: Artifact name, e.g. ``"structures"`` or ``"projection"``.
        source_dataset_id: Source dataset ID, or None if the source records none.
        source_md5: Hash of the source, from ``parquet.DatasetView.md5``.

    Returns:
        Stamps carrying the current library version and the artifact's lineage.
    """
    return Stamps(
        artifact=artifact,
        source_dataset_id=source_dataset_id,
        source_md5=source_md5,
        ord_schema_version=metadata.version("ord-schema"),
        artifact_lineage=lineage(artifact),
        rdkit_version=rdBase.rdkitVersion,
    )


def to_metadata(value: Stamps) -> dict[str, str]:
    """Returns ``value`` as Parquet footer key-value metadata."""
    result = {
        _META_ARTIFACT: value.artifact,
        _META_ARTIFACT_LINEAGE: value.artifact_lineage,
        _META_SOURCE_MD5: value.source_md5,
        _META_ORD_SCHEMA_VERSION: value.ord_schema_version,
    }
    if value.source_dataset_id is not None:
        result[_META_SOURCE_DATASET_ID] = value.source_dataset_id
    if value.rdkit_version is not None:
        result[_META_RDKIT_VERSION] = value.rdkit_version
    return result


def load_stamps(path: str | os.PathLike[str]) -> Stamps:
    """Reads an artifact's footer stamps.

    Args:
        path: Path to a derived artifact Parquet file.

    Returns:
        The stamps recorded when the artifact was written.

    Raises:
        ValueError: If ``path`` is missing a required stamp, e.g. because it is a source
            dataset rather than a derived artifact.
    """
    with pq.ParquetFile(path) as parquet_file:
        raw = parquet_file.schema_arrow.metadata or {}
    values = {}
    for key in _REQUIRED:
        value = raw.get(key.encode())
        if value is None:
            missing = [key for key in _REQUIRED if raw.get(key.encode()) is None]
            raise ValueError(
                f"not a derived artifact; missing footer keys: {sorted(missing)}"
            )
        values[key] = value.decode()
    source_dataset_id = raw.get(_META_SOURCE_DATASET_ID.encode())
    rdkit_version = raw.get(_META_RDKIT_VERSION.encode())
    return Stamps(
        artifact=values[_META_ARTIFACT],
        source_dataset_id=(
            source_dataset_id.decode() if source_dataset_id is not None else None
        ),
        source_md5=values[_META_SOURCE_MD5],
        ord_schema_version=values[_META_ORD_SCHEMA_VERSION],
        artifact_lineage=values[_META_ARTIFACT_LINEAGE],
        rdkit_version=rdkit_version.decode() if rdkit_version is not None else None,
    )


def missing_columns(path: str | os.PathLike[str], schema: pa.Schema) -> list[str]:
    """Returns the columns ``schema`` declares that ``path`` does not carry.

    Stamps say nothing about a file's columns, so this is what sees a column added to
    an artifact's definition: a file written before it stamps exactly as a file written
    after it, and the difference shows up only in the schema.

    Args:
        path: Path to check. Need not exist.
        schema: The columns this library's version of the artifact carries.

    Returns:
        The declared names the file lacks, in declaration order. A file that cannot be
        read at all lacks all of them.
    """
    try:
        with pq.ParquetFile(path) as parquet_file:
            names = set(parquet_file.schema_arrow.names)
    except (OSError, ValueError, pa.ArrowInvalid):
        return list(schema.names)
    return [name for name in schema.names if name not in names]


def is_current(
    path: str | os.PathLike[str],
    artifact: str,
    source_md5: str,
    schema: pa.Schema | None = None,
) -> bool:
    """Returns whether ``path`` holds ``artifact`` derived from ``source_md5`` by us.

    All of the source content, the library version, the artifact lineage, and the
    RDKit version must match; any of them changing can change a value. A missing,
    unreadable, or wrong-artifact file is not current.

    Args:
        path: Path to check. Need not exist.
        artifact: Artifact name the file is expected to hold.
        source_md5: Hash the file is expected to restate.
        schema: The columns the artifact carries, when the caller has them. Given, a
            file lacking one is not current, which is what makes a column added
            without a version bump derive again instead of being skipped.

    Returns:
        Whether the file is one this library would write today.
    """
    try:
        value = load_stamps(path)
    except (OSError, ValueError):
        # Broad on purpose: every way of failing to read stamps -- absent, truncated,
        # not Parquet at all -- answers "derive it again", which is safe.
        return False
    if schema is not None and missing_columns(path, schema):
        return False
    return value.source_md5 == source_md5 and stamps_are_current(value, artifact)


def stamps_are_current(value: Stamps, artifact: str) -> bool:
    """Returns whether ``value`` records ``artifact`` as this library defines it.

    The source content is deliberately not compared: a caller holding the stamps of a
    parent asks whether the file was written by the current definition, and the hash it
    carries is the answer to a different question.

    Args:
        value: Stamps read from an artifact.
        artifact: Artifact name the file is expected to hold.

    Returns:
        Whether the artifact name, the library version, the artifact lineage, and the
        RDKit version all match. A file with no RDKit stamp is not current -- and a
        lineage mismatch covers a redefinition of this artifact or of anything it was
        derived from.
    """
    return (
        value.artifact == artifact
        and value.ord_schema_version == metadata.version("ord-schema")
        and value.artifact_lineage == lineage(artifact)
        and value.rdkit_version == rdBase.rdkitVersion
    )


def is_artifact(path: str | os.PathLike[str]) -> bool:
    """Returns whether ``path`` carries derived-artifact stamps.

    Args:
        path: Path to a Parquet file.

    Returns:
        True if ``path`` records artifact stamps, False if it is readable Parquet that
        does not.

    Raises:
        ValueError: If ``path`` cannot be read as Parquet at all. This predicate decides
            whether to treat a file as a source, so answering False for a truncated one
            would feed it to a reader that fails later without naming it.
    """
    return read_stamps(path) is not None


def read_stamps(path: str | os.PathLike[str]) -> Stamps | None:
    """Returns ``path``'s artifact stamps, or None if it carries none.

    One footer read where a predicate followed by ``load_stamps`` would be two, which
    matters to ``derive_tree``: it needs the answer *and* the stamps behind it for every
    match of a corpus-wide glob.

    Args:
        path: Path to a Parquet file.

    Returns:
        The stamps, or None if ``path`` is readable Parquet that records none.

    Raises:
        ValueError: If ``path`` cannot be read as Parquet at all. This decides whether a
            glob match is a candidate parent, so answering None for a truncated one
            would either feed it to a reader that fails later without naming it, or drop
            it from the run without a word.
    """
    try:
        return load_stamps(path)
    except pa.ArrowInvalid as error:
        raise ValueError(f"{path} is not readable as Parquet") from error
    except (OSError, ValueError):
        return None


def glob_root(pattern: str) -> pathlib.PurePath:
    """Returns the leading directories of a glob pattern that contain no wildcards.

    Output paths are built relative to this, so ``data/*/x.parquet`` under
    ``output_dir=structures`` lands at ``structures/<subdir>/x.parquet``. A pattern
    naming a single file has no wildcard to stop at, so its own last component is the
    file rather than a directory and the root is the directory holding it.

    Args:
        pattern: Glob pattern, with ``*``, ``?``, or ``[`` marking the first wildcard.

    Returns:
        The wildcard-free leading directories, or ``.`` when the pattern names no
        directory or begins with a wildcard.
    """
    parts = pathlib.PurePath(pattern).parts
    fixed = []
    for part in parts:
        if any(char in part for char in "*?["):
            break
        fixed.append(part)
    else:
        fixed = fixed[:-1]
    return pathlib.PurePath(*fixed)


def output_path(
    source: str,
    pattern: str,
    output_dir: str,
    *,
    root: str | os.PathLike[str] | None = None,
) -> pathlib.Path:
    """Maps a source path to its artifact path under ``output_dir``.

    Args:
        source: Path the pattern matched.
        pattern: The glob it matched, whose wildcard-free prefix the output mirrors from
            when ``root`` is not given.
        output_dir: Directory to write beneath.
        root: The directory to mirror from, for a caller that has one. A caller building
            a pattern around a directory name has to escape it, and an escaped ``[`` is
            a character class as far as ``glob_root`` can tell -- so the prefix it finds
            stops short and the mirror carries the source layout twice.

    Returns:
        Where ``source``'s artifact goes.
    """
    relative = pathlib.PurePath(source).relative_to(
        glob_root(pattern) if root is None else pathlib.PurePath(root)
    )
    return pathlib.Path(output_dir) / relative


def _is_parent(match: str, parent_artifact: str | None, stamps: Stamps | None) -> bool:
    """Returns whether ``match`` is the kind of file this derivation reads.

    Args:
        match: Path a glob matched, named in the log line when it is ignored.
        parent_artifact: Artifact name the derivation reads, or None to read source
            datasets.
        stamps: ``match``'s stamps, or None if it carries none.

    Returns:
        True if ``match`` should be derived from. Anything else readable is ignored
        rather than refused, so a pattern reaching the output tree stays re-runnable.
    """
    if stamps is None:
        if parent_artifact is None:
            return True
        logger.info(
            "%s is a source dataset, not a %s; ignoring", match, parent_artifact
        )
        return False
    if parent_artifact is None:
        logger.info("%s is a derived artifact, not a source; ignoring", match)
        return False
    if stamps.artifact == parent_artifact:
        return True
    logger.info(
        "%s is a %s, not a %s; ignoring", match, stamps.artifact, parent_artifact
    )
    return False


def _is_irreplaceable(destination: pathlib.Path) -> bool:
    """Returns whether writing to ``destination`` would destroy something unrecoverable.

    A derived artifact is ours to rewrite; anything else at that path is not. Comparing
    destinations against the run's own parents is not enough once a derivation reads one
    tree and writes another, because the source datasets are then in neither set -- and
    a mistyped ``output_dir`` aimed at them would replace a corpus that cannot be
    regenerated with artifacts that can.

    The question is who wrote the file, not whether it is current, so ``ord.artifact``
    alone decides it: nothing but this library writes that key, and a file carrying it
    is one of ours whatever else its footer says. Requiring the whole stamp set would
    answer a different question -- whether the file matches *today's* format -- and so
    would read every artifact written before a key was added or renamed as a stranger's,
    refusing the very re-derive that would bring it up to date.

    Args:
        destination: Path an artifact would be written to.

    Returns:
        True if something is already there and this library did not write it.
    """
    if not destination.exists():
        return False
    try:
        with pq.ParquetFile(destination) as parquet_file:
            raw = parquet_file.schema_arrow.metadata or {}
    except (OSError, ValueError, pa.ArrowInvalid):
        # Not readable as Parquet at all, so certainly not an artifact this run wrote.
        return True
    return _META_ARTIFACT.encode() not in raw


def _parent_provenance(
    parent: str, parent_artifact: str | None, stamps: Stamps | None
) -> tuple[str, str | None]:
    """Returns the provenance an artifact reading ``parent`` should stamp.

    A derived parent passes its own stamps through rather than being hashed itself, so
    every artifact records the *source dataset* content it reflects however many
    derivations away it sits, and a consumer checks an artifact against the dataset in
    one hop.

    Passing the hash through carries the source content across the hop but not the
    version stamps, which describe whoever wrote each file. A stale parent is therefore
    refused rather than read: an artifact derived from one would stamp itself with the
    *current* versions while holding what an older definition produced, and since the
    dataset hash it inherits does not change when the parent is rebuilt, nothing would
    ever mark it stale again. Refusing keeps every term ``is_current`` compares true of
    a chained artifact as well as a directly derived one.

    Args:
        parent: Path to the file being derived from.
        parent_artifact: Artifact name ``parent`` holds, or None if it is a source
            dataset.
        stamps: ``parent``'s stamps, already read; None when it is a source dataset.

    Returns:
        The hash of the originating source dataset and its ID, if it records one.

    Raises:
        ValueError: If ``parent`` is not readable as the expected kind of file, or is a
            derived artifact that is itself out of date.
    """
    if parent_artifact is not None:
        # Selection already established that stamps holds this artifact.
        assert stamps is not None  # Type hint.
        if not stamps_are_current(stamps, parent_artifact):
            raise ValueError(
                f"{parent} is a stale {parent_artifact}; derive it again first, or "
                "what is written here would claim a provenance it does not have"
            )
        return stamps.source_md5, stamps.source_dataset_id
    try:
        view = parquet.DatasetView(parent)
    except ValueError as error:
        raise ValueError(f"{parent} is not a source dataset: {error}") from error
    return view.md5(), view.dataset_id or None


def derive_tree(
    input_pattern: str,
    output_dir: str,
    *,
    artifact: str,
    write: ArtifactWriter,
    schema: pa.Schema | None = None,
    force: bool = False,
    parent_artifact: str | None = None,
    root: str | os.PathLike[str] | None = None,
) -> tuple[int, int, int]:
    """Derives one artifact per file matching ``input_pattern``.

    Outputs mirror the inputs' directory layout beneath ``output_dir``. Matches that are
    not the kind of parent this derivation reads are ignored, so a recursive pattern
    reaching the output tree stays re-runnable. Artifacts whose footer already records
    the current source content and versions are skipped, so re-running is cheap.

    Args:
        input_pattern: Glob matching the files to derive from.
        output_dir: Directory to write artifacts beneath.
        artifact: Artifact name, used for the staleness check and the footer stamp.
        write: Writer taking ``(parent, output, source_md5=..., source_dataset_id=...)``
            and returning rows.
        schema: The columns ``write`` produces, when the caller has them. Given, an
            artifact lacking one is derived again rather than skipped.
        force: Rewrite artifacts that are already current.
        parent_artifact: Artifact the parents hold, for a derivation that reads another
            artifact rather than a source dataset. None means the parents are source
            datasets.
        root: The directory outputs mirror from, for a caller that has one rather than a
            pattern. Defaults to ``input_pattern``'s wildcard-free leading directories,
            which is what a caller who wrote the pattern means by it.

    Returns:
        ``(written, skipped, ignored)`` counts: artifacts derived, artifacts already
        current, and matches that were not the expected kind of parent.

    Raises:
        ValueError: If the input pattern matches nothing, if any match cannot be read as
            Parquet, if any destination would land on a parent, or if any destination
            already holds a file this library did not write.
    """
    # Deduplicated: a pattern with more than one ** reaches the same file by more than
    # one route, and a source derived twice is written twice -- the second pass reading
    # what the first just wrote, and counted again in what this returns.
    matches = sorted(set(glob.glob(input_pattern, recursive=True)))
    if not matches:
        raise ValueError(f"no datasets matched: {input_pattern}")
    sources = []
    parent_stamps: dict[str, Stamps | None] = {}
    ignored = 0
    for match in matches:
        stamps = read_stamps(match)
        if _is_parent(match, parent_artifact, stamps):
            sources.append(match)
            parent_stamps[match] = stamps
        else:
            ignored += 1
    logger.info(
        "Found %d matching %s files", len(sources), parent_artifact or "dataset"
    )
    # Every destination against every source, not against its own: an output_dir nested
    # in the source tree maps one dataset onto a *different* one, which a per-source
    # check passes, destroying a source the run had not reached yet.
    destinations = {
        source: output_path(source, input_pattern, output_dir, root=root)
        for source in sources
    }
    # Two different mistakes, reported apart: a destination that is one of this run's
    # own inputs, and one holding a file from outside this library. Naming both "inputs"
    # sends the reader looking for the input that landed there, and for the second kind
    # there is none -- the file was already sitting at the destination.
    resolved_sources = {pathlib.Path(source).resolve() for source in sources}
    over_inputs, over_strangers = [], []
    for destination in destinations.values():
        if destination.resolve() in resolved_sources:
            over_inputs.append(str(destination))
        elif _is_irreplaceable(destination):
            over_strangers.append(str(destination))
    refused = []
    if over_inputs:
        refused.append(f"its own inputs: {sorted(over_inputs)}")
    if over_strangers:
        refused.append(f"files this library did not write: {sorted(over_strangers)}")
    if refused:
        raise ValueError(
            f"output_dir {output_dir!r} would write over {' and '.join(refused)}"
        )
    written = skipped = 0
    for source in sources:
        destination = destinations[source]
        source_md5, source_dataset_id = _parent_provenance(
            source, parent_artifact, parent_stamps[source]
        )
        if not force and is_current(destination, artifact, source_md5, schema):
            logger.info("%s is current; skipping", destination)
            skipped += 1
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        rows = write(
            pathlib.Path(source),
            destination,
            source_md5=source_md5,
            source_dataset_id=source_dataset_id,
        )
        logger.info("Wrote %d rows to %s", rows, destination)
        written += 1
    logger.info(
        "Derived %d %s artifacts (%d already current, %d ignored)",
        written,
        artifact,
        skipped,
        ignored,
    )
    return written, skipped, ignored
