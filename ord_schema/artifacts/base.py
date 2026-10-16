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
* ``ord.artifact_schema_version`` -- the version of the schema the artifact was
  derived to: its columns and what each holds. One for every artifact, because they
  are published and fetched as one build.
* ``ord.source_md5`` -- the content of the source dataset it restates. An artifact
  derived from another artifact passes this through rather than hashing its parent, so
  every artifact names the dataset it reflects however many derivations away it sits,
  and one comparison answers "is this current for that dataset?"
* ``ord.rdkit_version`` -- the RDKit that derived it. Every artifact leans on RDKit --
  canonical SMILES in the projection, fingerprints in the structures artifact -- and
  RDKit releases have changed both canonicalization and pattern fingerprints, either of
  which silently breaks an artifact whose consumers run the new RDKit against files
  built by the old. Optional to *read* so that a file stamped before the key existed
  still loads; absent reads as stale, which rebuilds it.
* ``ord.source_dataset_id`` -- the source's ID, written only when the source records
  one, and so the only key not required to read stamps back.

``is_current`` requires the artifact name, the source content, the artifact schema
version, and the RDKit version to all match, since any of the four changing can
change a value or mean the file answers a different question.

No key names the ord-schema release. A release need not change any artifact, and a
development build reports the last release's number, so a release stamp says nothing
about whether a file is current -- ``ord.artifact_schema_version`` is what says the
columns moved -- while making a rebuild of an unchanged artifact write different
bytes, which an upload then carries as a change.
"""

import dataclasses
import glob
import os
import pathlib
import re
from collections.abc import Callable, Collection, Mapping
from typing import Protocol

import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import rdBase

from ord_schema import parquet
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# A bracket holding a single wildcard character, the way glob.escape writes one, matches
# that character alone.
_ESCAPED_WILDCARD = re.compile(r"\[([*?\[])\]")

# Bumped when a reader built for one schema would misread an artifact derived to
# another: a column renamed, removed, or given a different meaning. A change a reader of
# either reads correctly leaves it alone -- a column added, which ``missing_columns``
# sees, or values improved by a fixed normalization, which is derived again with --force
# and published in place. base_test pins what each artifact holds, so a change to it
# fails there until someone decides which kind it is. One schema for every artifact,
# because they are published and fetched as one build.
ARTIFACT_SCHEMA_VERSION = "1"


_META_ARTIFACT = "ord.artifact"
_META_ARTIFACT_SCHEMA_VERSION = "ord.artifact_schema_version"
_META_RDKIT_VERSION = "ord.rdkit_version"
_META_SOURCE_DATASET_ID = "ord.source_dataset_id"
_META_SOURCE_MD5 = "ord.source_md5"

_REQUIRED = (
    _META_ARTIFACT,
    _META_ARTIFACT_SCHEMA_VERSION,
    _META_SOURCE_MD5,
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
    artifact_schema_version: str
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
        Stamps carrying the artifact schema version and the RDKit version in use.
    """
    return Stamps(
        artifact=artifact,
        source_dataset_id=source_dataset_id,
        source_md5=source_md5,
        artifact_schema_version=ARTIFACT_SCHEMA_VERSION,
        rdkit_version=rdBase.rdkitVersion,
    )


def to_metadata(value: Stamps) -> dict[str, str]:
    """Returns ``value`` as Parquet footer key-value metadata."""
    result = {
        _META_ARTIFACT: value.artifact,
        _META_ARTIFACT_SCHEMA_VERSION: value.artifact_schema_version,
        _META_SOURCE_MD5: value.source_md5,
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
        return stamps_from_metadata(parquet_file.schema_arrow.metadata or {})


def stamps_from_metadata(raw: Mapping[bytes, bytes]) -> Stamps:
    """Reads stamps from Parquet footer key-value metadata, for a file read elsewhere.

    Args:
        raw: The footer's key-value metadata.

    Returns:
        The stamps recorded when the artifact was written.

    Raises:
        ValueError: If ``raw`` is missing a required stamp, e.g. because it is a source
            dataset's footer rather than a derived artifact's.
    """
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
        artifact_schema_version=values[_META_ARTIFACT_SCHEMA_VERSION],
        rdkit_version=rdkit_version.decode() if rdkit_version is not None else None,
    )


def _is_nested(dtype: pa.DataType) -> bool:
    """Returns whether ``dtype`` holds other fields rather than values of its own."""
    return (
        pa.types.is_struct(dtype)
        or pa.types.is_list(dtype)
        or pa.types.is_large_list(dtype)
        or pa.types.is_map(dtype)
    )


def _declared_paths(
    schema: pa.Schema | pa.DataType, prefix: str = ""
) -> list[tuple[str, pa.DataType | None]]:
    """Returns every field a schema declares, dotted and typed, in declaration order.

    A list's or map's own name is not a field a reader binds; what it holds is, so the
    descent passes through them without adding a segment of its own.

    Args:
        schema: A schema, or a nested type within one.
        prefix: What the walk has already named, for the recursion.

    Returns:
        Pairs of dotted path and leaf type, outermost first. The type is None for a
        field that holds other fields, whose own type is decided by the leaves beneath
        it and would otherwise report a change twice.
    """
    if isinstance(schema, pa.Schema):
        fields: list[pa.Field] = list(schema)
    elif pa.types.is_struct(schema):
        fields = [schema.field(index) for index in range(schema.num_fields)]
    elif pa.types.is_list(schema) or pa.types.is_large_list(schema):
        return _declared_paths(schema.value_type, prefix)
    elif pa.types.is_map(schema):
        return _declared_paths(schema.item_type, prefix)
    else:
        return []
    paths: list[tuple[str, pa.DataType | None]] = []
    for field in fields:
        name = f"{prefix}.{field.name}" if prefix else field.name
        paths.append((name, None if _is_nested(field.type) else field.type))
        paths.extend(_declared_paths(field.type, name))
    return paths


def missing_columns(path: str | os.PathLike[str], schema: pa.Schema) -> list[str]:
    """Returns the columns ``schema`` declares that ``path`` does not carry.

    Stamps say nothing about a file's columns, so this is what sees a column added to
    an artifact's definition: a file written before it stamps exactly as a file written
    after it, and the difference shows up only in the schema.

    Compared down through structs, lists and maps, not across the top level alone. The
    projection is one struct per message and nearly every field it has is nested, so a
    top-level comparison would call a file current that is missing everything added
    since it was written -- and a reader binding the new column would fail on a file
    nothing marked stale.

    A leaf whose type changed counts as missing too, because a path comparison alone
    cannot see it: widening a column or moving a timestamp from seconds to microseconds
    leaves every name where it was, so a file holding the old type would read as current
    and keep serving values the current definition does not describe.

    Args:
        path: Path to check. Need not exist.
        schema: The columns this library's version of the artifact carries.

    Returns:
        The declared paths the file lacks, dotted from the top level and in declaration
        order. A path present with the wrong type carries the type it should hold, so
        the message says what changed rather than naming a column the file visibly has.
        A file that cannot be read at all lacks all of them.
    """
    declared = _declared_paths(schema)
    try:
        with pq.ParquetFile(path) as parquet_file:
            found = dict(_declared_paths(parquet_file.schema_arrow))
    except (OSError, ValueError, pa.ArrowInvalid):
        return [name for name, _ in declared]
    missing = []
    for name, leaf in declared:
        if name not in found:
            missing.append(name)
        elif leaf is not None and found[name] != leaf:
            missing.append(f"{name}: {leaf}")
    return missing


def is_current(
    path: str | os.PathLike[str],
    artifact: str,
    source_md5: str,
    schema: pa.Schema | None = None,
) -> bool:
    """Returns whether ``path`` holds ``artifact`` derived from ``source_md5`` by us.

    All of the source content, the artifact schema version, and the RDKit version must
    match;
    any of them changing can change a value. A missing, unreadable, or wrong-artifact
    file is not current.

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
        Whether the artifact name, the artifact schema version, and the RDKit version
        all match. A file with no RDKit stamp is not current.
    """
    return (
        value.artifact == artifact
        and value.artifact_schema_version == ARTIFACT_SCHEMA_VERSION
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

    A bracket holding only ``*``, ``?``, or ``[`` matches that character alone, which
    is how ``glob.escape`` writes a directory name that contains one, so such a bracket
    is read as the character: ``glob.escape("data[1]") + "/*/x.parquet"`` has the root
    ``data[1]``.

    Args:
        pattern: Glob pattern, with ``*``, ``?``, or ``[`` marking the first wildcard.

    Returns:
        The wildcard-free leading directories, unescaped, or ``.`` when the pattern
        names no directory or begins with a wildcard.
    """
    parts = pathlib.PurePath(pattern).parts
    fixed = []
    for part in parts:
        if any(char in _ESCAPED_WILDCARD.sub("", part) for char in "*?["):
            break
        fixed.append(_ESCAPED_WILDCARD.sub(r"\1", part))
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
        root: The directory to mirror from, for a caller that has one; otherwise
            ``glob_root(pattern)``.

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
    # own inputs, and one holding a file from outside this library. One message for both
    # would send the reader looking for the input that landed there, and for the second
    # kind there is none -- the file was already sitting at the destination.
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


def left_behind(
    directory: str | os.PathLike[str],
    artifact: str,
    schemas: Mapping[str, pa.Schema],
    artifact_paths: Callable[[str | os.PathLike[str], str], list[pathlib.Path]],
    derived: Collection[str],
) -> dict[str, str]:
    """Returns the keys already under ``directory`` that a run left short of a column.

    An artifact split across a tree -- one subdirectory per pivot level, per indexed
    path -- is derived a subset at a time, which is the cheaper run and often the only
    one wanted. The tree that run writes into holds whatever earlier runs put there, and
    those files stamp current forever, so a schema that grew a column leaves every key
    nobody named short of it. A corpus reads the whole tree rather than the part the run
    touched.

    Columns only, matching what a reader refuses whatever its policy: an artifact short
    a column fails every corpus that opens it, while one whose stamps are merely old is
    refused where ``require_current`` says so and read otherwise. A build that stopped
    on the second would decide that question on the reader's behalf, and leave no way to
    finish a subset run against a tree an older library wrote.

    Args:
        directory: The directory the keys sit under.
        artifact: Artifact name the files are expected to hold; anything else is left
            alone, since a tree is not a build's to police.
        schemas: Every key this artifact covers, mapped to the columns it carries.
        artifact_paths: Lists one key's artifacts exactly as a reader finds them, taking
            ``directory`` and the key.
        derived: The keys the run covered, which are current by construction.

    Returns:
        One entry per uncovered key whose artifacts no longer match its schema, naming
        the first artifact that does not and what it lacks. Keys with no artifacts are
        not included: there is nothing there to be out of date.
    """
    behind: dict[str, str] = {}
    for key in sorted(set(schemas) - set(derived)):
        for path in artifact_paths(directory, key):
            try:
                stamps = load_stamps(path)
            except (OSError, ValueError):
                continue
            if stamps.artifact != artifact:
                continue
            missing = missing_columns(path, schemas[key])
            if missing:
                behind[key] = f"{path} is without {missing}"
                break
    return behind
