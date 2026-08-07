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

"""Shared mechanics for derived artifacts: footer stamps, and the derive-a-tree driver.

A derived artifact restates a source dataset in a queryable shape. It adds no
information, so it is wrong if it disagrees with its source and stale if its source has
moved on. These stamps make both conditions detectable without reading column data.

More than one artifact is expected -- a projection carrying every field, and pivoted
indexes over the paths that turn out to be hot -- so the parts that do not vary between
them live here: the stamps, the output-path mapping, and the driver that walks a glob,
refuses to write over its own sources, ignores matches that are already artifacts, and
skips what is current. An artifact module supplies only its schema, its projection, and
its name.

Every artifact carries these keys, in the ``ord.`` namespace that ``ord_schema.parquet``
already uses for source datasets:

* ``ord.artifact`` -- which artifact this is, so a reader can tell a view from a
  projection without inspecting the schema.
* ``ord.artifact_version`` -- one version across *all* artifacts, deliberately. A
  derivation change usually touches shared helpers, and a reader comparing artifacts to
  each other needs to know they were built by the same definition. Per-artifact versions
  would let a view and a projection of the same dataset disagree while both looked
  current.
* ``ord.source_md5`` -- what it was derived from.
* ``ord.ord_schema_version`` -- what derived it.
* ``ord.source_dataset_id`` -- the source's ID, written only when the source records
  one, and so the only key not required to read stamps back.

``is_current`` requires the artifact name, the source content, the library version, and
the artifact version to all match, since any of the four changing can change a value or
mean the file answers a different question.
"""

import dataclasses
import glob
import os
import pathlib
from importlib import metadata
from typing import Protocol

import pyarrow as pa
import pyarrow.parquet as pq

from ord_schema import parquet
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Bumped when any artifact's definition changes -- a column's meaning, how it is
# populated, or a shared helper that feeds one. Existing artifacts then read as stale.
ARTIFACT_VERSION = "1"

_META_ARTIFACT = "ord.artifact"
_META_ARTIFACT_VERSION = "ord.artifact_version"
_META_SOURCE_DATASET_ID = "ord.source_dataset_id"
_META_SOURCE_MD5 = "ord.source_md5"
_META_ORD_SCHEMA_VERSION = "ord.ord_schema_version"

_REQUIRED = (
    _META_ARTIFACT,
    _META_ARTIFACT_VERSION,
    _META_SOURCE_MD5,
    _META_ORD_SCHEMA_VERSION,
)


@dataclasses.dataclass(frozen=True)
class Stamps:
    """Footer stamps recording what an artifact was derived from, and by what."""

    artifact: str
    source_dataset_id: str | None
    source_md5: str
    ord_schema_version: str
    artifact_version: str


class ArtifactWriter(Protocol):
    """Writes one derived artifact from an open view of its source.

    The first two arguments are positional-only, so an implementation is free to name
    them for its own artifact.
    """

    def __call__(
        self,
        source: parquet.DatasetView,
        output: pathlib.Path,
        /,
        *,
        source_md5: str,
    ) -> int:
        """Returns the number of rows written to ``output``."""


def current_stamps(
    artifact: str, source_dataset_id: str | None, source_md5: str
) -> Stamps:
    """Returns the stamps to write for an artifact derived by the current library.

    Args:
        artifact: Artifact name, e.g. ``"view"`` or ``"projection"``.
        source_dataset_id: Source dataset ID, or None if the source records none.
        source_md5: Hash of the source, from ``parquet.DatasetView.md5``.

    Returns:
        Stamps carrying the current library and artifact versions.
    """
    return Stamps(
        artifact=artifact,
        source_dataset_id=source_dataset_id,
        source_md5=source_md5,
        ord_schema_version=metadata.version("ord-schema"),
        artifact_version=ARTIFACT_VERSION,
    )


def to_metadata(value: Stamps) -> dict[str, str]:
    """Returns ``value`` as Parquet footer key-value metadata."""
    result = {
        _META_ARTIFACT: value.artifact,
        _META_ARTIFACT_VERSION: value.artifact_version,
        _META_SOURCE_MD5: value.source_md5,
        _META_ORD_SCHEMA_VERSION: value.ord_schema_version,
    }
    if value.source_dataset_id is not None:
        result[_META_SOURCE_DATASET_ID] = value.source_dataset_id
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
    return Stamps(
        artifact=values[_META_ARTIFACT],
        source_dataset_id=(
            source_dataset_id.decode() if source_dataset_id is not None else None
        ),
        source_md5=values[_META_SOURCE_MD5],
        ord_schema_version=values[_META_ORD_SCHEMA_VERSION],
        artifact_version=values[_META_ARTIFACT_VERSION],
    )


def is_current(path: str | os.PathLike[str], artifact: str, source_md5: str) -> bool:
    """Returns whether ``path`` holds ``artifact`` derived from ``source_md5`` by us.

    All of the source content, the library version, and the artifact version must match;
    any of them changing can change a value. A missing, unreadable, or wrong-artifact
    file is not current.
    """
    try:
        value = load_stamps(path)
    except (OSError, ValueError):
        # Broad on purpose: every way of failing to read stamps -- absent, truncated,
        # not Parquet at all -- answers "derive it again", which is safe.
        return False
    return (
        value.artifact == artifact
        and value.source_md5 == source_md5
        and value.ord_schema_version == metadata.version("ord-schema")
        and value.artifact_version == ARTIFACT_VERSION
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
    try:
        load_stamps(path)
    except pa.ArrowInvalid as error:
        raise ValueError(f"{path} is not readable as Parquet") from error
    except (OSError, ValueError):
        return False
    return True


def glob_root(pattern: str) -> pathlib.PurePath:
    """Returns the leading directories of a glob pattern that contain no wildcards.

    Output paths are built relative to this, so ``data/*/x.parquet`` under
    ``output_dir=views`` lands at ``views/<subdir>/x.parquet``. A pattern naming a
    single file has no wildcard to stop at, so its own last component is the file
    rather than a directory and the root is the directory holding it.
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


def output_path(source: str, pattern: str, output_dir: str) -> pathlib.Path:
    """Maps a source path to its artifact path under ``output_dir``."""
    relative = pathlib.PurePath(source).relative_to(glob_root(pattern))
    return pathlib.Path(output_dir) / relative


def derive_tree(
    input_pattern: str,
    output_dir: str,
    *,
    artifact: str,
    write: ArtifactWriter,
    force: bool = False,
) -> tuple[int, int, int]:
    """Derives one artifact per dataset matching ``input_pattern``.

    Outputs mirror the inputs' directory layout beneath ``output_dir``. Matches that
    are themselves derived artifacts are ignored, so a recursive pattern reaching the
    output tree stays re-runnable. Artifacts whose footer already records the current
    source content, library version, and artifact version are skipped, so re-running is
    cheap.

    Args:
        input_pattern: Glob matching source Parquet datasets.
        output_dir: Directory to write artifacts beneath.
        artifact: Artifact name, used for the staleness check and the footer stamp.
        write: Writer taking ``(view, output, source_md5=...)`` and returning rows.
        force: Rewrite artifacts that are already current.

    Returns:
        ``(written, skipped, ignored)`` counts: artifacts derived, artifacts already
        current, and matches that were derived artifacts rather than sources.

    Raises:
        ValueError: If the input pattern matches nothing, if any match cannot be read as
            Parquet, or if any destination would land on a source.
    """
    matches = sorted(glob.glob(input_pattern, recursive=True))
    if not matches:
        raise ValueError(f"no datasets matched: {input_pattern}")
    sources = []
    ignored = 0
    for match in matches:
        if is_artifact(match):
            logger.info("%s is a derived artifact, not a source; ignoring", match)
            ignored += 1
        else:
            sources.append(match)
    logger.info("Found %d datasets", len(sources))
    # Every destination against every source, not against its own: an output_dir nested
    # in the source tree maps one dataset onto a *different* one, which a per-source
    # check passes, destroying a source the run had not reached yet.
    destinations = {
        source: output_path(source, input_pattern, output_dir) for source in sources
    }
    resolved_sources = {pathlib.Path(source).resolve() for source in sources}
    clobbered = sorted(
        str(destination)
        for destination in destinations.values()
        if destination.resolve() in resolved_sources
    )
    if clobbered:
        raise ValueError(
            f"output_dir {output_dir!r} would write over source datasets: {clobbered}"
        )
    written = skipped = 0
    for source in sources:
        destination = destinations[source]
        try:
            view = parquet.DatasetView(source)
        except ValueError as error:
            raise ValueError(f"{source} is not a source dataset: {error}") from error
        source_md5 = view.md5()
        if not force and is_current(destination, artifact, source_md5):
            logger.info("%s is current; skipping", destination)
            skipped += 1
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        rows = write(view, destination, source_md5=source_md5)
        logger.info("Wrote %d rows to %s", rows, destination)
        written += 1
    logger.info(
        "Derived %d %ss (%d already current, %d ignored)",
        written,
        artifact,
        skipped,
        ignored,
    )
    return written, skipped, ignored
