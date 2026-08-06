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
them live here: the stamps, the output-path mapping, and the driver that walks a glob
and skips what is already current. An artifact module supplies only its schema, its
projection, and its name.

Every artifact carries the same five keys, in the ``ord.`` namespace that
``ord_schema.parquet`` already uses for source datasets:

* ``ord.artifact`` -- which artifact this is, so a reader can tell a view from a
  projection without inspecting the schema.
* ``ord.artifact_version`` -- one version across *all* artifacts, deliberately. A
  derivation change usually touches shared helpers, and a reader comparing artifacts to
  each other needs to know they were built by the same definition. Per-artifact versions
  would let a view and a projection of the same dataset disagree while both looked
  current.
* ``ord.source_dataset_id``, ``ord.source_md5`` -- what it was derived from.
* ``ord.ord_schema_version`` -- what derived it.

``is_current`` requires content, library, and definition to all match, since any of the
three changing can change a value.
"""

import dataclasses
import glob
import os
import pathlib
from collections.abc import Callable
from importlib import metadata

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


def stamps(artifact: str, source_dataset_id: str | None, source_md5: str) -> Stamps:
    """Returns the stamps to write for an artifact derived now.

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
        return False
    return (
        value.artifact == artifact
        and value.source_md5 == source_md5
        and value.ord_schema_version == metadata.version("ord-schema")
        and value.artifact_version == ARTIFACT_VERSION
    )


def is_artifact(path: str | os.PathLike[str]) -> bool:
    """Returns whether ``path`` carries derived-artifact stamps.

    A derived artifact is never a source: a recursive glob that reaches the output tree
    would otherwise try to derive an artifact from an artifact.
    """
    try:
        load_stamps(path)
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
    write: Callable[..., int],
    force: bool = False,
) -> tuple[int, int]:
    """Derives one artifact per dataset matching ``input_pattern``.

    Outputs mirror the inputs' directory layout beneath ``output_dir``. Artifacts whose
    footer already records the current source content, library version, and artifact
    version are skipped, so re-running is cheap. Matches that are themselves derived
    artifacts are ignored, so a recursive pattern reaching the output tree stays
    re-runnable.

    Args:
        input_pattern: Glob matching source Parquet datasets.
        output_dir: Directory to write artifacts beneath.
        artifact: Artifact name, used for the staleness check and the footer stamp.
        write: Writer taking ``(view, output, source_md5=...)`` and returning rows.
        force: Rewrite artifacts that are already current.

    Returns:
        ``(written, skipped)`` counts, where ``skipped`` counts artifacts already
        current. Ignored matches are logged rather than counted.

    Raises:
        ValueError: If the input pattern matches nothing, or if a destination would
            land on its own source.
    """
    sources = sorted(glob.glob(input_pattern, recursive=True))
    if not sources:
        raise ValueError(f"no datasets matched: {input_pattern}")
    logger.info("Found %d datasets", len(sources))
    written = skipped = ignored = 0
    for source in sources:
        if is_artifact(source):
            logger.warning("%s is a derived artifact, not a source; ignoring", source)
            ignored += 1
            continue
        destination = output_path(source, input_pattern, output_dir)
        if destination.resolve() == pathlib.Path(source).resolve():
            # An output_dir that lands back on the source tree would publish the
            # artifact over the dataset it was derived from, losing the source.
            raise ValueError(
                f"{source} maps to itself under output_dir {output_dir!r}; "
                "an artifact cannot be written over its own source"
            )
        view = parquet.DatasetView(source)
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
    return written, skipped
