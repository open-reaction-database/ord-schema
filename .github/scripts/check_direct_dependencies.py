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

"""Rejects PEP 508 direct references in built distribution metadata.

PyPI refuses any upload whose ``Requires-Dist`` carries a direct reference -- ``name @ url``,
such as a ``git+https://`` pin -- with ``400 Can't have direct dependency``. ``twine check``
inspects only the description and core metadata fields, so it passes such a distribution and
the rejection lands mid-release, after the version bump has been committed.

Usage:
    python .github/scripts/check_direct_dependencies.py dist/*
"""

import sys
import tarfile
import zipfile
from email import message_from_string
from pathlib import Path

from packaging.requirements import Requirement


def read_metadata(path: Path) -> str:
    """Reads the core metadata out of a built distribution.

    Args:
        path: Path to a ``.whl`` or ``.tar.gz`` distribution.

    Returns:
        Contents of the distribution's METADATA (wheel) or PKG-INFO (sdist).

    Raises:
        ValueError: If ``path`` is not a wheel or sdist, or carries no metadata.
    """
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            for name in archive.namelist():
                if name.endswith(".dist-info/METADATA"):
                    return archive.read(name).decode()
        raise ValueError(f"no .dist-info/METADATA in {path}")
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path) as archive:
            for member in archive.getmembers():
                # Only the top-level PKG-INFO describes the distribution itself; the
                # sdist also carries one per .egg-info directory.
                if member.name.count("/") == 1 and member.name.endswith("/PKG-INFO"):
                    handle = archive.extractfile(member)
                    if handle is not None:
                        return handle.read().decode()
        raise ValueError(f"no top-level PKG-INFO in {path}")
    raise ValueError(f"not a wheel or sdist: {path}")


def direct_references(metadata: str) -> list[str]:
    """Returns the Requires-Dist entries that carry a direct URL reference.

    Args:
        metadata: Core metadata text, as returned by ``read_metadata``.

    Returns:
        The offending requirement strings, in declaration order.
    """
    message = message_from_string(metadata)
    return [dep for dep in message.get_all("Requires-Dist", []) if Requirement(dep).url]


def main(paths: list[str]) -> int:
    """Returns a process exit status; 1 if any distribution declares a direct dependency.

    Args:
        paths: Distribution files to check.

    Returns:
        0 when every distribution is free of direct references, else 1.

    Raises:
        ValueError: If ``paths`` is empty, which would otherwise pass vacuously.
    """
    if not paths:
        raise ValueError("no distributions to check")
    offenders = 0
    for path in paths:
        for dep in direct_references(read_metadata(Path(path))):
            print(f"{path}: PyPI rejects direct dependency: {dep}")
            offenders += 1
    if offenders:
        print(
            f"\n{offenders} direct dependency(ies) found. PyPI will reject this upload with "
            "400 Can't have direct dependency. Depend on a released version instead."
        )
        return 1
    print(f"No direct dependencies in {len(paths)} distribution(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
