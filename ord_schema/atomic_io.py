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
"""Atomic-write helpers.

Centralizes the write-to-tmp-then-rename pattern used throughout
``ord_schema`` so partial files never appear at a destination path on
crash, exception, or interrupt.
"""

import contextlib
import os
import tempfile
from collections.abc import Iterator


@contextlib.contextmanager
def atomic_path(path: str) -> Iterator[str]:
    """Yields a sibling temp path; renames it onto ``path`` on clean exit.

    Usage:
        with atomic_path(dest) as tmp, open(tmp, "wb") as f:
            f.write(...)

    On clean exit the temp path is renamed onto ``path`` via ``os.replace``,
    which is atomic on POSIX when source and destination live on the same
    filesystem (the temp is a sibling, so this holds). On any exception
    (including ``KeyboardInterrupt``) the temp is removed best-effort and
    the original ``path`` is left untouched.

    The temp is created via ``tempfile.mkstemp`` in the destination's
    directory with a unique random suffix, so concurrent writers to the
    same ``path`` do not collide on a fixed sentinel name. The temp file
    is created (empty) immediately on entry; callers typically open it for
    writing and overwrite it.
    """
    parent = os.path.dirname(path) or "."
    basename = os.path.basename(path)
    fd, tmp = tempfile.mkstemp(dir=parent, prefix=basename + ".", suffix=".tmp")
    os.close(fd)
    try:
        yield tmp
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise
    else:
        os.replace(tmp, path)
