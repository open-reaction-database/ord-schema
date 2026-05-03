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
from collections.abc import Iterator


@contextlib.contextmanager
def atomic_path(path: str) -> Iterator[str]:
    """Yields a sibling temp path; renames it onto ``path`` on clean exit.

    Usage:
        with atomic_path(dest) as tmp:
            with open(tmp, "wb") as f:
                f.write(...)

    On clean exit the temp path is renamed onto ``path`` via ``os.replace``,
    which is atomic on POSIX when source and destination live on the same
    filesystem (the temp is a sibling, so this holds). On any exception the
    temp is removed best-effort and the original ``path`` is left untouched.

    The temp uses a fixed ``.tmp`` suffix on the destination, matching the
    convention already used by ``scripts/process_dataset.py``. Concurrent
    writers to the same destination will collide; switch to
    ``tempfile.mkstemp`` if that ever becomes a real concern.
    """
    tmp = path + ".tmp"
    try:
        yield tmp
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise
    else:
        os.replace(tmp, path)
