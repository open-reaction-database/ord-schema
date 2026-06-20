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
"""Tests for ord_schema.atomic_io."""

import pathlib

import pytest

from ord_schema import atomic_io


def _read(path: str) -> bytes:
    with pathlib.Path(path).open("rb") as f:
        return f.read()


def _write(path: str, data: bytes) -> None:
    with pathlib.Path(path).open("wb") as f:
        f.write(data)


def test_atomic_path_publishes_on_clean_exit(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    with atomic_io.atomic_path(dest) as tmp:
        # The temp lives in the destination's directory with a unique random
        # suffix bracketed by the basename and ``.tmp``.
        assert pathlib.Path(tmp).parent == pathlib.Path(dest).parent
        assert tmp != dest
        assert pathlib.Path(tmp).name.startswith("out.bin.")
        assert tmp.endswith(".tmp")
        _write(tmp, b"hello")
        assert not pathlib.Path(dest).exists(), "destination must not exist before exit"
    assert _read(dest) == b"hello"
    # Only the destination remains; no orphan temp.
    assert sorted(q.name for q in pathlib.Path(tmp_path).iterdir()) == ["out.bin"]


def test_atomic_path_replaces_existing_destination(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    _write(dest, b"old")
    with atomic_io.atomic_path(dest) as tmp:
        _write(tmp, b"new")
    assert _read(dest) == b"new"
    assert sorted(q.name for q in pathlib.Path(tmp_path).iterdir()) == ["out.bin"]


def test_atomic_path_removes_temp_on_exception(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    with pytest.raises(RuntimeError, match="boom"), atomic_io.atomic_path(dest) as tmp:  # noqa: PT012  (abort-on-exception test: the multi-statement body in the context is under test)
        _write(tmp, b"partial")
        raise RuntimeError("boom")
    assert sorted(q.name for q in pathlib.Path(tmp_path).iterdir()) == []


def test_atomic_path_preserves_destination_on_exception(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    _write(dest, b"original")
    with pytest.raises(RuntimeError, match="boom"), atomic_io.atomic_path(dest) as tmp:  # noqa: PT012  (abort-on-exception test: the multi-statement body in the context is under test)
        _write(tmp, b"partial")
        raise RuntimeError("boom")
    assert _read(dest) == b"original"
    assert sorted(q.name for q in pathlib.Path(tmp_path).iterdir()) == ["out.bin"]


def test_atomic_path_removes_temp_on_keyboard_interrupt(tmp_path):
    """KeyboardInterrupt (a BaseException) also triggers cleanup; documents the contract."""
    dest = (tmp_path / "out.bin").as_posix()
    with pytest.raises(KeyboardInterrupt), atomic_io.atomic_path(dest) as tmp:  # noqa: PT012  (abort-on-exception test: the multi-statement body in the context is under test)
        _write(tmp, b"partial")
        raise KeyboardInterrupt
    assert sorted(q.name for q in pathlib.Path(tmp_path).iterdir()) == []


def test_atomic_path_swallows_missing_temp_on_exception(tmp_path):
    """Cleanup must succeed even if the temp was already removed by something else."""
    dest = (tmp_path / "out.bin").as_posix()
    with pytest.raises(RuntimeError, match="boom"), atomic_io.atomic_path(dest) as tmp:  # noqa: PT012  (abort-on-exception test: the multi-statement body in the context is under test)
        pathlib.Path(tmp).unlink()  # Remove the temp out from under us.
        raise RuntimeError("boom")
    assert sorted(q.name for q in pathlib.Path(tmp_path).iterdir()) == []


def test_atomic_path_unique_temps_for_concurrent_callers(tmp_path):
    """Two overlapping atomic_path contexts on the same destination get distinct temps."""
    dest = (tmp_path / "out.bin").as_posix()
    with atomic_io.atomic_path(dest) as tmp1, atomic_io.atomic_path(dest) as tmp2:
        assert tmp1 != tmp2
