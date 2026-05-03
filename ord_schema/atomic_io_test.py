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

import os

import pytest

from ord_schema import atomic_io


def _read(path: str) -> bytes:
    with open(path, "rb") as f:
        return f.read()


def _write(path: str, data: bytes) -> None:
    with open(path, "wb") as f:
        f.write(data)


def test_atomic_path_publishes_on_clean_exit(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    with atomic_io.atomic_path(dest) as tmp:
        assert tmp == dest + ".tmp"
        _write(tmp, b"hello")
        assert not os.path.exists(dest), "destination must not exist before exit"
    assert _read(dest) == b"hello"
    assert not os.path.exists(dest + ".tmp")


def test_atomic_path_replaces_existing_destination(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    _write(dest, b"old")
    with atomic_io.atomic_path(dest) as tmp:
        _write(tmp, b"new")
    assert _read(dest) == b"new"


def test_atomic_path_removes_temp_on_exception(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    with pytest.raises(RuntimeError, match="boom"), atomic_io.atomic_path(dest) as tmp:
        _write(tmp, b"partial")
        raise RuntimeError("boom")
    assert not os.path.exists(dest)
    assert not os.path.exists(dest + ".tmp")


def test_atomic_path_preserves_destination_on_exception(tmp_path):
    dest = (tmp_path / "out.bin").as_posix()
    _write(dest, b"original")
    with pytest.raises(RuntimeError, match="boom"), atomic_io.atomic_path(dest) as tmp:
        _write(tmp, b"partial")
        raise RuntimeError("boom")
    assert _read(dest) == b"original"
    assert not os.path.exists(dest + ".tmp")


def test_atomic_path_swallows_missing_temp_on_exception(tmp_path):
    """Exception before the writer creates the temp file: cleanup is a no-op, original error still propagates."""
    dest = (tmp_path / "out.bin").as_posix()
    with pytest.raises(RuntimeError, match="boom"), atomic_io.atomic_path(dest):
        raise RuntimeError("boom")
    assert not os.path.exists(dest)
    assert not os.path.exists(dest + ".tmp")
