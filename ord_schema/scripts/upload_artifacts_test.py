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

"""Tests for ord_schema.scripts.upload_artifacts."""

import pytest

from ord_schema import huggingface
from ord_schema.artifacts import occurrences, pivot
from ord_schema.scripts import upload_artifacts


def test_main_requires_every_level_and_path(monkeypatch):
    calls = []

    def fake_upload(root, **kwargs):
        calls.append((root, kwargs))
        return huggingface.ArtifactUpload(
            prefixes={}, uploaded=[], deleted=[], commit_url=None
        )

    monkeypatch.setattr(huggingface, "upload_artifacts", fake_upload)
    upload_artifacts.main(
        [
            "--root",
            "artifacts",
            "--message",
            "derive",
            "--revision",
            "next",
            "--dry_run",
        ]
    )
    assert calls == [
        (
            "artifacts",
            {
                "commit_message": "derive",
                "revision": "next",
                "pivot_levels": sorted(pivot.LEVELS),
                "occurrence_paths": sorted(occurrences.PATHS),
                "dry_run": True,
            },
        )
    ]


def test_message_is_required():
    with pytest.raises(SystemExit, match="2"):
        upload_artifacts.parse_args(["--root", "artifacts"])
