# Copyright 2022 Open Reaction Database Project Authors
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

"""Tests for ord_schema.orm.scripts.add_datasets.

Loading behavior is covered in ord_schema.orm.loading_test; these cover the CLI wiring
(argument parsing, stage validation) and that main delegates to loading.load_datasets.
"""

import pathlib

import pytest
from sqlalchemy.orm import Session

from ord_schema.orm.mappers import Mappers
from ord_schema.orm.scripts import add_datasets

_PBTXT_FIXTURE = str(
    pathlib.Path(__file__).parent / ".." / "testdata" / "ord-nielsen-example.pbtxt"
)


def test_main(prepared_engine):
    """main() parses args and runs both stages over the matching datasets."""
    argv = ["--dsn", str(prepared_engine.url), "--pattern", _PBTXT_FIXTURE]
    add_datasets.main(add_datasets.parse_args(argv))
    with Session(prepared_engine) as session:
        assert session.query(Mappers.Reaction).count() > 0


def test_unknown_stage_rejected():
    with pytest.raises(SystemExit):
        add_datasets.main(
            add_datasets.parse_args(["--pattern", _PBTXT_FIXTURE, "--stages", "bogus"])
        )


def test_classify_reactions_requires_extra(monkeypatch):
    """--classify_reactions exits early with a clear message when the extra is missing."""
    monkeypatch.setattr(add_datasets.database, "update_reaction_classes", None)
    with pytest.raises(SystemExit, match="reaction-class"):
        add_datasets.main(
            add_datasets.parse_args(
                ["--pattern", _PBTXT_FIXTURE, "--classify_reactions"]
            )
        )
