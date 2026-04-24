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

"""Tests for ord_schema.orm.scripts.add_datasets."""

import os

from sqlalchemy.orm import Session

from ord_schema import message_helpers, parquet_dataset
from ord_schema.orm import database
from ord_schema.orm.database import prepare_database
from ord_schema.orm.scripts import add_datasets
from ord_schema.proto import dataset_pb2

_PBTXT_FIXTURE = os.path.join(os.path.dirname(__file__), "..", "testdata", "ord-nielsen-example.pbtxt")


def test_main(test_engine):
    assert prepare_database(test_engine)
    argv = [
        "--dsn",
        str(test_engine.url),
        "--pattern",
        _PBTXT_FIXTURE,
    ]
    add_datasets.main(add_datasets.parse_args(argv))


def test_main_parquet(test_engine, tmp_path):
    rdkit_cartridge = prepare_database(test_engine)
    fixture = message_helpers.load_message(_PBTXT_FIXTURE, dataset_pb2.Dataset)
    parquet_path = (tmp_path / "dataset.parquet").as_posix()
    message_helpers.write_message(fixture, parquet_path)
    argv = ["--dsn", str(test_engine.url), "--pattern", parquet_path]
    add_datasets.main(add_datasets.parse_args(argv))
    # Verify the streaming path stored the streaming MD5 and the reaction count.
    expected_md5, expected_count = parquet_dataset.streaming_md5(parquet_path)
    with Session(test_engine) as session:
        assert database.get_dataset_md5(fixture.dataset_id, session) == expected_md5
        assert database.get_dataset_size(fixture.dataset_id, session) == expected_count
    # A second invocation without --overwrite is a no-op when the MD5 matches.
    add_datasets.main(add_datasets.parse_args(argv))
    del rdkit_cartridge  # Cartridge handling is exercised via main()'s add_rdkit step.
