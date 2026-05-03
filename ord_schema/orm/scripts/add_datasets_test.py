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

import pytest
from sqlalchemy import select
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session

from ord_schema import message_helpers, parquet_dataset
from ord_schema.orm import database
from ord_schema.orm.database import prepare_database
from ord_schema.orm.mappers import Mappers
from ord_schema.orm.scripts import add_datasets
from ord_schema.proto import dataset_pb2, reaction_pb2

_PBTXT_FIXTURE = os.path.join(os.path.dirname(__file__), "..", "testdata", "ord-nielsen-example.pbtxt")


@pytest.fixture(name="prepared_engine")
def prepared_engine_fixture(test_engine) -> Engine:
    """``test_engine`` with the ORM schema (and RDKit cartridge) installed."""
    assert prepare_database(test_engine)
    return test_engine


def _write_parquet_dataset(tmp_path) -> tuple[str, dataset_pb2.Dataset]:
    dataset = message_helpers.load_message(_PBTXT_FIXTURE, dataset_pb2.Dataset)
    parquet_path = (tmp_path / "dataset.parquet").as_posix()
    parquet_dataset.write_dataset(dataset, parquet_path)
    return parquet_path, dataset


def test_main(prepared_engine):
    argv = ["--dsn", str(prepared_engine.url), "--pattern", _PBTXT_FIXTURE]
    add_datasets.main(add_datasets.parse_args(argv))


def test_main_parquet(prepared_engine, tmp_path):
    """Streaming Parquet ingest stores the streaming MD5, reaction count, and Reaction rows."""
    parquet_path, dataset = _write_parquet_dataset(tmp_path)
    argv = ["--dsn", str(prepared_engine.url), "--pattern", parquet_path]
    add_datasets.main(add_datasets.parse_args(argv))
    expected_md5, expected_count = parquet_dataset.streaming_md5(parquet_path)
    with Session(prepared_engine) as session:
        # Dataset row's metadata columns reflect the streaming MD5/count.
        assert database.get_dataset_md5(dataset.dataset_id, session) == expected_md5
        assert database.get_dataset_size(dataset.dataset_id, session) == expected_count
        # Reaction rows were actually inserted (not just the num_reactions column).
        mapped_dataset = session.scalar(select(Mappers.Dataset).where(Mappers.Dataset.dataset_id == dataset.dataset_id))
        assert mapped_dataset is not None
        assert len(mapped_dataset.reactions) == expected_count
        # Spot-check that one reaction round-trips byte-identical through the proto column.
        original = dataset.reactions[0]
        stored = next(r for r in mapped_dataset.reactions if r.reaction_id == original.reaction_id)
        assert reaction_pb2.Reaction.FromString(stored.proto) == original


def test_main_parquet_skip_unchanged(prepared_engine, tmp_path):
    """A second invocation without --overwrite is a no-op when the MD5 matches."""
    parquet_path, _ = _write_parquet_dataset(tmp_path)
    argv = ["--dsn", str(prepared_engine.url), "--pattern", parquet_path]
    add_datasets.main(add_datasets.parse_args(argv))
    add_datasets.main(add_datasets.parse_args(argv))


def test_main_parquet_rejects_changed_without_overwrite(prepared_engine, tmp_path):
    """A re-ingest with mutated content but no --overwrite raises."""
    parquet_path, dataset = _write_parquet_dataset(tmp_path)
    argv = ["--dsn", str(prepared_engine.url), "--pattern", parquet_path]
    add_datasets.main(add_datasets.parse_args(argv))
    dataset.reactions[0].outcomes[0].conversion.value = 999.0
    parquet_dataset.write_dataset(dataset, parquet_path)
    with pytest.raises(RuntimeError):  # main() collects per-future failures and raises in aggregate.
        add_datasets.main(add_datasets.parse_args(argv))
