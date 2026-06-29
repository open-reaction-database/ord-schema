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

"""Tests for ord_schema.orm.loading."""

import pathlib

import pytest
from sqlalchemy import select, text
from sqlalchemy.orm import Session

from ord_schema import message_helpers, parquet
from ord_schema.orm import database, loading
from ord_schema.orm.mappers import Mappers
from ord_schema.proto import dataset_pb2, reaction_pb2

_PBTXT_FIXTURE = str(
    pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt"
)


def _write_parquet_dataset(tmp_path) -> tuple[str, dataset_pb2.Dataset]:
    dataset = message_helpers.load_message(_PBTXT_FIXTURE, dataset_pb2.Dataset)
    parquet_path = (tmp_path / "dataset.parquet").as_posix()
    parquet.save_dataset(dataset, parquet_path)
    return parquet_path, dataset


def _derived_counts(engine) -> dict[str, int]:
    with Session(engine) as session:
        return {
            table: session.execute(
                text(f"SELECT count(*) FROM derived.{table}")  # noqa: S608  (constant)
            ).scalar()
            for table in (
                "reaction_smiles",
                "compound_smiles",
                "product_compound_smiles",
            )
        }


def test_load_datasets(prepared_engine):
    loading.load_datasets(_PBTXT_FIXTURE, str(prepared_engine.url))
    with Session(prepared_engine) as session:
        assert session.query(Mappers.Reaction).count() > 0
    assert _derived_counts(prepared_engine)["reaction_smiles"] > 0


def test_load_datasets_parquet(prepared_engine, tmp_path):
    """Streaming Parquet ingest stores the streaming MD5, reaction count, and Reaction rows."""
    parquet_path, dataset = _write_parquet_dataset(tmp_path)
    loading.load_datasets(parquet_path, str(prepared_engine.url))
    expected_md5, expected_count = parquet.streaming_md5(parquet_path)
    with Session(prepared_engine) as session:
        # Dataset row's metadata columns reflect the streaming MD5/count.
        assert database.get_dataset_md5(dataset.dataset_id, session) == expected_md5
        assert database.get_dataset_size(dataset.dataset_id, session) == expected_count
        # Reaction rows were actually inserted (not just the num_reactions column).
        mapped_dataset = session.scalar(
            select(Mappers.Dataset).where(
                Mappers.Dataset.dataset_id == dataset.dataset_id
            )
        )
        assert mapped_dataset is not None
        assert len(mapped_dataset.reactions) == expected_count
        # Spot-check that one reaction round-trips byte-identical through public.reactions.
        original = dataset.reactions[0]
        stored = next(
            r for r in mapped_dataset.reactions if r.reaction_id == original.reaction_id
        )
        assert reaction_pb2.Reaction.FromString(stored.proto_row.proto) == original


def test_load_datasets_parquet_skip_unchanged(prepared_engine, tmp_path):
    """A second invocation without overwrite is a no-op when the MD5 matches."""
    parquet_path, dataset = _write_parquet_dataset(tmp_path)
    loading.load_datasets(parquet_path, str(prepared_engine.url))
    loading.load_datasets(parquet_path, str(prepared_engine.url))
    # Skip path must not re-insert: reaction count is unchanged after the second call.
    _, expected_count = parquet.streaming_md5(parquet_path)
    with Session(prepared_engine) as session:
        mapped_dataset = session.scalar(
            select(Mappers.Dataset).where(
                Mappers.Dataset.dataset_id == dataset.dataset_id
            )
        )
        assert mapped_dataset is not None
        assert len(mapped_dataset.reactions) == expected_count


def test_load_datasets_parquet_rejects_changed_without_overwrite(
    prepared_engine, tmp_path
):
    """A re-ingest with mutated content but no overwrite raises."""
    parquet_path, dataset = _write_parquet_dataset(tmp_path)
    loading.load_datasets(parquet_path, str(prepared_engine.url))
    dataset.reactions[0].outcomes[0].conversion.value = 999.0
    parquet.save_dataset(dataset, parquet_path)
    with pytest.raises(
        RuntimeError
    ):  # load_datasets collects per-future failures and raises in aggregate.
        loading.load_datasets(parquet_path, str(prepared_engine.url))


def test_stages_ingest_then_derived(prepared_engine):
    """stages=ingest writes ord.*/public.* with no derived rows; a later stages=derived backfills them."""
    url = str(prepared_engine.url)
    loading.load_datasets(_PBTXT_FIXTURE, url, stages=["ingest"])
    with Session(prepared_engine) as session:
        assert session.query(Mappers.Reaction).count() > 0
    assert _derived_counts(prepared_engine) == {
        "reaction_smiles": 0,
        "compound_smiles": 0,
        "product_compound_smiles": 0,
    }
    loading.load_datasets(_PBTXT_FIXTURE, url, stages=["derived"])
    assert _derived_counts(prepared_engine)["reaction_smiles"] > 0
    with Session(prepared_engine) as session:
        assert (
            session.execute(text("SELECT count(*) FROM rdkit.reactions")).scalar() > 0
        )


def test_unknown_stage_raises():
    # Stage validation happens before any database access, so no engine is needed.
    with pytest.raises(ValueError, match="unknown stages"):
        loading.load_datasets(_PBTXT_FIXTURE, "postgresql://", stages=["bogus"])
