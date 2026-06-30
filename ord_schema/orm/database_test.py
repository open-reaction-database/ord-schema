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

"""Tests for ord_schema.orm.database."""

import datetime

import pytest
from sqlalchemy import select, text

from ord_schema.orm.database import (
    backfill_submission_times,
    delete_dataset,
    get_dataset_md5,
    get_dataset_size,
    update_derived_tables,
    update_rdkit_tables,
)
from ord_schema.orm.mappers import Mappers
from ord_schema.orm.public_mappers import DatasetMetadata
from ord_schema.proto import reaction_pb2


def test_orm(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionOutcome)
        .join(Mappers.ProductCompound)
        .join(Mappers.ProductMeasurement)
        .join(Mappers.Percentage)
        .where(
            Mappers.ProductMeasurement.type == "YIELD", Mappers.Percentage.value >= 70
        )
    )
    results = test_session.execute(query)
    reactions = [
        reaction_pb2.Reaction.FromString(result[0].proto_row.proto)
        for result in results
    ]
    assert len(reactions) == 12


def test_delete_dataset(test_session):
    assert test_session.query(Mappers.Reaction).count() == 80
    delete_dataset("test_dataset", test_session)
    assert test_session.query(Mappers.Reaction).count() == 0


def test_update_rdkit_tables_idempotent(test_session):
    """Re-running update_rdkit_tables inserts no duplicate rows.

    Guards the EXCEPT->NOT EXISTS rewrite: the fixture has a reaction with no
    SMILES (NULL reaction_smiles), which a missing IS NOT NULL guard would
    re-insert as a junk row on every run.
    """
    before_reactions = test_session.execute(
        text("SELECT count(*) FROM rdkit.reactions")
    ).scalar()
    before_mols = test_session.execute(text("SELECT count(*) FROM rdkit.mols")).scalar()
    assert before_reactions > 0
    assert before_mols > 0
    update_rdkit_tables("test_dataset", test_session)
    assert (
        test_session.execute(text("SELECT count(*) FROM rdkit.reactions")).scalar()
        == before_reactions
    )
    assert (
        test_session.execute(text("SELECT count(*) FROM rdkit.mols")).scalar()
        == before_mols
    )
    # Invariant: the IS NOT NULL guards keep NULL mol/reaction rows out of the cartridge tables.
    assert (
        test_session.execute(
            text("SELECT count(*) FROM rdkit.mols WHERE mol IS NULL")
        ).scalar()
        == 0
    )
    assert (
        test_session.execute(
            text("SELECT count(*) FROM rdkit.reactions WHERE reaction IS NULL")
        ).scalar()
        == 0
    )


def test_update_derived_tables_idempotent(test_session):
    """Re-running update_derived_tables inserts no duplicate derived rows (NOT EXISTS guard)."""
    tables = ("reaction_smiles", "compound_smiles", "product_compound_smiles")
    before = {
        table: test_session.execute(
            text(f"SELECT count(*) FROM derived.{table}")  # noqa: S608  (constant)
        ).scalar()
        for table in tables
    }
    assert before["reaction_smiles"] > 0
    update_derived_tables("test_dataset", test_session)
    for table, count in before.items():
        assert (
            test_session.execute(
                text(f"SELECT count(*) FROM derived.{table}")  # noqa: S608  (constant)
            ).scalar()
            == count
        )


def test_update_derived_tables_batched(test_session, monkeypatch):
    """A small batch size reproduces the full result, exercising the multi-batch path."""
    tables = ("reaction_smiles", "compound_smiles", "product_compound_smiles")
    full = {
        table: test_session.execute(
            text(f"SELECT count(*) FROM derived.{table}")  # noqa: S608  (constant)
        ).scalar()
        for table in tables
    }
    # Every derived table must be populated up front, otherwise the re-derivation below would
    # pass trivially (0 == 0) without exercising the reaction or compound batch paths.
    assert all(count > 0 for count in full.values()), full
    # Clear the derived rows, then re-derive in batches far smaller than the dataset (80
    # reactions) so the result is built across many batches rather than one.
    for table in tables:
        test_session.execute(text(f"DELETE FROM derived.{table}"))  # noqa: S608  (constant)
    monkeypatch.setattr("ord_schema.orm.database._DERIVED_BATCH", 7)
    update_derived_tables("test_dataset", test_session)
    for table, count in full.items():
        assert (
            test_session.execute(
                text(f"SELECT count(*) FROM derived.{table}")  # noqa: S608  (constant)
            ).scalar()
            == count
        )


def test_unlinked_partial_indexes(test_session):
    """prepare_database creates partial indexes over unlinked rows to keep incremental linking cheap."""
    indexes = dict(
        test_session.execute(
            text(
                "SELECT indexname, indexdef FROM pg_indexes "
                "WHERE schemaname IN ('ord', 'derived')"
            )
        ).all()
    )
    for name, predicate in (
        ("reaction_smiles_unlinked_index", "rdkit_reaction_id IS NULL"),
        ("compound_smiles_unlinked_index", "rdkit_mol_id IS NULL"),
        ("product_compound_smiles_unlinked_index", "rdkit_mol_id IS NULL"),
    ):
        assert name in indexes, f"missing index {name}"
        assert predicate in indexes[name], (
            f"{name} is not partial on {predicate!r}: {indexes[name]}"
        )


def test_enum_types_in_ord_schema(test_session):
    """Mapped enum types live in the ord schema so they are rendered schema-qualified.

    Without inherit_schema the enum type is created wherever the create-time
    search_path points and referenced unqualified; that breaks ingest when the
    default search_path is pinned to public and the connecting role owns an
    eponymous ord schema (the prod 'ord' role), so the unqualified cast cannot find
    the type. Pinning the types to ord makes resolution search_path-independent.
    """
    schemas = test_session.execute(
        text(
            "SELECT t.typname, n.nspname FROM pg_type t "
            "JOIN pg_namespace n ON n.oid = t.typnamespace "
            "WHERE t.typtype = 'e'"
        )
    ).all()
    assert schemas, "no mapped enum types found"
    for typname, schema in schemas:
        assert schema == "ord", f"{typname} is in {schema}, expected ord"


def test_default_search_path_is_public(test_session):
    """prepare_database pins the database default search_path to public, not the role's ord schema."""
    setting = test_session.execute(
        text(
            "SELECT array_to_string(s.setconfig, ',') "
            "FROM pg_db_role_setting s JOIN pg_database d ON d.oid = s.setdatabase "
            "WHERE d.datname = current_database() AND s.setrole = 0"
        )
    ).scalar()
    assert setting is not None
    assert "search_path=public" in setting


def test_get_dataset_md5(test_session):
    assert (
        get_dataset_md5("test_dataset", test_session)
        == "0343d39a98d38eb39abd69d899af2bdf"
    )
    assert get_dataset_md5("other_dataset", test_session) is None


def test_get_dataset_size(test_session):
    assert get_dataset_size("test_dataset", test_session) == 80
    with pytest.raises(ValueError, match="other_dataset"):
        get_dataset_size("other_dataset", test_session)


def test_submitted_at(test_session):
    # add_dataset() populates submitted_at from a reaction's last record_modified
    # entry; the fixture's reactions were modified on 2021-02-25 and created on
    # 2020-11-28, so record_modified must win.
    submitted_at = test_session.execute(
        select(DatasetMetadata.submitted_at)
    ).scalar_one()
    assert submitted_at == datetime.date(2021, 2, 25)


def test_backfill_submission_times(test_session):
    test_session.execute(text("UPDATE public.datasets SET submitted_at = NULL"))
    backfill_submission_times(test_session)
    submitted_at = test_session.execute(
        select(DatasetMetadata.submitted_at)
    ).scalar_one()
    assert submitted_at is not None
