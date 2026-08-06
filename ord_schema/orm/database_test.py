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
import pathlib
from typing import cast
from unittest.mock import Mock

import pytest
from rdkit import Chem
from sqlalchemy import select, text
from sqlalchemy.orm import Session

from ord_schema import message_helpers
from ord_schema.datasets import load_dataset
from ord_schema.orm import database as _orm_database
from ord_schema.orm.database import (
    add_dataset,
    backfill_submission_times,
    classify_dataset,
    delete_dataset,
    get_dataset_md5,
    get_dataset_size,
    update_derived_data,
    update_derived_tables,
    update_rdkit_ids,
    update_rdkit_tables,
)
from ord_schema.orm.mappers import Mappers, to_proto
from ord_schema.orm.public_mappers import DatasetMetadata
from ord_schema.proto import reaction_pb2

_HAS_REACTION_CLASS = _orm_database.update_reaction_classes is not None


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

    Guards the EXCEPT->NOT EXISTS rewrite: the fixture has a reaction with no SMILES
    (NULL reaction_smiles), which a missing IS NOT NULL guard would re-insert as a junk
    row on every run.
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
    # Invariant: the IS NOT NULL guards keep NULL mol/reaction rows out of the cartridge
    # tables.
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
    """Re-running update_derived_tables inserts no duplicate rows (NOT EXISTS guard)."""
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
    """A small batch size reproduces the full result, exercising multi-batch path."""
    tables = ("reaction_smiles", "compound_smiles", "product_compound_smiles")
    full = {
        table: test_session.execute(
            text(f"SELECT count(*) FROM derived.{table}")  # noqa: S608  (constant)
        ).scalar()
        for table in tables
    }
    # Every derived table must be populated up front, otherwise the re-derivation below
    # would pass trivially (0 == 0) without exercising the reaction or compound batch
    # paths.
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


def test_update_derived_tables_sharded_covers_all(prepared_engine):
    """Sharded SMILES derivation partitions every derived table and covers every row.

    Each shard derives a disjoint hash-partition of the reaction/compound ids (the
    derived-stage analog of row-group sharding for ingest). Guards two properties: the
    predicate actually splits each table's rows across more than one shard, and all
    shards together derive exactly what the unsharded pass would -- a following whole-
    dataset pass finds nothing new.
    """
    dataset = load_dataset(
        pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
        as_dataset=True,
    )
    with Session(prepared_engine) as session, session.begin():
        add_dataset(dataset, session)
    tables = ("reaction_smiles", "compound_smiles", "product_compound_smiles")

    def counts() -> dict[str, int]:
        with Session(prepared_engine) as session:
            return {
                table: session.execute(
                    text(f"SELECT count(*) FROM derived.{table}")  # noqa: S608  (constant)
                ).scalar()
                for table in tables
            }

    # Derive one shard at a time, recording how many rows each shard adds per table
    # (idempotent inserts, so a shard only adds its own partition's rows).
    num_shards = 4
    prev = dict.fromkeys(tables, 0)
    per_shard: list[dict[str, int]] = []
    with Session(prepared_engine) as session:
        for shard_index in range(num_shards):
            with session.begin():
                update_derived_tables(
                    dataset.dataset_id, session, shard=(shard_index, num_shards)
                )
            current = counts()
            per_shard.append({table: current[table] - prev[table] for table in tables})
            prev = current
    sharded = counts()
    for table in tables:
        assert sharded[table] > 0, (table, sharded)
        # Rows land in more than one shard, so the hash predicate really partitions each
        # table (guarded by a size floor so a genuinely tiny table can't flake the
        # check).
        if sharded[table] >= 2 * num_shards:
            non_empty = sum(1 for shard in per_shard if shard[table] > 0)
            assert non_empty >= 2, (table, [shard[table] for shard in per_shard])
    # A whole-dataset pass adds nothing: the shards already covered every row.
    with Session(prepared_engine) as session, session.begin():
        update_derived_tables(dataset.dataset_id, session)
    assert counts() == sharded


def test_rdkit_sharded_matches_serial(prepared_engine):
    """Sharded RDKit population inserts and links the same rows as the serial pass.

    Every RDKit sub-step partitions by the same SMILES hash, so a shard inserts exactly
    the structures its links reference. Guards that a single shard is a strict subset
    (partitioning is real) and that all shards together produce what the unsharded pass
    would -- a following serial pass inserts and links nothing new.
    """
    dataset = load_dataset(
        pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
        as_dataset=True,
    )
    with Session(prepared_engine) as session, session.begin():
        add_dataset(dataset, session)
    with Session(prepared_engine) as session, session.begin():
        update_derived_tables(
            dataset.dataset_id, session
        )  # SMILES first; RDKit reads them.

    def counts() -> dict[str, int]:
        with Session(prepared_engine) as session:
            queries = {
                "mols": "SELECT count(*) FROM rdkit.mols",
                "reactions": "SELECT count(*) FROM rdkit.reactions",
                "linked_reactions": (
                    "SELECT count(*) FROM derived.reaction_smiles "
                    "WHERE rdkit_reaction_id IS NOT NULL"
                ),
                "linked_compounds": (
                    "SELECT count(*) FROM derived.compound_smiles "
                    "WHERE rdkit_mol_id IS NOT NULL"
                ),
                "linked_products": (
                    "SELECT count(*) FROM derived.product_compound_smiles "
                    "WHERE rdkit_mol_id IS NOT NULL"
                ),
            }
            return {
                key: session.execute(text(query)).scalar()
                for key, query in queries.items()
            }

    def run_shard(shard: tuple[int, int] | None) -> None:
        with Session(prepared_engine) as session:
            with session.begin():
                update_rdkit_tables(dataset.dataset_id, session, shard=shard)
            with session.begin():
                update_rdkit_ids(dataset.dataset_id, session, shard=shard)

    num_shards = 4
    run_shard((0, num_shards))
    first = counts()
    for shard_index in range(1, num_shards):
        run_shard((shard_index, num_shards))
    sharded = counts()
    # Shard 0 alone was a strict, non-empty subset of the inserted structures.
    assert 0 < first["mols"] < sharded["mols"], (first, sharded)
    for key in (
        "mols",
        "reactions",
        "linked_reactions",
        "linked_compounds",
        "linked_products",
    ):
        assert sharded[key] > 0, (key, sharded)
    # A serial (whole-dataset) pass inserts and links nothing new: the shards
    # already covered it all.
    run_shard(None)
    assert counts() == sharded


@pytest.mark.skipif(
    not _HAS_REACTION_CLASS, reason="reaction-class extra not installed"
)
def test_classify_sharded_covers_all(prepared_engine):
    """Sharded classification partitions reactions and together classifies every one.

    Analogous to the SMILES/RDKit sharding tests; skipped unless the reaction-class
    extra is installed (it loads a transformer model per call). Shard 0 is a strict,
    non-empty subset and all shards together classify exactly what the unsharded pass
    would -- a following whole-dataset pass adds nothing.
    """
    dataset = load_dataset(
        pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
        as_dataset=True,
    )
    with Session(prepared_engine) as session, session.begin():
        add_dataset(dataset, session)
    with Session(prepared_engine) as session, session.begin():
        update_derived_tables(
            dataset.dataset_id, session
        )  # SMILES; classification reads them.

    def count() -> int:
        with Session(prepared_engine) as session:
            return session.execute(
                text("SELECT count(*) FROM derived.reaction_classes")
            ).scalar()

    num_shards = 4
    with Session(prepared_engine) as session, session.begin():
        classify_dataset(dataset.dataset_id, session, shard=(0, num_shards))
    first = count()
    with Session(prepared_engine) as session:
        for shard_index in range(1, num_shards):
            with session.begin():
                classify_dataset(
                    dataset.dataset_id, session, shard=(shard_index, num_shards)
                )
    sharded = count()
    # Shard 0 is a strict, non-empty subset: the fixture's ~80 reactions reliably hash
    # into more than one of the 4 buckets, so this holds for it (a 1-reaction fixture
    # would not).
    assert 0 < first < sharded, (first, sharded)
    # A whole-dataset pass adds nothing: the shards already classified every reaction.
    with Session(prepared_engine) as session, session.begin():
        classify_dataset(dataset.dataset_id, session)
    assert count() == sharded


@pytest.mark.parametrize(
    ("derived_table", "id_column", "mapper"),
    [
        ("derived.compound_smiles", "compound_id", Mappers.Compound),
        (
            "derived.product_compound_smiles",
            "product_compound_id",
            Mappers.ProductCompound,
        ),
    ],
)
def test_compound_smiles_match_reference(
    test_session, derived_table, id_column, mapper
):
    """The set-based SMILES equal the per-compound smiles_from_compound reference.

    Guards value parity of the bulk SMILES-identifier path against the
    message-reconstruction reference. The fallback (compounds without a stored
    SMILES, which this fixture does not contain) is covered by
    test_compound_smiles_fallback_without_stored_smiles.
    """
    rows = test_session.execute(
        text(f"SELECT {id_column}, smiles FROM {derived_table}")  # noqa: S608  (constant)
    ).all()
    assert rows, derived_table
    for compound_id, smiles in rows:
        compound = test_session.get(mapper, compound_id)
        expected = message_helpers.smiles_from_compound(
            cast(
                "reaction_pb2.Compound | reaction_pb2.ProductCompound",
                to_proto(compound),
            )
        )
        assert smiles == expected, (derived_table, compound_id, smiles, expected)


def test_compound_smiles_fallback_without_stored_smiles(prepared_engine):
    """A compound with no stored SMILES is derived via message-reconstruction fallback.

    The set-based fast path only covers compounds with a SMILES identifier; this
    exercises the session.get + to_proto + smiles_from_compound branch, which the
    standard fixture (every compound carries a SMILES) never reaches.
    """
    dataset = load_dataset(
        pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
        as_dataset=True,
    )
    # Replace one compound's SMILES identifier with the equivalent InChI so only the
    # fallback (structure reconstruction from a non-SMILES identifier) can derive its
    # SMILES.
    compound = next(
        component
        for reaction in dataset.reactions
        for reaction_input in reaction.inputs.values()
        for component in reaction_input.components
    )
    inchi = Chem.MolToInchi(
        Chem.MolFromSmiles(message_helpers.get_compound_smiles(compound))
    )
    del compound.identifiers[:]
    compound.identifiers.add(type=reaction_pb2.CompoundIdentifier.INCHI, value=inchi)
    expected = Chem.MolToSmiles(Chem.MolFromInchi(inchi))
    with Session(prepared_engine) as session:
        with session.begin():
            add_dataset(dataset, session)
        with session.begin():
            update_derived_tables(dataset.dataset_id, session)
        # Exactly the modified compound lacks a SMILES identifier, so its derived row
        # proves the fallback ran and produced the InChI-derived canonical SMILES.
        fallback_smiles = (
            session.execute(
                text(
                    "SELECT smiles FROM derived.compound_smiles cs "
                    "WHERE NOT EXISTS (SELECT 1 FROM ord.compound_identifier ci "
                    "WHERE ci.compound_id = cs.compound_id AND ci.type = 'SMILES')"
                )
            )
            .scalars()
            .all()
        )
    assert fallback_smiles == [expected], fallback_smiles


def test_underivable_compound_skips_reconstruction(prepared_engine, monkeypatch):
    """A compound with only non-structural identifiers is not reconstructed each run.

    Such a compound never yields a derived row, so the NOT EXISTS guard re-selects it on
    every derived pass. The set-based structural pre-filter must skip it before the
    expensive session.get + to_proto fallback (which would raise ValueError anyway);
    otherwise a database of underivable compounds pays a full reconstruction pass for
    zero new rows on every re-run.
    """
    dataset = load_dataset(
        pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
        as_dataset=True,
    )
    # Strip one compound down to a NAME-only identifier so it has no structural
    # identifier to derive a SMILES from; every other compound keeps its SMILES fast
    # path.
    compound = next(
        component
        for reaction in dataset.reactions
        for reaction_input in reaction.inputs.values()
        for component in reaction_input.components
    )
    del compound.identifiers[:]
    compound.identifiers.add(
        type=reaction_pb2.CompoundIdentifier.NAME, value="mystery reagent"
    )
    to_proto_calls = []
    real_to_proto = _orm_database.to_proto

    def _spy_to_proto(mapper):
        to_proto_calls.append(mapper)
        return real_to_proto(mapper)

    monkeypatch.setattr(_orm_database, "to_proto", _spy_to_proto)
    with Session(prepared_engine) as session:
        with session.begin():
            add_dataset(dataset, session)
        with session.begin():
            update_derived_tables(dataset.dataset_id, session)
        underivable = (
            session.execute(
                text(
                    "SELECT count(*) FROM ord.compound c "
                    "WHERE NOT EXISTS (SELECT 1 FROM ord.compound_identifier ci "
                    "WHERE ci.compound_id = c.id AND ci.type != 'NAME') "
                    "AND NOT EXISTS (SELECT 1 FROM derived.compound_smiles cs "
                    "WHERE cs.compound_id = c.id)"
                )
            )
            .scalars()
            .one()
        )
    # The NAME-only compound derives nothing (no derived row) ...
    assert underivable == 1, underivable
    # ... and the every-compound-has-SMILES fixture plus the pre-filter mean the
    # reconstruction fallback never runs, so to_proto is not called at all.
    assert to_proto_calls == [], to_proto_calls


# Counts every ord.compound by how it reaches its reaction, alongside how many of
# those rows made it into derived.compound_smiles and how many of those are still
# missing an rdkit.mols link. [Ti+5] structures are excluded from the unlinked tally
# because _update_rdkit_mols deliberately keeps them out of rdkit.mols (see the NOT
# LIKE guard there, and issue #672).
_COMPOUND_ATTACHMENT_COVERAGE = text(
    """
    SELECT CASE
             WHEN ord.reaction_input.reaction_id IS NOT NULL THEN 'reaction_input'
             WHEN ord.compound.reaction_input_id IS NOT NULL THEN 'workup_input'
             ELSE 'product_measurement'
           END AS attachment,
           count(*) AS compounds,
           count(derived.compound_smiles.compound_id) AS derived,
           count(*) FILTER (
               WHERE derived.compound_smiles.compound_id IS NOT NULL
                 AND derived.compound_smiles.rdkit_mol_id IS NULL
                 AND derived.compound_smiles.smiles NOT LIKE '%[Ti+5]%'
           ) AS unlinked
    FROM ord.compound
    LEFT JOIN ord.reaction_input
        ON ord.compound.reaction_input_id = ord.reaction_input.id
    LEFT JOIN derived.compound_smiles
        ON derived.compound_smiles.compound_id = ord.compound.id
    GROUP BY 1
"""
)


def test_derived_compound_smiles_covers_every_attachment(prepared_engine):
    """Every ord.compound is derived and RDKit-linked, whichever parent it hangs off.

    A compound reaches its reaction as a reaction input, as a workup's input (whose
    ord.reaction_input row has reaction_id NULL), or as a product measurement's
    authentic standard. A dataset-scoping join that follows only
    reaction_input.reaction_id silently drops the latter two, so assert all three
    attachments appear and none loses rows to the derived or RDKit passes.
    """
    dataset = load_dataset(
        pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
        as_dataset=True,
    )
    # The fixture has workup inputs but no authentic standard; attach one so the
    # product_measurement path is populated.
    measurement = next(
        measurement
        for reaction in dataset.reactions
        for outcome in reaction.outcomes
        for product in outcome.products
        for measurement in product.measurements
    )
    measurement.uses_authentic_standard = True
    measurement.authentic_standard.identifiers.add(
        type=reaction_pb2.CompoundIdentifier.SMILES, value="c1ccccc1"
    )
    with Session(prepared_engine) as session:
        with session.begin():
            add_dataset(dataset, session)
            update_derived_data(dataset.dataset_id, session)
        coverage = {
            row.attachment: row
            for row in session.execute(_COMPOUND_ATTACHMENT_COVERAGE)
        }
    # Guard against the assertions below passing vacuously on an attachment the fixture
    # lacks.
    assert set(coverage) == {
        "reaction_input",
        "workup_input",
        "product_measurement",
    }, coverage
    for attachment, row in coverage.items():
        assert row.compounds == row.derived, (attachment, row.compounds, row.derived)
        assert row.unlinked == 0, (attachment, row.unlinked)


def test_unlinked_partial_indexes(test_session):
    """prepare_database creates partial indexes on unlinked rows for cheap linking."""
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


def test_polymorphic_fk_indexes_are_partial(test_session):
    """Each polymorphic foreign-key index is partial (WHERE ...

    IS NOT NULL).     Under single-table inheritance a message's FK column is NULL for
    every row of a sibling     subclass, so a full index would store a NULL entry for
    most rows and every insert would touch     all of them. Indexing only the non-NULL
    rows keeps the index and the per-insert write     proportional to the rows that use
    each parent. ord.time is the widest such table (one FK per     possible parent), so
    it is a good representative.
    """
    indexes = dict(
        test_session.execute(
            text(
                "SELECT indexname, indexdef FROM pg_indexes "
                "WHERE schemaname = 'ord' AND tablename = 'time'"
            )
        ).all()
    )
    fk_indexes = {name: d for name, d in indexes.items() if name.startswith("ix_")}
    assert len(fk_indexes) >= 3, indexes  # ord.time has many parent FK columns.
    for name, indexdef in fk_indexes.items():
        assert "IS NOT NULL" in indexdef, f"{name} is not partial: {indexdef}"


def test_enum_types_in_ord_schema(test_session):
    """Mapped enum types live in the ord schema so they are rendered schema-qualified.

    Without inherit_schema the enum type is created wherever the create-time search_path
    points and referenced unqualified; that breaks ingest when the default search_path
    is pinned to public and the connecting role owns an eponymous ord schema (the prod
    'ord' role), so the unqualified cast cannot find the type. Pinning the types to ord
    makes resolution search_path-independent.
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


def test_surrogate_keys_are_uuidv7(test_session):
    """ord.* surrogate keys are uuid columns with version-7 (time-sortable) values."""
    types = dict(
        test_session.execute(
            text(
                "SELECT table_name || '.' || column_name, data_type "
                "FROM information_schema.columns "
                "WHERE table_schema = 'ord' "
                "AND ((table_name = 'reaction' AND column_name = 'id') "
                "  OR (table_name = 'compound' "
                "AND column_name IN ('id', 'reaction_input_id')))"
            )
        ).all()
    )
    assert types == {
        "reaction.id": "uuid",
        "compound.id": "uuid",
        "compound.reaction_input_id": "uuid",
    }, types
    # psycopg returns uuid columns as uuid.UUID; ingest mints version-7 ids.
    reaction_id = test_session.execute(
        text("SELECT id FROM ord.reaction LIMIT 1")
    ).scalar_one()
    assert reaction_id.version == 7


@pytest.mark.parametrize(
    ("version_num", "raises"),
    [("110000", True), ("120000", False), ("160005", False)],
)
def test_check_server_version(version_num, raises):
    """The version guard rejects pre-12 servers and accepts 12+."""
    connection = Mock()
    connection.execute.return_value.scalar_one.return_value = version_num
    if raises:
        with pytest.raises(RuntimeError, match="PostgreSQL 12"):
            _orm_database._check_server_version(connection)
    else:
        _orm_database._check_server_version(connection)


def test_check_server_version_memoized_per_connection():
    """The check queries once per connection, then reuses the cached result on info."""
    connection = Mock()
    connection.info = {}
    connection.execute.return_value.scalar_one.return_value = "160005"
    _orm_database._check_server_version(connection)
    _orm_database._check_server_version(connection)
    connection.execute.assert_called_once()


def test_default_search_path_is_public(test_session):
    """prepare_database pins the database default search_path to public, not ord."""
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
