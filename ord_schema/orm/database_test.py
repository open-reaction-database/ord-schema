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

import pytest
from sqlalchemy import select, text

from ord_schema.orm.database import delete_dataset, get_dataset_md5, get_dataset_size, update_rdkit_tables
from ord_schema.orm.mappers import Mappers
from ord_schema.proto import reaction_pb2


def test_orm(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionOutcome)
        .join(Mappers.ProductCompound)
        .join(Mappers.ProductMeasurement)
        .join(Mappers.Percentage)
        .where(Mappers.ProductMeasurement.type == "YIELD", Mappers.Percentage.value >= 70)
    )
    results = test_session.execute(query)
    reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]
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
    before_reactions = test_session.execute(text("SELECT count(*) FROM rdkit.reactions")).scalar()
    before_mols = test_session.execute(text("SELECT count(*) FROM rdkit.mols")).scalar()
    assert before_reactions > 0 and before_mols > 0  # Sanity: the fixture populated the cartridge tables.
    update_rdkit_tables("test_dataset", test_session)
    assert test_session.execute(text("SELECT count(*) FROM rdkit.reactions")).scalar() == before_reactions
    assert test_session.execute(text("SELECT count(*) FROM rdkit.mols")).scalar() == before_mols
    # Invariant: the IS NOT NULL guards keep NULL mol/reaction rows out of the cartridge tables.
    assert test_session.execute(text("SELECT count(*) FROM rdkit.mols WHERE mol IS NULL")).scalar() == 0
    assert test_session.execute(text("SELECT count(*) FROM rdkit.reactions WHERE reaction IS NULL")).scalar() == 0


def test_get_dataset_md5(test_session):
    assert get_dataset_md5("test_dataset", test_session) == "0343d39a98d38eb39abd69d899af2bdf"
    assert get_dataset_md5("other_dataset", test_session) is None


def test_get_dataset_size(test_session):
    assert get_dataset_size("test_dataset", test_session) == 80
    with pytest.raises(ValueError):
        get_dataset_size("other_dataset", test_session)
