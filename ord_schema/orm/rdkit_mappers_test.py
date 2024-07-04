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

"""Tests for ord_schema.orm.rdkit_mappers."""
import platform

import pytest
from sqlalchemy import cast, func, select

from ord_schema.orm.mappers import Mappers
from ord_schema.orm.rdkit_mappers import CString, FingerprintType, RDKitMol, RDKitMols, RDKitReactions

pytestmark = pytest.mark.skipif(platform.machine() != "x86_64", reason="RDKit cartridge is required")


def test_tanimoto_operator(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionInput)
        .join(Mappers.Compound)
        .join(RDKitMols)
        .where(RDKitMols.morgan_bfp % FingerprintType.MORGAN_BFP(cast("c1ccccc1CCC(O)C", RDKitMol)))
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20


@pytest.mark.parametrize("fp_type", list(FingerprintType))
def test_tanimoto(test_session, fp_type):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionInput)
        .join(Mappers.Compound)
        .join(RDKitMols)
        .where(RDKitMols.tanimoto("c1ccccc1CCC(O)C", fp_type=fp_type) > 0.5)
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20


def test_substructure_operator(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionInput)
        .join(Mappers.Compound)
        .join(RDKitMols)
        .where(RDKitMols.mol.op("@>")("c1ccccc1CCC(O)C"))
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20


def test_contains_substructure(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionInput)
        .join(Mappers.Compound)
        .join(RDKitMols)
        .where(RDKitMols.contains_substructure("c1ccccc1CCC(O)C"))
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20


def test_smarts_operator(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionInput)
        .join(Mappers.Compound)
        .join(RDKitMols)
        .where(RDKitMols.mol.op("@>")(func.qmol_from_smarts(cast("c1ccccc1CCC(O)[#6]", CString))))
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20


def test_matches_smarts(test_session):
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionInput)
        .join(Mappers.Compound)
        .join(RDKitMols)
        .where(RDKitMols.matches_smarts("c1ccccc1CCC(O)[#6]"))
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20


def test_reaction_smarts_operator(test_session):
    query = (
        select(Mappers.Reaction)
        .join(RDKitReactions)
        .where(
            RDKitReactions.reaction.op("@>")(func.reaction_from_smarts(cast("[#6:1].[#9:2]>>[#6:1][#9:2]", CString)))
        )
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 79  # One reaction has a BYPRODUCT outcome, so no reaction SMILES.


def test_reaction_matches_smarts(test_session):
    query = (
        select(Mappers.Reaction)
        .join(RDKitReactions)
        .where(RDKitReactions.matches_smarts("[#6:1].[#9:2]>>[#6:1][#9:2]"))
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 79  # One reaction has a BYPRODUCT outcome, so no reaction SMILES.
