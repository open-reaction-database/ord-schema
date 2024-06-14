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
import pytest
from sqlalchemy import select
from sqlalchemy.exc import ProgrammingError

from ord_schema.orm.mappers import Mappers
from ord_schema.orm.rdkit_mappers import FingerprintType, RDKitMol, RDKitReaction


def test_tanimoto_operator(test_session):
    try:
        query = (
            select(Mappers.Reaction)
            .join(Mappers.ReactionInput)
            .join(Mappers.Compound)
            .join(RDKitMol)
            .where(RDKitMol.morgan_bfp % FingerprintType.MORGAN_BFP("c1ccccc1CCC(O)C"))
        )
        results = test_session.execute(query)
        assert len(results.fetchall()) == 20
    except ProgrammingError as error:
        pytest.skip(f"RDKit cartridge is required: {error}")


@pytest.mark.parametrize("fp_type", list(FingerprintType))
def test_tanimoto(test_session, fp_type):
    try:
        query = (
            select(Mappers.Reaction)
            .join(Mappers.ReactionInput)
            .join(Mappers.Compound)
            .join(RDKitMol)
            .where(RDKitMol.tanimoto("c1ccccc1CCC(O)C", fp_type=fp_type) > 0.5)
        )
        results = test_session.execute(query)
        assert len(results.fetchall()) == 20
    except ProgrammingError as error:
        pytest.skip(f"RDKit cartridge is required: {error}")


def test_substructure_operator(test_session):
    try:
        query = (
            select(Mappers.Reaction)
            .join(Mappers.ReactionInput)
            .join(Mappers.Compound)
            .join(RDKitMol)
            .where(RDKitMol.mol.op("@>")("c1ccccc1CCC(O)C"))
        )
        results = test_session.execute(query)
        assert len(results.fetchall()) == 20
    except ProgrammingError as error:
        pytest.skip(f"RDKit cartridge is required: {error}")


def test_substructure(test_session):
    try:
        query = (
            select(Mappers.Reaction)
            .join(Mappers.ReactionInput)
            .join(Mappers.Compound)
            .join(RDKitMol)
            .where(RDKitMol.check_substructure("c1ccccc1CCC(O)C"))
        )
        results = test_session.execute(query)
        assert len(results.fetchall()) == 20
    except ProgrammingError as error:
        pytest.skip(f"RDKit cartridge is required: {error}")


def test_smarts(test_session):
    try:
        query = (
            select(Mappers.Reaction)
            .join(Mappers.ReactionInput)
            .join(Mappers.Compound)
            .join(RDKitMol)
            .where(RDKitMol.check_smarts("c1ccccc1CCC(O)[#6]"))
        )
        results = test_session.execute(query)
        assert len(results.fetchall()) == 20
    except ProgrammingError as error:
        pytest.skip(f"RDKit cartridge is required: {error}")


def test_reaction_smarts(test_session):
    try:
        query = select(Mappers.Reaction).join(RDKitReaction).where(RDKitReaction.check_smarts("[#6]>>[#6]"))
        results = test_session.execute(query)
        assert len(results.fetchall()) == 20
    except ProgrammingError as error:
        pytest.skip(f"RDKit cartridge is required: {error}")
