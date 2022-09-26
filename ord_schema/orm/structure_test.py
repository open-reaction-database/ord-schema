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

"""Tests for ord_schema.orm.structure."""
import pytest
from sqlalchemy import select
from ord_schema.orm.mappers import Compound, Reaction, ReactionInput
from ord_schema.orm.structure import FingerprintType, Structure


def test_tanimoto_operator(test_session):
    query = (
        select(Reaction)
        .join(ReactionInput)
        .join(Compound)
        .join(Structure)
        .where(Structure.morgan_bfp % FingerprintType.MORGAN_BFP("c1ccccc1CCC(O)C"))
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20


@pytest.mark.parametrize("fp_type", list(FingerprintType))
def test_tanimoto(test_session, fp_type):
    query = (
        select(Reaction)
        .join(ReactionInput)
        .join(Compound)
        .join(Structure)
        .where(Structure.tanimoto("c1ccccc1CCC(O)C", fp_type=fp_type) > 0.5)
    )
    results = test_session.execute(query)
    assert len(results.fetchall()) == 20
