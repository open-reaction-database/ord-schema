# Copyright 2020 Open Reaction Database Project Authors
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
"""Tests for ord_schema.units."""

import os
import pytest
import xmltodict

from ord_schema import units
from ord_schema.proto import reaction_pb2
from ord_schema.translators import scilligence


def test_scilligence_translator():
    with open(os.path.join(os.path.dirname(__file__),
            'scilligence_example.xml')) as fid:
        reaction = scilligence.convert_reaction(fid.read())
    assert len(reaction.inputs) == 5
    assert all(get_compound_name(_input.compounds[0])
        for _input in reaction.inputs)
    assert reaction.provenance.record_reacted.time.value != ''
    # TODO: write additional tests


def test_scilligence_translator_should_fail():
    with pytest.raises(ValueError, match=r"invalid xml"):
        scilligence.convert_reaction('garbage string')
    garbage_xml = xmltodict.unparse({'garbage': 'xml'})
    with pytest.raises(ValueError, match=r"scilligence eln"):
        scilligence.convert_reaction(garbage_xml)
