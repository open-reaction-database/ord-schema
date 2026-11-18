# Copyright 2026 Open Reaction Database Project Authors
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
"""Tests for ord_schema.macros.workups."""

import textwrap

import pytest

from ord_schema import validations
from ord_schema.macros import solutions, workups
from ord_schema.proto import reaction_pb2


def test_add_solution_default_type():
    solution = solutions.simple_solution(solvent_smiles="O", volume="100mL")
    workup = workups.add_solution(solution)
    assert workup.type == reaction_pb2.ReactionWorkup.ADDITION
    assert len(workup.input.components) == len(solution)
    for component in workup.input.components:
        assert component.reaction_role == reaction_pb2.ReactionRole.WORKUP
    validations.validate_message(workup)


@pytest.mark.parametrize("workup_type", ["ADDITION", "EXTRACTION", "WASH"])
def test_add_solution_type_options(workup_type):
    solution = solutions.simple_solution(
        solvent_smiles="O",
        solute_smiles="[Na+].[Cl-]",
        concentration="1M",
        volume="100mL",
    )
    workup = workups.add_solution(solution, workup_type=workup_type)
    assert workup.type == getattr(reaction_pb2.ReactionWorkup, workup_type)
    assert len(workup.input.components) == 2
    for component in workup.input.components:
        assert component.reaction_role == reaction_pb2.ReactionRole.WORKUP
    validations.validate_message(workup)


@pytest.mark.parametrize("keep_phase", ["aqueous", "organic"])
def test_separate_phases(keep_phase):
    workup = workups.separate_phases(keep_phase=keep_phase)
    assert workup.type == reaction_pb2.ReactionWorkup.EXTRACTION
    assert workup.keep_phase == keep_phase
    validations.validate_message(workup)


def test_drying_agent():
    workup = workups.drying_agent("[Mg+2].[O-]S([O-])(=O)=O")
    assert workup.type == reaction_pb2.ReactionWorkup.DRY_WITH_MATERIAL
    assert len(workup.input.components) == 1
    component = workup.input.components[0]
    assert component.reaction_role == reaction_pb2.ReactionRole.WORKUP
    assert len(component.identifiers) == 1
    identifier = component.identifiers[0]
    assert identifier.type == reaction_pb2.CompoundIdentifier.SMILES
    assert identifier.value == "[Mg+2].[O-]S([O-])(=O)=O"
    validations.validate_message(workup)


def test_drying_agent_records_cxsmiles_as_cxsmiles():
    workup = workups.drying_agent("C[C@H](N)C(=O)O |&1:1|")
    identifier = workup.input.components[0].identifiers[0]
    assert identifier.type == reaction_pb2.CompoundIdentifier.CXSMILES
    assert not validations.validate_message(workup).warnings


@pytest.mark.parametrize("keep_phase", ["filtrate", "solid"])
def test_filtration(keep_phase):
    workup = workups.filtration(keep_phase=keep_phase)
    assert workup.type == reaction_pb2.ReactionWorkup.FILTRATION
    assert workup.keep_phase == keep_phase
    validations.validate_message(workup)


def test_rotovap():
    workup = workups.rotovap()
    assert workup.type == reaction_pb2.ReactionWorkup.CONCENTRATION
    validations.validate_message(workup)


def test_module_docstring_example():
    assert workups.__doc__ is not None
    _, marker, example = workups.__doc__.partition("is recorded as::")
    assert marker, "the module docstring's example follows 'is recorded as::'"
    reaction = reaction_pb2.Reaction()
    # The code is this module's own docstring.
    exec(textwrap.dedent(example), {"reaction": reaction})  # noqa: S102
    assert len(reaction.workups) == 11
    for workup in reaction.workups:
        assert not validations.validate_message(workup, raise_on_error=False).errors
