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
"""Macros for quickly creating workup steps.

Example usage:

from ord_schema.macros import workups
from ord_schema.macros import solutions

The reaction was quenched with 300 ml of 0.5M NaHCO3 solution and the phases were separated.
The aqueous layer was back-extracted twice with 100 ml of EtOAc. The organic layers were
combined and washed with 200 ml of saturated NaCl solution, dried over MgSO4 and filtered.
Upon concentration to an oil, the crude product was purified by normal phase column
chromatography

reaction.workups.MergeFrom([
    workups.add_solution(solutions.simple_solution(
        solvent_smiles='O', solute_smiles='[Na+][HCO3-]', volume='300mL', concentration='0.5M')),
    workups.separate_phases(keep_phase='organic'),
    workups.add_solution(solutions.simple_solution(
        solvent_smiles='CC(=O)OCC', volume='100mL'), type='EXTRACTION'),
    workups.separate_phases(keep_phase='organic'),
    workups.add_solution(solutions.simple_solution(
        solvent_smiles='CC(=O)OCC', volume='100mL'), type='EXTRACTION'),
    workups.separate_phases(keep_phase='organic'),
    workups.add_solution(solutions.brine('200mL'), type='WASH'),
    workups.drying_agent('[Mg+2].[O-]S([O-])(=O)=O'),
    workups.filter(keep_phase='filtrate'),
    workups.rotovap(),
    reaction_pb2.ReactionWorkup(type='OTHER_CHROMATOGRAPHY'),
])
"""
from collections.abc import Iterable

from ord_schema.proto import reaction_pb2
from ord_schema import units

UNITS_RESOLVER = units.UnitResolver()
CONCENTRATION_RESOLVER = units.UnitResolver(units.CONCENTRATION_UNIT_SYNONYMS)


def add_solution(
    solution: Iterable[reaction_pb2.Compound], workup_type: str = "ADDITION"
) -> reaction_pb2.ReactionWorkup:
    """Create a workup representing addition of a solution.

    type is commonly one of 'ADDITION', 'EXTRACTION', or 'WASH'; see
    ReactionWorkup.WorkupType enum for full list of possible values.
    """
    workup = reaction_pb2.ReactionWorkup(type=workup_type)
    workup.input.components.MergeFrom(solution)
    for component in workup.input.components:
        component.reaction_role = reaction_pb2.ReactionRole.WORKUP
    return workup


def separate_phases(keep_phase: str) -> reaction_pb2.ReactionWorkup:
    """Create a workup representing a phase separation.

    keep_phase is commonly either 'aqueous' or 'organic'.
    """
    return reaction_pb2.ReactionWorkup(type="EXTRACTION", keep_phase=keep_phase)


def drying_agent(agent_smiles: str) -> reaction_pb2.ReactionWorkup:
    """Create a workup representing addition of a drying agent."""
    workup = reaction_pb2.ReactionWorkup(type="DRY_WITH_MATERIAL")
    component = workup.input.components.add()
    component.identifiers.add(value=agent_smiles, type="SMILES")
    component.reaction_role = reaction_pb2.ReactionRole.WORKUP
    return workup


def filtration(keep_phase: str) -> reaction_pb2.ReactionWorkup:
    """Create a workup representing a filtration step.

    keep_phase should be one of 'filtrate' or 'solid'.
    """
    return reaction_pb2.ReactionWorkup(type="FILTRATION", keep_phase=keep_phase)


def rotovap() -> reaction_pb2.ReactionWorkup:
    """Create a workup representing a rotary evaporation step."""
    return reaction_pb2.ReactionWorkup(type="CONCENTRATION")
