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
"""Macros for programmatic message creation."""
from typing import Optional, Union

from ord_schema.proto import reaction_pb2
from ord_schema import units

UNITS_RESOLVER = units.UnitResolver()
CONCENTRATION_RESOLVER = units.UnitResolver(units.CONCENTRATION_UNIT_SYNONYMS, forbidden_units={})


def simple_solution(
    solvent_smiles: str,
    solute_smiles: Optional[str] = None,
    volume: Union[None, str, reaction_pb2.Volume] = None,
    concentration: Union[None, str, reaction_pb2.Concentration] = None,
    saturated: bool = False,
) -> list[reaction_pb2.Compound]:
    """Creates a solution with at most one solvent and one solute.

    Args:
        solvent_smiles: SMILES of the solvent.
        solute_smiles: SMILES of the solute. If not specified, a pure solvent
            is generated.
        volume: Volume of solution (total volume including the solute).
        concentration: Concentration of solution. If both concentration and
            volume are specified, a solute quantity is computed. The schema
            does not currently have a way to represent an unknown volume of
            solution with a known concentration.
        saturated: Whether the solution is saturated. `Saturated` and
            `concentrated` cannot both be specified.
    Returns: A list of solvent/solute Compounds.
    """
    if saturated:
        if concentration is not None:
            raise ValueError("Cannot specify both `saturated=True` and a concentration.")
        if solute_smiles is None:
            raise ValueError("Must specify a solute if `saturated=True`")
    if isinstance(volume, str):
        volume = UNITS_RESOLVER.resolve(volume)
    if isinstance(concentration, str):
        concentration = CONCENTRATION_RESOLVER.resolve(concentration)

    if volume is not None and concentration is not None:
        solute_amount = units.compute_solute_quantity(volume, concentration)
    elif saturated:
        solute_amount = reaction_pb2.Amount(unmeasured=reaction_pb2.UnmeasuredAmount(type="SATURATED"))
    else:
        # This case covers 1. pure solvents 2. solution with unknown volume.
        solute_amount = None

    output_compounds = []
    solvent_pb = reaction_pb2.Compound()
    solvent_pb.identifiers.add(value=solvent_smiles, type="SMILES")
    if volume is not None:
        solvent_pb.amount.volume.MergeFrom(volume)
        if solute_smiles is not None:
            solvent_pb.amount.volume_includes_solutes = True
    output_compounds.append(solvent_pb)

    if solute_smiles:
        solute_pb = reaction_pb2.Compound()
        solute_pb.identifiers.add(value=solute_smiles, type="SMILES")
        if solute_amount is not None:
            solute_pb.amount.MergeFrom(solute_amount)
        output_compounds.append(solute_pb)
    return output_compounds


def brine(volume=None):
    return simple_solution(solute_smiles="[Na+].[Cl-]", solvent_smiles="O", volume=volume, saturated=True)
