# Copyright 2020 The Open Reaction Database Authors
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

from typing import List, Union

from ord_schema.proto import reaction_pb2
from ord_schema import units

UNITS_RESOLVER = units.UnitResolver()
CONCENTRATION_RESOLVER = units.UnitResolver(units.CONCENTRATION_UNIT_SYNONYMS,
                                            forbidden_units={})


def simple_solution(
    solvent_smiles=None,
    solute_smiles=None,
    volume: Union[None, str, reaction_pb2.Volume] = None,
    concentration: Union[None, str, reaction_pb2.Concentration] = None
) -> List[reaction_pb2.Compound]:
    if isinstance(volume, str):
        volume = UNITS_RESOLVER.resolve(volume)
    if isinstance(concentration, str):
        concentration = CONCENTRATION_RESOLVER.resolve(concentration)

    if volume is not None and concentration is not None:
        solute_amount = units.compute_solute_quantity(volume, concentration)
    else:
        solute_amount = None

    output_compounds = []
    solvent_pb = reaction_pb2.Compound()
    solvent_pb.identifiers.add(value=solvent_smiles, type='SMILES')
    if volume is not None:
        solvent_pb.amount.volume.MergeFrom(volume)
        if solute_smiles is not None:
            solvent_pb.amount.volume_includes_solutes = True
    output_compounds.append(solvent_pb)

    if solute_smiles:
        solute_pb = reaction_pb2.Compound()
        solute_pb.identifiers.add(value=solute_smiles, type='SMILES')
        if solute_amount is not None:
            solute_pb.amount.MergeFrom(solute_amount)
        output_compounds.append(solute_pb)
    return output_compounds


def BRINE(volume=None):
    return simple_solution(solute_smiles='[Na+].[Cl-]',
                           solvent_smiles='O',
                           volume=volume,
                           concentration='6.14 M')
