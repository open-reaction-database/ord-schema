from typing import List, Union

from ord_schema.proto import reaction_pb2
from ord_schema import units

UNITS_RESOLVER = units.UnitResolver()
CONCENTRATION_RESOLVER = units.UnitResolver(units.CONCENTRATION_UNIT_SYNONYMS, forbidden_units={})

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

    output_compounds = []
    if solvent_smiles:
        component = reaction_pb2.Compound()
        component.identifiers.add(value=solvent_smiles, type='SMILES')
        if solute_smiles:
            component.amount.volume_includes_solutes = True
        component.amount.volume.MergeFrom(volume)
        output_compounds.append(component)

    if solute_smiles:
        component = reaction_pb2.Compound()
        component.identifiers.add(value=solute_smiles, type='SMILES')
        component.amount.MergeFrom(solute_amount)
        output_compounds.append(component)
    return output_compounds

def BRINE(volume=None):
    return simple_solution(solute_smiles='[Na+][Cl-]', solvent_smiles='O',
        volume=volume, concentration='6.14 M')
