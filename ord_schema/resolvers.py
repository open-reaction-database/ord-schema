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
"""Name/string resolution to structured messages or identifiers."""
import email.message
import re
import urllib.parse
import urllib.request
import urllib.error

from rdkit import Chem

from ord_schema.logging import get_logger
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)

_COMPOUND_STRUCTURAL_IDENTIFIERS = [
    reaction_pb2.CompoundIdentifier.SMILES,
    reaction_pb2.CompoundIdentifier.INCHI,
    reaction_pb2.CompoundIdentifier.MOLBLOCK,
    reaction_pb2.CompoundIdentifier.CXSMILES,
    reaction_pb2.CompoundIdentifier.XYZ,
]

_USERNAME = "github-actions"
_EMAIL = "github-actions@github.com"


def canonicalize_smiles(smiles: str) -> str:
    """Canonicalizes a SMILES string.

    Args:
        smiles: SMILES string.

    Returns:
        Canonicalized SMILES string.

    Raises:
        ValueError: If the SMILES cannot be parsed by RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    return Chem.MolToSmiles(mol)


def name_resolve(value_type: str, value: str) -> tuple[str, str]:
    """Resolves compound identifiers to SMILES via multiple APIs."""
    for resolver, resolver_func in _NAME_RESOLVERS.items():
        try:
            smiles = resolver_func(value_type, value)
            if smiles is not None:
                return smiles, resolver
        except urllib.error.HTTPError as error:
            logger.info(f"{resolver} could not resolve {value_type} {value}: {error}")
    raise ValueError(f"Could not resolve {value_type} {value} to SMILES")


def resolve_names(message: reaction_pb2.Reaction) -> bool:
    """Attempts to resolve compound NAME identifiers to SMILES.

    When a NAME identifier is resolved, a SMILES identifier is added to the list
    of identifiers for that compound. Note that this function moves on to the
    next Compound after the first successful name resolution.

    Args:
        message: Reaction proto.

    Returns:
        Boolean whether `message` was modified.
    """
    modified = False
    compounds = message_helpers.find_submessages(message, reaction_pb2.Compound)
    for compound in compounds:
        if any(identifier.type in _COMPOUND_STRUCTURAL_IDENTIFIERS for identifier in compound.identifiers):
            continue  # Compound already has a structural identifier.
        for identifier in compound.identifiers:
            if identifier.type == identifier.NAME:
                try:
                    smiles, resolver = name_resolve("name", identifier.value)
                    new_identifier = compound.identifiers.add()
                    new_identifier.type = new_identifier.SMILES
                    new_identifier.value = smiles
                    new_identifier.details = f"NAME resolved by the {resolver}"
                    modified = True
                    break
                except ValueError:
                    pass
    return modified


def _pubchem_resolve(value_type: str, value: str) -> str:
    """Resolves compound identifiers to SMILES via the PubChem REST API."""
    with urllib.request.urlopen(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{value_type}/"
        f"{urllib.parse.quote(value)}/property/IsomericSMILES/txt"
    ) as response:
        return response.read().decode().strip()


def _cactus_resolve(value_type: str, value: str) -> str:
    """Resolves compound identifiers to SMILES via the CACTUS API."""
    del value_type  # Unused.
    with urllib.request.urlopen(
        "https://cactus.nci.nih.gov/chemical/structure/" f"{urllib.parse.quote(value)}/smiles"
    ) as response:
        return response.read().decode().strip()


def _emolecules_resolve(value_type: str, value: str) -> str:
    """Resolves compound identifiers to SMILES via the eMolecules API."""
    del value_type  # Unused.
    with urllib.request.urlopen("https://www.emolecules.com/" f"lookup?q={urllib.parse.quote(value)}") as response:
        response_text = response.read().decode().strip()
    if response_text == "__END__":
        raise urllib.error.HTTPError("", 404, "eMolecules lookup unsuccessful", email.message.Message(), None)
    return response_text.split("\t")[0]


def resolve_input(input_string: str) -> reaction_pb2.ReactionInput:
    """Resolve a text-based description of an input in one of the following
    formats:
        (1) [AMOUNT] of [NAME]
        (2) [AMOUNT] of [CONCENTRATION] [SOLUTE] in [SOLVENT]

    Args:
        input_string: String describing the input.

    Returns:
        ReactionInput message.

    Raises:
        ValueError: if the string cannot be parsed properly.
    """
    reaction_input = reaction_pb2.ReactionInput()
    if " of " not in input_string:
        raise ValueError("String does not match template!")
    amount_string, description = input_string.split(" of ")
    if " in " not in description:
        component_name = description
        component = reaction_input.components.add()
        component.CopyFrom(message_helpers.build_compound(name=component_name.strip(), amount=amount_string))
        resolve_names(reaction_input)
        return reaction_input
    pattern = re.compile(r"(\d+.?\d*)\s?(\w+)\s(.+)\sin\s(.+)")
    match = pattern.fullmatch(description.strip())
    if not match:
        raise ValueError("String did not match template!")
    conc_value, conc_units, solute_name, solvent_name = match.groups()
    assert solute_name is not None  # Type hint.
    assert solvent_name is not None  # Type hint.
    solute = reaction_input.components.add()
    solvent = reaction_input.components.add()
    solute.CopyFrom(message_helpers.build_compound(name=solute_name.strip()))
    solvent.CopyFrom(message_helpers.build_compound(name=solvent_name.strip(), amount=amount_string))
    if solvent.amount.WhichOneof("kind") != "volume":
        raise ValueError("Total amount of solution must be a volume!")
    solvent.amount.volume_includes_solutes = True
    message_helpers.set_solute_moles(solute, [solvent], f"{conc_value} {conc_units}")
    resolve_names(reaction_input)
    return reaction_input


# Standard name resolvers.
_NAME_RESOLVERS = {
    "PubChem API": _pubchem_resolve,
    "NCI/CADD Chemical Identifier Resolver": _cactus_resolve,
    "eMolecules Lookup Service": _emolecules_resolve,
}
