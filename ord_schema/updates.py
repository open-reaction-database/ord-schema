"""Automated updates for Reaction messages."""

import datetime
import urllib
import uuid

from rdkit import Chem

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

_COMPOUND_STRUCTURAL_IDENTIFIERS = [
    reaction_pb2.CompoundIdentifier.SMILES,
    reaction_pb2.CompoundIdentifier.INCHI,
    reaction_pb2.CompoundIdentifier.MOLBLOCK,
    reaction_pb2.CompoundIdentifier.CXSMILES,
    reaction_pb2.CompoundIdentifier.XYZ,
]


def _pubchem_resolve(value_type, value):
    """Resolves compound identifiers to SMILES via the PubChem REST API."""
    response = urllib.request.urlopen(
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'
        f'{value_type}/{value}/property/IsomericSMILES/txt')
    return response.read().decode().strip()


def resolve_names(message):
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
        if any(identifier.type in _COMPOUND_STRUCTURAL_IDENTIFIERS
               for identifier in compound.identifiers):
            continue  # Compound already has a structural identifier.
        for identifier in compound.identifiers:
            if identifier.type == identifier.NAME:
                try:
                    smiles = _pubchem_resolve('name', identifier.value)
                    new_identifier = compound.identifiers.add()
                    new_identifier.type = new_identifier.SMILES
                    new_identifier.value = smiles
                    new_identifier.details = 'NAME resolved by PubChem'
                    modified = True
                    break
                except urllib.error.HTTPError:
                    pass
    return modified


def add_binary_identifiers(message):
    """Adds RDKIT_BINARY identifiers for compounds with valid structures.

    Note that the RDKIT_BINARY representations are mostly useful in the context
    of searching the database. Accordingly, this function is not included in the
    standard set of Reaction updates in update_reaction().

    Args:
        message: Reaction proto.

    Returns:
        Boolean whether `message` was modified.
    """
    modified = False
    compounds = message_helpers.find_submessages(message, reaction_pb2.Compound)
    for compound in compounds:
        if any(identifier.type == identifier.RDKIT_BINARY
               for identifier in message.identifiers):
            continue
        for identifier in compound.identifiers:
            mol = None
            if Chem and identifier.type == identifier.SMILES:
                mol = Chem.MolFromSmiles(identifier.value)
            elif identifier.type == identifier.INCHI:
                mol = Chem.MolFromInchi(identifier.value)
            elif identifier.type == identifier.MOLBLOCK:
                mol = Chem.MolFromMolBlock(identifier.value)
            if mol is not None:
                source = reaction_pb2.CompoundIdentifier.IdentifierType.Name(
                    identifier.type)
                compound.identifiers.add(
                    bytes_value=mol.ToBinary(),
                    type='RDKIT_BINARY',
                    details=f'Generated from {source}')
                modified = True
                break  # Only add one RDKIT_BINARY per Compound.
    return modified


def update_reaction(reaction):
    """Updates a Reaction message.

    Current updates:
      * Sets record_id if not already set.

    Args:
        reaction: reaction_pb2.Reaction message.
    """
    # NOTE(kearnes): Provenance updates do not trigger "modified" events.
    if not reaction.provenance.HasField('record_created'):
        reaction.provenance.record_created.time.value = (
            datetime.datetime.utcnow().ctime())
    # TODO(kearnes): There's more complexity here regarding record_id values
    # that are set outside of the submission pipeline. It's not as simple as
    # requiring them to be empty, since a submission may include edits to
    # existing reactions.
    if not reaction.provenance.record_id:
        record_id = f'ord-{uuid.uuid4().hex}'
        reaction.provenance.record_id = record_id
    modified = False
    for func in _UPDATES:
        modified = modified or func(reaction)
    if modified:
        reaction.provenance.record_modified.add().time.value = (
            datetime.datetime.utcnow().ctime())


# Standard updates.
_UPDATES = [
    resolve_names,
]
