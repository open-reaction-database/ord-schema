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
"""Automated updates for Reaction messages."""

import datetime
import urllib
import uuid
import re

from absl import logging
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
    compounds = message_helpers.find_submessages(message,
                                                 reaction_pb2.Compound)
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
                except urllib.error.HTTPError as error:
                    logging.info('PubChem could not resolve NAME %s: %s',
                                 identifier.value, error)
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
    compounds = message_helpers.find_submessages(message,
                                                 reaction_pb2.Compound)
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
                compound.identifiers.add(bytes_value=mol.ToBinary(),
                                         type='RDKIT_BINARY',
                                         details=f'Generated from {source}')
                modified = True
                break  # Only add one RDKIT_BINARY per Compound.
    return modified


def update_reaction(reaction):
    """Updates a Reaction message.

    Current updates:
      * Sets reaction_id if not already set.
      * Adds a record modification event to the provenance.
      * Resolves compound identifier names if no structural identifiers are
        defined for a given compound.

    Args:
        reaction: reaction_pb2.Reaction message.

    Returns:
        A dictionary mapping placeholder reaction_ids to newly-assigned
            reaction_ids.
    """
    modified = False
    id_substitutions = {}
    if not reaction.provenance.HasField('record_created'):
        reaction.provenance.record_created.time.value = (
            datetime.datetime.utcnow().ctime())
        modified = True
    if not re.fullmatch('^ord-[0-9a-f]{32}$', reaction.reaction_id):
        # NOTE(kearnes): This does not check for the case where a Dataset is
        # edited and reaction_id values are changed inappropriately. This will
        # need to be either (1) caught in review or (2) found by a complex
        # check of the diff.
        reaction_id = f'ord-{uuid.uuid4().hex}'
        if reaction.reaction_id:
            id_substitutions[reaction.reaction_id] = reaction_id
        reaction.reaction_id = reaction_id
        modified = True
    for func in _UPDATES:
        # NOTE(kearnes): Order is important here; if you write
        # `modified or func(reaction)` and modified is True, the interpreter
        # will skip the evaluation of func(reaction).
        modified = func(reaction) or modified
    if modified:
        event = reaction.provenance.record_modified.add()
        event.time.value = datetime.datetime.utcnow().ctime()
        event.details = 'Automatic updates from the submission pipeline.'
    return id_substitutions


def update_dataset(dataset):
    """Updates a Dataset message.

    Current updates:
      * All reaction-level updates in update_reaction.
      * reaction_id cross-references between Reactions in the dataset.

    Args:
        dataset: dataset_pb2.Dataset message.
    """
    # Reaction-level updates
    id_substitutions = {}
    for reaction in dataset.reactions:
        id_substitutions.update(update_reaction(reaction))
    # Dataset-level updates of cross-references
    # Note: if the Dataset has been validated, then we know that all
    # cross-referened reaction_ids correspond to placeholder reaction_ids of
    # another Reaction in the Dataset. If not validated, this update might
    # raise a KeyError exception.
    for reaction in dataset.reactions:
        for reaction_input in reaction.inputs:
            for component in reaction_input.components:
                for preparation in component.preparations:
                    if preparation.reaction_id:
                        preparation.reaction_id = id_substitutions[
                            preparation.reaction_id]
            for crude_component in reaction_input.crude_components:
                crude_component.reaction_id = id_substitutions[
                    crude_component.reaction_id]


# Standard updates.
_UPDATES = [
    resolve_names,
]
