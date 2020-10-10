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
import re
import urllib.parse
import urllib.request
import urllib.error
import uuid

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

_USERNAME = 'github-actions'
_EMAIL = 'github-actions@github.com'


def canonicalize_smiles(smiles):
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
        raise ValueError(f'Could not parse SMILES: {smiles}')
    return Chem.MolToSmiles(mol)


def name_resolve(value_type, value):
    """Resolves compound identifiers to SMILES via multiple APIs."""
    for resolver, resolver_func in _NAME_RESOLVERS.items():
        try:
            smiles = resolver_func(value_type, value)
            if smiles is not None:

                return smiles, resolver
        except urllib.error.HTTPError as error:
            logging.info('%s could not resolve %s %s: %s', resolver, value_type,
                         value, error)
    raise ValueError(f'Could not resolve {value_type} {value} to SMILES')


def _pubchem_resolve(value_type, value):
    """Resolves compound identifiers to SMILES via the PubChem REST API."""
    response = urllib.request.urlopen(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{value_type}/'
        f'{urllib.parse.quote(value)}/property/IsomericSMILES/txt')
    return response.read().decode().strip()


def _cactus_resolve(value_type, value):
    """Resolves compound identifiers to SMILES via the CACTUS API."""
    del value_type
    response = urllib.request.urlopen(
        'https://cactus.nci.nih.gov/chemical/structure/'
        f'{urllib.parse.quote(value)}/smiles')
    return response.read().decode().strip()


def _emolecules_resolve(value_type, value):
    """Resolves compound identifiers to SMILES via the eMolecules API."""
    del value_type
    response = urllib.request.urlopen(
        'https://www.emolecules.com/lookup?q={}'.format(
            urllib.parse.quote(value)))
    response_text = response.read().decode().strip()
    if response_text == '__END__':
        raise urllib.error.HTTPError(None, None,
                                     'eMolecules lookup unsuccessful', None,
                                     None)
    return response_text.split('\t')[0]


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
                    smiles, resolver = name_resolve('name', identifier.value)
                    new_identifier = compound.identifiers.add()
                    new_identifier.type = new_identifier.SMILES
                    new_identifier.value = smiles
                    new_identifier.details = f'NAME resolved by the {resolver}'
                    modified = True
                    break
                except ValueError:
                    pass
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
        event.person.username = _USERNAME
        event.person.email = _EMAIL
        event.details = 'Automatic updates from the submission pipeline.'
    return id_substitutions


def update_dataset(dataset):
    """Updates a Dataset message.

    Current updates:
      * All reaction-level updates in update_reaction.
      * reaction_id cross-references between Reactions in the dataset.

    Args:
        dataset: dataset_pb2.Dataset message.

    Raises:
        KeyError: if the dataset has not been validated and there exists a
            cross-referenced reaction_id in any Reaction that is not defined
            elsewhere in the Dataset.
    """
    # Reaction-level updates
    id_substitutions = {}
    for reaction in dataset.reactions:
        id_substitutions.update(update_reaction(reaction))
    # Dataset-level updates of cross-references
    for reaction in dataset.reactions:
        for reaction_input in reaction.inputs.values():
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

# Standard name resolvers.
_NAME_RESOLVERS = {
    'PubChem API': _pubchem_resolve,
    'NCI/CADD Chemical Identifier Resolver': _cactus_resolve,
    'eMolecules Lookup Service': _emolecules_resolve,
}
