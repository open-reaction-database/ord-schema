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
"""Helper functions for constructing Protocol Buffer messages."""

import enum
import os

import flask
from google import protobuf
from google.protobuf import json_format
from google.protobuf import text_format
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from ord_schema import units
from ord_schema.proto import reaction_pb2

_COMPOUND_IDENTIFIER_LOADERS = {
    reaction_pb2.CompoundIdentifier.SMILES: Chem.MolFromSmiles,
    reaction_pb2.CompoundIdentifier.INCHI: Chem.MolFromInchi,
    reaction_pb2.CompoundIdentifier.MOLBLOCK: Chem.MolFromMolBlock,
}

# pylint: disable=too-many-arguments
# pylint: disable=too-many-branches


def build_compound(smiles=None,
                   name=None,
                   amount=None,
                   role=None,
                   is_limiting=None,
                   prep=None,
                   prep_details=None,
                   vendor=None):
    """Builds a Compound message with the most common fields.

    Args:
        smiles: Text compound SMILES.
        name: Text compound name.
        amount: Text amount string, e.g. '1.25 g'.
        role: Text reaction role. Must match a value in ReactionRoleType.
        is_limiting: Boolean whether this compound is limiting for the reaction.
        prep: Text compound preparation type. Must match a value in
            PreparationType.
        prep_details: Text compound preparation details. If provided, `prep` is
            required.
        vendor: Text compound vendor/supplier.

    Returns:
        Compound message.

    Raises:
        KeyError: if `role` or `prep` does not match a supported enum value.
        TypeError: if `amount` units are not supported.
        ValueError: if `prep_details` is provided and `prep` is None.
    """
    compound = reaction_pb2.Compound()
    if smiles:
        compound.identifiers.add(value=smiles, type='SMILES')
    if name:
        compound.identifiers.add(value=name, type='NAME')
    if amount:
        resolver = units.UnitResolver()
        amount_pb = resolver.resolve(amount)
        if isinstance(amount_pb, reaction_pb2.Mass):
            compound.mass.CopyFrom(amount_pb)
        elif isinstance(amount_pb, reaction_pb2.Moles):
            compound.moles.CopyFrom(amount_pb)
        elif isinstance(amount_pb, reaction_pb2.Volume):
            compound.volume.CopyFrom(amount_pb)
        else:
            raise TypeError(f'unsupported units for amount: {amount_pb}')
    if role:
        field = reaction_pb2.Compound.DESCRIPTOR.fields_by_name['reaction_role']
        values_dict = field.enum_type.values_by_name
        try:
            compound.reaction_role = values_dict[role.upper()].number
        except KeyError as error:
            raise KeyError(
                f'{role} is not a supported type: {values_dict.keys()}'
            ) from error
    if is_limiting is not None:
        if not (is_limiting is True or is_limiting is False):
            raise TypeError(
                f'is_limiting must be a boolean value: {is_limiting}')
        compound.is_limiting = is_limiting
    if prep:
        field = reaction_pb2.CompoundPreparation.DESCRIPTOR.fields_by_name[
            'type']
        values_dict = field.enum_type.values_by_name
        try:
            compound.preparations.add().type = values_dict[prep.upper()].number
        except KeyError as error:
            raise KeyError(
                f'{prep} is not a supported type: {values_dict.keys()}'
            ) from error
        if (compound.preparations[0].type
                == reaction_pb2.CompoundPreparation.CUSTOM and
                not prep_details):
            raise ValueError(
                'prep_details must be provided when CUSTOM prep is used')
    if prep_details:
        if not prep:
            raise ValueError('prep must be provided when prep_details is used')
        compound.preparations[0].details = prep_details
    if vendor:
        compound.vendor_source = vendor
    return compound


def set_solute_moles(solute, solvents, concentration, overwrite=False):
    """Helps define components for stock solution inputs with a single solute
    and a one or more solvent compounds.

    Args:
        solute: Compound with identifiers, roles, etc.; this argument is
            modified in place to define an amount in moles.
        solvents: list of Compounds each with defined volume.
        concentration: string defining solute concentration.
        overwrite: whether to overwrite an existing solute amount if defined.
            Defaults to False

    Raises:
        ValueError: if any solvent does not have a defined volume.
        ValueError: if the solute has an existing amount field and overwrite
            is set to False.

    Returns:
        List of Compounds to assign to a repeated components field.
    """
    # Check solute definition
    if solute.WhichOneof('amount') and not overwrite:
        raise ValueError('solute has defined amount and overwrite is False')

    # Get total solvent volume in liters.
    volume_liter = 0
    for solvent in solvents:
        if not solvent.HasField('volume') or not solvent.volume.value:
            raise ValueError('solvent must have defined volume')
        if solvent.volume.units == solvent.volume.LITER:
            volume_liter += solvent.volume.value
        elif solvent.volume.units == solvent.volume.MILLILITER:
            volume_liter += solvent.volume.value * 1e-3
        elif solvent.volume.units == solvent.volume.MICROLITER:
            volume_liter += solvent.volume.value * 1e-6
        elif solvent.volume.units == solvent.volume.NANOLITER:
            volume_liter += solvent.volume.value * 1e-9
        else:
            raise ValueError('solvent units not recognized by set_solute_moles')
    # Get solute concentration in molar.
    resolver = units.UnitResolver(
        unit_synonyms=units.CONCENTRATION_UNIT_SYNONYMS, forbidden_units={})
    concentration_pb = resolver.resolve(concentration)
    if concentration_pb.units == concentration_pb.MOLAR:
        concentration_molar = concentration_pb.value
    elif concentration_pb.units == concentration_pb.MILLIMOLAR:
        concentration_molar = concentration_pb.value * 1e-3
    elif concentration_pb.units == concentration_pb.MICROMOLAR:
        concentration_molar = concentration_pb.value * 1e-6
    else:
        raise ValueError(f'unsupported units: {concentration_pb.units}')
    # Assign moles amount and return.
    moles = volume_liter * concentration_molar
    if moles < 1e-6:
        value = moles * 1e9
        unit = reaction_pb2.Moles.NANOMOLE
    elif moles < 1e-3:
        value = moles * 1e6
        unit = reaction_pb2.Moles.MICROMOLE
    elif moles < 1:
        value = moles * 1e3
        unit = reaction_pb2.Moles.MILLIMOLE
    else:
        value = moles
        unit = reaction_pb2.Moles.MOLE
    solute.moles.value = value
    solute.moles.units = unit
    return [solute] + solvents


def build_data(filename, description):
    """Reads raw data from a file and creates a Data message.

    Args:
        filename: Text filename.
        description: Text description of the data.

    Returns:
        Data message.
    """
    _, extension = os.path.splitext(filename)
    if not extension.startswith('.'):
        raise ValueError(f'cannot deduce the file format for {filename}')
    data = reaction_pb2.Data()
    data.format = extension[1:]
    with open(filename, 'rb') as f:
        data.bytes_value = f.read()
    data.description = description
    return data


def find_submessages(message, submessage_type):
    """Recursively finds all submessages of a specified type.

    Args:
        message: Protocol buffer.
        submessage_type: Protocol buffer type.

    Returns:
        List of messages.

    Raises:
        TypeError: if `submessage_type` is not a protocol buffer type.
    """
    if not issubclass(submessage_type, protobuf.message.Message):
        raise TypeError('submessage_type must be a Protocol Buffer type')
    submessage_name = submessage_type.DESCRIPTOR.full_name
    submessages = []
    for field, value in message.ListFields():
        if field.type != field.TYPE_MESSAGE:
            continue
        if field.message_type.full_name == submessage_name:
            if field.label == field.LABEL_REPEATED:
                submessages.extend(value)
            else:
                submessages.append(value)
        elif field.message_type.GetOptions().map_entry:
            # Map field.
            field_value = field.message_type.fields_by_name['value']
            if field_value.type != field_value.TYPE_MESSAGE:
                continue
            if field_value.message_type.full_name == submessage_name:
                submessages.extend(value.values())
            else:
                for submessage in value.values():
                    submessages.extend(
                        find_submessages(submessage, submessage_type))
        elif field.label == field.LABEL_REPEATED:
            # Standard repeated field.
            for submessage in value:
                submessages.extend(find_submessages(submessage,
                                                    submessage_type))
        else:
            submessages.extend(find_submessages(value, submessage_type))
    return submessages


def smiles_from_compound(compound):
    """Fetches or generates a SMILES identifier for a compound.

    If a SMILES identifier already exists, it is simply returned.

    Args:
        compound: reaction_pb2.Compound message.

    Returns:
        Text SMILES.
    """
    for identifier in compound.identifiers:
        if identifier.type == reaction_pb2.CompoundIdentifier.SMILES:
            return identifier.value
    return Chem.MolToSmiles(mol_from_compound(compound))


def molblock_from_compound(compound):
    """Fetches or generates a MolBlock identifier for a compound.

    Args:
        compound: reaction_pb2.Compound message.

    Returns:
        molblock: MolBlock identifier.
    """
    for identifier in compound.identifiers:
        if identifier.type == reaction_pb2.CompoundIdentifier.MOLBLOCK:
            return identifier.value
    return Chem.MolToMolBlock(mol_from_compound(compound))


# pylint: disable=inconsistent-return-statements
def mol_from_compound(compound, return_identifier=False):
    """Creates an RDKit Mol from a Compound message.

    Args:
        compound: reaction_pb2.Compound message.
        return_identifier: If True, return the CompoundIdentifier used to
            create the Mol.

    Returns:
        mol: RDKit Mol.
        identifier: The identifier that was used to create `mol`. Only returned
            if `return_identifier` is True.

    Raises:
        ValueError: If no structural identifier is available, or if the
            resulting Mol object is invalid.
    """
    for identifier in compound.identifiers:
        if identifier.type in _COMPOUND_IDENTIFIER_LOADERS:
            mol = _COMPOUND_IDENTIFIER_LOADERS[identifier.type](
                identifier.value)
            if not mol:
                raise ValueError(
                    f'invalid structural identifier for Compound: {identifier}')
            if return_identifier:
                return mol, identifier
            return mol
    raise ValueError(f'no valid structural identifier for Compound: {compound}')


# pylint: enable=inconsistent-return-statements


def check_compound_identifiers(compound):
    """Verifies that structural compound identifiers are consistent.

    Args:
        compound: reaction_pb2.Compound message.

    Raises:
        ValueError: If structural identifiers are not consistent or are invalid.
    """
    smiles = set()
    for identifier in compound.identifiers:
        if identifier.type in _COMPOUND_IDENTIFIER_LOADERS:
            mol = _COMPOUND_IDENTIFIER_LOADERS[identifier.type](
                identifier.value)
        else:
            continue
        if not mol:
            raise ValueError(
                f'invalid structural identifier for Compound: {identifier}')
        smiles.add(Chem.MolToSmiles(mol))
    if len(smiles) > 1:
        raise ValueError(f'structural identifiers are inconsistent: {smiles}')


def get_reaction_smiles(message, allow_incomplete=True, validate=True):
    """Fetches or generates a reaction SMILES.

    Args:
        message: reaction_pb2.Reaction message.
        allow_incomplete: Boolean whether to allow "incomplete" reaction SMILES
            that do not include all components (e.g. if a component does not
            have a structural identifier).
        validate: Boolean whether to validate the reaction SMILES with rdkit.
            Only used if allow_incomplete is False.

    Returns:
        Text reaction SMILES.

    Raises:
        ValueError: If the reaction contains errors.
    """
    for identifier in message.identifiers:
        if identifier.type == reaction_pb2.ReactionIdentifier.REACTION_SMILES:
            return identifier.value
    reactants, agents, products = set(), set(), set()
    roles = reaction_pb2.Compound.ReactionRole()
    for key in sorted(message.inputs):
        for compound in message.inputs[key].components:
            try:
                smiles = smiles_from_compound(compound)
            except ValueError as error:
                if allow_incomplete:
                    continue
                raise error
            if compound.reaction_role in [
                    roles.REAGENT, roles.SOLVENT, roles.CATALYST
            ]:
                agents.add(smiles)
            elif compound.reaction_role == roles.INTERNAL_STANDARD:
                continue
            else:
                reactants.add(smiles)
    for outcome in message.outcomes:
        for product in outcome.products:
            try:
                smiles = smiles_from_compound(product.compound)
            except ValueError as error:
                if allow_incomplete:
                    continue
                raise error
            products.add(smiles)
    if not allow_incomplete and (not reactants or not products):
        raise ValueError(
            'reaction must contain at least one reactant and one product')
    if not reactants and not products:
        raise ValueError('reaction contains no valid reactants or products')
    reaction_smiles = '%s>%s>%s' % ('.'.join(sorted(reactants)), '.'.join(
        sorted(agents)), '.'.join(sorted(products)))
    if validate and not allow_incomplete:
        reaction_smiles = validate_reaction_smiles(reaction_smiles)
    return reaction_smiles


def validate_reaction_smiles(reaction_smiles):
    """Validates reaction SMILES.

    Args:
        reaction_smiles: Text reaction SMILES.

    Returns:
        Updated reaction SMILES.

    Raises:
        ValueError: If the reaction contains errors.
    """
    try:
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles,
                                                      useSmiles=True)
        rdChemReactions.SanitizeRxn(reaction)
    except ValueError as error:
        raise ValueError(
            f'reaction contains errors: {reaction_smiles}') from error
    _, num_errors = reaction.Validate()
    if num_errors:
        raise ValueError(f'reaction contains errors: {reaction_smiles}')
    return rdChemReactions.ReactionToSmiles(reaction)


class MessageFormat(enum.Enum):
    """Input/output types for protocol buffer messages."""
    BINARY = '.pb'
    JSON = '.json'
    PBTXT = '.pbtxt'


# pylint: disable=inconsistent-return-statements
def load_message(filename, message_type):
    """Loads a protocol buffer message from a file.

    Args:
        filename: Text filename containing a serialized protocol buffer message.
        message_type: google.protobuf.message.Message subclass.

    Returns:
        Message object.

    Raises:
        ValueError: if the message cannot be parsed, or if `input_format` is not
            supported.
    """
    _, extension = os.path.splitext(filename)
    input_format = MessageFormat(extension)
    if input_format == MessageFormat.BINARY:
        mode = 'rb'
    else:
        mode = 'r'
    with open(filename, mode) as f:
        try:
            if input_format == MessageFormat.JSON:
                return json_format.Parse(f.read(), message_type())
            if input_format == MessageFormat.PBTXT:
                return text_format.Parse(f.read(), message_type())
            if input_format == MessageFormat.BINARY:
                return message_type.FromString(f.read())
        except (json_format.ParseError, protobuf.message.DecodeError,
                text_format.ParseError) as error:
            raise ValueError(f'error parsing {filename}: {error}') from error


# pylint: enable=inconsistent-return-statements


def write_message(message, filename):
    """Writes a protocol buffer message to disk.

    Args:
        message: Protocol buffer message.
        filename: Text output filename.

    Raises:
        ValueError: if `filename` does not have the expected suffix.
    """
    _, extension = os.path.splitext(filename)
    output_format = MessageFormat(extension)
    if output_format == MessageFormat.BINARY:
        mode = 'wb'
    else:
        mode = 'w'
    with open(filename, mode) as f:
        if output_format == MessageFormat.JSON:
            f.write(json_format.MessageToJson(message))
        elif output_format == MessageFormat.PBTXT:
            f.write(text_format.MessageToString(message))
        elif output_format == MessageFormat.BINARY:
            f.write(message.SerializeToString(deterministic=True))


def id_filename(filename):
    """Converts a filename into a relative path for the repository.

    Args:
        filename: Text basename including an ID.

    Returns:
        Text filename relative to the root of the repository.
    """
    basename = os.path.basename(filename)
    prefix, suffix = basename.split('-')
    if not prefix.startswith('ord'):
        raise ValueError(
            'basename does not have the required "ord" prefix: {basename}')
    return flask.safe_join('data', suffix[:2], basename)


def create_message(message_name):
    """Converts a message name into an instantiation of that class, where
    the message belongs to the reaction_pb2 module.

    Args:
        message_name: Text name of a message field. For example, "Reaction" or
            "TemperatureConditions.Measurement".

    Returns:
        Initialized message of the requested type.

    Raises:
        ValueError if the name cannot be resolved.
    """
    message_class = reaction_pb2
    try:
        for name in message_name.split('.'):
            message_class = getattr(message_class, name)
        return message_class()
    except (AttributeError, TypeError) as error:
        raise ValueError(
            f'Cannot resolve message name {message_name}') from error
