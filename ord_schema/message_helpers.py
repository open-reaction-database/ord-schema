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

from google import protobuf
from google.protobuf import json_format
from google.protobuf import text_format

from ord_schema import units
from ord_schema.proto import reaction_pb2

try:
    from rdkit import Chem
    from rdkit import __version__ as RDKIT_VERSION
except ImportError:
    Chem = None
    RDKIT_VERSION = None


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
        field = reaction_pb2.Compound.DESCRIPTOR.fields_by_name[
            'reaction_role']
        values_dict = field.enum_type.values_by_name
        try:
            compound.reaction_role = values_dict[role.upper()].number
        except KeyError:
            raise KeyError(
                f'{role} is not a supported type: {values_dict.keys()}')
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
        except KeyError:
            raise KeyError(
                f'{prep} is not a supported type: {values_dict.keys()}')
        if (compound.preparations[0].type
                == reaction_pb2.CompoundPreparation.CUSTOM
                and not prep_details):
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
        else:
            raise ValueError(
                'solvent units not recognized by set_solute_moles')
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
                submessages.extend(
                    find_submessages(submessage, submessage_type))
        else:
            submessages.extend(find_submessages(value, submessage_type))
    return submessages


def get_compound_mol(compound):
    """Retrieves an RDKit molecule object from a Compound.

    Args:
        compound: reaction_pb2.Compound message.

    Returns:
        RDKit molecule.
    """
    if not Chem:
        raise ImportError('Missing RDKit dependency')
    for identifier in compound.identifiers:
        if identifier.type == reaction_pb2.CompoundIdentifier.RDKIT_BINARY:
            return Chem.Mol(identifier.bytes_value)
    return None


def get_compound_smiles(compound):
    """Retrieves a SMILES string identifier from a Compound.

    Args:
        compound: reaction_pb2.Compound message.

    Returns:
        The SMILES string identifier for the compound or None.
    """
    for identifier in compound.identifiers:
        if identifier.type == reaction_pb2.CompoundIdentifier.SMILES:
            return identifier.value
    return None


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
            raise ValueError(f'error parsing {filename}: {error}')


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
            f.write(message.SerializeToString())


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
    return os.path.join('data', suffix[:2], basename)
