"""Helper functions for constructing Protocol Buffer messages."""

import hashlib
import os
import sys

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

# Prefix for filenames that store reaction_pb2.Data values.
DATA_PREFIX = 'ord_data'


# pylint: disable=too-many-arguments
# pylint: disable=too-many-branches


def build_compound(smiles=None, name=None, amount=None, role=None,
                   is_limiting=None, prep=None, prep_details=None, vendor=None):
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
        except KeyError:
            raise KeyError(
                f'{role} is not a supported type: {values_dict.keys()}')
    if is_limiting is not None:
        compound.is_limiting = is_limiting
    if prep:
        field = reaction_pb2.CompoundPreparation.DESCRIPTOR.fields_by_name[
            'type']
        values_dict = field.enum_type.values_by_name
        try:
            compound.preparation.type = values_dict[prep.upper()].number
        except KeyError:
            raise KeyError(
                f'{prep} is not a supported type: {values_dict.keys()}')
        if (compound.preparation.type == reaction_pb2.CompoundPreparation.CUSTOM
                and not prep_details):
            raise ValueError(
                'prep_details must be provided when CUSTOM prep is used')
    if prep_details:
        if not prep:
            raise ValueError('prep must be provided when prep_details is used')
        compound.preparation.details = prep_details
    if vendor:
        compound.vendor_source = vendor
    return compound


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
        elif field.label == field.LABEL_REPEATED:
            # Standard repeated field.
            for submessage in value:
                submessages.extend(
                    find_submessages(submessage, submessage_type))
        else:
            submessages.extend(find_submessages(value, submessage_type))
    return submessages


def write_data(message, dirname, min_size=0.0, max_size=1.0):
    """Writes a Data value to a file.

    If a value is a URL or is smaller than `min_size`, it is left unchanged.

    Args:
        message: Data message.
        dirname: Text output directory.
        min_size: Float minimum size of data before it will be written (in MB).
        max_size: Float maximum size of data to write (in MB).

    Returns:
        Text filename containing the written data, or None if the value is a
        URL or is smaller than `min_size`.

    Raises:
        ValueError: if there is no value defined in `message` or if the value is
            larger than `max_size`.
    """
    kind = message.WhichOneof('kind')
    if kind == 'value':
        value = message.value.encode()  # Convert to bytes.
    elif kind == 'bytes_value':
        value = message.bytes_value
    elif kind == 'url':
        return None
    else:
        raise ValueError('no value to write')
    value_size = sys.getsizeof(value) / 1e6
    if value_size < min_size:
        return None
    if value_size > max_size:
        raise ValueError(
            f'value is larger than max_size ({value_size} vs {max_size}')
    value_hash = hashlib.sha256(value).hexdigest()
    suffix = message.format or 'txt'
    basename = f'{DATA_PREFIX}-{value_hash}.{suffix}'
    filename = os.path.join(dirname, basename)
    with open(filename, 'wb') as f:
        f.write(value)
    return filename


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


def load_message(filename, message_type, input_format):
    """Loads a Reaction proto from a file.

    Args:
        filename: Text filename containing a serialized protocol buffer message.
        message_type: google.protobuf.message.Message subclass.
        input_format: Text input format. Supported options are
            ['binary', 'json', 'pbtxt'].

    Returns:
        Message object.

    Raises:
        ValueError: if the message cannot be parsed, or if `input_format` is not
            supported.
    """
    if input_format == 'binary':
        mode = 'rb'
    else:
        mode = 'r'
    with open(filename, mode) as f:
        if input_format == 'json':
            return json_format.Parse(f.read(), message_type())
        if input_format == 'pbtxt':
            return text_format.Parse(f.read(), message_type())
        if input_format == 'binary':
            message = message_type()
            message.ParseFromString(f.read())
            return message
    raise ValueError(f'unsupported input_format: {input_format}')
