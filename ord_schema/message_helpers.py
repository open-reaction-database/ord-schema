"""Helper functions for constructing Protocol Buffer messages."""

import os

from ord_schema import units
from ord_schema.proto import reaction_pb2


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


def build_binary_data(filename, description=None):
    """Reads raw data from a file and creates a BinaryData message.

    Args:
        filename: Text filename.
        description: Text description of the data.

    Returns:
        BinaryData message.
    """
    _, extension = os.path.splitext(filename)
    if not extension.startswith('.'):
        raise ValueError(f'cannot deduce the file format for {filename}')
    data = reaction_pb2.BinaryData()
    data.format = extension[1:]
    with open(filename, 'rb') as f:
        data.value = f.read()
    if description:
        data.description = description
    return data
