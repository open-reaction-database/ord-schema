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

import contextlib
import enum
import functools
import gzip
import os
import pathlib
import posixpath
import re
import warnings
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any, TypeVar, cast

import pandas as pd
from google.protobuf import (
    json_format,
    text_format,
)
from google.protobuf.descriptor import Descriptor
from google.protobuf.message import DecodeError
from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdmolfiles

import ord_schema
from ord_schema import atomic_io, units
from ord_schema.proto import reaction_pb2

# Enhanced stereochemistry is structural; atom labels and coordinates are not.
_CX_ENHANCED_STEREO = rdmolfiles.CXSmilesFields.CX_ENHANCEDSTEREO

_COMPOUND_IDENTIFIER_LOADERS = {
    reaction_pb2.CompoundIdentifier.SMILES: Chem.MolFromSmiles,
    # MolFromSmiles reads CXSMILES, so it is as structural as SMILES.
    reaction_pb2.CompoundIdentifier.CXSMILES: Chem.MolFromSmiles,
    reaction_pb2.CompoundIdentifier.INCHI: Chem.MolFromInchi,
    reaction_pb2.CompoundIdentifier.MOLBLOCK: Chem.MolFromMolBlock,
}
# The identifier types a Mol can be built from, keyed off the loaders above so
# every caller shares this module's single source of truth for what "structural"
# means: a type absent from the loaders yields no structure anywhere downstream,
# so treating it as one only hides the compound from the passes that could fill
# the gap. STRUCTURAL_IDENTIFIER_TYPE_NAMES matches the values stored in
# ord.compound_identifier.type, for callers querying the database.
STRUCTURAL_IDENTIFIER_TYPES = frozenset(_COMPOUND_IDENTIFIER_LOADERS)
STRUCTURAL_IDENTIFIER_TYPE_NAMES = frozenset(
    reaction_pb2.CompoundIdentifier.CompoundIdentifierType.Name(identifier_type)
    for identifier_type in _COMPOUND_IDENTIFIER_LOADERS
)
# Structural identifiers repeat heavily -- shared solvents, reagents and
# catalysts recur across the reactions of a dataset -- and validation parses
# each one twice, once to check that it parses and once to canonicalize it for
# the consistency check. Memoizing the derived SMILES collapses both.
#
# The cached value is a string, never the Mol: an RDKit Mol is mutable, so
# handing the same one to every caller would let any of them corrupt the entry.
# Bounded because a large dataset's distinct set does not saturate (uspto grows
# by ~1.3 new molecules per reaction); this holds strings, so the ceiling is
# tens of MB.
#
# Entries assume RDKit parses and writes the same way for the life of the process.
# A caller that flips a global mid-run (Chem.SetUseLegacyStereoPerception and the
# like) has to call cache_clear() on both functions below.
_IDENTIFIER_CACHE_SIZE = 100_000

MessageType = TypeVar("MessageType")  # Generic for setting return types


@functools.lru_cache(maxsize=_IDENTIFIER_CACHE_SIZE)
def canonical_smiles_for_identifier(identifier_type: int, value: str) -> str | None:
    """Canonicalizes a structural identifier, or returns None if it will not parse.

    Args:
        identifier_type: CompoundIdentifier type enum value.
        value: The identifier value.

    Returns:
        Canonical SMILES, or None if ``value`` is not a valid identifier of that
        type. Returns None for types no Mol can be built from.
    """
    loader = _COMPOUND_IDENTIFIER_LOADERS.get(identifier_type)
    if loader is None:
        return None
    mol = loader(value)
    if mol is None:
        return None
    return canonical_smiles(mol)


@functools.lru_cache(maxsize=_IDENTIFIER_CACHE_SIZE)
def identifier_parses_unsanitized(identifier_type: int, value: str) -> bool:
    """Tests whether an identifier parses at all once sanitization is skipped.

    Separate from :func:`canonical_smiles_for_identifier` so the extra parse
    happens only for the identifiers that failed, which are rare.

    Args:
        identifier_type: CompoundIdentifier type enum value.
        value: The identifier value.

    Returns:
        True if ``value`` parses with ``sanitize=False``.
    """
    loader = _COMPOUND_IDENTIFIER_LOADERS.get(identifier_type)
    if loader is None:
        return False
    return loader(value, sanitize=False) is not None


def _resolve_enum_value(descriptor: Descriptor, field_name: str, key: str) -> int:
    """Returns the numeric value of ``key`` in ``descriptor``'s named enum field."""
    assert descriptor is not None  # Type hint.
    field = descriptor.fields_by_name[field_name]
    assert field.enum_type is not None  # Type hint.
    values = field.enum_type.values_by_name
    try:
        return values[key.upper()].number
    except KeyError as error:
        raise KeyError(f"{key} is not a supported type: {values.keys()}") from error


def build_compound(
    smiles: str | None = None,
    name: str | None = None,
    amount: str | None = None,
    role: str | None = None,
    is_limiting: bool | None = None,
    prep: str | None = None,
    prep_details: str | None = None,
    vendor: str | None = None,
) -> reaction_pb2.Compound:
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
        compound.identifiers.add(value=smiles, type="SMILES")
    if name:
        compound.identifiers.add(value=name, type="NAME")
    if amount:
        if amount.lower() in ("saturated", "catalytic", "titrated"):
            compound.amount.unmeasured.CopyFrom(
                reaction_pb2.UnmeasuredAmount(type=amount.upper())
            )
        else:
            resolver = units.UnitResolver()
            amount_pb = resolver.resolve(amount)
            if isinstance(amount_pb, reaction_pb2.Mass):
                compound.amount.mass.CopyFrom(amount_pb)
            elif isinstance(amount_pb, reaction_pb2.Moles):
                compound.amount.moles.CopyFrom(amount_pb)
            elif isinstance(amount_pb, reaction_pb2.Volume):
                compound.amount.volume.CopyFrom(amount_pb)
            else:
                raise TypeError(f"unsupported units for amount: {amount_pb}")
    if role:
        compound_desc = reaction_pb2.Compound.DESCRIPTOR
        assert compound_desc is not None  # Type hint.
        compound.reaction_role = _resolve_enum_value(
            compound_desc, "reaction_role", role
        )
    if is_limiting is not None:
        if not (is_limiting is True or is_limiting is False):
            raise TypeError(f"is_limiting must be a boolean value: {is_limiting}")
        compound.is_limiting = is_limiting
    if prep:
        prep_desc = reaction_pb2.CompoundPreparation.DESCRIPTOR
        assert prep_desc is not None  # Type hint.
        compound.preparations.add().type = _resolve_enum_value(prep_desc, "type", prep)
        if (
            compound.preparations[0].type == reaction_pb2.CompoundPreparation.CUSTOM
            and not prep_details
        ):
            raise ValueError("prep_details must be provided when CUSTOM prep is used")
    if prep_details:
        if not prep:
            raise ValueError("prep must be provided when prep_details is used")
        compound.preparations[0].details = prep_details
    if vendor:
        compound.source.vendor = vendor
    return compound


def set_solute_moles(
    solute: reaction_pb2.Compound,
    solvents: Sequence[reaction_pb2.Compound],
    concentration: str,
    overwrite: bool = False,
) -> list[reaction_pb2.Compound]:
    """Helps define components for stock solution inputs.

    Handles a single solute and one or more solvent compounds.

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
    if solute.amount.WhichOneof("kind") and not overwrite:
        raise ValueError("solute has defined amount and overwrite is False")

    # Get total solvent volume in liters.
    volume_liter = 0.0
    for solvent in solvents:
        amount = solvent.amount
        if not amount.HasField("volume") or not amount.volume.value:
            raise ValueError("solvent must have defined volume")
        if amount.volume.units == amount.volume.LITER:
            volume_liter += amount.volume.value
        elif amount.volume.units == amount.volume.MILLILITER:
            volume_liter += amount.volume.value * 1e-3
        elif amount.volume.units == amount.volume.MICROLITER:
            volume_liter += amount.volume.value * 1e-6
        elif amount.volume.units == amount.volume.NANOLITER:
            volume_liter += amount.volume.value * 1e-9
        else:
            raise ValueError("solvent units not recognized by set_solute_moles")
    # Get solute concentration in molar.
    resolver = units.UnitResolver(
        unit_synonyms=units.CONCENTRATION_UNIT_SYNONYMS, forbidden_units={}
    )
    concentration_pb = resolver.resolve(concentration)
    assert isinstance(concentration_pb, reaction_pb2.Concentration)  # Type hint.

    solute_moles = units.compute_solute_quantity(
        reaction_pb2.Volume(value=volume_liter, units=reaction_pb2.Volume.LITER),
        concentration_pb,
    )
    solute.amount.MergeFrom(solute_moles)
    return [solute, *list(solvents)]


def build_data(filename: str, description: str) -> reaction_pb2.Data:
    """Reads raw data from a file and creates a Data message.

    Args:
        filename: Text filename.
        description: Text description of the data.

    Returns:
        Data message.
    """
    extension = pathlib.Path(filename).suffix
    if not extension.startswith("."):
        raise ValueError(f"cannot deduce the file format for {filename}")
    data = reaction_pb2.Data()
    data.format = extension[1:]
    with pathlib.Path(filename).open("rb") as f:
        data.bytes_value = f.read()
    data.description = description
    return data


def find_submessages(
    message: ord_schema.Message, submessage_type: type[MessageType]
) -> list[MessageType]:
    """Recursively finds all submessages of a specified type.

    Args:
        message: Protocol buffer.
        submessage_type: Protocol buffer type.

    Returns:
        List of messages.

    Raises:
        TypeError: if `submessage_type` is not a protocol buffer type.
    """
    if not issubclass(submessage_type, ord_schema.Message):
        raise TypeError("submessage_type must be a Protocol Buffer type")
    sub_desc = submessage_type.DESCRIPTOR
    assert sub_desc is not None  # Type hint.
    submessage_name = sub_desc.full_name
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
            map_msg_type = field.message_type
            assert map_msg_type is not None  # Type hint.
            field_value = map_msg_type.fields_by_name["value"]
            if field_value.type != field_value.TYPE_MESSAGE:
                continue
            if field_value.message_type.full_name == submessage_name:
                submessages.extend(value.values())
            else:
                for submessage in value.values():
                    submessages.extend(find_submessages(submessage, submessage_type))
        elif field.label == field.LABEL_REPEATED:
            # Standard repeated field.
            for submessage in value:
                submessages.extend(find_submessages(submessage, submessage_type))
        else:
            submessages.extend(find_submessages(value, submessage_type))
    return submessages


def structural_identifiers(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
) -> Iterator[reaction_pb2.CompoundIdentifier]:
    """Yields a compound's structural identifiers, best first.

    CXSMILES comes before SMILES because it is a superset and RDKit reads either, so
    preferring it keeps enhanced stereochemistry: a plain SMILES would assert one
    configuration where the source recorded a group. Everything else follows in message
    order.

    Yielding rather than picking one lets a caller keep trying, so a malformed value
    does not hide a good identifier behind it. Nothing here parses the values, so a
    yielded identifier is a candidate rather than a promise.

    Args:
        compound: Compound or ProductCompound message.

    Yields:
        Each identifier with a non-empty value whose type a Mol can be built from.
    """
    preferred = (
        reaction_pb2.CompoundIdentifier.CXSMILES,
        reaction_pb2.CompoundIdentifier.SMILES,
    )
    for identifier_type in preferred:
        for identifier in compound.identifiers:
            if identifier.type == identifier_type and identifier.value:
                yield identifier
    for identifier in compound.identifiers:
        if (
            identifier.type not in preferred
            and identifier.type in _COMPOUND_IDENTIFIER_LOADERS
            and identifier.value
        ):
            yield identifier


def smiles_from_compound(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
) -> str | None:
    """Returns canonical SMILES for a compound, preferring its CXSMILES form.

    Identifiers are tried in :func:`structural_identifiers` order, so this and
    :func:`mol_from_compound` agree about which one describes the compound.

    A compound with nothing readable is an ordinary state in ORD -- ligands and reagents
    are routinely recorded by name alone -- so it reads as None rather than raising.

    Args:
        compound: Compound or ProductCompound message.

    Returns:
        Canonical SMILES, carrying an enhanced-stereochemistry block where the structure
        has one (see :func:`canonical_smiles`), or None if the compound records no
        structure any loader can read.
    """
    for identifier in structural_identifiers(compound):
        canonical = canonical_smiles_for_identifier(identifier.type, identifier.value)
        if canonical:
            return canonical
    return None


def reaction_smiles_without_agents(reaction_smiles: str) -> str | None:
    """Returns a reaction SMILES with its agent block removed, or None if unreadable.

    Splitting the string on ``>`` would corrupt a CXSMILES extension, which indexes
    atoms positionally: every index in the block would point at the wrong atom once the
    agents are gone. Rebuilding through RDKit recomputes them.

    The surviving templates are added in canonical order because RDKit numbers the block
    against the order they were added and then writes them sorted; adding them in any
    other order leaves a stereo group marking whichever atom lands at that index.

    Atom mapping survives, being part of the atoms themselves, and so does enhanced
    stereochemistry. Fragment grouping does not, and canonical atom ordering replaces
    whatever the source used.

    Args:
        reaction_smiles: A reaction SMILES or CXSMILES, with or without agents.

    Returns:
        Canonical ``reactants>>products`` CXSMILES. None if RDKit cannot read the input,
        reports errors in it, or it has no reactant or no product once agents are gone:
        a half reaction is not a restatement of the reaction the source recorded.
    """
    try:
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        if reaction is None:
            return None
        without_agents = rdChemReactions.ChemicalReaction()
        for mol in sorted(reaction.GetReactants(), key=Chem.MolToSmiles):
            without_agents.AddReactantTemplate(mol)
        for mol in sorted(reaction.GetProducts(), key=Chem.MolToSmiles):
            without_agents.AddProductTemplate(mol)
        if not without_agents.GetNumReactantTemplates():
            return None
        if not without_agents.GetNumProductTemplates():
            return None
        # Recorded values are held to the same bar as generated ones; being deposited
        # is not evidence of being well formed.
        _validate_reaction(without_agents)
        return rdChemReactions.ReactionToCXSmiles(without_agents) or None
    except (ValueError, RuntimeError):
        # RDKit raises RuntimeError, not ValueError, for a malformed extension block.
        return None


def derived_reaction_smiles(reaction: reaction_pb2.Reaction) -> str | None:
    """Returns the reaction SMILES that derived artifacts store, or None.

    A recorded ``REACTION_CXSMILES`` or ``REACTION_SMILES`` wins over generating one,
    because it carries atom mapping, which generation cannot reconstruct. It is
    normalized, not copied: agents removed and the result canonicalized, so agent
    placement and atom ordering stop deciding whether two reactions look alike. A
    recorded value RDKit cannot read, or that survives as a half reaction, falls through
    to generation.

    Either way the result carries no agents. An empty agent block is idiomatic, and
    excluding agents keeps a reagent, solvent, or catalyst recorded only by name --
    routine for ligands -- from deciding whether the reaction gets a SMILES at all.
    Generation is otherwise strict: every reactant and product must be readable, since a
    SMILES silently missing one describes a different reaction and nothing marks it.

    Args:
        reaction: Reaction message.

    Returns:
        Canonical ``reactants>>products`` SMILES, or None if nothing recorded can be
        read and nothing complete can be generated.
    """
    for identifier_type in (
        reaction_pb2.ReactionIdentifier.REACTION_CXSMILES,
        reaction_pb2.ReactionIdentifier.REACTION_SMILES,
    ):
        for identifier in reaction.identifiers:
            if identifier.type != identifier_type or not identifier.value:
                continue
            without_agents = reaction_smiles_without_agents(identifier.value)
            if without_agents is not None:
                return without_agents
    try:
        return (
            generate_reaction_smiles(
                reaction, allow_incomplete=False, include_agents=False
            )
            or None
        )
    except ValueError:
        return None


def molblock_from_compound(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
) -> str:
    """Fetches or generates a MolBlock identifier for a compound.

    Args:
        compound: reaction_pb2.Compound or ProductCompound message.

    Returns:
        molblock: MolBlock identifier.

    Raises:
        ValueError: if no structural identifiers are defined.
    """
    return get_compound_molblock(compound) or Chem.MolToMolBlock(
        mol_from_compound(compound)
    )


def canonical_smiles(mol: Chem.Mol) -> str:
    """Returns a canonical SMILES, keeping enhanced stereochemistry if present.

    Plain SMILES cannot express AND/OR stereo groups, so writing one would silently
    assert a more specific structure than the source recorded. Only that field is
    emitted; atom labels and coordinates are presentation, not structure.

    Args:
        mol: RDKit Mol.

    Returns:
        Canonical SMILES, with a ``|a:...|``/``|o...|`` block only for molecules that
        have enhanced stereochemistry; others match ``Chem.MolToSmiles`` exactly.
    """
    return Chem.MolToCXSmiles(mol, Chem.SmilesWriteParams(), _CX_ENHANCED_STEREO)


def split_cxsmiles_extension(value: str) -> tuple[str, str | None]:
    """Splits a CXSMILES extension block off a value, if it has one.

    The block follows whitespace and is introduced by ``|``. Whitespace alone does not
    identify one -- records exist whose SMILES carry trailing non-breaking spaces -- so
    anything else is returned intact rather than truncated at the space.

    Args:
        value: A SMILES or CXSMILES string.

    Returns:
        ``(smiles, extension)``, where extension is None when there is no block.
    """
    tokens = value.split(maxsplit=1)
    if len(tokens) == 2 and tokens[1].lstrip().startswith("|"):
        return tokens[0], tokens[1]
    return value, None


def mol_from_compound(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
    return_identifier: bool = False,
) -> Chem.Mol | tuple[Chem.Mol, reaction_pb2.CompoundIdentifier]:
    """Creates an RDKit Mol from a Compound message.

    Identifiers are tried in :func:`structural_identifiers` order, so this and
    :func:`smiles_from_compound` agree about which one describes the compound, and one
    malformed identifier does not mask a readable one recorded alongside it.

    Args:
        compound: reaction_pb2.Compound message.
        return_identifier: If True, return the CompoundIdentifier used to
            create the Mol.

    Returns:
        mol: RDKit Mol.
        identifier: The identifier that was used to create `mol`. Only returned
            if `return_identifier` is True.

    Raises:
        ValueError: If no structural identifier reads as a Mol. Unlike
            :func:`smiles_from_compound`, which reports that as None, callers here want
            a Mol in hand and have nothing to do with its absence.
    """
    for identifier in structural_identifiers(compound):
        mol = _COMPOUND_IDENTIFIER_LOADERS[identifier.type](identifier.value)
        if mol is not None:
            if return_identifier:
                return mol, identifier
            return mol
    raise ValueError(f"no valid structural identifier for Compound: {compound}")


def check_compound_identifiers(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
) -> None:
    """Verifies that structural compound identifiers are consistent.

    Compared through :func:`canonical_smiles`, so identifiers that disagree about
    enhanced stereochemistry are inconsistent: they assert different things about the
    molecule even though they share a skeleton.

    Args:
        compound: reaction_pb2.Compound message.

    Raises:
        ValueError: If structural identifiers are not consistent or are invalid.
    """
    smiles = set()
    for identifier in compound.identifiers:
        if identifier.type not in _COMPOUND_IDENTIFIER_LOADERS:
            continue
        canonical = canonical_smiles_for_identifier(identifier.type, identifier.value)
        if canonical is None:
            raise ValueError(
                f"invalid structural identifier for Compound: {identifier}"
            )
        smiles.add(canonical)
    if len(smiles) > 1:
        # Sorted because a set renders in hash order, which varies between runs:
        # the same dataset would otherwise produce messages that cannot be diffed
        # against a previous run's. Kept as a list literal so the quoting stays
        # and SMILES with trailing whitespace remain visible.
        raise ValueError(f"structural identifiers are inconsistent: {sorted(smiles)}")


def get_reaction_smiles(
    message: reaction_pb2.Reaction,
    generate_if_missing: bool = False,
    allow_incomplete: bool = True,
    allow_unspecified_roles: bool = True,
) -> str | None:
    """Fetches or generates a reaction SMILES.

    A stored REACTION_CXSMILES identifier is returned whole, extension block and all.
    The block is chemistry the source recorded -- enhanced stereochemistry, fragment
    grouping -- and RDKit reads it, so there is nothing to gain by dropping it.
    ``split_cxsmiles_extension`` separates the two parts for a caller that needs the
    bare SMILES, e.g. to look for atom maps.

    Args:
        message: reaction_pb2.Reaction message.
        generate_if_missing: Whether to generate a reaction SMILES from the
            inputs and outputs if one is not defined explicitly.
        allow_incomplete: Boolean whether to allow "incomplete" reaction SMILES
            that do not include all components (e.g. if a component does not
            have a structural identifier).
        allow_unspecified_roles: If True, reactants and products with the
            UNSPECIFIED reaction role will be included when generating a reaction
            SMILES.

    Returns:
        Text reaction SMILES, or None.

    Raises:
        ValueError: If the reaction contains errors.
    """
    types = [
        reaction_pb2.ReactionIdentifier.REACTION_SMILES,
        reaction_pb2.ReactionIdentifier.REACTION_CXSMILES,
    ]
    for identifier in message.identifiers:
        if identifier.type in types:
            return identifier.value
    if not generate_if_missing:
        return None
    return generate_reaction_smiles(
        message,
        allow_incomplete=allow_incomplete,
        allow_unspecified_roles=allow_unspecified_roles,
    )


def generate_reaction_smiles(
    message: reaction_pb2.Reaction,
    *,
    allow_incomplete: bool = True,
    allow_unspecified_roles: bool = True,
    include_agents: bool = True,
) -> str:
    """Builds a reaction SMILES from a Reaction's components, ignoring its identifiers.

    Args:
        message: reaction_pb2.Reaction message.
        allow_incomplete: Whether to omit components with no readable structure rather
            than raising. Only components the result would carry are considered, so an
            unreadable agent is irrelevant when ``include_agents`` is False. Tolerates
            missing components but not a missing half: RDKit calls a reaction with no
            reactants or no products an error however it got that way.
        allow_unspecified_roles: If True, components with the UNSPECIFIED reaction role
            are treated as reactants and products.
        include_agents: Whether to emit the middle ``>agents>`` block. Reagents,
            solvents, and catalysts are dropped when False, which leaves the SMILES
            describing the transformation alone.

    Returns:
        Canonical reaction CXSMILES, so enhanced stereochemistry recorded on a component
        survives into the reaction.

    Raises:
        ValueError: If a component the result would carry has no readable structure and
            ``allow_incomplete`` is False, if there is no reactant or no product, or if
            RDKit reports errors in the assembled reaction.
    """
    reactants, agents, products = set(), set(), set()
    roles = reaction_pb2.ReactionRole
    agent_roles = [roles.REAGENT, roles.SOLVENT, roles.CATALYST]
    reactant_roles = [roles.REACTANT]
    product_roles = [roles.PRODUCT]
    if allow_unspecified_roles:
        reactant_roles.append(roles.UNSPECIFIED)
        product_roles.append(roles.UNSPECIFIED)

    def _add(
        compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
        target: set[str] | None,
    ) -> None:
        """Adds a compound's SMILES to ``target``, or skips it if nothing holds it."""
        # Checked before parsing, so a component the result would not carry cannot fail
        # the whole reaction: a NAME-only ligand is silent when agents are excluded.
        if target is None:
            return
        smiles = smiles_from_compound(compound)
        if smiles is None:
            if allow_incomplete:
                return
            raise ValueError(
                f"no readable structure for a component: {compound.identifiers}"
            )
        target.add(smiles)

    for key in sorted(message.inputs):
        for compound in message.inputs[key].components:
            if compound.reaction_role in agent_roles:
                _add(compound, agents if include_agents else None)
            elif compound.reaction_role in reactant_roles:
                _add(compound, reactants)
    for outcome in message.outcomes:
        for product in outcome.products:
            _add(product, products if product.reaction_role in product_roles else None)

    if not allow_incomplete and (not reactants or not products):
        raise ValueError("reaction must contain at least one reactant and one product")
    if not reactants and not products:
        raise ValueError("reaction contains no valid reactants or products")
    # Assembled from parsed components rather than by joining their SMILES: a component
    # carrying a CXSMILES extension puts a `|...|` block mid-string, which is not a
    # reaction SMILES at all. Building the reaction lets RDKit re-emit one block for the
    # whole reaction, with the atom indices it refers to recomputed.
    built = rdChemReactions.ChemicalReaction()
    for add, group in (
        (built.AddReactantTemplate, reactants),
        (built.AddAgentTemplate, agents),
        (built.AddProductTemplate, products),
    ):
        for smiles in sorted(group):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                # RDKit accepts a null template and then dies with SIGSEGV writing it,
                # which no caller can catch. A ValueError is what callers handle.
                raise ValueError(f"cannot parse component SMILES: {smiles!r}")
            add(mol)
    try:
        reaction_smiles = rdChemReactions.ReactionToCXSmiles(built)
        # Checked on the reaction in hand; sanitizing mutates it, so the string is
        # written first.
        _validate_reaction(built)
    except (RuntimeError, ValueError) as error:
        raise ValueError(f"cannot write reaction SMILES: {error!s}") from error
    return reaction_smiles


def _validate_reaction(reaction: rdChemReactions.ChemicalReaction) -> None:
    """Sanitizes a reaction in place and raises if RDKit reports any error.

    Args:
        reaction: Reaction to check. Sanitized in place, so pass a reaction whose
            SMILES has already been written if the unsanitized form is wanted.

    Raises:
        ValueError: If validation reports errors.
        RuntimeError: If RDKit cannot sanitize the reaction.
    """
    rdChemReactions.SanitizeRxn(reaction)
    _, num_errors = reaction.Validate()
    if num_errors:
        raise ValueError("reaction SMILES contains errors")


def validate_reaction_smiles(reaction_smiles: str) -> None:
    """Validates reaction SMILES.

    Prefer :func:`_validate_reaction` where a parsed reaction is already at hand; this
    is for callers holding only the string, e.g. checking a recorded identifier.

    Args:
        reaction_smiles: Text reaction SMILES.

    Raises:
        ValueError: If the reaction contains errors.
    """
    try:
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        if not reaction:
            raise ValueError("reaction SMILES could not be parsed")
        _validate_reaction(reaction)
    except (RuntimeError, ValueError) as error:
        raise ValueError(
            f"bad reaction SMILES ({error!s}): {reaction_smiles}"
        ) from error


def reaction_from_smiles(reaction_smiles: str) -> reaction_pb2.Reaction:
    """Builds a Reaction by splitting a reaction SMILES."""
    reaction = rdChemReactions.ReactionFromSmarts(reaction_smiles, useSmiles=True)
    rdChemReactions.RemoveMappingNumbersFromReactions(reaction)
    message = reaction_pb2.Reaction()
    message.identifiers.add(value=reaction_smiles, type="REACTION_SMILES")
    reaction_input = message.inputs["from_reaction_smiles"]
    for mol in reaction.GetReactants():
        component = reaction_input.components.add()
        component.identifiers.add(
            value=Chem.MolToSmiles(mol),
            type="SMILES",
            details="Extracted from reaction SMILES",
        )
        component.reaction_role = reaction_pb2.ReactionRole.REACTANT
        component.amount.unmeasured.type = component.amount.unmeasured.CUSTOM
        component.amount.unmeasured.details = "Extracted from reaction SMILES"
    for smiles in reaction_smiles.split(">")[1].split("."):
        if not smiles:
            continue
        component = reaction_input.components.add()
        component.identifiers.add(
            value=smiles, type="SMILES", details="Extracted from reaction SMILES"
        )
        component.reaction_role = reaction_pb2.ReactionRole.REAGENT
        component.amount.unmeasured.type = component.amount.unmeasured.CUSTOM
        component.amount.unmeasured.details = "Extracted from reaction SMILES"
    outcome = message.outcomes.add()
    for mol in reaction.GetProducts():
        component = outcome.products.add()
        component.identifiers.add(
            value=Chem.MolToSmiles(mol),
            type="SMILES",
            details="Extracted from reaction SMILES",
        )
        component.reaction_role = reaction_pb2.ReactionRole.PRODUCT
    return message


def get_product_yield(
    product: reaction_pb2.ProductCompound, as_measurement: bool = False
) -> reaction_pb2.ProductMeasurement | float | None:
    """Returns the value of a product's yield if it is defined.

    If multiple measurements of type YIELD exist, only the first is returned.

    A yield recorded as ``float_value``, ``amount``, or text has no percentage to
    report: reading one as a percentage would assume a scale the source did not state.
    Pass ``as_measurement`` to reach those. ``views`` applies the same rule across every
    outcome and counts what it drops.

    Args:
        product: ProductCompound message.
        as_measurement: Whether to return the full ProductMeasurement that
            corresponds to the yield measurement. Defaults to False.

    Returns:
        The ProductMeasurement when ``as_measurement`` is set. Otherwise the yield as a
        percentage, or None when the product records no yield or records one that is not
        a percentage.
    """
    for measurement in product.measurements:
        if measurement.type == measurement.YIELD:
            if as_measurement:
                return measurement
            # An unset percentage, and an unset value within one, both read as 0.0
            # through protobuf; a fabricated zero yield is worse than a null.
            if not measurement.HasField(
                "percentage"
            ) or not measurement.percentage.HasField("value"):
                return None
            return measurement.percentage.value
    return None


def get_compound_identifier(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
    identifier_type: reaction_pb2.CompoundIdentifier.CompoundIdentifierType,
) -> str | None:
    """Returns the value of a compound identifier if it exists.

    If multiple identifiers of that type exist, only the first is returned.

    Args:
        compound: Compound message.
        identifier_type: The CompoundIdentifier type to retrieve the value of.

    Returns:
        Identifier value or None if the identifier is not defined.
    """
    for identifier in compound.identifiers:
        if identifier.type == identifier_type:
            return identifier.value
    return None


def set_compound_identifier(
    compound: reaction_pb2.Compound,
    identifier_type: reaction_pb2.CompoundIdentifier.CompoundIdentifierType,
    value: str,
) -> reaction_pb2.CompoundIdentifier:
    """Sets the value of a compound identifier if it exists or creates one.

    If multiple identifiers of that type exist, only the first is overwritten.

    Args:
        compound: Compound message.
        identifier_type: The CompoundIdentifier type to retrieve the value of.
        value: The value to set.

    Returns:
        The compound identifier that was modified or created.
    """
    for identifier in compound.identifiers:
        if identifier.type == identifier_type:
            identifier.value = value
            return identifier
    return compound.identifiers.add(type=identifier_type, value=value)


def get_compound_smiles(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
) -> str | None:
    """Returns the value of the compound's SMILES identifier if it exists.

    Args:
        compound: Compound message.

    Returns:
        SMILES string or None if the compound has no SMILES identifier.
    """
    return get_compound_identifier(compound, reaction_pb2.CompoundIdentifier.SMILES)


def set_compound_smiles(
    compound: reaction_pb2.Compound, value: str
) -> reaction_pb2.CompoundIdentifier:
    """Sets the value of the compound's SMILES identifier if it exists or creates one.

    Args:
        compound: Compound message.
        value: The value to set.

    Returns:
        The compound identifier that was modified or created.
    """
    return set_compound_identifier(
        compound, reaction_pb2.CompoundIdentifier.SMILES, value
    )


def is_transition_metal(atom: Chem.Atom) -> bool:
    """Determines if an atom is a transition metal.

    Args:
        atom: The atom in question. Should be of type rdkit.Chem.rdchem.Atom

    Returns:
        Boolean for whether the atom is a transition metal.
    """
    atom_n = atom.GetAtomicNum()
    return 22 <= atom_n <= 29 or 40 <= atom_n <= 47 or 72 <= atom_n <= 79


def has_transition_metal(mol: Chem.Mol) -> bool:
    """Determines if a molecule contains a transition metal.

    Args:
        mol: The molecule in question. Should be of type rdkit.Chem.rdchem.Mol

    Returns:
        Boolean for whether the molecule has a transition metal.
    """
    return any(is_transition_metal(atom) for atom in mol.GetAtoms())


def set_dative_bonds(
    mol: Chem.Mol, from_atoms: tuple[str, ...] = ("N", "P")
) -> Chem.Mol:
    """Converts metal-ligand bonds to dative.

    Replaces some single bonds between metals and atoms with atomic numbers
    in fromAtoms with dative bonds. For all atoms except carbon, the
    replacement is only done if the atom has "too many" bonds. To handle
    metal-carbene complexes, metal-carbon bonds are converted to dative
    if the sum of the explicit and implicit valence of the carbon atom
    does not equal its default valence, 4.

    Args:
        mol: The molecule to be converted.
        from_atoms: tuple of atomic symbols corresponding to atom types that
            should have atom-metal bonds converted to dative. Default is N and P

    Returns:
        The modified molecule.
    """
    p_table = Chem.GetPeriodicTable()
    edit_mol = Chem.RWMol(mol)
    edit_mol.UpdatePropertyCache(strict=False)
    metals = [atom for atom in edit_mol.GetAtoms() if is_transition_metal(atom)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            nbr_atom = nbr.GetSymbol()
            # Handles carbon-bound (e.g., NHC-type or CO) ligands
            # Converts carbon-metal bond to dative if carbon's total valence +
            # formal charge does not equal 4
            if nbr_atom in from_atoms and nbr_atom == "C":
                if nbr.GetFormalCharge() > 0:
                    warnings.warn(
                        f"A positively charged C atom bound to "
                        f"{metal.GetSymbol()} was found in the compound "
                        f"with SMILES {Chem.MolToSmiles(mol)}. If this is "
                        f"a datively bound metal-carbene complex, "
                        f"the positive charge should be removed from "
                        f"the SMILES string before setting dative bonds"
                    )
                if (
                    nbr.GetTotalValence() + nbr.GetFormalCharge()
                    != p_table.GetDefaultValence(nbr_atom)
                    and edit_mol.GetBondBetweenAtoms(
                        nbr.GetIdx(), metal.GetIdx()
                    ).GetBondType()
                    == Chem.BondType.SINGLE
                ):
                    edit_mol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                    edit_mol.AddBond(nbr.GetIdx(), metal.GetIdx(), Chem.BondType.DATIVE)

            # Handles atoms other than carbon (P, N, O, S, etc.)
            # Converts atom-metal bond to dative if bonds to atom
            # excedes its default valence
            elif nbr_atom in from_atoms and nbr_atom != "C":
                if (
                    nbr.GetExplicitValence() > p_table.GetDefaultValence(nbr_atom)
                    and edit_mol.GetBondBetweenAtoms(
                        nbr.GetIdx(), metal.GetIdx()
                    ).GetBondType()
                    == Chem.BondType.SINGLE
                ):
                    edit_mol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                    edit_mol.AddBond(nbr.GetIdx(), metal.GetIdx(), Chem.BondType.DATIVE)

    return edit_mol.GetMol()


def get_compound_name(compound: reaction_pb2.Compound) -> str | None:
    """Returns the value of the compound's NAME identifier if it exists.

    Args:
        compound: Compound message.

    Returns:
        NAME string or None if the compound has no NAME identifier.
    """
    return get_compound_identifier(compound, reaction_pb2.CompoundIdentifier.NAME)


def set_compound_name(
    compound: reaction_pb2.Compound, value: str
) -> reaction_pb2.CompoundIdentifier:
    """Sets the value of the compound's NAME identifier if it exists or creates one.

    Args:
        compound: Compound message.
        value: The value to set.

    Returns:
        The compound identifier that was modified or created.
    """
    return set_compound_identifier(
        compound, reaction_pb2.CompoundIdentifier.NAME, value
    )


def get_compound_molblock(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
) -> str | None:
    """Returns the value of the compound's MOLBLOCK identifier if it exists.

    Args:
        compound: Compound message.

    Returns:
        MOLBLOCK string or None if the compound has no MOLBLOCK identifier.
    """
    return get_compound_identifier(compound, reaction_pb2.CompoundIdentifier.MOLBLOCK)


def set_compound_molblock(
    compound: reaction_pb2.Compound, value: str
) -> reaction_pb2.CompoundIdentifier:
    """Sets the value of the compound's MOLBLOCK identifier if it exists or creates one.

    Args:
        compound: Compound message.
        value: The value to set.

    Returns:
        The compound identifier that was modified or created.
    """
    return set_compound_identifier(
        compound, reaction_pb2.CompoundIdentifier.MOLBLOCK, value
    )


class MessageFormat(enum.Enum):
    """Input/output types for protocol buffer messages.

    BINARY/BINPB and PBTXT/TXTPB pairs use the same wire format; the second of each pair
    is the newer canonical suffix recommended by protobuf.dev.
    """

    BINARY = ".pb"
    BINPB = ".binpb"
    JSON = ".json"
    PBTXT = ".pbtxt"
    TXTPB = ".txtpb"


_BINARY_FORMATS = {MessageFormat.BINARY, MessageFormat.BINPB}
_TEXT_FORMATS = {MessageFormat.PBTXT, MessageFormat.TXTPB}


def _message_format(path: pathlib.Path) -> MessageFormat:
    """Returns the MessageFormat implied by a filename's suffix, ignoring ``.gz``."""
    suffixes = path.suffixes
    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]
    return MessageFormat(suffixes[-1] if suffixes else "")


def load_message(
    filename: str | os.PathLike[str], message_type: type[MessageType]
) -> MessageType:
    """Loads a protocol buffer message from a file.

    Args:
        filename: Filename (``str`` or path-like) containing a serialized message.
        message_type: Message subclass.

    Returns:
        Message object.

    Raises:
        ValueError: if the message cannot be parsed, if ``filename`` is a Parquet
            dataset (use ``datasets.load_dataset``), or if the format is unsupported.
    """
    path = pathlib.Path(filename)
    if path.suffix == ".parquet":
        raise ValueError(
            f"{path} is a Parquet dataset; load_message reads a single serialized "
            "message. Use ord_schema.datasets.load_dataset for the whole Dataset, or "
            "ord_schema.parquet.DatasetView to stream large ones."
        )
    this_open = gzip.open if path.suffix == ".gz" else open
    input_format = _message_format(path)
    mode = "rb" if input_format in _BINARY_FORMATS else "rt"
    with this_open(path, mode) as f:
        try:
            if input_format == MessageFormat.JSON:
                return json_format.Parse(f.read(), message_type())
            if input_format in _TEXT_FORMATS:
                return text_format.Parse(f.read(), message_type())
            if input_format in _BINARY_FORMATS:
                return message_type.FromString(f.read())  # ty: ignore[unresolved-attribute]
            raise ValueError(f"unsupported MessageFormat: {input_format}")
        except (
            json_format.ParseError,
            DecodeError,
            text_format.ParseError,
        ) as error:
            raise ValueError(f"error parsing {path}: {error}") from error


def save_message(message: ord_schema.Message, filename: str | os.PathLike[str]) -> None:
    """Writes a protocol buffer message to disk.

    Args:
        message: Protocol buffer message.
        filename: Output filename (``str`` or path-like).

    Raises:
        ValueError: if `filename` does not have the expected suffix.
    """
    path = pathlib.Path(filename)
    output_format = _message_format(path)
    with (
        atomic_io.atomic_path(path) as tmp_filename,
        _open_for_write(tmp_filename, dest=path) as f,
    ):
        if output_format == MessageFormat.JSON:
            f.write(json_format.MessageToJson(message).encode())
        elif output_format in _TEXT_FORMATS:
            # Protobuf 5+ MessageToBytes() defaults to ASCII; non-ASCII string fields
            # (common in chemistry text) then raise UnicodeEncodeError. MessageToString
            # returns a Unicode str; UTF-8 bytes are the portable wire form for .pbtxt.
            f.write(text_format.MessageToString(message).encode("utf-8"))
        elif output_format in _BINARY_FORMATS:
            f.write(message.SerializeToString(deterministic=True))


@contextlib.contextmanager
def _open_for_write(
    tmp_path: str | os.PathLike[str], *, dest: str | os.PathLike[str]
) -> Iterator[Any]:
    """Opens ``tmp_path`` for binary writing.

    For ``.gz`` destinations, wraps the file in a ``GzipFile`` with a fixed mtime and
    pins the gzip header's filename to the destination's basename so the encoded bytes
    are deterministic regardless of the random temp name used during writing.
    """
    dest = pathlib.Path(dest)
    with pathlib.Path(tmp_path).open("wb") as raw:
        if dest.suffix == ".gz":
            with gzip.GzipFile(
                filename=dest.name, mode="wb", mtime=1, fileobj=raw
            ) as f:
                yield f
        else:
            yield raw


def id_filename(filename: str) -> str:
    """Converts a filename into a relative path for the repository.

    Args:
        filename: Text basename including an ID.

    Returns:
        Text filename relative to the root of the repository.
    """
    basename = pathlib.Path(filename).name
    prefix, suffix = basename.split("-")
    if not prefix.startswith("ord"):
        raise ValueError(
            f'basename does not have the required "ord" prefix: {basename}'
        )
    shard = suffix[:2]
    # Reject anything that could let the shard escape the "data/" root
    # (e.g. "..", "/x").
    if not shard.isalnum():
        raise ValueError(f"basename shard must be alphanumeric: {basename}")
    result = posixpath.join("data", shard, basename)
    # Defense in depth: basename carries no separator and shard is alphanumeric, so
    # this is redundant -- but it is what pins the result inside the "data/" root.
    if posixpath.normpath(result) != result or not result.startswith("data/"):
        raise ValueError(f"unsafe path from basename: {basename}")
    return result


def create_message(message_name: str) -> ord_schema.Message:
    """Converts a message name into an instantiation of that class.

    The message belongs to the reaction_pb2 module.

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
        for name in message_name.split("."):
            message_class = getattr(message_class, name)
        ctor = cast(type[ord_schema.Message], message_class)
        return ctor()
    except (AttributeError, TypeError) as error:
        raise ValueError(f"Cannot resolve message name {message_name}") from error


def messages_to_dataframe(
    messages: Iterable[ord_schema.Message], drop_constant_columns: bool = False
) -> pd.DataFrame:
    """Converts a list of protos to a pandas DataFrame.

    Args:
        messages: List of protos.
        drop_constant_columns: Whether to drop columns that have the same value
            for all rows.

    Returns:
        DataFrame.
    """
    rows = [message_to_row(message) for message in messages]
    df = pd.DataFrame(rows)
    if drop_constant_columns:
        drop = [column for column in df.columns if len(df[column].unique()) == 1]
        for column in drop:
            del df[column]
    return df


def message_to_row(
    message: ord_schema.Message, trace: tuple[str, ...] | None = None
) -> dict[str, ord_schema.ScalarType]:
    """Converts a proto into a flat dictionary mapping fields to values.

    The keys indicate any nesting; for instance a proto that looks like this:

    value: {
      subvalue: 5
    }

    will show up as {'value.subvalue': 5} in the dict.

    Args:
        message: Proto to convert.
        trace: Tuple of strings; the trace of nested field names.

    Returns:
        Dict mapping string field names to scalar value types.
    """
    if trace is None:
        trace = ()
    row = {}
    for field, value in message.ListFields():
        if field.label == field.LABEL_REPEATED:
            if (
                field.type == field.TYPE_MESSAGE
                and field.message_type.GetOptions().map_entry
            ):
                value_field = field.message_type.fields_by_name["value"]
                for key, subvalue in value.items():
                    this_trace = (*trace, f'{field.name}["{key}"]')
                    safe_update(
                        row,
                        _message_to_row(
                            field=value_field, value=subvalue, trace=this_trace
                        ),
                    )
            else:
                for i, subvalue in enumerate(value):
                    this_trace = (*trace, f"{field.name}[{i}]")
                    safe_update(
                        row,
                        _message_to_row(field=field, value=subvalue, trace=this_trace),
                    )
        else:
            this_trace = (*trace, field.name)
            safe_update(
                row, _message_to_row(field=field, value=value, trace=this_trace)
            )
    return row


def safe_update(target: dict, update: Mapping) -> None:
    """Checks that `update` will not clobber any keys in `target`."""
    for key in update:
        if key in target:
            raise KeyError(f"key already exists: {key}")
    target.update(update)


def _message_to_row(
    field: ord_schema.FieldDescriptor,
    value: ord_schema.Message | ord_schema.ScalarType,
    trace: tuple[str, ...],
) -> dict[str, ord_schema.ScalarType]:
    """Recursively creates a dict for a single value.

    Args:
        field: FieldDescriptor for this field.
        value: Value for this field. Should be either a proto or a scalar.
        trace: Tuple of strings; the trace of nested field names.

    Returns:
        Dict mapping string field names to scalar values.
    """
    if field.type == field.TYPE_MESSAGE:
        assert isinstance(value, ord_schema.Message)  # Type hint.
        return message_to_row(message=value, trace=trace)
    if field.type == field.TYPE_ENUM:
        enum_value = field.enum_type.values_by_number[value].name
        return {".".join(trace): enum_value}
    assert isinstance(value, (str, bytes, float, int, bool))  # Type hint.
    return {".".join(trace): value}


def parse_doi(doi: str) -> str:
    """Parses a DOI from e.g. a URL.

    Args:
        doi: DOI string.

    Returns:
        The (possibly trimmed) DOI.

    Raises:
        ValueError: if the DOI cannot be parsed.
    """
    # See https://www.doi.org/doi_handbook/2_Numbering.html#2.2.
    match = re.search(r"(10\.[\d.]+/[a-zA-Z\d.-]+)", doi)
    if not match:
        raise ValueError(f"could not parse DOI: {doi}")
    return match.group(1)
