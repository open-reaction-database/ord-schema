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

import http.client
import re
import urllib.error
import urllib.parse
import urllib.request

from rdkit import Chem

import ord_schema
from ord_schema import message_helpers
from ord_schema.logging import get_logger
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)

# Applied to each socket operation. Without it a resolver that accepts the connection
# and then goes quiet blocks its caller forever, which costs a worker per hung lookup
# when resolution is dispatched to a thread pool.
_TIMEOUT_SECONDS = 10

# Resolvers answer with a single SMILES, so a response this large means the upstream is
# misbehaving (an error page, a redirect loop) and is not worth buffering.
_MAX_RESPONSE_BYTES = 1 << 20


class _ResolverError(Exception):
    """A resolver did not produce a structure.

    Resolvers signal failure with this alone, so the fallback chain does not have to
    enumerate the transport library's exception types on their behalf.
    """


def canonicalize_smiles(smiles: str) -> str:
    """Canonicalizes a SMILES string, raising if it will not parse.

    A thin wrapper over :func:`message_helpers.canonical_smiles`, so a resolved
    structure is written the same way as one derived from a Compound. Callers holding a
    Mol should use that directly; this exists for the string-in, string-out case.

    Args:
        smiles: SMILES string, which may carry a CXSMILES extension block.

    Returns:
        Canonical SMILES, carrying a block where the structure has enhanced
        stereochemistry or coordinate bonds.

    Raises:
        ValueError: If the SMILES cannot be parsed by RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    return message_helpers.canonical_smiles(mol)


def resolve_name(value_type: str, value: str) -> tuple[str, str]:
    """Resolves compound identifiers to SMILES via multiple APIs.

    Resolvers are tried in order until one answers. Any of them failing falls through
    to the next, so only an exhausted chain is a failure.

    Args:
        value_type: The kind of identifier being resolved, e.g. "name".
        value: The identifier to resolve.

    Returns:
        A tuple of SMILES and the name of the resolver that produced it.

    Raises:
        ValueError: If no resolver returns a structure.
    """
    for resolver, resolver_func in _NAME_RESOLVERS.items():
        try:
            smiles = resolver_func(value_type, value)
            # A resolver that answers with something RDKit cannot read is no more use
            # than one that answers with nothing, and worse to keep: written out as a
            # SMILES identifier it counts as structural, which masks the compound from
            # every later resolution pass, so a working resolver never gets to correct
            # it. An empty body arrives here as "" and is the same problem. Both fall
            # through to the next resolver.
            if smiles and Chem.MolFromSmiles(smiles) is not None:
                return smiles, resolver
            logger.info(
                f"{resolver} returned no usable structure for {value_type} {value}: "
                f"{smiles!r}"
            )
        except _ResolverError as error:
            logger.info(f"{resolver} could not resolve {value_type} {value}: {error}")
    raise ValueError(f"Could not resolve {value_type} {value} to SMILES")


def resolve_names(message: ord_schema.Message) -> bool:
    """Attempts to resolve compound NAME identifiers to SMILES.

    When a NAME identifier is resolved, a SMILES identifier is added to the list
    of identifiers for that compound. The first success ends work on that
    Compound; any remaining NAME identifiers on it are left unresolved.

    A compound already carrying an identifier a Mol can be built from is skipped.
    Coordinates alone do not count, since nothing in this library reads them into a
    structure, so a compound recorded as XYZ plus a name is worth resolving.

    Args:
        message: Protocol buffer tree containing Compound submessages (e.g. Reaction
            or ReactionInput).

    Returns:
        Boolean whether `message` was modified.
    """
    modified = False
    compounds = message_helpers.find_submessages(message, reaction_pb2.Compound)
    for compound in compounds:
        if any(
            identifier.type in message_helpers.STRUCTURAL_IDENTIFIER_TYPES
            for identifier in compound.identifiers
        ):
            continue  # Compound already has a structural identifier.
        for identifier in compound.identifiers:
            if identifier.type == identifier.NAME:
                try:
                    smiles, resolver = resolve_name("name", identifier.value)
                    new_identifier = compound.identifiers.add()
                    new_identifier.type = new_identifier.SMILES
                    new_identifier.value = smiles
                    new_identifier.details = f"NAME resolved by the {resolver}"
                    modified = True
                    break
                except ValueError:
                    pass
    return modified


def _fetch_text(url: str) -> str:
    """Fetches ``url`` as stripped text, bounded by a timeout and a size cap.

    Every way a request can fail arrives as :class:`_ResolverError`, so the fallback
    chain treats an unreachable service, an HTTP status, and an unreadable body alike:
    all of them mean this resolver has no answer, and the next one should be asked.
    ``OSError`` covers the transport failures, ``URLError`` and ``TimeoutError`` among
    them; ``HTTPException`` covers a response that is malformed or ends early.

    Args:
        url: The URL to fetch.

    Returns:
        The decoded response body with surrounding whitespace removed. Undecodable
        bytes become replacement characters rather than failing the lookup, since the
        answer is a SMILES and will not parse as a structure either way.

    Raises:
        _ResolverError: If the request fails, or the body exceeds
            ``_MAX_RESPONSE_BYTES``.
    """
    try:
        # S310: the scheme is a literal in each caller's f-string, so no caller-supplied
        # value can reach urlopen as a scheme.
        with urllib.request.urlopen(url, timeout=_TIMEOUT_SECONDS) as response:  # noqa: S310
            body = response.read(_MAX_RESPONSE_BYTES + 1)
            # Passing a size means EOF ends the read quietly, where the no-argument
            # form would raise: a connection dropped mid-body leaves a short answer
            # that is not detectably broken, since a truncated SMILES usually parses
            # as a different molecule. ``length`` is what a declared body still owes,
            # and is None for a chunked one, which raises IncompleteRead instead.
            if response.length:
                raise _ResolverError(f"{url} ended early")
    except (OSError, http.client.HTTPException) as error:
        raise _ResolverError(f"{url} could not be read: {error}") from error
    if len(body) > _MAX_RESPONSE_BYTES:
        raise _ResolverError(f"{url} sent too much data")
    return body.decode(errors="replace").strip()


def _pubchem_resolve(value_type: str, value: str) -> str:
    """Resolves compound identifiers to SMILES via the PubChem REST API.

    A 503 means genuine service load or an IP-level block; it is not retried, and
    ``resolve_name`` falls through to the next resolver.
    """
    # Both segments are quoted with safe="" so that a name containing slashes or dot
    # segments cannot rewrite the request path and resolve a different compound.
    return _fetch_text(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"{urllib.parse.quote(value_type, safe='')}/"
        f"{urllib.parse.quote(value, safe='')}/property/IsomericSMILES/txt"
    )


def _cactus_resolve(value_type: str, value: str) -> str:
    """Resolves compound identifiers to SMILES via the NCI/CADD CIR web service.

    CIR resolves trade names, trivial names, abbreviations, and CAS numbers, so it
    serves as a fast fallback when PubChem is unavailable. Unknown names come back as
    HTTP 500 (not 404), which is reported like any other failed lookup and falls
    through to the next resolver.
    """
    del value_type  # CIR infers the identifier type from the string.
    body = _fetch_text(
        f"https://cactus.nci.nih.gov/chemical/structure/"
        f"{urllib.parse.quote(value, safe='')}/smiles"
    )
    # CIR can return several representations newline-separated; take the first.
    return body.split("\n", 1)[0]


def _opsin_resolve(value_type: str, value: str) -> str:
    """Resolves systematic IUPAC names to SMILES via the OPSIN web service.

    OPSIN is a parser for systematic IUPAC nomenclature, not a lookup service: it will
    refuse trade names, trivial names, and abbreviations (e.g. "aspirin", "THF") with an
    HTTP 404. Treat it strictly as a complement to PubChem.
    """
    del value_type  # OPSIN only supports names.
    return _fetch_text(
        f"https://www.ebi.ac.uk/opsin/ws/{urllib.parse.quote(value, safe='')}.smi"
    )


def resolve_input(input_string: str) -> reaction_pb2.ReactionInput:
    """Resolves a text-based description of an input.

    Supported formats:
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
        component.CopyFrom(
            message_helpers.build_compound(
                name=component_name.strip(), amount=amount_string
            )
        )
        resolve_names(reaction_input)
        return reaction_input
    pattern = re.compile(r"(\d+.?\d*)\s?(\w+)\s(.+)\sin\s(.+)")
    match = pattern.fullmatch(description.strip())
    if not match:
        raise ValueError("String does not match template!")
    conc_value, conc_units, solute_name, solvent_name = match.groups()
    assert solute_name is not None  # Type hint.
    assert solvent_name is not None  # Type hint.
    solute = reaction_input.components.add()
    solvent = reaction_input.components.add()
    solute.CopyFrom(message_helpers.build_compound(name=solute_name.strip()))
    solvent.CopyFrom(
        message_helpers.build_compound(name=solvent_name.strip(), amount=amount_string)
    )
    if solvent.amount.WhichOneof("kind") != "volume":
        raise ValueError("Total amount of solution must be a volume!")
    solvent.amount.volume_includes_solutes = True
    message_helpers.set_solute_moles(solute, [solvent], f"{conc_value} {conc_units}")
    resolve_names(reaction_input)
    return reaction_input


# Standard name resolvers, tried in order until one returns a SMILES.
#
# PubChem handles the bulk of inputs (trade names, trivial names, abbreviations,
# CAS numbers). NCI/CADD CIR (cactus) covers the same broad space and is a fast
# fallback when PubChem is rate-limited. OPSIN handles systematic IUPAC names
# that lookup services may not index.
#
# eMolecules' public /lookup?q= endpoint is unusable (it always replies "__END__";
# their current API requires authentication), and ChemSpider requires a registered
# RSC API key, so neither is included.
_NAME_RESOLVERS = {
    "PubChem API": _pubchem_resolve,
    "NCI/CADD CIR": _cactus_resolve,
    "OPSIN": _opsin_resolve,
}
