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

import re
import time
import urllib.error
import urllib.parse
import urllib.request

from rdkit import Chem

import ord_schema
from ord_schema import message_helpers
from ord_schema.logging import get_logger
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)

# Hard wall-clock ceiling on one resolver request, covering connect and body together.
# An upstream that accepts the connection and then stalls -- or drip-feeds a byte at a
# time to stay under a per-socket timeout -- would otherwise block the calling thread
# for good; a caller dispatching to a worker pool loses a worker per hung lookup.
_REQUEST_TIMEOUT_SECONDS = 10

# How long the socket may go silent, during connect or mid-body. Must stay below
# _REQUEST_TIMEOUT_SECONDS: the read loop only begins a read it can finish inside the
# request budget, so a stall tolerance at or above it leaves no room to read at all.
_STALL_TIMEOUT_SECONDS = 5

# Resolvers answer with a single SMILES, so a response this large means the upstream is
# misbehaving (an error page, a redirect loop) and is not worth buffering.
_MAX_RESPONSE_BYTES = 1 << 20

_COMPOUND_STRUCTURAL_IDENTIFIERS = [
    reaction_pb2.CompoundIdentifier.SMILES,
    reaction_pb2.CompoundIdentifier.INCHI,
    reaction_pb2.CompoundIdentifier.MOLBLOCK,
    reaction_pb2.CompoundIdentifier.CXSMILES,
    reaction_pb2.CompoundIdentifier.XYZ,
]


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


def resolve_name(value_type: str, value: str) -> tuple[str, str]:
    """Resolves compound identifiers to SMILES via multiple APIs.

    Any single resolver being unreachable, erroring, or timing out falls through to the
    next one; only an exhausted chain is a failure. ``URLError`` covers ``HTTPError``
    and connect-time timeouts, while a stall part-way through a response body surfaces
    as a bare ``TimeoutError``.
    """
    for resolver, resolver_func in _NAME_RESOLVERS.items():
        try:
            smiles = resolver_func(value_type, value)
            if smiles is not None:
                return smiles, resolver
        except (urllib.error.URLError, TimeoutError) as error:
            logger.info(f"{resolver} could not resolve {value_type} {value}: {error}")
    raise ValueError(f"Could not resolve {value_type} {value} to SMILES")


def resolve_names(message: ord_schema.Message) -> bool:
    """Attempts to resolve compound NAME identifiers to SMILES.

    When a NAME identifier is resolved, a SMILES identifier is added to the list
    of identifiers for that compound. Note that this function moves on to the
    next Compound after the first successful name resolution.

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
            identifier.type in _COMPOUND_STRUCTURAL_IDENTIFIERS
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
    """Fetches ``url`` as stripped text under a wall-clock deadline and a size cap.

    The socket timeout bounds a single read, not the request: an upstream sending a
    little data just often enough to reset it would stream forever. Two rules together
    make the ceiling hard. The body is consumed with ``read1``, which returns as soon as
    any bytes arrive, so control comes back between chunks rather than being held inside
    one read that never returns; and a read is only begun when the stall tolerance still
    fits inside the deadline, so the last read cannot run past it. Total time is
    therefore at most ``_REQUEST_TIMEOUT_SECONDS``, not that plus a trailing read.

    Args:
        url: The URL to fetch.

    Returns:
        The decoded response body with surrounding whitespace removed.

    Raises:
        TimeoutError: If the body cannot be finished within the request budget.
        urllib.error.URLError: If the response exceeds ``_MAX_RESPONSE_BYTES``.
    """
    deadline = time.monotonic() + _REQUEST_TIMEOUT_SECONDS
    chunks: list[bytes] = []
    size = 0
    # S310: every caller builds an https literal here, interpolating only a quoted path
    # segment, so no caller-supplied scheme can reach urlopen.
    with urllib.request.urlopen(  # noqa: S310
        url, timeout=_STALL_TIMEOUT_SECONDS
    ) as response:
        while True:
            if time.monotonic() + _STALL_TIMEOUT_SECONDS > deadline:
                raise TimeoutError(f"Timed out reading a response from {url}")
            chunk = response.read1(_MAX_RESPONSE_BYTES)
            if not chunk:
                break
            size += len(chunk)
            if size > _MAX_RESPONSE_BYTES:
                raise urllib.error.URLError(f"Response from {url} is too large")
            chunks.append(chunk)
    return b"".join(chunks).decode().strip()


def _pubchem_resolve(value_type: str, value: str) -> str:
    """Resolves compound identifiers to SMILES via the PubChem REST API.

    A 503 (genuine service load or an IP-level block) propagates to ``resolve_name``,
    which falls through to the next resolver instead of retrying.
    """
    return _fetch_text(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{value_type}/"
        f"{urllib.parse.quote(value)}/property/IsomericSMILES/txt"
    )


def _cactus_resolve(value_type: str, value: str) -> str:
    """Resolves compound identifiers to SMILES via the NCI/CADD CIR web service.

    CIR resolves trade names, trivial names, abbreviations, and CAS numbers, so it
    serves as a fast fallback when PubChem is unavailable. Unknown names come back as
    HTTP 500 (not 404), which ``resolve_name`` treats like any other HTTPError and falls
    through.
    """
    del value_type  # CIR infers the identifier type from the string.
    body = _fetch_text(
        f"https://cactus.nci.nih.gov/chemical/structure/{urllib.parse.quote(value)}/smiles"
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
        f"https://www.ebi.ac.uk/opsin/ws/{urllib.parse.quote(value)}.smi"
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
# their current API requires authentication), and ChemSpider now requires a
# registered RSC API key, so neither is included.
_NAME_RESOLVERS = {
    "PubChem API": _pubchem_resolve,
    "NCI/CADD CIR": _cactus_resolve,
    "OPSIN": _opsin_resolve,
}
