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
# A caller dispatching these to a worker pool loses a worker for as long as a hung
# lookup runs, so the ceiling is what stops an upstream outage from draining the pool.
_REQUEST_TIMEOUT_SECONDS = 10

# How long the socket may go silent, during connect or mid-body. Must stay below
# _REQUEST_TIMEOUT_SECONDS; see _fetch_text for how the two combine.
_STALL_TIMEOUT_SECONDS = 5

# Resolvers answer with a single SMILES, so a response this large means the upstream is
# misbehaving (an error page, a redirect loop) and is not worth buffering.
_MAX_RESPONSE_BYTES = 1 << 20

# Read-size hint, kept well below _MAX_RESPONSE_BYTES: read1 clamps the request to the
# framing only when the length is known, so a close-delimited response would otherwise
# allocate the whole cap to receive a few dozen bytes.
_READ_CHUNK_BYTES = 1 << 16

# Statuses that say the service, not the identifier, is the problem: PubChem answers 503
# both under load and when it has blocked the caller's IP. Repeating the request for
# every remaining compound would deepen a rate limit rather than resolve anything.
_BACK_OFF_STATUSES = frozenset({429, 503})


class _ResolverError(Exception):
    """A resolver did not produce a structure.

    Resolvers signal failure with this hierarchy alone, so the fallback chain does not
    have to enumerate the transport library's exception types on their behalf.
    """


class _ResolverUnavailableError(_ResolverError):
    """A resolver could not be reached at all.

    Separate from the base class because no response arrived, which distinguishes the
    failures that cost a full timeout from an HTTP status that comes back promptly.
    A resolver that is unreachable for one lookup is unreachable for the next, so
    :func:`resolve_names` stops asking it for the rest of a traversal.
    """


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


def _resolve_name(
    value_type: str, value: str, unavailable: set[str]
) -> tuple[str, str]:
    """Resolves one identifier, skipping and recording resolvers that are unreachable.

    Args:
        value_type: The kind of identifier being resolved, e.g. "name".
        value: The identifier to resolve.
        unavailable: Names of resolvers already found unreachable, added to in place.
            Callers resolving a single identifier pass an empty set and discard it.

    Returns:
        A tuple of SMILES and the name of the resolver that produced it.

    Raises:
        ValueError: If every resolver is skipped or fails.
    """
    for resolver, resolver_func in _NAME_RESOLVERS.items():
        if resolver in unavailable:
            continue
        try:
            smiles = resolver_func(value_type, value)
            if smiles is not None:
                return smiles, resolver
        except _ResolverError as error:
            if isinstance(error, _ResolverUnavailableError):
                unavailable.add(resolver)
            logger.info(f"{resolver} could not resolve {value_type} {value}: {error}")
    raise ValueError(f"Could not resolve {value_type} {value} to SMILES")


def resolve_name(value_type: str, value: str) -> tuple[str, str]:
    """Resolves compound identifiers to SMILES via multiple APIs.

    Any single resolver being unreachable, erroring, or timing out falls through to the
    next one, so only an exhausted chain is a failure. Each resolver is bounded
    separately, which puts the worst case at the length of the chain times
    ``_REQUEST_TIMEOUT_SECONDS``.
    """
    return _resolve_name(value_type, value, set())


def resolve_names(message: ord_schema.Message) -> bool:
    """Attempts to resolve compound NAME identifiers to SMILES.

    When a NAME identifier is resolved, a SMILES identifier is added to the list
    of identifiers for that compound. Note that this function moves on to the
    next Compound after the first successful name resolution.

    A resolver found unreachable is not tried again for the rest of this message.
    Without that, an outage costs the full request budget once per compound, so a
    message carrying dozens of names would spend minutes re-confirming the same
    dead service. The cost is that a resolver which recovers mid-traversal stays
    unused until the next call.

    Args:
        message: Protocol buffer tree containing Compound submessages (e.g. Reaction
            or ReactionInput).

    Returns:
        Boolean whether `message` was modified.
    """
    modified = False
    unavailable: set[str] = set()
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
                    smiles, resolver = _resolve_name(
                        "name", identifier.value, unavailable
                    )
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
    fits inside the deadline, so the last read cannot run past it.

    Args:
        url: The URL to fetch.

    Returns:
        The decoded response body with surrounding whitespace removed.

    Raises:
        _ResolverError: If the resolver answers with an HTTP error status other than a
            back-off. The service is up and the status describes this identifier, so it
            is asked again for the next one.
        _ResolverUnavailableError: If no response arrives, the status asks us to back
            off, the body arrives incomplete, cannot be finished within the request
            budget, or exceeds ``_MAX_RESPONSE_BYTES``.
    """
    deadline = time.monotonic() + _REQUEST_TIMEOUT_SECONDS
    body = bytearray()
    try:
        # S310: every caller builds an https literal here, interpolating only a quoted
        # path segment, so no caller-supplied scheme can reach urlopen.
        with urllib.request.urlopen(  # noqa: S310
            url, timeout=_STALL_TIMEOUT_SECONDS
        ) as response:
            while True:
                if time.monotonic() + _STALL_TIMEOUT_SECONDS > deadline:
                    raise _ResolverUnavailableError(f"{url} timed out")
                chunk = response.read1(_READ_CHUNK_BYTES)
                if not chunk:
                    break
                body += chunk
                if len(body) > _MAX_RESPONSE_BYTES:
                    raise _ResolverUnavailableError(f"{url} sent too much data")
            # read1 reports a connection that dies mid-body as EOF rather than raising,
            # so a short read has to be caught here: a truncated SMILES can still parse,
            # as a different molecule. read() would have raised, but cannot be used
            # without giving up the guarantees above.
            declared = response.getheader("Content-Length")
            if (
                declared is not None
                and declared.isdigit()
                and len(body) != int(declared)
            ):
                raise _ResolverUnavailableError(
                    f"{url} sent {len(body)} bytes of a declared {declared}"
                )
    except http.client.HTTPException as error:
        # A chunked body that ends early raises IncompleteRead, which is not an OSError
        # and so is not covered by the URLError clause below.
        raise _ResolverUnavailableError(
            f"{url} sent a bad response: {error}"
        ) from error
    except urllib.error.HTTPError as error:
        if error.code in _BACK_OFF_STATUSES:
            raise _ResolverUnavailableError(
                f"{url} returned {error.code}; backing off"
            ) from error
        # Any other status describes the identifier rather than the service -- notably
        # CIR, which answers 500 for a name it does not know -- so keep asking.
        raise _ResolverError(f"{url} returned {error.code}") from error
    except (urllib.error.URLError, TimeoutError) as error:
        raise _ResolverUnavailableError(
            f"{url} could not be reached: {error}"
        ) from error
    return body.decode().strip()


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
