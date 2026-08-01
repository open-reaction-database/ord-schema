# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Translating a natural-language question into a structured ORD query.

The model only ever emits a *structured* query, through a forced tool call -- it never
writes SQL and never invents SMILES. Compound names are resolved deterministically by
``ord_schema.resolvers`` (PubChem/CIR/OPSIN), so the chemistry is grounded in a resolver
rather than in model recall.

This module is the translation half only: text in, a resolved interpretation out. It
holds no HTTP, no database, and no cache of its own, so the same translation can serve a
web API, a notebook, and an evaluation harness without any of them inheriting the
others' dependencies. Serving concerns stay with the server:

* **Errors** are library exceptions (:class:`NLQueryError` and its subclasses), which a
  web layer maps to status codes. Nothing here knows what a 503 is.
* **Caching** is injected. Pass any :class:`Cache` implementation to reuse name
  resolutions across requests; omit it and every lookup is live. The cache is always an
  optimization, never a dependency -- an implementation that fails or times out should
  return None rather than raise, so translation falls back to a live call.

Requires the ``agent`` extra; see :mod:`ord_schema.agent`.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
from importlib import resources
from typing import Any, Literal, Protocol, cast

import anthropic
from pydantic import BaseModel, Field, ValidationError
from rdkit import Chem

from ord_schema.logging import get_logger
from ord_schema.resolvers import canonicalize_smiles, resolve_name

logger = get_logger(__name__)

# Haiku is fast and cheap and the translation task is tightly constrained; override with
# ORD_NL_QUERY_MODEL (e.g. a Sonnet snapshot) if disambiguation needs more reasoning.
DEFAULT_MODEL = "claude-haiku-4-5"
MAX_TOKENS = 1024

# Bump when the prompt or the NLQuery schema changes, so a caller's cache does not serve
# interpretations produced by a definition that no longer exists.
TRANSLATION_CACHE_VERSION = "v3"
TRANSLATION_CACHE_TTL_SECONDS = 60 * 60

# Name -> SMILES resolutions are cached separately and for much longer: they are stable
# (a name maps to the same structure) and shared across different questions that mention
# the same compound, which spares the name resolvers repeated lookups.
RESOLVE_CACHE_VERSION = "v1"
RESOLVE_CACHE_TTL_SECONDS = 60 * 60 * 24 * 30

Target = Literal["INPUT", "OUTPUT"]
MatchMode = Literal["EXACT", "SIMILAR", "SUBSTRUCTURE", "SMARTS"]


class NLQueryError(Exception):
    """Base class for failures translating a natural-language question."""


class ModelUnavailableError(NLQueryError):
    """The language model could not be reached, or is not configured."""


class ModelRateLimitedError(NLQueryError):
    """The language model is rate limited; the same request may succeed later."""


class MalformedQueryError(NLQueryError):
    """The model returned no structured query, or one that failed schema validation."""


class Cache(Protocol):
    """An optional key-value cache for name resolutions.

    Implementations must treat every operation as best-effort: return None on a miss
    *or* on any failure, and swallow errors on write. Translation is correct without a
    cache, so an unreachable one must never turn into a failed request.
    """

    async def get(self, key: str) -> str | None:
        """Returns the cached value for ``key``, or None on a miss or any failure."""

    async def set(self, key: str, value: str, ttl_seconds: int) -> None:
        """Stores ``value`` under ``key``, ignoring failures."""


class NLComponent(BaseModel):
    """A single component constraint extracted from the question."""

    identifier: str = Field(
        description=(
            "The compound as the user named it: a common/trade/IUPAC name (e.g. "
            "'ibuprofen', 'benzene'), a SMILES, or -- only when mode is SMARTS -- a "
            "SMARTS pattern. Prefer the plain name; do not translate names to SMILES."
        )
    )
    target: Target = Field(
        description=(
            "INPUT if the compound is a reactant/reagent/solvent consumed by the "
            "reaction ('using X', 'from X', 'with X'); OUTPUT if it is produced by "
            "the reaction ('synthesizing X', 'to make X', 'yields X')."
        )
    )
    mode: MatchMode = Field(
        description=(
            "EXACT for a specific named molecule (the default). SUBSTRUCTURE when the "
            "user wants molecules that merely contain a group/scaffold ('containing a "
            "pyridine ring'). SIMILAR for 'similar to'/'like X'. SMARTS only when the "
            "user supplies or describes an explicit query pattern."
        )
    )


class NLQuery(BaseModel):
    """Structured query produced by the language model from natural language."""

    components: list[NLComponent] = Field(
        default_factory=list,
        description="Per-compound constraints; AND-combined with the other fields.",
    )
    min_yield: float | None = Field(
        default=None, description="Minimum percent yield (0-100), if requested."
    )
    max_yield: float | None = Field(
        default=None, description="Maximum percent yield (0-100), if requested."
    )
    min_conversion: float | None = Field(
        default=None, description="Minimum percent conversion (0-100), if requested."
    )
    max_conversion: float | None = Field(
        default=None, description="Maximum percent conversion (0-100), if requested."
    )
    reaction_smarts: str | None = Field(
        default=None,
        description=(
            "A reaction SMARTS (reactants>>products) when the user describes a "
            "transformation rather than individual components. Usually omitted."
        ),
    )
    similarity_threshold: float | None = Field(
        default=None,
        description="Tanimoto threshold (0-1) for SIMILAR components; default 0.5.",
    )
    use_stereochemistry: bool | None = Field(
        default=None,
        description="True only if the user asks to respect stereochemistry/chirality.",
    )
    limit: int | None = Field(
        default=None, description="Maximum number of reactions to return, if stated."
    )


class ResolvedComponent(BaseModel):
    """A component after name resolution, surfaced for transparency."""

    identifier: str
    smiles: str
    resolver: str
    target: Target
    mode: MatchMode


# The system prompt lives in nl_query_prompt.md (alongside this module) so it reads as
# plain markdown and can be edited without touching Python string escaping.
SYSTEM_PROMPT = (
    (resources.files("ord_schema.agent") / "nl_query_prompt.md")
    .read_text(encoding="utf-8")
    .strip()
)

_TOOL: dict[str, Any] = {
    "name": "build_query",
    "description": "Build a structured ORD search query from the user's question.",
    "input_schema": NLQuery.model_json_schema(),
}


def get_client() -> anthropic.AsyncAnthropic:
    """Returns an Anthropic async client.

    Returns:
        A client reading ANTHROPIC_API_KEY from the environment.

    Raises:
        ModelUnavailableError: If ANTHROPIC_API_KEY is not configured.
    """
    if not os.getenv("ANTHROPIC_API_KEY"):
        raise ModelUnavailableError(
            "Natural-language search is unavailable: ANTHROPIC_API_KEY is not set."
        )
    return anthropic.AsyncAnthropic()


def model_name() -> str:
    """Returns the model to translate with, honoring ORD_NL_QUERY_MODEL."""
    return os.getenv("ORD_NL_QUERY_MODEL", DEFAULT_MODEL)


async def translate(query: str, client: anthropic.AsyncAnthropic) -> NLQuery:
    """Translates a natural-language question into a structured query.

    Args:
        query: The user's free-text question.
        client: Anthropic async client.

    Returns:
        The structured NLQuery produced by the model.

    Raises:
        ModelRateLimitedError: If the model is rate limited.
        ModelUnavailableError: If the model is unreachable or erroring.
        MalformedQueryError: If the model does not return a usable structured query.
    """
    try:
        response = await client.messages.create(
            model=model_name(),
            max_tokens=MAX_TOKENS,
            system=SYSTEM_PROMPT,
            tools=[cast(anthropic.types.ToolParam, _TOOL)],
            tool_choice={"type": "tool", "name": "build_query"},
            messages=[{"role": "user", "content": query}],
        )
    except anthropic.RateLimitError as error:
        raise ModelRateLimitedError(
            "Natural-language search is busy right now; please retry shortly."
        ) from error
    except anthropic.APIError as error:
        # Connection failures, server errors, overloads, and auth problems all degrade
        # to a graceful "temporarily unavailable" rather than an unhandled error.
        logger.warning(f"Anthropic API error during NL translation: {error}")
        raise ModelUnavailableError(
            "Natural-language search is temporarily unavailable."
        ) from error
    for block in response.content:
        if (
            isinstance(block, anthropic.types.ToolUseBlock)
            and block.name == "build_query"
        ):
            try:
                return NLQuery.model_validate(block.input)
            except ValidationError as error:
                # The forced tool schema makes this unlikely, but a payload that slips
                # through is reported as malformed rather than crashing the caller.
                logger.warning(f"Model tool call failed schema validation: {error}")
                raise MalformedQueryError(
                    "Language model returned a malformed structured query."
                ) from error
    raise MalformedQueryError("Language model did not return a structured query.")


def translation_cache_key(query: str) -> str:
    """Returns the cache key for a translation, scoped to the prompt version and model.

    The question is *not* case-folded, unlike a name lookup: a question can carry a
    SMILES, and ``CCO`` and ``cco`` are different molecules. Two questions differing
    only in case therefore get separate entries rather than one wrong answer.
    """
    digest = hashlib.sha256(query.strip().encode()).hexdigest()
    return f"nl_translate:{TRANSLATION_CACHE_VERSION}:{model_name()}:{digest}"


def resolve_cache_key(name: str) -> str:
    """Returns the cache key for a name -> SMILES resolution."""
    digest = hashlib.sha256(name.strip().lower().encode()).hexdigest()
    return f"nl_resolve:{RESOLVE_CACHE_VERSION}:{digest}"


async def _resolve_name_cached(name: str, cache: Cache | None) -> tuple[str, str]:
    """Resolves a compound name to (SMILES, resolver), caching successful lookups.

    The blocking PubChem/CIR/OPSIN call runs in a worker thread. Failures are not cached
    so a transient PubChem outage does not poison the cache with a permanent miss.

    Args:
        name: The compound name to resolve.
        cache: Optional cache; every lookup is live when omitted.

    Returns:
        A tuple of canonical SMILES and the resolver that produced it (e.g. "PubChem
        API"); the resolver is suffixed with " (cached)" on a cache hit.

    Raises:
        ValueError: If the name cannot be resolved to a structure.
    """
    key = resolve_cache_key(name)
    if cache is not None:
        raw = await cache.get(key)
        if raw is not None:
            try:
                smiles, resolver = json.loads(raw)
            except (ValueError, TypeError) as error:
                logger.warning(
                    f"Discarding bad cached resolution for {name!r}: {error}"
                )
            else:
                return smiles, f"{resolver} (cached)"
    smiles, resolver = await asyncio.to_thread(resolve_name, "name", name)
    if cache is not None:
        await cache.set(key, json.dumps([smiles, resolver]), RESOLVE_CACHE_TTL_SECONDS)
    return smiles, resolver


async def resolve_component(
    component: NLComponent, cache: Cache | None = None
) -> ResolvedComponent:
    """Resolves a component's identifier to canonical SMILES.

    SMARTS patterns pass through untouched once validated. Otherwise the identifier
    is treated as a SMILES if RDKit can parse it, and falls back to (cached) name
    resolution via PubChem/CIR/OPSIN.

    Args:
        component: The component to resolve.
        cache: Optional cache for name resolutions.

    Returns:
        The resolved component.

    Raises:
        ValueError: If a SMARTS pattern is unparseable, or a non-SMARTS identifier
            resolves to neither SMILES nor a name.
    """
    if component.mode == "SMARTS":
        # The model authors SMARTS directly; validate up front so a bad pattern fails
        # here rather than deep in query execution (and is caught on dry runs).
        if Chem.MolFromSmarts(component.identifier) is None:
            raise ValueError(f"Invalid SMARTS pattern: {component.identifier!r}")
        return ResolvedComponent(
            identifier=component.identifier,
            smiles=component.identifier,
            resolver="SMARTS (verbatim)",
            target=component.target,
            mode=component.mode,
        )
    try:
        smiles = canonicalize_smiles(component.identifier)
        resolver = "SMILES (verbatim)"
    except ValueError:
        smiles, resolver = await _resolve_name_cached(component.identifier, cache)
    return ResolvedComponent(
        identifier=component.identifier,
        smiles=smiles,
        resolver=resolver,
        target=component.target,
        mode=component.mode,
    )
