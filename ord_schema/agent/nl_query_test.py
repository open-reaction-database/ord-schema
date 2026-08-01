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

"""Tests for ord_schema.nl_query."""

import json
from types import SimpleNamespace
from unittest import mock

import anthropic
import httpx
import pytest

from ord_schema.agent import nl_query


class _FakeCache:
    """An in-memory Cache implementation; records calls for assertions."""

    def __init__(self, initial=None):
        self.store = dict(initial or {})
        self.writes = []

    async def get(self, key):
        return self.store.get(key)

    async def set(self, key, value, ttl_seconds):
        self.writes.append((key, value, ttl_seconds))
        self.store[key] = value


def _tool_response(payload):
    block = anthropic.types.ToolUseBlock(
        type="tool_use", id="toolu_test", name="build_query", input=payload
    )
    return SimpleNamespace(content=[block])


def _client(response=None, error=None):
    client = mock.AsyncMock()
    if error is not None:
        client.messages.create.side_effect = error
    else:
        client.messages.create.return_value = response
    return client


# Translation


@pytest.mark.asyncio
async def test_translate_parses_a_forced_tool_call():
    client = _client(
        _tool_response(
            {
                "components": [
                    {"identifier": "ibuprofen", "target": "OUTPUT", "mode": "EXACT"}
                ],
                "min_yield": 70,
            }
        )
    )
    result = await nl_query.translate("reactions making ibuprofen over 70%", client)
    assert result.components[0].identifier == "ibuprofen"
    assert result.components[0].target == "OUTPUT"
    assert result.min_yield == 70


@pytest.mark.asyncio
async def test_translate_without_a_tool_call_is_malformed():
    client = _client(SimpleNamespace(content=[SimpleNamespace(type="text")]))
    with pytest.raises(nl_query.MalformedQueryError):
        await nl_query.translate("hello", client)


@pytest.mark.asyncio
async def test_translate_with_an_invalid_payload_is_malformed():
    client = _client(_tool_response({"components": [{"identifier": "x"}]}))
    with pytest.raises(nl_query.MalformedQueryError):
        await nl_query.translate("hello", client)


@pytest.mark.asyncio
async def test_translate_maps_a_rate_limit_to_its_own_error():
    request = httpx.Request("POST", "https://api.anthropic.com")
    error = anthropic.RateLimitError(
        "slow down", response=httpx.Response(429, request=request), body=None
    )
    with pytest.raises(nl_query.ModelRateLimitedError):
        await nl_query.translate("hello", _client(error=error))


@pytest.mark.asyncio
async def test_translate_maps_an_api_error_to_unavailable():
    error = anthropic.APIConnectionError(
        request=httpx.Request("POST", "https://api.anthropic.com")
    )
    with pytest.raises(nl_query.ModelUnavailableError):
        await nl_query.translate("hello", _client(error=error))


def test_translate_errors_share_a_base_class():
    # A caller that does not care which failure it was can catch one exception.
    for error in (
        nl_query.ModelUnavailableError,
        nl_query.ModelRateLimitedError,
        nl_query.MalformedQueryError,
    ):
        assert issubclass(error, nl_query.NLQueryError)


def test_get_client_without_an_api_key_is_unavailable(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    with pytest.raises(nl_query.ModelUnavailableError, match="ANTHROPIC_API_KEY"):
        nl_query.get_client()


def test_model_name_honors_the_environment(monkeypatch):
    monkeypatch.delenv("ORD_NL_QUERY_MODEL", raising=False)
    assert nl_query.model_name() == nl_query.DEFAULT_MODEL
    monkeypatch.setenv("ORD_NL_QUERY_MODEL", "claude-sonnet-5")
    assert nl_query.model_name() == "claude-sonnet-5"


def test_the_tool_schema_is_the_query_model():
    assert nl_query._TOOL["name"] == "build_query"
    properties = nl_query._TOOL["input_schema"]["properties"]
    assert "components" in properties
    assert "min_yield" in properties


def test_the_prompt_ships_with_the_package():
    assert "build_query" in nl_query.SYSTEM_PROMPT


# Resolution


@pytest.mark.asyncio
async def test_a_verbatim_smiles_is_not_sent_to_a_name_resolver(monkeypatch):
    monkeypatch.setattr(
        nl_query, "resolve_name", mock.Mock(side_effect=AssertionError("resolved"))
    )
    component = nl_query.NLComponent(
        identifier="c1ccccc1", target="INPUT", mode="EXACT"
    )
    resolved = await nl_query.resolve_component(component)
    assert resolved.resolver == "SMILES (verbatim)"
    assert resolved.smiles == "c1ccccc1"


@pytest.mark.asyncio
async def test_a_name_falls_back_to_the_resolver(monkeypatch):
    monkeypatch.setattr(
        nl_query, "resolve_name", mock.Mock(return_value=("CCO", "PubChem API"))
    )
    component = nl_query.NLComponent(identifier="ethanol", target="INPUT", mode="EXACT")
    resolved = await nl_query.resolve_component(component)
    assert (resolved.smiles, resolved.resolver) == ("CCO", "PubChem API")


@pytest.mark.asyncio
async def test_smarts_passes_through_untouched():
    component = nl_query.NLComponent(
        identifier="[#6][Br]", target="OUTPUT", mode="SMARTS"
    )
    resolved = await nl_query.resolve_component(component)
    assert resolved.smiles == "[#6][Br]"
    assert resolved.resolver == "SMARTS (verbatim)"


@pytest.mark.asyncio
async def test_an_unparseable_smarts_is_rejected():
    component = nl_query.NLComponent(
        identifier="not a smarts", target="OUTPUT", mode="SMARTS"
    )
    with pytest.raises(ValueError, match="Invalid SMARTS"):
        await nl_query.resolve_component(component)


@pytest.mark.asyncio
async def test_an_unresolvable_name_raises(monkeypatch):
    monkeypatch.setattr(
        nl_query, "resolve_name", mock.Mock(side_effect=ValueError("no such compound"))
    )
    component = nl_query.NLComponent(
        identifier="unobtainium", target="INPUT", mode="EXACT"
    )
    with pytest.raises(ValueError, match="no such compound"):
        await nl_query.resolve_component(component)


# Caching is injected, and optional


@pytest.mark.asyncio
async def test_resolution_works_with_no_cache_at_all(monkeypatch):
    resolver = mock.Mock(return_value=("CCO", "PubChem API"))
    monkeypatch.setattr(nl_query, "resolve_name", resolver)
    component = nl_query.NLComponent(identifier="ethanol", target="INPUT", mode="EXACT")
    assert (await nl_query.resolve_component(component)).smiles == "CCO"
    assert resolver.call_count == 1


@pytest.mark.asyncio
async def test_a_cache_hit_skips_the_resolver(monkeypatch):
    monkeypatch.setattr(
        nl_query, "resolve_name", mock.Mock(side_effect=AssertionError("resolved"))
    )
    key = nl_query.resolve_cache_key("ethanol")
    cache = _FakeCache({key: json.dumps(["CCO", "PubChem API"])})
    component = nl_query.NLComponent(identifier="ethanol", target="INPUT", mode="EXACT")
    resolved = await nl_query.resolve_component(component, cache)
    assert resolved.smiles == "CCO"
    assert resolved.resolver == "PubChem API (cached)"


@pytest.mark.asyncio
async def test_a_cache_miss_writes_the_resolution_back(monkeypatch):
    monkeypatch.setattr(
        nl_query, "resolve_name", mock.Mock(return_value=("CCO", "PubChem API"))
    )
    cache = _FakeCache()
    component = nl_query.NLComponent(identifier="ethanol", target="INPUT", mode="EXACT")
    await nl_query.resolve_component(component, cache)
    assert len(cache.writes) == 1
    key, value, ttl = cache.writes[0]
    assert key == nl_query.resolve_cache_key("ethanol")
    assert json.loads(value) == ["CCO", "PubChem API"]
    assert ttl == nl_query.RESOLVE_CACHE_TTL_SECONDS


@pytest.mark.asyncio
async def test_a_corrupt_cache_entry_falls_back_to_the_resolver(monkeypatch):
    monkeypatch.setattr(
        nl_query, "resolve_name", mock.Mock(return_value=("CCO", "PubChem API"))
    )
    key = nl_query.resolve_cache_key("ethanol")
    cache = _FakeCache({key: "not json"})
    component = nl_query.NLComponent(identifier="ethanol", target="INPUT", mode="EXACT")
    resolved = await nl_query.resolve_component(component, cache)
    assert resolved.smiles == "CCO"


@pytest.mark.asyncio
async def test_a_failing_resolution_is_not_cached(monkeypatch):
    monkeypatch.setattr(
        nl_query, "resolve_name", mock.Mock(side_effect=ValueError("transient"))
    )
    cache = _FakeCache()
    component = nl_query.NLComponent(
        identifier="unobtainium", target="INPUT", mode="EXACT"
    )
    with pytest.raises(ValueError, match="transient"):
        await nl_query.resolve_component(component, cache)
    assert cache.writes == []


def test_cache_keys_are_scoped_by_version_and_model(monkeypatch):
    assert nl_query.resolve_cache_key("THF").startswith(
        f"nl_resolve:{nl_query.RESOLVE_CACHE_VERSION}:"
    )
    monkeypatch.setenv("ORD_NL_QUERY_MODEL", "model-a")
    first = nl_query.translation_cache_key("q")
    monkeypatch.setenv("ORD_NL_QUERY_MODEL", "model-b")
    # A different model must not read another model's interpretation.
    assert nl_query.translation_cache_key("q") != first


def test_resolution_cache_keys_ignore_case_and_padding():
    assert nl_query.resolve_cache_key("  ThF ") == nl_query.resolve_cache_key("thf")


def test_translation_cache_keys_are_case_sensitive():
    # A question can contain a SMILES, and CCO is not cco.
    assert nl_query.translation_cache_key("with CCO") != nl_query.translation_cache_key(
        "with cco"
    )
    assert nl_query.translation_cache_key(" q ") == nl_query.translation_cache_key("q")
