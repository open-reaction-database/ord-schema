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

"""Tests for ord_schema.search.nl.

Nothing here reaches the network: a stub client returns canned tool calls, which is
enough to pin every behavior that is this module's own rather than the model's.
"""

import dataclasses
import json
import types
from typing import Any

import anthropic
import httpx
import pyarrow as pa
import pytest
from anthropic.types import TextBlock, ToolUseBlock

from ord_schema import parquet
from ord_schema.artifacts import projection, structures
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.search import execute, nl


@pytest.fixture(scope="module")
def corpus(tmp_path_factory):
    """A two-reaction corpus, one of them using pyridine as a solvent."""
    root = tmp_path_factory.mktemp("nl")
    reaction = reaction_pb2.Reaction(reaction_id="ord-nl01")
    component = reaction.inputs["in"].components.add()
    component.identifiers.add(type="SMILES", value="c1ccncc1")
    component.reaction_role = reaction_pb2.ReactionRole.SOLVENT
    other = reaction_pb2.Reaction(reaction_id="ord-nl02")
    component = other.inputs["in"].components.add()
    component.identifiers.add(type="SMILES", value="CCO")
    component.reaction_role = reaction_pb2.ReactionRole.REACTANT
    source = root / "data" / "ord_dataset-nl.parquet"
    source.parent.mkdir(parents=True, exist_ok=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-nl",
            name="test",
            description="test",
            reactions=[reaction, other],
        ),
        str(source),
    )
    projected = root / "projections" / source.name
    projected.parent.mkdir(parents=True, exist_ok=True)
    projection.write_projection(source, projected)
    structured = root / "structures" / source.name
    structured.parent.mkdir(parents=True, exist_ok=True)
    structures.write_structures(projected, structured)
    with execute.Corpus(
        str(projected), str(structured), resolver={}.__getitem__
    ) as value:
        yield value


_SOLVENT = {
    "op": "exists",
    "path": "inputs.components",
    "where": {"op": "eq", "path": "reaction_role", "value": {"literal": "SOLVENT"}},
}
# The syntax a model invents when it wants "any identifier", and the compiler's answer
# to it names a real field.
_BAD_PATH = {"op": "eq", "path": "identifiers[*].value", "value": {"literal": "x"}}


@dataclasses.dataclass(frozen=True)
class _Refusal:
    """A canned ``cannot_answer`` call."""

    reason: str


class _StubClient:
    """Returns canned responses in order, and records the requests it was given.

    A string entry becomes a text response, a dict becomes a ``build_query`` tool call,
    and an exception is raised as though the API had raised it.
    """

    def __init__(self, *responses):
        self._responses = list(responses)
        self.requests: list[dict] = []
        self.messages = self

    def create(self, **kwargs):
        self.requests.append(kwargs)
        payload = self._responses.pop(0)
        if isinstance(payload, Exception):
            raise payload
        usage = types.SimpleNamespace(
            input_tokens=1,
            output_tokens=1,
            cache_creation_input_tokens=0,
            cache_read_input_tokens=0,
        )
        if isinstance(payload, str):
            block: TextBlock | ToolUseBlock = TextBlock(type="text", text=payload)
            stop = "end_turn"
        elif isinstance(payload, _Refusal):
            block = ToolUseBlock(
                type="tool_use",
                id="toolu_stub",
                name="cannot_answer",
                input={"reason": payload.reason},
            )
            stop = "tool_use"
        else:
            block = ToolUseBlock(
                type="tool_use", id="toolu_stub", name="build_query", input=payload
            )
            stop = "tool_use"
        return types.SimpleNamespace(content=[block], usage=usage, stop_reason=stop)


def _stub(*responses) -> Any:
    """Returns a stub client, typed loosely so it stands in for the real one.

    The module's signatures name ``anthropic.Anthropic`` because that is what a caller
    passes; only this file substitutes something else.
    """
    return _StubClient(*responses)


def _where(result) -> Any:
    """Returns a translated query's predicate, which a caller has to have."""
    assert result.where is not None
    return result.where


def _status_error(status_code: int) -> anthropic.APIStatusError:
    """Returns an APIStatusError the SDK would have raised for a status code."""
    request = httpx.Request("POST", "https://api.anthropic.com/v1/messages")
    response = httpx.Response(status_code, request=request, json={"error": {}})
    if status_code == 429:
        return anthropic.RateLimitError("slow down", response=response, body=None)
    return anthropic.InternalServerError("boom", response=response, body=None)


def test_the_client_prefers_the_dedicated_key(monkeypatch):
    # A shell holding ANTHROPIC_API_KEY for whatever else its owner is doing must not
    # end up paying for a corpus search.
    monkeypatch.setenv(nl.API_KEY_VARIABLE, "sk-ord")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-general")
    assert nl.get_client().api_key == "sk-ord"


def test_the_client_falls_back_to_the_sdk_lookup(monkeypatch):
    monkeypatch.delenv(nl.API_KEY_VARIABLE, raising=False)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-general")
    assert nl.get_client().api_key == "sk-general"


def test_the_errors_share_one_base():
    # A caller catches the base and still tells a rate limit from a query that could
    # not be built; ord-interface maps these onto status codes.
    assert issubclass(nl.ModelRateLimitedError, nl.NLQueryError)
    assert issubclass(nl.ModelUnavailableError, nl.NLQueryError)
    assert issubclass(nl.MalformedQueryError, nl.NLQueryError)


def test_a_tree_returned_as_a_json_string_is_still_understood():
    # What both measured models do most of the time: the predicate arrives encoded in a
    # string rather than as an object.
    client = _stub({"where": json.dumps(_SOLVENT)})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"


def test_a_tree_returned_as_an_object_is_understood_too():
    client = _stub({"where": _SOLVENT})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"


def test_the_prefix_is_marked_cacheable():
    # Cache reads are most of what a query costs; an uncached prefix is a tenfold bill.
    client = _stub({"where": _SOLVENT})
    nl.translate("solvent reactions", client=client)
    assert client.requests[0]["system"][0]["cache_control"] == {
        "type": "ephemeral",
        "ttl": "1h",
    }


def test_the_schema_and_the_rules_reach_the_prompt():
    client = _stub({"where": _SOLVENT})
    nl.translate("solvent reactions", client=client)
    system = client.requests[0]["system"][0]["text"]
    assert "reaction_role" in system
    assert "identifiers[*]" in system


def test_a_bad_path_is_handed_back_once_and_recovered():
    client = _stub({"where": _BAD_PATH}, {"where": _SOLVENT})
    assert _where(nl.translate("aspirin reactions", client=client)).op == "exists"
    assert len(client.requests) == 2


def test_the_repair_carries_the_compiler_s_suggestion():
    # The suggestion is what makes one repair worth paying for.
    client = _stub({"where": _BAD_PATH}, {"where": _SOLVENT})
    nl.translate("aspirin reactions", client=client)
    sent = client.requests[1]["messages"][-1]["content"][0]
    assert sent["is_error"]
    assert "did you mean" in sent["content"]


def test_a_second_failure_raises_rather_than_looping():
    client = _stub({"where": _BAD_PATH}, {"where": _BAD_PATH})
    with pytest.raises(nl.MalformedQueryError, match="identifiers"):
        nl.translate("aspirin reactions", client=client)
    assert len(client.requests) == 2


def test_repair_can_be_turned_off_for_measurement():
    client = _stub({"where": _BAD_PATH})
    with pytest.raises(nl.MalformedQueryError):
        nl.translate("aspirin reactions", client=client, repair=False)
    assert len(client.requests) == 1


def test_a_query_that_does_not_validate_is_repaired_too():
    # Failing pydantic and failing the compiler are the same problem to a caller, and
    # both carry a message worth handing back.
    client = _stub({"where": {"op": "nonsense"}}, {"where": _SOLVENT})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"


def test_a_rate_limit_is_named_as_one():
    client = _stub(_status_error(429))
    with pytest.raises(nl.ModelRateLimitedError):
        nl.translate("solvent reactions", client=client)


def test_an_unreachable_model_is_named_as_one():
    client = _stub(_status_error(500))
    with pytest.raises(nl.ModelUnavailableError):
        nl.translate("solvent reactions", client=client)


def test_a_response_with_no_tool_call_is_refused():
    client = _stub("I would rather chat about the weather.")
    with pytest.raises(nl.MalformedQueryError, match="no query"):
        nl.translate("solvent reactions", client=client)


def test_the_summary_states_the_row_count():
    table = pa.table({"reaction_id": [str(value) for value in range(100_000)]})
    assert "100000 rows" in nl.summarize(table)


def test_the_summary_is_bounded_by_the_row_cap_not_the_table():
    # The prompt cost of describing a result must not grow with the result: a hundred
    # thousand rows and a million differ only by the digits in the count.
    hundred_thousand = nl.summarize(
        pa.table({"reaction_id": [str(value) for value in range(100_000)]})
    )
    million = nl.summarize(
        pa.table({"reaction_id": [str(value) for value in range(1_000_000)]})
    )
    assert abs(len(million) - len(hundred_thousand)) < 20


def test_the_summary_says_how_much_it_left_out():
    table = pa.table({"reaction_id": [str(value) for value in range(9)]})
    assert "4 more rows not shown" in nl.summarize(table)


def test_the_answer_call_sees_the_summary_and_not_the_rows():
    table = pa.table({"reaction_id": [str(value) for value in range(100_000)]})
    client = _stub("A hundred thousand reactions.")
    text = nl.answer("how many?", table, client=client)
    sent = client.requests[0]["messages"][0]["content"]
    assert text == "A hundred thousand reactions."
    assert "100000 rows" in sent
    # Only the first few rows travel; a row from the middle would mean the table did.
    assert '"reaction_id": "50000"' not in sent


def test_ask_returns_the_query_it_ran(corpus):
    # A caller shows the user what was searched, and offers to run it again.
    client = _stub({"where": _SOLVENT}, "Two reactions, both with pyridine.")
    answer = nl.ask("solvent reactions", corpus, client=client)
    assert answer.question == "solvent reactions"
    assert _where(answer.query).op == "exists"
    assert answer.table.num_rows == corpus.search(answer.query).num_rows
    assert answer.text == "Two reactions, both with pyridine."


def test_ask_passes_the_timeout_through(corpus, monkeypatch):
    # A slow question has to fail rather than hold the caller, which it can only do if
    # the bound actually reaches the search.
    seen = {}
    real = corpus.search

    def recording(translated, **kwargs):
        seen.update(kwargs)
        return real(translated)

    monkeypatch.setattr(corpus, "search", recording)
    client = _stub({"where": _SOLVENT}, "Some reactions.")
    nl.ask("solvent reactions", corpus, client=client, timeout_seconds=12.5)
    assert seen == {"timeout_seconds": 12.5}


def test_a_declined_question_is_named_as_unanswerable():
    # Forcing build_query would leave the model no way out, and a model with no way out
    # invents a query rather than refusing.
    client = _stub(_Refusal("a value is a literal, never another column"))
    with pytest.raises(nl.UnanswerableError, match="another column"):
        nl.translate("reactions that ran longer than their workup took", client=client)


def test_a_refusal_after_a_failed_query_is_marked_as_such():
    # Declining on the repair turn still answers the caller, but the model got there by
    # way of a query the compiler refused rather than by reading the question.
    client = _stub({"where": _BAD_PATH}, _Refusal("no column comparison"))
    with pytest.raises(nl.UnanswerableError) as caught:
        nl.translate("reactions that ran longer than their workup took", client=client)
    assert caught.value.attempted
    assert len(client.requests) == 2


def test_an_outright_refusal_is_not_marked_as_attempted():
    client = _stub(_Refusal("no column comparison"))
    with pytest.raises(nl.UnanswerableError) as caught:
        nl.translate("anything unanswerable", client=client)
    assert not caught.value.attempted


def test_a_refusal_is_not_repaired():
    # Nothing was wrong with the model's reasoning, so asking again only costs money.
    client = _stub(_Refusal("the schema does not hold that"))
    with pytest.raises(nl.UnanswerableError):
        nl.translate("anything unanswerable", client=client)
    assert len(client.requests) == 1


def test_both_tools_are_offered_and_one_is_required():
    client = _stub({"where": _SOLVENT})
    nl.translate("solvent reactions", client=client)
    request = client.requests[0]
    assert [tool["name"] for tool in request["tools"]] == [
        "build_query",
        "cannot_answer",
    ]
    assert request["tool_choice"] == {"type": "any"}


def test_declining_is_told_apart_from_failing():
    # A caller shows these differently: one is "ask me another way", the other is a bug
    # report, and both are NLQueryError to anything that only wants to catch one thing.
    assert issubclass(nl.UnanswerableError, nl.NLQueryError)
    assert not issubclass(nl.UnanswerableError, nl.MalformedQueryError)
