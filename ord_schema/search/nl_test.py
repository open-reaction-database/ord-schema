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
from anthropic.types import TextBlock, ThinkingBlock, ToolUseBlock

from ord_schema import parquet
from ord_schema.artifacts import projection, structures
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.search import execute, nl, nl_log, placeholders


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
        str(projected),
        str(structured),
        resolver={"pyridine": "C1=CC=NC=C1"}.__getitem__,
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


@dataclasses.dataclass(frozen=True)
class _SafetyRefusal:
    """A canned decline from the model's safety classifiers."""

    category: str


class _StubClient:
    """Returns canned responses in order, and records the requests it was given.

    A string entry becomes a text response, a dict becomes a ``build_query`` tool call,
    a list is the response's content as given, and an exception is raised as though the
    API had raised it.
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
        if isinstance(payload, _SafetyRefusal):
            return types.SimpleNamespace(
                content=[],
                usage=usage,
                stop_reason="refusal",
                stop_details=types.SimpleNamespace(category=payload.category),
            )
        if isinstance(payload, list):
            content = payload
        elif isinstance(payload, str):
            content = [TextBlock(type="text", text=payload)]
        elif isinstance(payload, _Refusal):
            content = [
                ToolUseBlock(
                    type="tool_use",
                    id="toolu_stub",
                    name="cannot_answer",
                    input={"reason": payload.reason},
                )
            ]
        else:
            content = [
                ToolUseBlock(
                    type="tool_use", id="toolu_stub", name="build_query", input=payload
                )
            ]
        stop = (
            "tool_use"
            if any(isinstance(block, ToolUseBlock) for block in content)
            else "end_turn"
        )
        return types.SimpleNamespace(content=content, usage=usage, stop_reason=stop)


def _stub(*responses) -> Any:
    """Returns a stub client, typed loosely so it stands in for the real one.

    The module's signatures name ``anthropic.Anthropic`` because that is what a caller
    passes; only this file substitutes something else.
    """
    return _StubClient(*responses)


def _predicate(translated) -> Any:
    """Returns a query's predicate, which a caller has to have."""
    assert translated.where is not None
    return translated.where


def _where(translation) -> Any:
    """Returns a translation's predicate."""
    return _predicate(translation.query)


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


def test_a_tree_closed_with_a_brace_too_many_is_still_understood():
    # Haiku 5.5 often closes an encoded predicate with the outer query's brace too.
    client = _stub({"where": json.dumps(_SOLVENT) + "}"})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"
    assert len(client.requests) == 1


def test_a_string_with_more_than_closers_after_it_stays_a_string():
    # Only stray closers are dropped; anything else after the value means the string
    # was not an encoded value at all.
    assert nl._coerce('{"a": 1} trailing') == '{"a": 1} trailing'
    assert nl._coerce("aspirin") == "aspirin"
    assert nl._coerce('{"a": "[Pd]"}]}') == {"a": "[Pd]"}


def test_a_tree_returned_as_an_object_is_understood_too():
    client = _stub({"where": _SOLVENT})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"


def test_a_query_nested_under_the_tool_s_own_name_is_unwrapped():
    # Sonnet nests its arguments under the tool's parameter name on most questions, and
    # repeats it when the failure is handed back. Query declares no "query" field and
    # refuses undeclared ones, so this payload can mean nothing else.
    client = _stub({"query": {"where": _SOLVENT}})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"


def test_a_nested_query_encoded_as_a_string_is_unwrapped_too():
    # The two model habits compose: the envelope arrives with the tree inside it still
    # JSON-encoded.
    client = _stub({"query": json.dumps({"where": _SOLVENT})})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"


def test_a_key_named_query_beside_others_is_not_an_envelope():
    # Only a payload that is exactly the envelope is unwrapped. Anything else is the
    # model having gone wrong in a way this cannot read, and belongs in the repair turn
    # where it is told what failed.
    client = _stub({"query": {"where": _SOLVENT}, "limit": 5})
    with pytest.raises(nl.MalformedQueryError):
        nl.translate("solvent reactions", client=client, repair=False)


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


def test_a_response_with_no_tool_call_is_handed_back_once():
    # No tool is forced, so a turn in prose is possible and gets the repair turn.
    client = _stub("Pyridine is a common solvent.", {"where": _SOLVENT})
    assert _where(nl.translate("solvent reactions", client=client)).op == "exists"
    reminder = client.requests[1]["messages"][-1]["content"]
    assert "build_query" in reminder
    assert "cannot_answer" in reminder


def test_a_repair_after_a_reply_in_prose_counts_both_turns():
    # One attempt, since the first turn wrote no query; the turns are what say a repair
    # ran, and the log reads "repaired" from them.
    client = _stub("Pyridine is a common solvent.", {"where": _SOLVENT})
    translated = nl.translate("solvent reactions", client=client)
    assert len(translated.attempts) == 1
    assert translated.turns == 2


def test_a_first_try_is_one_turn():
    assert (
        nl.translate("solvent reactions", client=_stub({"where": _SOLVENT})).turns == 1
    )


def test_a_failure_carries_the_turns_it_took():
    client = _stub("I would rather chat.", _Refusal("no"))
    with pytest.raises(nl.UnanswerableError) as raised:
        nl.translate("solvent reactions", client=client)
    assert raised.value.turns == 2
    assert not raised.value.attempted


def test_a_second_response_with_no_tool_call_raises():
    client = _stub("I would rather chat.", "Still chatting.")
    with pytest.raises(nl.MalformedQueryError, match="no query"):
        nl.translate("solvent reactions", client=client)
    assert len(client.requests) == 2


def test_a_response_with_no_tool_call_raises_when_repair_is_off():
    client = _stub("I would rather chat about the weather.")
    with pytest.raises(nl.MalformedQueryError, match="no query") as raised:
        nl.translate("solvent reactions", client=client, repair=False)
    assert len(client.requests) == 1
    # Paid for, though it built nothing to record as an attempt.
    assert raised.value.usage.output == 1
    assert not raised.value.attempts


def test_the_repair_replays_the_whole_turn():
    # A model that thinks needs its thinking back beside the tool call it led to.
    thinking = ThinkingBlock(type="thinking", thinking="", signature="sig")
    rejected = ToolUseBlock(
        type="tool_use", id="toolu_stub", name="build_query", input={"where": _BAD_PATH}
    )
    client = _stub([thinking, rejected], {"where": _SOLVENT})
    nl.translate("aspirin reactions", client=client)
    assert client.requests[1]["messages"][1] == {
        "role": "assistant",
        "content": [thinking, rejected],
    }


def test_effort_is_sent_only_when_asked_for():
    # Models without effort levels reject the parameter, so the default omits it.
    client = _stub({"where": _SOLVENT}, {"where": _SOLVENT})
    nl.translate("solvent reactions", client=client)
    nl.translate("solvent reactions", client=client, effort="low")
    assert "output_config" not in client.requests[0]
    assert client.requests[1]["output_config"] == {"effort": "low"}


def test_the_answer_call_takes_an_effort_too():
    client = _stub("Two reactions.")
    nl.answer(
        "how many?", pa.table({"reaction_id": ["a", "b"]}), client=client, effort="low"
    )
    assert client.requests[0]["output_config"] == {"effort": "low"}


def test_a_safety_refusal_is_named_as_one():
    client = _stub(_SafetyRefusal("bio"))
    with pytest.raises(nl.ModelRefusedError, match="bio") as raised:
        nl.translate("solvent reactions", client=client)
    # Nothing a repair could fix, so it is not handed back.
    assert len(client.requests) == 1
    assert nl.outcome_of(raised.value) == nl_log.Outcome.REFUSED
    assert not issubclass(nl.ModelRefusedError, nl.UnanswerableError)


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
    described = nl.answer("how many?", table, client=client)
    sent = client.requests[0]["messages"][0]["content"]
    assert described.text == "A hundred thousand reactions."
    assert "100000 rows" in sent
    # Only the first few rows travel; a row from the middle would mean the table did.
    assert '"reaction_id": "50000"' not in sent


def test_ask_returns_the_query_it_ran(corpus):
    # A caller shows the user what was searched, and offers to run it again.
    client = _stub({"where": _SOLVENT}, "Two reactions, both with pyridine.")
    answer = nl.ask("solvent reactions", corpus, client=client)
    assert answer.question == "solvent reactions"
    assert _predicate(answer.query).op == "exists"
    assert answer.table.num_rows == corpus.search(answer.query).num_rows
    assert answer.text == "Two reactions, both with pyridine."
    assert answer.resolved == {}


def test_ask_says_what_each_compound_resolved_to(corpus):
    # The query keeps the name the model wrote; the structure searched for sits beside
    # it, so a reader can see what the name became.
    where = {
        "op": "exists",
        "path": "inputs.components",
        "where": {"op": "same_compound", "path": "smiles", "compound": "pyridine"},
    }
    client = _stub({"where": where}, "One reaction used pyridine.")
    answer = nl.ask("reactions using pyridine", corpus, client=client)
    assert _predicate(answer.query).where.compound == "pyridine"
    assert answer.resolved == {"pyridine": "c1ccncc1"}
    assert answer.table.num_rows == 1


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


def test_both_tools_are_offered_and_neither_is_forced():
    # Newer models refuse a forced tool_choice; one call at most leaves a repair one
    # tool_use to answer.
    client = _stub({"where": _SOLVENT})
    nl.translate("solvent reactions", client=client)
    request = client.requests[0]
    assert [tool["name"] for tool in request["tools"]] == [
        "build_query",
        "cannot_answer",
    ]
    assert request["tool_choice"] == {"type": "auto", "disable_parallel_tool_use": True}


def test_declining_is_told_apart_from_failing():
    # A caller shows these differently: one is "ask me another way", the other is a bug
    # report, and both are NLQueryError to anything that only wants to catch one thing.
    assert issubclass(nl.UnanswerableError, nl.NLQueryError)
    assert not issubclass(nl.UnanswerableError, nl.MalformedQueryError)


def test_a_query_that_asks_nothing_is_handed_back():
    # An empty query compiles and returns the corpus, which reads as an answer rather
    # than as the failure it is. Opus often writes one for "reactions with no solvent
    # at all".
    client = _stub({}, {"where": _SOLVENT})
    result = nl.translate("reactions with no solvent at all", client=client)
    assert _where(result) is not None
    assert len(client.requests) == 2


def test_a_query_that_asks_nothing_fails_when_repair_is_off():
    client = _stub({})
    with pytest.raises(nl.MalformedQueryError, match="asks nothing"):
        nl.translate("anything", client=client, repair=False)


def test_a_bare_limit_is_a_question_someone_asks():
    # "Show me ten reactions" needs no predicate, and the corpus is not what comes back.
    client = _stub({"limit": 10})
    assert nl.translate("show me ten reactions", client=client).query.limit == 10


class _CollectingSink:
    """Keeps every record written to it, so a test can read what was logged."""

    def __init__(self):
        self.events: list[dict] = []

    def write(self, event):
        self.events.append(event)


def test_a_translation_carries_what_it_cost():
    # Every turn's price rides on the translation, which is what makes the cost of a
    # question observable and the cheap-model-plus-repair question answerable.
    client = _stub({"where": _SOLVENT})
    translation = nl.translate("solvent reactions", client=client)
    assert _where(translation).op == "exists"
    assert len(translation.attempts) == 1
    assert translation.attempts[0].usage.input == 1
    assert translation.attempts[0].error is None


def test_a_repaired_translation_keeps_the_query_the_compiler_refused():
    # The rejected query and the error rejecting it are what prompt work runs on.
    client = _stub({"where": _BAD_PATH}, {"where": _SOLVENT})
    translation = nl.translate("solvent reactions", client=client)
    assert len(translation.attempts) == 2
    assert translation.attempts[0].translation == {"where": _BAD_PATH}
    assert "identifiers" in str(translation.attempts[0].error)
    assert translation.attempts[1].error is None


def test_a_malformed_failure_reports_the_attempts_it_paid_for():
    # An exception carries no return value, so without this the two turns a failed
    # translation burned would go unrecorded -- and failures are the records worth
    # having.
    client = _stub({"where": _BAD_PATH}, {"where": _BAD_PATH})
    with pytest.raises(nl.MalformedQueryError) as caught:
        nl.translate("solvent reactions", client=client)
    assert len(caught.value.attempts) == 2
    assert all(attempt.error for attempt in caught.value.attempts)


def test_ask_records_the_effort_and_the_turns(corpus):
    sink = _CollectingSink()
    nl.ask(
        "solvent reactions",
        corpus,
        client=_stub("Prose first.", {"where": _SOLVENT}, "One reaction."),
        effort="low",
        sink=sink,
    )
    (event,) = sink.events
    assert event["effort"] == "low"
    assert event["turns"] == 2


def test_ask_records_the_question_and_what_it_became(corpus):
    sink = _CollectingSink()
    client = _stub({"where": _SOLVENT}, "Two reactions, both with pyridine.")
    answer = nl.ask(
        "solvent reactions",
        corpus,
        client=client,
        sink=sink,
        source="web",
        session_id="a91b",
    )
    (event,) = sink.events
    assert event["question"] == "solvent reactions"
    assert event["outcome"] == nl_log.Outcome.ANSWERED
    assert event["source"] == "web"
    assert event["session_id"] == "a91b"
    assert event["corpus_fingerprint"] == corpus.fingerprint
    assert event["row_count"] == answer.table.num_rows
    # The caller needs the identifier back, or a thumb has nothing to reference.
    assert event["record_id"] == answer.record_id


def test_the_prompt_fingerprint_tracks_the_cached_prefix(monkeypatch):
    # It answers "which translator wrote this" when an old record is compared against
    # today's behavior, and it moves exactly when the cached prefix does -- so a run
    # where it changes unexpectedly is a cache miss already paid for.
    before = nl.prompt_fingerprint()
    assert before == nl.prompt_fingerprint()
    monkeypatch.setattr(nl, "SYSTEM_PROMPT", nl.SYSTEM_PROMPT + " and one more rule")
    # Cached against the module constants, so anything replacing one has to say so;
    # a served process digests the whole prefix once rather than once per question.
    assert nl.prompt_fingerprint() == before
    nl.prompt_fingerprint.cache_clear()
    assert nl.prompt_fingerprint() != before
    nl.prompt_fingerprint.cache_clear()


def test_a_corpus_whose_artifacts_disagree_is_not_an_unresolved_compound(corpus):
    # PairingError is a ValueError, and the catch-all for ValueError is the resolver's.
    # An unresolvable compound is the question's fault and a corpus that does not pair
    # is the deployment's, which is the distinction the taxonomy exists to keep.
    assert (
        nl.outcome_of(execute.PairingError("projections and structures do not pair up"))
        is nl_log.Outcome.SEARCH_FAILED
    )
    assert (
        nl.outcome_of(ValueError("cannot resolve 'unobtainium'"))
        is nl_log.Outcome.UNRESOLVED_COMPOUND
    )


def test_a_question_nobody_logs_never_opens_an_artifact(corpus, monkeypatch):
    # Reading the stamps opens a footer per artifact, half a second on a real corpus.
    # A caller that passes no sink is not logging and must not pay for it.
    def _refuse(self):
        raise AssertionError("the fingerprint was read for a question nobody logs")

    monkeypatch.setattr(type(corpus), "fingerprint", property(_refuse))
    answer = nl.ask(
        "solvent reactions", corpus, client=_stub({"where": _SOLVENT}, "Two.")
    )
    # And nothing hands back an id for a record that was never written: a thumb built
    # against it would reference an ask no row carries, and read() would drop it.
    assert answer.record_id == ""


def test_an_unreadable_artifact_costs_the_record_not_the_answer(
    corpus, monkeypatch, caplog
):
    def _fail(self):
        raise OSError("artifact footer unreadable")

    monkeypatch.setattr(type(corpus), "fingerprint", property(_fail))
    sink = _CollectingSink()
    answer = nl.ask(
        "solvent reactions",
        corpus,
        client=_stub({"where": _SOLVENT}, "Two."),
        sink=sink,
    )
    assert answer.text == "Two."
    (event,) = sink.events
    assert event["corpus_fingerprint"] == ""
    assert "artifact footer unreadable" in caplog.text


def test_a_record_names_the_translator_that_wrote_it(corpus):
    sink = _CollectingSink()
    client = _stub({"where": _SOLVENT}, "Two reactions.")
    nl.ask("solvent reactions", corpus, client=client, sink=sink)
    (event,) = sink.events
    assert event["prompt_fingerprint"] == nl.prompt_fingerprint()


def test_what_a_question_cost_counts_every_call_it_made(corpus):
    # Three calls answer this one: a translation the compiler refuses, its repair, and
    # the sentence describing the result. A total summed from the attempts alone would
    # miss the third, and a decline would go unpriced entirely.
    sink = _CollectingSink()
    client = _stub({"where": _BAD_PATH}, {"where": _SOLVENT}, "Two reactions.")
    nl.ask("solvent reactions", corpus, client=client, sink=sink)
    (event,) = sink.events
    assert len(json.loads(event["attempts"])) == 2
    assert event["usage"]["input"] == 3


def test_a_declined_question_is_priced_even_though_it_built_nothing(corpus):
    sink = _CollectingSink()
    client = _stub(_Refusal("values are literals, never other columns"))
    with pytest.raises(nl.UnanswerableError):
        nl.ask("compare two columns", corpus, client=client, sink=sink)
    (event,) = sink.events
    assert json.loads(event["attempts"]) == []
    assert event["usage"]["input"] == 1


def test_a_record_says_what_each_compound_resolved_to(corpus):
    # The record keeps the query as the model wrote it, naming the compound, and a
    # resolver can answer differently later; the structure searched is what reproduces
    # the answer.
    sink = _CollectingSink()
    where = {
        "op": "exists",
        "path": "inputs.components",
        "where": {"op": "same_compound", "path": "smiles", "compound": "pyridine"},
    }
    client = _stub({"where": where}, "One reaction used pyridine.")
    nl.ask("reactions using pyridine", corpus, client=client, sink=sink)
    nl.ask(
        "solvent reactions",
        corpus,
        client=_stub({"where": _SOLVENT}, "Two."),
        sink=sink,
    )
    named, unnamed = sink.events
    assert json.loads(named["resolved"]) == {"pyridine": "c1ccncc1"}
    assert unnamed["resolved"] == "{}"


def test_a_search_whose_description_failed_still_records_its_names(corpus):
    sink = _CollectingSink()
    where = {
        "op": "exists",
        "path": "inputs.components",
        "where": {"op": "same_compound", "path": "smiles", "compound": "pyridine"},
    }
    client = _stub({"where": where}, _status_error(429))
    with pytest.raises(nl.ModelRateLimitedError):
        nl.ask("reactions using pyridine", corpus, client=client, sink=sink)
    (event,) = sink.events
    assert event["row_count"] == 1
    assert json.loads(event["resolved"]) == {"pyridine": "c1ccncc1"}


def test_a_question_that_matched_nothing_is_recorded_as_empty(corpus):
    sink = _CollectingSink()
    where = {"op": "eq", "path": "reaction_id", "value": {"literal": "ord-nope"}}
    client = _stub({"where": where}, "No reactions matched.")
    nl.ask("reactions that do not exist", corpus, client=client, sink=sink)
    (event,) = sink.events
    assert event["outcome"] == nl_log.Outcome.EMPTY
    assert event["row_count"] == 0


def test_a_failed_ask_is_still_recorded(corpus):
    # The failures are the records worth having, so the emit cannot sit on the path
    # that only runs when everything worked.
    sink = _CollectingSink()
    client = _stub({"where": _BAD_PATH}, {"where": _BAD_PATH})
    with pytest.raises(nl.MalformedQueryError):
        nl.ask("solvent reactions", corpus, client=client, sink=sink)
    (event,) = sink.events
    assert event["outcome"] == nl_log.Outcome.MALFORMED
    assert len(json.loads(event["attempts"])) == 2
    assert event["row_count"] is None
    assert event["resolved"] is None


def test_a_declined_question_records_whether_a_query_came_first(corpus):
    sink = _CollectingSink()
    client = _stub(_Refusal("values are literals, never other columns"))
    with pytest.raises(nl.UnanswerableError):
        nl.ask("compare two columns", corpus, client=client, sink=sink)
    (event,) = sink.events
    assert event["outcome"] == nl_log.Outcome.DECLINED
    assert event["declined_reason"] == "values are literals, never other columns"
    # No attempt: the model read the question rather than backing off after a refusal.
    assert json.loads(event["attempts"]) == []


# Placeholders

_PYRIDINE = placeholders.Placeholder(
    format=placeholders.StructureFormat.COMPOUND_SMILES, text="c1ccncc1"
)


def _uses(slot: dict) -> dict:
    """Returns a query for reactions with an input matching a structure predicate."""
    return {"where": {"op": "exists", "path": "inputs.components", "where": slot}}


def test_placeholders_follow_the_question():
    client = _stub(_uses({"op": "same_compound", "path": "smiles", "compound": "@A"}))
    nl.translate("reactions using @A", structures={"@A": _PYRIDINE}, client=client)
    assert client.requests[0]["messages"][0]["content"] == (
        "reactions using @A\n\nPlaceholders in this question:\n@A is a compound SMILES."
    )


def test_a_question_without_placeholders_is_sent_alone():
    client = _stub({"where": _SOLVENT})
    nl.translate("reactions run in a solvent", client=client)
    assert client.requests[0]["messages"][0]["content"] == "reactions run in a solvent"


def test_a_placeholder_is_substituted_and_the_attempt_keeps_it():
    # The query names the structure outright, so it runs again as it stands; the
    # attempt keeps what the model wrote, which is what prompt work reads.
    client = _stub(_uses({"op": "same_compound", "path": "smiles", "compound": "@A"}))
    translated = nl.translate(
        "reactions using @A", structures={"@A": _PYRIDINE}, client=client
    )
    clause = _where(translated).where
    assert (clause.smiles, clause.compound) == ("c1ccncc1", None)
    assert translated.attempts[0].translation["where"]["where"]["compound"] == "@A"


def test_a_format_the_predicate_refuses_is_handed_back():
    pattern = placeholders.Placeholder(
        format=placeholders.StructureFormat.COMPOUND_SMARTS, text="c[Br,I]"
    )
    client = _stub(
        _uses(
            {"op": "similarity", "path": "smiles", "compound": "@B", "threshold": 0.6}
        ),
        _uses({"op": "substructure", "path": "smiles", "compound": "@B"}),
    )
    translated = nl.translate(
        "inputs like @B", structures={"@B": pattern}, client=client
    )
    error = translated.attempts[0].error
    assert error is not None
    assert "@B is a compound SMARTS, which similarity does not take" in error
    assert translated.turns == 2


def test_an_unbound_placeholder_never_reaches_the_resolver(corpus):
    # The corpus's resolver knows only pyridine, so a placeholder looked up as a name
    # would fail differently: as an unresolved compound, after the search began.
    unbound = _uses({"op": "same_compound", "path": "smiles", "compound": "@D"})
    client = _stub(unbound, unbound)
    with pytest.raises(nl.MalformedQueryError, match="@D is not a placeholder"):
        nl.ask(
            "reactions using @A", corpus, structures={"@A": _PYRIDINE}, client=client
        )


def test_a_question_that_fails_the_check_makes_no_call():
    client = _stub()
    with pytest.raises(
        placeholders.PlaceholderError,
        match="@A is in the question but has no structure",
    ):
        nl.translate("reactions using @A", client=client)
    assert client.requests == []


def test_ask_searches_with_the_placeholder(corpus):
    client = _stub(
        _uses({"op": "same_compound", "path": "smiles", "compound": "@A"}),
        "One reaction uses @A.",
    )
    answer = nl.ask(
        "reactions using @A", corpus, structures={"@A": _PYRIDINE}, client=client
    )
    assert answer.table.column("reaction_id").to_pylist() == ["ord-nl01"]
    assert _predicate(answer.query).where.smiles == "c1ccncc1"
