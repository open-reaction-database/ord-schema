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

"""Tests for ord_schema.search.nl_eval.

These cover the scoring, the shipped cases, and the outcomes ``run_case`` maps a
translation onto, with a stub client standing in for the model. Measuring against a real
model and a real corpus is a command rather than a test.
"""

import json
import types
from typing import Any, cast

import pyarrow as pa
import pytest
from anthropic.types import ToolUseBlock

from ord_schema.search import execute, nl_eval, nl_log, query

_CASE = nl_eval.EvalCase(
    question="which reactions use pyridine as a solvent?",
    why="two conditions on one element, which a wrong translation splits in two",
    reference={"op": "not_null", "path": "reaction_id"},
)
_INEXPRESSIBLE = nl_eval.EvalCase(
    question="reactions where the temperature exceeds the pressure",
    why="comparing two columns, which the grammar cannot say",
    compiles=False,
)
# The syntax a model invents when it wants "any identifier"; the compiler refuses it.
_BAD_PATH = {
    "where": {"op": "eq", "path": "identifiers[*].value", "value": {"literal": "x"}}
}


class _StubClient:
    """Returns one canned tool call, standing in for the model.

    Attributes:
        messages: Itself, so ``client.messages.create`` reaches ``create``.
    """

    def __init__(self, name: str, payload: dict) -> None:
        """Initializes the stub.

        Args:
            name: Which tool the canned response calls.
            payload: The tool call's input.
        """
        self._name = name
        self._payload = payload
        self.messages = self

    def create(self, **kwargs: Any) -> Any:
        """Returns the canned response.

        Args:
            **kwargs: The request, which the stub ignores.

        Returns:
            A response carrying one ``tool_use`` block.
        """
        del kwargs
        block = ToolUseBlock(
            type="tool_use", id="toolu_stub", name=self._name, input=self._payload
        )
        usage = types.SimpleNamespace(
            input_tokens=1,
            output_tokens=1,
            cache_creation_input_tokens=0,
            cache_read_input_tokens=0,
        )
        return types.SimpleNamespace(
            content=[block], usage=usage, stop_reason="tool_use"
        )


class _SequenceClient:
    """Answers each call with the next canned tool call, for a repair turn.

    Attributes:
        messages: Itself, so ``client.messages.create`` reaches ``create``.
    """

    def __init__(self, *calls: tuple[str, dict]) -> None:
        """Initializes the stub.

        Args:
            *calls: ``(tool name, input)`` pairs, returned in order.
        """
        self._calls = list(calls)
        self.messages = self

    def create(self, **kwargs: Any) -> Any:
        """Returns the next canned response.

        Args:
            **kwargs: The request, which the stub ignores.

        Returns:
            A response carrying one ``tool_use`` block.
        """
        del kwargs
        name, payload = self._calls.pop(0)
        block = ToolUseBlock(type="tool_use", id="toolu_stub", name=name, input=payload)
        usage = types.SimpleNamespace(
            input_tokens=1,
            output_tokens=1,
            cache_creation_input_tokens=0,
            cache_read_input_tokens=0,
        )
        return types.SimpleNamespace(
            content=[block], usage=usage, stop_reason="tool_use"
        )


def _run(case, client, *, repair: bool = False) -> nl_eval.CaseResult:
    """Runs a case whose translation fails before any corpus is reached.

    Args:
        case: The case to run.
        client: The stub standing in for the model.
        repair: Whether a failed query gets handed back once.

    Returns:
        What ``run_case`` decided.
    """
    return nl_eval.run_case(
        case,
        cast(execute.Corpus, _StubCorpus()),
        client=cast(Any, client),
        model="stub",
        repair=repair,
    )


# Scoring


def test_the_same_reactions_in_another_order_pass():
    # Several queries are right, and they need not return the same rows in the same
    # order; the set is what a case holds a translation to.
    assert nl_eval.score(_CASE, ["ord-bb", "ord-aa"], ["ord-aa", "ord-bb"]).passed


def test_a_translation_that_missed_reactions_fails():
    result = nl_eval.score(_CASE, ["ord-aa"], ["ord-aa", "ord-bb"])
    assert not result.passed
    assert "missed 1" in result.detail
    assert "ord-bb" in result.detail


def test_a_translation_that_returned_too_many_fails():
    # The over-broad translation, which a case pinning only required reactions passes:
    # everything required is trivially among everything.
    result = nl_eval.score(_CASE, ["ord-aa", "ord-zz"], ["ord-aa"])
    assert not result.passed
    assert "wrongly returned 1" in result.detail
    assert "ord-zz" in result.detail


def test_a_failure_names_both_directions_at_once():
    result = nl_eval.score(_CASE, ["ord-zz"], ["ord-aa"])
    assert "missed 1" in result.detail
    assert "wrongly returned 1" in result.detail


def test_a_long_difference_is_summarized_rather_than_listed():
    result = nl_eval.score(_CASE, [], [f"ord-{index:02d}" for index in range(50)])
    assert "and 47 more" in result.detail


def test_the_report_names_the_question_and_why_the_case_exists():
    failure = nl_eval.score(_CASE, [], ["ord-aa"])
    text = nl_eval.report([failure])
    assert "0/1 passed" in text
    assert _CASE.question in text
    assert "splits in two" in text


# The shipped cases


def test_the_shipped_cases_load_and_say_why_they_exist():
    cases = nl_eval.load_cases()
    assert cases
    for case in cases:
        assert case.why


def test_every_shipped_reference_compiles():
    # A hand-written reference is as able to name a path the schema lacks as a model
    # is, and a case that cannot run is worse than no case.
    for case in nl_eval.load_cases():
        if case.reference is not None:
            query.compile_query(query.Query.model_validate({"where": case.reference}))


def test_a_case_the_grammar_cannot_express_is_marked_as_such():
    assert any(not case.compiles for case in nl_eval.load_cases())


def test_a_case_needs_a_reason_to_exist():
    with pytest.raises(ValueError, match="why"):
        nl_eval.EvalCase.model_validate(
            {"question": "anything?", "reference": {"op": "not_null", "path": "x"}}
        )


def test_a_case_that_compiles_needs_a_reference():
    # Without one there is nothing to compare a translation against, so the case would
    # pass whatever the model wrote.
    with pytest.raises(ValueError, match="needs a reference"):
        nl_eval.EvalCase.model_validate(
            {"question": "anything?", "why": "a placeholder"}
        )


def test_an_inexpressible_case_cannot_carry_a_reference():
    with pytest.raises(ValueError, match="no answer to reference"):
        nl_eval.EvalCase.model_validate(
            {
                "question": "anything?",
                "why": "a placeholder",
                "compiles": False,
                "reference": {"op": "not_null", "path": "reaction_id"},
            }
        )


def test_load_cases_reads_the_file_it_is_given(tmp_path):
    path = tmp_path / "cases.json"
    path.write_text(
        json.dumps(
            [{"question": "anything?", "why": "a placeholder", "compiles": False}]
        ),
        encoding="utf-8",
    )
    assert nl_eval.load_cases(path)[0].question == "anything?"


# How a translation can land


def test_declining_an_inexpressible_question_passes():
    result = _run(
        _INEXPRESSIBLE, _StubClient("cannot_answer", {"reason": "no column comparison"})
    )
    assert result.passed
    assert "declined" in result.detail


def test_a_query_that_does_not_compile_never_passes():
    # The layer's answer to an inexpressible question is cannot_answer. A build_query
    # call that happens not to compile reached the same verdict by the wrong road, and
    # scoring it as a refusal would hide the day the model starts answering instead.
    result = _run(_INEXPRESSIBLE, _StubClient("build_query", _BAD_PATH))
    assert not result.passed
    assert "did not compile" in result.detail


def test_a_refusal_reached_only_after_a_failed_query_fails():
    # The caller is served either way, but this case exists to measure whether the
    # model reads the question, not whether it eventually stops.
    # What marks it is the attempt translate recorded before the refusal, so the two
    # turns have to be run rather than an error handed in.
    client = _SequenceClient(
        ("build_query", _BAD_PATH),
        ("cannot_answer", {"reason": "no column comparison"}),
    )
    result = _run(_INEXPRESSIBLE, client, repair=True)
    assert not result.passed
    assert "built a query, then declined" in result.detail


def test_an_expressible_question_that_does_not_compile_fails():
    assert not _run(_CASE, _StubClient("build_query", _BAD_PATH)).passed


def test_declining_an_expressible_question_fails():
    assert not _run(_CASE, _StubClient("cannot_answer", {"reason": "gave up"})).passed


def test_a_run_records_its_cases(tmp_path):
    # The measurement that exists today runs from a laptop, and recording it is the
    # reason the log is a file and a bucket rather than a database behind a tunnel.
    sink = nl_log.JsonlSink(tmp_path / "run.jsonl")
    nl_eval.run_case(
        _INEXPRESSIBLE,
        cast(execute.Corpus, _StubCorpus()),
        client=cast(Any, _StubClient("cannot_answer", {"reason": "no comparison"})),
        model="stub",
        repair=False,
        sink=sink,
        session_id="run-1",
    )
    (line,) = (tmp_path / "run.jsonl").read_text(encoding="utf-8").splitlines()
    event = json.loads(line)
    assert event["outcome"] == nl_log.Outcome.DECLINED
    assert event["question"] == _INEXPRESSIBLE.question
    # A run is a session, so one run's cases group the way one visitor's questions do.
    assert event["session_id"] == "run-1"


class _StubCorpus:
    """Answers every query with one reaction, standing in for a real corpus.

    Attributes:
        fingerprint: What a record names the corpus by.
    """

    fingerprint = "stubcorpus"

    def search(self, translated: Any, **kwargs: Any) -> pa.Table:
        """Returns the canned result.

        Args:
            translated: The query, which the stub ignores.
            **kwargs: Bounds the stub ignores.

        Returns:
            A one-row table.
        """
        del translated, kwargs
        return pa.table({"reaction_id": ["ord-aa"]})


class _FailingCorpus(_StubCorpus):
    """Refuses every query the way a resolver refuses an unobtainable compound."""

    def search(self, translated: Any, **kwargs: Any) -> pa.Table:
        """Raises.

        Args:
            translated: The query, which the stub ignores.
            **kwargs: Bounds the stub ignores.

        Raises:
            ValueError: Always.
        """
        del translated, kwargs
        raise ValueError("cannot resolve 'unobtainium'")


def test_a_case_whose_search_fails_costs_that_case_and_not_the_run(tmp_path):
    # main runs the cases in a comprehension, so an exception escaping run_case costs
    # every other case's result along with this one.
    sink = nl_log.JsonlSink(tmp_path / "run.jsonl")
    result = nl_eval.run_case(
        _INEXPRESSIBLE,
        cast(execute.Corpus, _FailingCorpus()),
        client=cast(
            Any,
            _StubClient("build_query", {"where": {"op": "not_null", "path": "smiles"}}),
        ),
        model="stub",
        repair=False,
        sink=sink,
    )
    assert not result.passed
    assert "unobtainium" in result.detail
    (line,) = (tmp_path / "run.jsonl").read_text(encoding="utf-8").splitlines()
    assert json.loads(line)["outcome"] == nl_log.Outcome.UNRESOLVED_COMPOUND


def test_a_query_built_for_an_inexpressible_question_is_recorded(tmp_path):
    # The case fails, but the record is the point: a model answering a question the
    # grammar supposedly cannot express means either the model is wrong or the case
    # is, and the query it wrote is the evidence for deciding which.
    sink = nl_log.JsonlSink(tmp_path / "run.jsonl")
    result = nl_eval.run_case(
        _INEXPRESSIBLE,
        cast(execute.Corpus, _StubCorpus()),
        client=cast(
            Any,
            _StubClient("build_query", {"where": {"op": "not_null", "path": "smiles"}}),
        ),
        model="stub",
        repair=False,
        sink=sink,
    )
    assert not result.passed
    assert "cannot express" in result.detail
    (line,) = (tmp_path / "run.jsonl").read_text(encoding="utf-8").splitlines()
    event = json.loads(line)
    assert event["outcome"] == nl_log.Outcome.ANSWERED
    (attempt,) = json.loads(event["attempts"])
    assert json.loads(attempt["translation"])["where"]["op"] == "not_null"
    # A case file marched through end to end is not evidence about what people ask, and
    # it lands in the same place as the questions that are.
    assert event["source"] == "eval"


def test_a_run_without_a_sink_records_nothing(tmp_path):
    result = _run(_CASE, _StubClient("build_query", _BAD_PATH))
    assert not result.passed
    assert not list(tmp_path.iterdir())


def test_a_case_is_scored_over_every_match(monkeypatch):
    # A case is scored by comparing the reactions a translation returns against the
    # reference's. Bounded, both sides are arbitrary pages of their match sets -- SQL
    # leaves an unordered LIMIT free to pick any rows -- so two queries that agree would
    # be scored as disagreeing. The corpora in this file are far under any bound, so
    # only the argument itself can say this.
    seen = {}

    class _RecordedError(Exception):
        pass

    def _corpus(*args, **kwargs):
        seen.update(kwargs)
        raise _RecordedError

    monkeypatch.setattr(nl_eval.nl, "get_client", lambda: None)
    monkeypatch.setattr(nl_eval.execute, "Corpus", _corpus)
    with pytest.raises(_RecordedError):
        nl_eval.main(["--projections=x", "--structures=y"])
    assert seen["max_rows"] is None
