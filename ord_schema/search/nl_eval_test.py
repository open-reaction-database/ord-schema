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

These cover the scoring, and the outcomes ``run_case`` maps a translation onto, with a
stub client standing in for the model. Measuring against a real model and a real corpus
is a command rather than a test.
"""

import json
import types
from typing import Any, cast

import pytest
from anthropic.types import ToolUseBlock

from ord_schema.search import execute, nl_eval, nl_log

_CASE = nl_eval.EvalCase(
    question="which reactions use pyridine as a solvent?",
    why="two conditions on one element, which a wrong translation splits in two",
    must_return=["ord-aa", "ord-bb"],
    must_not_return=["ord-zz"],
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


class _RaisingClient:
    """Raises a canned error instead of answering.

    Attributes:
        messages: Itself, so ``client.messages.create`` reaches ``create``.
    """

    def __init__(self, error: Exception) -> None:
        """Initializes the stub.

        Args:
            error: What every request raises.
        """
        self._error = error
        self.messages = self

    def create(self, **kwargs: Any) -> Any:
        """Raises the canned error.

        Args:
            **kwargs: The request, which the stub ignores.

        Raises:
            Exception: Always, with whatever this stub was built around.
        """
        del kwargs
        raise self._error


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
        cast(execute.Corpus, None),
        client=cast(Any, client),
        model="stub",
        repair=repair,
    )


def test_a_query_returning_what_it_must_passes():
    assert nl_eval.score(_CASE, ["ord-aa", "ord-bb", "ord-cc"]).passed


def test_a_query_missing_a_required_reaction_fails():
    result = nl_eval.score(_CASE, ["ord-aa"])
    assert not result.passed
    assert "ord-bb" in result.detail


def test_a_query_returning_a_forbidden_reaction_fails():
    # The near-miss reaction is the whole point: a translation that finds pyridine in
    # one component and a solvent in another is wrong rather than differently spelled.
    result = nl_eval.score(_CASE, ["ord-aa", "ord-bb", "ord-zz"])
    assert not result.passed
    assert "ord-zz" in result.detail


def test_scoring_does_not_care_about_order_or_extras():
    # Several queries are right, and they need not return the same rows in the same
    # order; only the reactions named either way are pinned.
    assert nl_eval.score(_CASE, ["ord-cc", "ord-bb", "ord-dd", "ord-aa"]).passed


def test_the_report_names_the_question_and_why_the_case_exists():
    failure = nl_eval.score(_CASE, [])
    text = nl_eval.report([failure])
    assert "0/1 passed" in text
    assert _CASE.question in text
    assert "splits in two" in text


def test_the_shipped_cases_load_and_say_why_they_exist():
    cases = nl_eval.load_cases()
    assert cases
    for case in cases:
        assert case.why


def test_every_expressible_case_pins_both_sides():
    # must_return alone is passed by a translation that drops the predicate and hands
    # back the corpus, since everything required is trivially among everything. The
    # counterexamples are what make an over-broad query fail.
    for case in nl_eval.load_cases():
        if case.compiles:
            assert case.must_return
            assert case.must_not_return


def test_a_case_the_grammar_cannot_express_is_marked_as_such():
    cases = nl_eval.load_cases()
    assert any(not case.compiles for case in cases)


def test_load_cases_reads_the_file_it_is_given(tmp_path):
    path = tmp_path / "cases.json"
    path.write_text(
        json.dumps(
            [{"question": "anything?", "why": "a placeholder", "compiles": False}]
        ),
        encoding="utf-8",
    )
    assert nl_eval.load_cases(path)[0].question == "anything?"


def test_a_case_needs_a_reason_to_exist():
    with pytest.raises(ValueError, match="why"):
        nl_eval.EvalCase.model_validate({"question": "anything?"})


def test_declining_an_inexpressible_question_passes():
    result = _run(
        _INEXPRESSIBLE, _StubClient("cannot_answer", {"reason": "no column comparison"})
    )
    assert result.passed
    assert "declined" in result.detail


def test_a_refusal_reached_only_after_a_failed_query_fails():
    # The caller is served either way, but this case exists to measure whether the
    # model reads the question, not whether it eventually stops.
    # The model builds a query, the compiler refuses it, and only on the repair turn
    # does the model decline. What marks that is the attempt translate recorded before
    # the refusal, so the sequence has to be run rather than an error handed in.
    client = _SequenceClient(
        ("build_query", _BAD_PATH),
        ("cannot_answer", {"reason": "no column comparison"}),
    )
    result = _run(_INEXPRESSIBLE, client, repair=True)
    assert not result.passed
    assert "built a query, then declined" in result.detail


def test_a_query_that_does_not_compile_never_passes():
    # The layer's answer to an inexpressible question is cannot_answer. A build_query
    # call that happens not to compile reached the same verdict by the wrong road, and
    # scoring it as a refusal would hide the day the model starts answering instead.
    result = _run(_INEXPRESSIBLE, _StubClient("build_query", _BAD_PATH))
    assert not result.passed
    assert "did not compile" in result.detail


def test_an_expressible_question_that_does_not_compile_fails():
    result = _run(_CASE, _StubClient("build_query", _BAD_PATH))
    assert not result.passed


def test_declining_an_expressible_question_fails():
    result = _run(_CASE, _StubClient("cannot_answer", {"reason": "gave up"}))
    assert not result.passed


def test_a_run_records_its_cases(tmp_path):
    # The measurement that exists today runs from a laptop, and recording it is the
    # reason the log is a file and a bucket rather than a database behind a tunnel.
    sink = nl_log.JsonlSink(tmp_path / "run.jsonl")
    nl_eval.run_case(
        _INEXPRESSIBLE,
        cast(execute.Corpus, None),
        client=cast(Any, _StubClient("cannot_answer", {"reason": "no comparison"})),
        model="stub",
        repair=False,
        sink=sink,
        session_id="run-1",
    )
    (line,) = (tmp_path / "run.jsonl").read_text(encoding="utf-8").splitlines()
    event = json.loads(line)
    assert event["outcome"] == nl_log.DECLINED
    assert event["question"] == _INEXPRESSIBLE.question
    # A run is a session, so one run's cases group the way one visitor's questions do.
    assert event["session_id"] == "run-1"


def test_a_run_without_a_sink_records_nothing(tmp_path):
    result = _run(_CASE, _StubClient("build_query", _BAD_PATH))
    assert not result.passed
    assert not list(tmp_path.iterdir())
