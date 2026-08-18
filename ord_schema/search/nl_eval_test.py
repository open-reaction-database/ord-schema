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

The scoring is what these cover. Running a case reaches a model and a corpus, which is a
command rather than a test.
"""

import json

import pytest

from ord_schema.search import nl_eval

_CASE = nl_eval.EvalCase(
    question="which reactions use pyridine as a solvent?",
    why="two conditions on one element, which a wrong translation splits in two",
    must_return=["ord-aa", "ord-bb"],
    must_not_return=["ord-zz"],
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
        # A case with neither an expectation nor a refusal to make would pass whatever
        # the model wrote.
        assert case.must_return or case.must_not_return or not case.compiles


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
