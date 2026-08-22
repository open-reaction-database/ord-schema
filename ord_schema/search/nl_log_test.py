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

"""Tests for ord_schema.search.nl_log."""

import dataclasses
import json

import pytest

from ord_schema.search import nl_log

_QUERY = {"op": "eq", "path": "reaction_id", "value": {"literal": "ord-1"}}
_ANSWERED = nl_log.Ask(
    question="which reactions use pyridine as a solvent?",
    model="claude-haiku-4-5",
    corpus_fingerprint="b42f19ac",
    outcome=nl_log.ANSWERED,
    attempts=(nl_log.Attempt(translation=_QUERY, error=None),),
    row_count=3,
)


def _ask(**kwargs) -> nl_log.Ask:
    """Returns an answered ask with the fields a test cares about replaced."""
    return dataclasses.replace(_ANSWERED, **kwargs)


def test_a_translation_is_written_as_a_string_not_a_struct():
    # Compacting a month of records into Parquet infers a struct from whatever query
    # shapes that month held, and a later month's fields then read back as missing.
    event = nl_log.event(_ask())
    translation = event["attempts"][0]["translation"]
    assert isinstance(translation, str)
    assert json.loads(translation) == _QUERY


def test_usage_totals_the_repair_turn_onto_the_first_try():
    # What a question cost is both turns, which is the number deciding whether a cheap
    # model plus a repair is worth its accuracy.
    event = nl_log.event(
        _ask(
            attempts=(
                nl_log.Attempt(
                    translation={"op": "eq"},
                    error="no such path",
                    usage=nl_log.Usage(input=412, output=233, cache_read=15104),
                ),
                nl_log.Attempt(
                    translation=_QUERY,
                    error=None,
                    usage=nl_log.Usage(input=689, output=241, cache_read=15104),
                ),
            )
        )
    )
    assert event["usage"] == {
        "input": 1101,
        "output": 474,
        "cache_read": 30208,
        "cache_creation": 0,
    }


def test_the_prose_call_counts_toward_what_the_question_cost():
    # Two model calls answer a question: the translation and the sentence describing
    # the result. The total is what the question cost, and it is the one number here
    # that the attempts alone cannot give.
    event = nl_log.event(
        _ask(
            attempts=(
                nl_log.Attempt(
                    translation=_QUERY, error=None, usage=nl_log.Usage(input=10)
                ),
            ),
            answer_usage=nl_log.Usage(input=4, output=30),
        )
    )
    assert event["usage"]["input"] == 14
    assert event["usage"]["output"] == 30


def test_nothing_derivable_from_the_attempts_is_stored():
    # A stored copy of what the attempts already say is a second source of truth, and
    # the reader derives both of these without one.
    event = nl_log.event(_ask())
    assert "repaired" not in event
    assert "translation" not in event
    assert "declined_attempted" not in event


def test_every_outcome_writes_the_same_keys():
    # The envelope is what compaction infers a Parquet schema from, so a field that
    # appears only on some outcomes would give two months two schemas -- the failure
    # the string payload exists to avoid, reintroduced by the envelope.
    answered = nl_log.event(_ask())
    failed = nl_log.event(
        _ask(
            outcome=nl_log.MALFORMED,
            attempts=(nl_log.Attempt(translation={}, error="no such path"),),
            row_count=None,
            error="did not compile",
            translate_ms=1840.0,
        )
    )
    assert answered.keys() == failed.keys()


def test_a_stage_that_never_ran_is_a_null_rather_than_a_missing_key(tmp_path):
    event = nl_log.event(_ask(outcome=nl_log.MALFORMED, translate_ms=12.5))
    assert event["translate_ms"] == 12.5
    assert event["search_ms"] is None
    assert event["answer_ms"] is None


def test_a_file_sink_appends_one_line_per_event(tmp_path):
    path = tmp_path / "log.jsonl"
    sink = nl_log.JsonlSink(path)
    sink.write(nl_log.event(_ask(question="first")))
    sink.write(nl_log.event(_ask(question="second")))
    lines = path.read_text(encoding="utf-8").splitlines()
    assert [json.loads(line)["question"] for line in lines] == ["first", "second"]


def test_an_object_sink_keys_by_date_and_record(fake_store):
    sink = nl_log.S3Sink(fake_store, bucket="ord-internal", prefix="nl-log")
    sink.write(
        nl_log.event(_ask(record_id="0f3c", timestamp="2026-08-22T18:04:11+00:00"))
    )
    (key,) = fake_store.objects
    assert key == "nl-log/dt=2026-08-22/0f3c.json"


def test_a_failing_sink_never_fails_the_question(caplog):
    # The log is an optimization, never a dependency: an unreachable bucket must cost
    # the record rather than the answer the reader is waiting for.
    class _Broken:
        def write(self, event):
            raise RuntimeError("bucket unreachable")

    nl_log.emit(_Broken(), _ask())
    assert "bucket unreachable" in caplog.text


def test_emitting_to_no_sink_at_all_is_allowed():
    assert nl_log.emit(None, _ask()) is None


@pytest.fixture(name="fake_store")
def _fake_store():
    """An object store recording what was put, standing in for a boto3 S3 client."""

    class _Store:
        def __init__(self):
            self.objects: dict[str, bytes] = {}

        def put_object(self, *, Bucket: str, Key: str, Body: bytes) -> None:
            self.bucket = Bucket
            self.objects[Key] = Body

    return _Store()
