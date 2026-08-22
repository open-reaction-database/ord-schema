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
# Two predicates the grammar really produces, agreeing on no part of their shape: one
# nests a list of clauses, the other a quantifier over a structure predicate.
_CONJUNCTION = nl_log.Attempt(
    translation={
        "op": "and",
        "clauses": [{"op": "gt", "path": "y", "value": {"literal": 1}}],
    },
    error=None,
)
_QUANTIFIER = nl_log.Attempt(
    translation={
        "op": "exists",
        "path": "inputs.components",
        "where": {"op": "similarity", "smiles": "c1ccccc1", "threshold": 0.4},
    },
    error=None,
)
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


def test_each_turn_is_priced_on_its_own_attempt():
    # Per-attempt usage is what prices the repair turn separately from the first try,
    # which is the open question about choosing a cheap model and repairing it.
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
    assert [attempt["usage"]["input"] for attempt in event["attempts"]] == [412, 689]


def test_the_total_is_recorded_rather_than_summed_from_the_attempts():
    # A turn can spend tokens without producing an attempt -- declining costs a call,
    # and so does writing the sentence -- so the total is stated, not derived.
    event = nl_log.event(
        _ask(
            attempts=(
                nl_log.Attempt(
                    translation=_QUERY, error=None, usage=nl_log.Usage(input=10)
                ),
            ),
            usage=nl_log.Usage(input=14, output=30),
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


def test_reading_folds_a_thumb_and_a_label_onto_their_ask(tmp_path):
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c")))
    sink.write(nl_log.event(_ask(record_id="7b21", outcome=nl_log.EMPTY)))
    sink.write(nl_log.thumb("0f3c", "down", comment="wrong solvent"))
    sink.write(nl_log.label("0f3c", "wrong", note="two quantifiers"))
    found = nl_log.read(
        str(tmp_path / "*.jsonl"), statement="SELECT record_id FROM log"
    ).fetchall()
    # Two asks in, two rows out: a thumb and a label fold onto an ask rather than
    # arriving as rows of their own.
    assert {row[0] for row in found} == {"0f3c", "7b21"}
    verdicts = nl_log.read(
        str(tmp_path / "*.jsonl"),
        statement="SELECT thumb, verdict FROM log WHERE record_id = '0f3c'",
    ).fetchall()
    assert verdicts == [("down", "wrong")]


def test_reading_derives_the_translation_and_whether_a_repair_fired(tmp_path):
    # Neither is stored, so a reader that cannot produce them makes the log unusable
    # for the question it exists to answer.
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(
        nl_log.event(
            _ask(
                record_id="0f3c",
                attempts=(
                    nl_log.Attempt(translation={"op": "bad"}, error="no such path"),
                    nl_log.Attempt(translation=_QUERY, error=None),
                ),
            )
        )
    )
    sink.write(nl_log.event(_ask(record_id="7b21")))
    rows = nl_log.read(
        str(tmp_path / "*.jsonl"),
        statement="SELECT record_id, repaired, translation FROM log ORDER BY record_id",
    ).fetchall()
    assert [(row[0], row[1]) for row in rows] == [("0f3c", True), ("7b21", False)]
    assert json.loads(rows[0][2]) == _QUERY


def test_a_translation_that_never_compiled_reads_as_no_translation(tmp_path):
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(
        nl_log.event(
            _ask(
                outcome=nl_log.MALFORMED,
                attempts=(
                    nl_log.Attempt(translation={"op": "bad"}, error="no such path"),
                ),
            )
        )
    )
    ((translation,),) = nl_log.read(
        str(tmp_path / "*.jsonl"), statement="SELECT translation FROM log"
    ).fetchall()
    assert translation is None


def test_compaction_leaves_a_month_readable_the_same_way(tmp_path):
    # Two months whose queries have different shapes must still read as one table,
    # which is what the string payload buys and what a struct would lose.
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c")))
    sink.write(
        nl_log.event(
            _ask(
                record_id="7b21",
                attempts=(
                    nl_log.Attempt(
                        translation={"op": "exists", "where": {"op": "not_null"}},
                        error=None,
                    ),
                ),
            )
        )
    )
    compacted = tmp_path / "2026-08.parquet"
    nl_log.compact(str(tmp_path / "*.jsonl"), compacted)
    rows = nl_log.read(
        str(compacted),
        statement="SELECT record_id, translation FROM log ORDER BY record_id",
    ).fetchall()
    assert [row[0] for row in rows] == ["0f3c", "7b21"]
    assert json.loads(rows[1][1])["where"] == {"op": "not_null"}


def _month(tmp_path, name, attempt, record_id):
    """Writes one record, compacts it as a month of its own, and returns the file."""
    directory = tmp_path / name
    nl_log.JsonlSink(directory / "log.jsonl").write(
        nl_log.event(_ask(record_id=record_id, attempts=(attempt,)))
    )
    compacted = tmp_path / f"{name}.parquet"
    nl_log.compact(str(directory / "*.jsonl"), compacted)
    return compacted


def test_two_compacted_months_read_as_one_table(tmp_path):
    # Each month's Parquet bakes in a schema, and these two predicates share no shape
    # -- one nests a list of clauses, the other a quantifier over a structure
    # predicate. As strings the two months reconcile by construction, so every field of
    # both is reachable. Stored as structs it depends on what DuckDB can coerce, which
    # is measured rather than guessed: see the ord-logbook entry's probe for shapes
    # where the read either loses a field or refuses.
    july = _month(tmp_path, "2026-07", _CONJUNCTION, "0f3c")
    august = _month(tmp_path, "2026-08", _QUANTIFIER, "7b21")
    rows = nl_log.read(
        str(july),
        str(august),
        statement="SELECT record_id, translation FROM log ORDER BY record_id",
    ).fetchall()
    assert [row[0] for row in rows] == ["0f3c", "7b21"]
    assert json.loads(rows[0][1])["clauses"][0]["path"] == "y"
    assert json.loads(rows[1][1])["where"]["threshold"] == 0.4


def test_a_compacted_month_and_a_loose_day_read_as_one_table(tmp_path):
    # Compaction has to be able to arrive late, so a reader that could see only one
    # format would force a migration the day it did.
    compacted = _month(tmp_path, "2026-07", _CONJUNCTION, "0f3c")
    nl_log.JsonlSink(tmp_path / "today" / "log.jsonl").write(
        nl_log.event(_ask(record_id="7b21", attempts=(_QUANTIFIER,)))
    )
    rows = nl_log.read(
        str(compacted),
        str(tmp_path / "today" / "*.jsonl"),
        statement="SELECT record_id, translation FROM log ORDER BY record_id",
    ).fetchall()
    assert [row[0] for row in rows] == ["0f3c", "7b21"]
    assert json.loads(rows[1][1])["where"]["threshold"] == 0.4


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
