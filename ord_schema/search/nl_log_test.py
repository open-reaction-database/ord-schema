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
import uuid

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
    outcome=nl_log.Outcome.ANSWERED,
    attempts=(nl_log.Attempt(translation=_QUERY, error=None),),
    row_count=3,
)


def _ask(**kwargs) -> nl_log.Ask:
    """Returns an answered ask with the fields a test cares about replaced."""
    return dataclasses.replace(_ANSWERED, **kwargs)


def test_every_model_authored_payload_is_written_as_a_string():
    # Compacting a month of records into Parquet infers a struct from whatever shapes
    # that month held, and a field of a later month then reads back as missing. That
    # covers the attempts list as much as the query inside it: a month where everything
    # compiled has "error": null throughout and infers the field as JSON.
    event = nl_log.event(_ask())
    assert isinstance(event["attempts"], str)
    (attempt,) = json.loads(event["attempts"])
    assert isinstance(attempt["translation"], str)
    assert json.loads(attempt["translation"]) == _QUERY


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
    attempts = json.loads(event["attempts"])
    assert [attempt["usage"]["input"] for attempt in attempts] == [412, 689]


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


def test_an_outcome_survives_the_round_trip_as_a_plain_string():
    # The taxonomy is an enumeration so the members and the list of them cannot drift
    # apart, but what reaches the log has to stay a JSON string: a reader querying it
    # holds a string, and so does every record written before a member was added.
    event = nl_log.event(_ask(outcome=nl_log.Outcome.EMPTY))
    assert json.loads(json.dumps(event))["outcome"] == "empty"
    assert event["outcome"] == nl_log.Outcome.EMPTY


def test_every_outcome_is_reachable_by_iterating_the_enumeration():
    assert nl_log.Outcome.ANSWERED in list(nl_log.Outcome)
    assert len(set(nl_log.Outcome)) == len(list(nl_log.Outcome))


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
            outcome=nl_log.Outcome.MALFORMED,
            attempts=(nl_log.Attempt(translation={}, error="no such path"),),
            row_count=None,
            error="did not compile",
            translate_ms=1840.0,
        )
    )
    assert answered.keys() == failed.keys()


def test_a_stage_that_never_ran_is_a_null_rather_than_a_missing_key(tmp_path):
    event = nl_log.event(_ask(outcome=nl_log.Outcome.MALFORMED, translate_ms=12.5))
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


def test_an_object_sink_partitions_by_the_events_own_date(fake_store):
    # The partition comes from the record rather than the clock, so an event written
    # late still lands in the day it happened.
    sink = nl_log.S3Sink(fake_store, bucket="ord-internal", prefix="nl-log")
    written = nl_log.event(
        _ask(record_id="0f3c", timestamp="2026-08-22T18:04:11+00:00")
    )
    sink.write(written)
    (key,) = fake_store.objects
    assert key == f"nl-log/dt=2026-08-22/{written['event_id']}.json"
    assert fake_store.bucket == "ord-internal"


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
    sink.write(nl_log.event(_ask(record_id="7b21", outcome=nl_log.Outcome.EMPTY)))
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


def test_a_thumb_does_not_overwrite_the_ask_it_judges(fake_store):
    # Every event on one question shares its record_id, so a key built from that alone
    # collides: the thumb lands on the ask's object and the question is gone.
    sink = nl_log.S3Sink(fake_store, bucket="ord-internal")
    sink.write(nl_log.event(_ask(record_id="0f3c")))
    sink.write(nl_log.thumb("0f3c", "down"))
    assert len(fake_store.objects) == 2
    kinds = {json.loads(body)["event"] for body in fake_store.objects.values()}
    assert kinds == {"ask", "thumb"}


def test_two_thumbs_on_one_question_are_both_kept(fake_store):
    # An append-only log keeps both and lets the reader decide; overwriting one with
    # the other loses the history the log exists to hold.
    sink = nl_log.S3Sink(fake_store, bucket="ord-internal")
    sink.write(nl_log.thumb("0f3c", "up"))
    sink.write(nl_log.thumb("0f3c", "down"))
    assert len(fake_store.objects) == 2


def test_the_latest_verdict_wins_rather_than_the_lexically_largest(tmp_path):
    # "up" sorts after "down", so a reader taking the maximum reports the wrong verdict
    # whenever someone changes their mind from up to down.
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c")))
    sink.write({**nl_log.thumb("0f3c", "up"), "timestamp": "2026-08-22T09:00:00+00:00"})
    sink.write(
        {**nl_log.thumb("0f3c", "down"), "timestamp": "2026-08-22T10:00:00+00:00"}
    )
    sink.write(
        {
            **nl_log.label("0f3c", "wrong", note="first look", promoted=True),
            "timestamp": "2026-08-22T09:00:00+00:00",
        }
    )
    sink.write(
        {
            **nl_log.label("0f3c", "correct", note="second look"),
            "timestamp": "2026-08-22T10:00:00+00:00",
        }
    )
    ((thumb, verdict, note, promoted),) = nl_log.read(
        str(tmp_path / "*.jsonl"),
        statement="SELECT thumb, verdict, note, promoted FROM log",
    ).fetchall()
    assert thumb == "down"
    assert (verdict, note) == ("correct", "second look")
    # The retraction has to win too, or a promotion can never be taken back.
    assert promoted is False


def test_a_later_verdict_does_not_inherit_the_field_it_left_empty(tmp_path):
    # arg_max skips a row whose value is NULL, so per-field aggregates answer from
    # different events the moment one omits an optional field: the comment left beside
    # an earlier thumb would outlive the thumb itself, and the row would report a
    # verdict nobody gave.
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c")))
    sink.write(
        {
            **nl_log.thumb("0f3c", "up", comment="exactly right"),
            "timestamp": "2026-08-22T09:00:00+00:00",
        }
    )
    sink.write(
        {**nl_log.thumb("0f3c", "down"), "timestamp": "2026-08-22T10:00:00+00:00"}
    )
    sink.write(
        {
            **nl_log.label("0f3c", "wrong", note="two quantifiers", reference=_QUERY),
            "timestamp": "2026-08-22T09:00:00+00:00",
        }
    )
    sink.write(
        {**nl_log.label("0f3c", "correct"), "timestamp": "2026-08-22T10:00:00+00:00"}
    )
    ((thumb, comment, verdict, note, reference),) = nl_log.read(
        str(tmp_path / "*.jsonl"),
        statement="SELECT thumb, thumb_comment, verdict, note, reference FROM log",
    ).fetchall()
    assert (thumb, comment) == ("down", None)
    assert (verdict, note, reference) == ("correct", None, None)


def test_a_month_of_clean_attempts_reads_beside_one_holding_a_rejection(tmp_path):
    # The whole point of the string payload: a file where every attempt compiled writes
    # "error": null throughout, and an inferred schema calls that field JSON where a
    # file holding one rejection calls it text. The two then refuse to union, and a
    # compacted month has the wrong type in it for good.
    clean = nl_log.JsonlSink(tmp_path / "clean.jsonl")
    clean.write(nl_log.event(_ask(record_id="0f3c")))
    rejected = nl_log.JsonlSink(tmp_path / "rejected.jsonl")
    rejected.write(
        nl_log.event(
            _ask(
                record_id="7b21",
                attempts=(
                    nl_log.Attempt(translation={"op": "eq"}, error="no such path"),
                    nl_log.Attempt(translation=_QUERY, error=None),
                ),
            )
        )
    )
    found = nl_log.read(
        str(tmp_path / "clean.jsonl"),
        str(tmp_path / "rejected.jsonl"),
        statement="SELECT record_id, repaired FROM log ORDER BY record_id",
    ).fetchall()
    assert found == [("0f3c", False), ("7b21", True)]


def test_a_log_of_uuids_reads_beside_one_of_opaque_tokens(tmp_path):
    # A client mints its own session_id. All-dashed-UUID infers as UUID and anything
    # else infers as VARCHAR, and a read spanning both refuses to union them.
    uuids = nl_log.JsonlSink(tmp_path / "uuids.jsonl")
    uuids.write(nl_log.event(_ask(session_id=str(uuid.uuid4()))))
    tokens = nl_log.JsonlSink(tmp_path / "tokens.jsonl")
    tokens.write(nl_log.event(_ask(session_id="sess-abc123")))
    sessions = nl_log.read(
        str(tmp_path / "uuids.jsonl"),
        str(tmp_path / "tokens.jsonl"),
        statement="SELECT session_id FROM log",
    ).fetchall()
    assert "sess-abc123" in {row[0] for row in sessions}


def test_an_event_read_through_two_sources_counts_once(tmp_path):
    # Compaction leaves the objects it folded in place, so the documented glob pair
    # reads a compacted month beside the raw days inside it. Counting an ask twice
    # doubles every rate the log exists to report.
    raw = tmp_path / "raw"
    sink = nl_log.JsonlSink(raw / "day.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c")))
    folded = tmp_path / "2026-08.parquet"
    assert nl_log.compact(str(raw / "*.jsonl"), folded) == 1
    rows = nl_log.read(
        str(folded), str(raw / "*.jsonl"), statement="SELECT record_id FROM log"
    ).fetchall()
    assert rows == [("0f3c",)]


def test_redaction_empties_what_a_person_typed_not_only_the_question(tmp_path):
    # A thumb comment and a label note are free text the same way a question is, and
    # emptying only the columns an ask names would leave both of them.
    sink = nl_log.JsonlSink(tmp_path / "day.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c")))
    nl_log.write(sink, nl_log.thumb("0f3c", "down", comment="reader comment text"))
    nl_log.write(sink, nl_log.label("0f3c", "wrong", note="maintainer note text"))
    folded = tmp_path / "2026-08.parquet"
    nl_log.compact(str(tmp_path / "*.jsonl"), folded, redact=True)
    written = folded.read_bytes()
    assert b"reader comment text" not in written
    assert b"maintainer note text" not in written
    # The outcome and the timings still have to be there, or redaction costs the
    # measurement it was meant to keep.
    outcomes = nl_log.read(str(folded), statement="SELECT outcome FROM log").fetchall()
    assert outcomes == [("answered",)]


def test_a_day_holding_only_feedback_still_compacts(tmp_path):
    # A real partition for a deployment where somebody left a thumb on a day nobody
    # asked anything, so the ask columns redaction empties are not in it.
    sink = nl_log.JsonlSink(tmp_path / "day.jsonl")
    nl_log.write(sink, nl_log.thumb("0f3c", "down", comment="wrong solvent"))
    folded = tmp_path / "2026-08.parquet"
    assert nl_log.compact(str(tmp_path / "*.jsonl"), folded, redact=True) == 1
    assert b"wrong solvent" not in folded.read_bytes()


def test_reading_nothing_at_all_says_so(tmp_path):
    with pytest.raises(ValueError, match="at least one glob"):
        nl_log.read(statement="SELECT * FROM log")


@pytest.mark.parametrize(
    ("build", "value"),
    [(nl_log.thumb, "Down"), (nl_log.label, "incorrect")],
)
def test_a_verdict_outside_the_vocabulary_is_refused(build, value):
    # The log is append-only, so a typo can be superseded but never corrected, and a
    # verdict spelled "Down" is invisible to every query that asks for "down".
    with pytest.raises(ValueError, match=value):
        build("0f3c", value)


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
                outcome=nl_log.Outcome.MALFORMED,
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


def test_compaction_can_drop_the_free_text_and_keep_everything_else(tmp_path):
    # The free text goes and everything the analysis runs on stays, so a redacted file
    # is still worth reading.
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(
        nl_log.event(
            _ask(record_id="0f3c", question="who used pyridine", answer_text="Two.")
        )
    )
    compacted = tmp_path / "2026-08.parquet"
    nl_log.compact(str(tmp_path / "*.jsonl"), compacted, redact=True)
    ((question, answer, outcome, row_count, translation),) = nl_log.read(
        str(compacted),
        statement=(
            "SELECT question, answer_text, outcome, row_count, translation FROM log"
        ),
    ).fetchall()
    assert question is None
    assert answer is None
    assert outcome == nl_log.Outcome.ANSWERED
    assert row_count == 3
    # The translation stays: it is what reproduces the rows the log does not store.
    assert json.loads(translation) == _QUERY
    assert "pyridine" not in compacted.read_bytes().decode("utf-8", "replace")


def test_a_redacted_month_reads_beside_an_unredacted_one(tmp_path):
    # Redaction empties the columns rather than dropping them, so the two tiers still
    # have one schema and a query spanning them binds.
    sink = nl_log.JsonlSink(tmp_path / "old" / "log.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c", question="an old question")))
    old = tmp_path / "2026-07.parquet"
    nl_log.compact(str(tmp_path / "old" / "*.jsonl"), old, redact=True)
    recent = nl_log.JsonlSink(tmp_path / "new" / "log.jsonl")
    recent.write(nl_log.event(_ask(record_id="7b21", question="a recent question")))
    new = tmp_path / "2026-08.parquet"
    nl_log.compact(str(tmp_path / "new" / "*.jsonl"), new)
    rows = nl_log.read(
        str(old),
        str(new),
        statement="SELECT record_id, question FROM log ORDER BY record_id",
    ).fetchall()
    assert rows == [("0f3c", None), ("7b21", "a recent question")]


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


def test_an_ask_records_which_surface_it_arrived_from():
    # A served endpoint, an eval run, and a script all write to one prefix, and the
    # published PubMed log shows what it costs to leave them mixed: hand-typed queries,
    # bare accession lookups, and tool-generated ones sit side by side there with no
    # way to tell them apart afterwards.
    assert nl_log.event(_ask(source="web"))["source"] == "web"


def test_a_run_and_a_reader_are_told_apart_by_source(tmp_path):
    sink = nl_log.JsonlSink(tmp_path / "log.jsonl")
    sink.write(nl_log.event(_ask(record_id="0f3c", source="web")))
    sink.write(nl_log.event(_ask(record_id="7b21", source="eval")))
    # Nothing states it where nothing supplies it, rather than a default that would
    # read as a claim about where the question came from.
    sink.write(nl_log.event(_ask(record_id="c40a")))
    found = nl_log.read(
        str(tmp_path / "*.jsonl"), statement="SELECT record_id, source FROM log"
    ).fetchall()
    assert dict(found) == {"0f3c": "web", "7b21": "eval", "c40a": None}
