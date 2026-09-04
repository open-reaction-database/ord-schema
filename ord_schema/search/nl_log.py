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

"""A durable record of every question, the query it became, and how that went.

The questions people ask are the only material that can refine a translator, and a
prompt, a grammar, and an eval set are guesses until real ones land against them. This
module is where a question goes after it has been answered.

The log is an append-only stream of events rather than a table of rows. Three kinds
share a ``record_id``: the **ask**, written here; the **thumb** a reader leaves; the
**label** a maintainer applies. They are three assertions, by three parties, at three
times, and folding them onto one mutable row would lose which was which on the day a
label contradicts a thumb.

A record is a *typed envelope around an opaque payload*. Identifiers, outcome, usage,
timings, and fingerprints are fields; the question and the whole ``attempts`` list are
strings. They stay strings because a model-authored ``Query`` is recursive and
differently shaped record to record, and because a month in which everything compiled
carries ``"error": null`` throughout, where a month holding one rejection carries text.
Either way, compacting into Parquet would leave DuckDB to infer a struct from whichever
shapes a month happened to hold, and whether two such months then read back together
depends on what it can coerce: measured cases lose a field silently -- the rows still
count in every aggregate while the part worth recording is unreachable -- or refuse the
read outright. As strings the months reconcile by construction, and no inference is
involved. The measurements are in the ord-logbook entry "Where the question log lives".

Results are not recorded. A translation and a corpus fingerprint reproduce them, and a
question that returned a hundred thousand reactions should not cost a hundred thousand
rows of log.
"""

import contextlib
import dataclasses
import enum
import json
import os
import pathlib
import sys
import uuid
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from typing import Any, Protocol

import duckdb

from ord_schema.logging import get_logger
from ord_schema.search.execute import sql_string

logger = get_logger(__name__)


def new_identifier() -> str:
    """Returns a fresh identifier, as hex rather than as a dashed UUID.

    DuckDB infers a column of dashed UUIDs as ``UUID`` and one of anything else as
    ``VARCHAR``, and a read spanning both refuses to union them. Hex is 32 characters
    of the same randomness that no source can infer as anything but text.
    """
    return uuid.uuid4().hex


class Thumb(enum.StrEnum):
    """A reader's verdict on an answer."""

    UP = "up"
    DOWN = "down"


class Verdict(enum.StrEnum):
    """A maintainer's verdict on a translation."""

    CORRECT = "correct"
    WRONG = "wrong"
    UNCLEAR = "unclear"


class Outcome(enum.StrEnum):
    """What became of a question.

    The first two are the ordinary endings and the rest are the error taxonomy ``nl``
    raises, one outcome each, so a failure is queryable rather than a string someone
    has to grep. Members are strings, so a record serializes to ordinary JSON and a
    reader compares against ordinary text.
    """

    ANSWERED = "answered"
    EMPTY = "empty"
    DECLINED = "declined"
    MALFORMED = "malformed"
    RATE_LIMITED = "rate_limited"
    UNAVAILABLE = "unavailable"
    UNRESOLVED_COMPOUND = "unresolved_compound"
    TIMEOUT = "timeout"
    SEARCH_FAILED = "search_failed"


@dataclasses.dataclass(frozen=True)
class Usage:
    """What one call to the model cost, in tokens.

    Attributes:
        input: Uncached input tokens.
        output: Generated tokens.
        cache_read: Input tokens served from the cached prefix, which is most of a
            query's input and a fraction of its price.
        cache_creation: Input tokens written into the cache.
    """

    input: int = 0
    output: int = 0
    cache_read: int = 0
    cache_creation: int = 0

    def __add__(self, other: "Usage") -> "Usage":
        """Returns the total of two calls, for summing a repair turn onto its first."""
        return Usage(
            input=self.input + other.input,
            output=self.output + other.output,
            cache_read=self.cache_read + other.cache_read,
            cache_creation=self.cache_creation + other.cache_creation,
        )


@dataclasses.dataclass(frozen=True)
class Attempt:
    """One ``build_query`` call, and what the compiler said about it.

    Attributes:
        translation: The tool input, coerced but not necessarily valid. It is held as
            raw JSON rather than as a ``Query`` because the attempt worth recording is
            usually the one that failed to validate, where there is no ``Query`` to
            hold.
        error: Why it was rejected, or None where it was not.
        usage: What this turn cost.
    """

    translation: dict[str, Any]
    error: str | None
    usage: Usage = dataclasses.field(default_factory=Usage)


@dataclasses.dataclass(frozen=True)
class Ask:
    """One question, and what answering it did.

    Attributes:
        question: The question, in English.
        model: Which model translated.
        corpus_fingerprint: The corpus the query ran against, so the record stays
            reproducible without holding the rows.
        prompt_fingerprint: Which translator wrote this, so a record from months ago
            can be told apart from one written by today's prompt.
        outcome: How the question ended.
        attempts: Every ``build_query`` call, in order. Empty where the model declined
            without writing a query, which is how "read the question and said no" is
            told apart from "backed off once the compiler refused a guess".
        row_count: How many reactions came back, where the search ran.
        declined_reason: What the model said the grammar lacks.
        error: The failure, where one ended the ask.
        answer_text: The prose the reader saw, which is what a thumb is a verdict on.
        usage: What the whole question cost, every model call included. More than the
            attempts add up to, because a turn can spend tokens without producing one:
            declining costs a call, and so does writing the sentence describing the
            result.
        translate_ms: Wall time turning the question into a query.
        search_ms: Wall time running it, or None where it never ran.
        answer_ms: Wall time writing the prose, or None where none was written. A
            stage that did not run is None rather than an absent field, so every
            record carries the same keys whatever became of the question -- which is
            what compaction infers a Parquet schema from.
        source: Which surface the question arrived from -- ``"web"``, ``"eval"``,
            ``"api"``. A served endpoint, an eval run, and a script all write to one
            prefix, and nothing else in the record tells them apart: a rephrasing after
            a bad answer and a harness marching through a case file look alike once
            they are mixed. Low cardinality on purpose, so it groups.
        session_id: Groups a visit, so a rephrasing reads as a correction of what came
            before it. Minted by the client, which decides what it spans.
        record_id: What a thumb and a label reference.
        timestamp: When the question was asked, in UTC.
    """

    question: str
    model: str
    corpus_fingerprint: str
    outcome: Outcome
    prompt_fingerprint: str = ""
    attempts: Sequence[Attempt] = ()
    row_count: int | None = None
    declined_reason: str | None = None
    error: str | None = None
    answer_text: str | None = None
    usage: Usage = dataclasses.field(default_factory=Usage)
    translate_ms: float | None = None
    search_ms: float | None = None
    answer_ms: float | None = None
    source: str | None = None
    session_id: str | None = None
    record_id: str = dataclasses.field(default_factory=new_identifier)
    timestamp: str = dataclasses.field(
        default_factory=lambda: datetime.now(UTC).isoformat()
    )


def event(ask: Ask) -> dict[str, Any]:
    """Returns the JSON-ready ask event.

    Args:
        ask: What happened.

    Returns:
        A dict of JSON primitives, the model-authored payloads among them serialized
        to strings. It holds nothing ``read`` can derive: the translation is the last
        attempt no error rejected, and a repair fired when there is more than one.

        ``attempts`` is one string rather than a list of structs for the same reason a
        translation is: a month in which every attempt compiled writes ``"error": null``
        throughout, and DuckDB infers the field as JSON where a month holding one
        rejection infers it as text. The two then refuse to read together, and a
        compacted month has the wrong type baked in for good.
    """
    return {
        "event": "ask",
        "event_id": new_identifier(),
        "record_id": ask.record_id,
        "source": ask.source,
        "session_id": ask.session_id,
        "timestamp": ask.timestamp,
        "question": ask.question,
        "attempts": json.dumps(
            [
                {
                    "translation": json.dumps(attempt.translation, sort_keys=True),
                    "error": attempt.error,
                    "usage": dataclasses.asdict(attempt.usage),
                }
                for attempt in ask.attempts
            ],
            sort_keys=True,
        ),
        "outcome": str(ask.outcome),
        "row_count": ask.row_count,
        "declined_reason": ask.declined_reason,
        "error": ask.error,
        "answer_text": ask.answer_text,
        "model": ask.model,
        "corpus_fingerprint": ask.corpus_fingerprint,
        "prompt_fingerprint": ask.prompt_fingerprint,
        "usage": dataclasses.asdict(ask.usage),
        "translate_ms": ask.translate_ms,
        "search_ms": ask.search_ms,
        "answer_ms": ask.answer_ms,
    }


class Sink(Protocol):
    """Somewhere a record goes."""

    def write(self, event: Mapping[str, Any]) -> None:
        """Records one event."""


class ObjectStore(Protocol):
    """The one operation an object sink needs, which a boto3 S3 client satisfies.

    Taking the client rather than building one keeps this package free of an AWS
    dependency: a deployment passes ``boto3.client("s3")``, and a test passes something
    that records what it was handed.
    """

    def put_object(self, *, Bucket: str, Key: str, Body: bytes) -> Any:
        """Stores one object."""


class JsonlSink:
    """Appends records to a local file, one JSON object per line.

    Needs no credentials and no network, so an eval run on a laptop is recorded, and
    the file reads back through the same query as the bucket.
    """

    def __init__(self, path: str | os.PathLike[str]) -> None:
        """Initializes the sink.

        Args:
            path: File to append to. Its parent directory is created if missing.
        """
        self._path = pathlib.Path(path)
        self._path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, event: Mapping[str, Any]) -> None:
        """Appends one event."""
        with self._path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event, sort_keys=True) + "\n")


class StdoutSink:
    """Writes records to stdout, for a container whose logs already go somewhere."""

    def write(self, event: Mapping[str, Any]) -> None:
        """Prints one event as a single line."""
        print(json.dumps(event, sort_keys=True), file=sys.stdout, flush=True)


class S3Sink:
    """Writes one object per record, partitioned by date.

    The caller supplies both the store and the bucket, so this holds no opinion about
    where a record lands and this package carries no AWS dependency.
    """

    def __init__(
        self, store: ObjectStore, *, bucket: str, prefix: str = "nl-log/raw"
    ) -> None:
        """Initializes the sink.

        Args:
            store: Something with boto3's ``put_object``.
            bucket: Bucket to write into.
            prefix: Key prefix; the date partition goes under it. Keep it disjoint from
                where compacted months are written, so a lifecycle rule can expire the
                raw days without touching the month that folded them.
        """
        self._store = store
        self._bucket = bucket
        self._prefix = prefix.rstrip("/")

    def write(self, event: Mapping[str, Any]) -> None:
        """Stores one event under ``<prefix>/dt=<date>/<event_id>.json``.

        Keyed by the event rather than by the question it concerns: an ask, the thumb
        left on it, and a maintainer's label all carry one ``record_id``, so a key
        built from that would have each overwrite the last and the question itself
        would be the one lost.

        The date comes from the event's own timestamp rather than from the clock, so a
        record lands in the partition it belongs to even when it is written late.
        """
        day = str(event["timestamp"])[:10]
        self._store.put_object(
            Bucket=self._bucket,
            Key=f"{self._prefix}/dt={day}/{event['event_id']}.json",
            Body=json.dumps(event, sort_keys=True).encode("utf-8"),
        )


def write(sink: Sink | None, event: Mapping[str, Any]) -> None:
    """Records one event, and never lets the log cost the answer.

    Every event goes through here -- the ask this module builds and the thumb or label a
    caller builds alike. A sink reaches a file or a network, and a reader is waiting on
    the question rather than on the record of it, so a failure is warned about and
    dropped: an optimization, never a dependency, which is the contract the search cache
    holds too.

    Args:
        sink: Where the record goes; None records nothing.
        event: The event to record.
    """
    if sink is None:
        return
    try:
        sink.write(event)
    except Exception as error:  # noqa: BLE001
        logger.warning("dropping a question log record: %s", error)


def emit(sink: Sink | None, ask: Ask) -> None:
    """Records an ask.

    Args:
        sink: Where the record goes; None records nothing.
        ask: What happened.
    """
    write(sink, event(ask))


def thumb(
    record_id: str,
    value: Thumb | str,
    *,
    comment: str | None = None,
    session_id: str | None = None,
) -> dict[str, Any]:
    """Returns a reader's verdict on an answer.

    Args:
        record_id: The ask this judges.
        value: A ``Thumb``, or the string spelling one.
        comment: Whatever the reader wanted to add.
        session_id: The visit it came from.

    Returns:
        A JSON-ready thumb event.

    Raises:
        ValueError: If ``value`` is not one of the members. The log is append-only, so
            a typo can be superseded but never corrected, and a verdict spelled
            ``"Down"`` is invisible to every query that asks for ``"down"``.
    """
    return {
        "event": "thumb",
        "event_id": new_identifier(),
        "record_id": record_id,
        "session_id": session_id,
        "timestamp": datetime.now(UTC).isoformat(),
        "value": str(Thumb(value)),
        "comment": comment,
    }


def label(
    record_id: str,
    verdict: Verdict | str,
    *,
    reference: Mapping[str, Any] | None = None,
    note: str | None = None,
    promoted: bool = False,
) -> dict[str, Any]:
    """Returns a maintainer's verdict on a translation.

    Args:
        record_id: The ask this judges.
        verdict: A ``Verdict``, or the string spelling one.
        reference: A query that answers the question, where the translation did not.
            This is what an eval case needs and what a thumb cannot supply.
        note: Why.
        promoted: Whether this became an eval case.

    Returns:
        A JSON-ready label event.

    Raises:
        ValueError: If ``verdict`` is not one of the members; see ``thumb``.
    """
    return {
        "event": "label",
        "event_id": new_identifier(),
        "record_id": record_id,
        "timestamp": datetime.now(UTC).isoformat(),
        "verdict": str(Verdict(verdict)),
        "reference": (
            None if reference is None else json.dumps(reference, sort_keys=True)
        ),
        "note": note,
        "promoted": promoted,
    }


# One ask per row, with the later events folded onto it and the fields the record does
# not store derived here. A thumb or a label can be left more than once, so these
# aggregate rather than join, and they aggregate by timestamp: the newest verdict is the
# one that stands, where the largest would make "up" outrank a later "down" and would
# let a promotion outlive its retraction.
#
# The whole event is carried through one aggregate as a struct rather than each of its
# fields through its own. arg_max skips a row whose value is NULL, so per-field
# aggregates answer from different events the moment one of them omits an optional
# field: a reader who thumbs down without a comment would keep the comment left beside
# an earlier thumb up, and the row would be a verdict nobody gave.
#
# The empty row carrying no rows is what lets this bind against a log nobody has left
# feedback on yet: the columns a thumb and a label contribute do not exist until one has
# been written, and a reader that only works once someone has replied is no reader.
_SKELETON = """
SELECT
    NULL::VARCHAR AS event,
    NULL::VARCHAR AS event_id,
    NULL::VARCHAR AS record_id,
    NULL::VARCHAR AS session_id,
    NULL::VARCHAR AS attempts,
    NULL::VARCHAR AS question,
    NULL::VARCHAR AS answer_text,
    NULL::VARCHAR AS value,
    NULL::VARCHAR AS comment,
    NULL::VARCHAR AS verdict,
    NULL::VARCHAR AS reference,
    NULL::VARCHAR AS note,
    NULL::BOOLEAN AS promoted
WHERE FALSE
"""

# The identifiers a source can infer a type for. A file whose record_ids are all dashed
# UUIDs reads as UUID and one holding anything else reads as VARCHAR, and the union of
# the two refuses. This module mints hex so its own files never diverge; the cast is
# what covers a client that minted its own session_id as a dashed UUID.
_TEXT_COLUMNS = ("event_id", "record_id", "session_id")

_LOG = """
WITH events AS (
    SELECT * FROM ({sources})
    -- An event is immutable and carries its own identifier, so one arriving twice is
    -- the same event twice: a compacted month read beside the raw days it folded, or
    -- two globs that overlap. Counting it twice silently doubles every aggregate.
    QUALIFY row_number() OVER (PARTITION BY event_id ORDER BY timestamp) = 1
),
asks AS (SELECT * FROM events WHERE event = 'ask'),
verdicts AS (
    SELECT
        record_id,
        arg_max(
            {{'value': value, 'comment': comment}}, timestamp
        ) FILTER (event = 'thumb') AS thumb,
        arg_max(
            {{
                'verdict': verdict,
                'reference': reference,
                'note': note,
                'promoted': promoted
            }},
            timestamp
        ) FILTER (event = 'label') AS label
    FROM events
    WHERE event IN ('thumb', 'label')
    GROUP BY record_id
)
SELECT
    asks.* EXCLUDE (event, value, comment, verdict, reference, note, promoted),
    coalesce(json_array_length(asks.attempts), 0) > 1 AS repaired,
    list_last(
        list_transform(
            list_filter(
                from_json(asks.attempts, '["JSON"]'),
                attempt -> json_extract(attempt, '$.error') = 'null'::JSON
            ),
            attempt -> json_extract_string(attempt, '$.translation')
        )
    ) AS translation,
    verdicts.thumb['value'] AS thumb,
    verdicts.thumb['comment'] AS thumb_comment,
    verdicts.label['verdict'] AS verdict,
    verdicts.label['reference'] AS reference,
    verdicts.label['note'] AS note,
    coalesce(verdicts.label['promoted'], FALSE) AS promoted
FROM asks LEFT JOIN verdicts USING (record_id)
"""


def _source(pattern: str) -> str:
    """Returns the reader clause for one glob, chosen by what the glob names.

    The skeleton is unioned in per source rather than once over all of them, so every
    column exists whatever that source happens to hold -- which is what lets the
    identifier casts bind against a file of nothing but thumbs.
    """
    reader = "read_parquet" if pattern.endswith(".parquet") else "read_json_auto"
    casts = ", ".join(f"{name}::VARCHAR AS {name}" for name in _TEXT_COLUMNS)
    # S608: the reader is one of two constants, the path is quoted as a literal, and
    # the select list is built from this module's own constants.
    return (
        f"SELECT * REPLACE ({casts}) FROM ("  # noqa: S608
        f"SELECT * FROM {reader}({sql_string(pattern)}, union_by_name=true) "
        f"UNION ALL BY NAME {_SKELETON})"
    )


def read(
    *sources: str,
    statement: str = "SELECT * FROM log",
    connection: duckdb.DuckDBPyConnection | None = None,
) -> duckdb.DuckDBPyRelation:
    """Returns the log as one relation, thumbs and labels folded onto their asks.

    Args:
        *sources: Globs of JSON or Parquet, mixed freely, so a compacted month and the
            days written since read as one table.
        statement: SQL over the ``log`` view this defines.
        connection: DuckDB connection; an in-process one is opened if omitted.

    Returns:
        One row per question, carrying ``translation`` and ``repaired`` derived from
        the attempts, and whatever verdicts were left on it. An event read twice --
        through overlapping globs, or a compacted month beside the days it folded --
        counts once.

    Raises:
        ValueError: If no source was given. A caller splatting a computed list has an
            empty one to answer for, and the alternative is a parser error naming SQL
            it never wrote.
    """
    if not sources:
        raise ValueError("read needs at least one glob of JSON or Parquet to read")
    connection = connection if connection is not None else duckdb.connect()
    union = "\nUNION ALL BY NAME\n".join(_source(pattern) for pattern in sources)
    connection.execute(
        "CREATE OR REPLACE TEMP VIEW log AS " + _LOG.format(sources=union)
    )
    return connection.sql(statement)


# What ``redact`` empties: every column holding free text. A reader's thumb comment and
# a maintainer's label note are as free-form as the question is, and a column named on
# an ask is not the only place text arrives. The translation is deliberately not here:
# it is what reproduces the rows the log declines to store, and it is bounded and
# structured where free text is neither.
_FREE_TEXT = ("question", "answer_text", "comment", "note")


def compact(
    pattern: str, output: str | os.PathLike[str], *, redact: bool = False
) -> int:
    """Rewrites a stretch of the log into one Parquet file.

    One object per question keeps writing simple and reading slow: a year of them is
    tens of thousands of files to open. Folding a month into one is what keeps the log
    queryable, and it is safe only because every model-authored payload is a string --
    nested ones would give each month its own inferred struct, and a field of one month
    would then be unreachable from a read spanning both.

    This writes the folded copy and leaves the objects it read where they are. ``read``
    counts an event once however many sources carry it, so the two can be read together
    until whoever owns the prefix decides the raw days can go.

    Args:
        pattern: Glob of the JSON objects to fold.
        output: Parquet file to write.
        redact: Empty the free-text fields, keeping the outcomes, timings, usage, and
            fingerprints. The columns are emptied rather than dropped, so a redacted
            month and an unredacted one share one schema and a query spanning both
            binds.

    Returns:
        How many events were written.
    """
    selected = "*"
    if redact:
        emptied = ", ".join(f"NULL::VARCHAR AS {name}" for name in _FREE_TEXT)
        selected = f"* EXCLUDE ({', '.join(_FREE_TEXT)}), {emptied}"
    with contextlib.closing(duckdb.connect()) as connection:
        # S608: a quoted literal around this module's own reader clause, and a select
        # list built from this module's own constants.
        connection.execute(
            f"COPY (SELECT {selected} FROM ({_source(pattern)})) "  # noqa: S608
            f"TO {sql_string(str(output))} (FORMAT parquet)"
        )
        row = connection.execute(
            f"SELECT count(*) FROM read_parquet({sql_string(str(output))})"  # noqa: S608
        ).fetchone()
    return 0 if row is None else int(row[0])
