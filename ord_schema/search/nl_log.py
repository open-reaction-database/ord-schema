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
timings, and fingerprints are fields; the question, the translations, and the compiler's
errors are strings. The payload has to stay a string because a model-authored ``Query``
is recursive and differently shaped record to record: compacting a month of them into
Parquet makes DuckDB infer a struct from whichever shapes that month held, and a
cross-month read then takes the first file's struct, so a later month's rows still count
in every aggregate while the field that made them worth recording is unreachable. The
measurement is in the ord-logbook entry "Where the question log lives".

Results are not recorded. A translation and a corpus fingerprint reproduce them, and a
question that returned a hundred thousand reactions should not cost a hundred thousand
rows of log.
"""

import dataclasses
import json
import os
import pathlib
import sys
import uuid
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from typing import Any, Protocol

from ord_schema.logging import get_logger

logger = get_logger(__name__)

# What became of a question. The first two are the ordinary endings and the rest are
# the error taxonomy ``nl`` raises, one outcome each, so a failure is queryable rather
# than a string someone has to grep.
ANSWERED = "answered"
EMPTY = "empty"
DECLINED = "declined"
MALFORMED = "malformed"
RATE_LIMITED = "rate_limited"
UNAVAILABLE = "unavailable"
UNRESOLVED_COMPOUND = "unresolved_compound"
TIMEOUT = "timeout"
SEARCH_FAILED = "search_failed"

OUTCOMES = (
    ANSWERED,
    EMPTY,
    DECLINED,
    MALFORMED,
    RATE_LIMITED,
    UNAVAILABLE,
    UNRESOLVED_COMPOUND,
    TIMEOUT,
    SEARCH_FAILED,
)


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
        outcome: One of ``OUTCOMES``.
        attempts: Every ``build_query`` call, in order. Empty where the model declined
            without writing a query, which is how "read the question and said no" is
            told apart from "backed off once the compiler refused a guess".
        row_count: How many reactions came back, where the search ran.
        declined_reason: What the model said the grammar lacks.
        error: The failure, where one ended the ask.
        answer_text: The prose the reader saw, which is what a thumb is a verdict on.
        answer_usage: What writing that sentence cost. Kept apart from the attempts so
            translation can be priced on its own, and summed into the record's total so
            what the question cost is one field.
        translate_ms: Wall time turning the question into a query.
        search_ms: Wall time running it, or None where it never ran.
        answer_ms: Wall time writing the prose, or None where none was written. Each
            stage is its own field rather than a key in a timings mapping, because a
            mapping whose keys depend on how far the ask got is a struct whose shape
            depends on the same thing -- and compaction would infer one schema per
            month from it.
        session_id: Groups a visit, so a rephrasing reads as a correction of what came
            before it. Minted by the client; there are no accounts.
        record_id: What a thumb and a label reference.
        timestamp: When the question was asked, in UTC.
    """

    question: str
    model: str
    corpus_fingerprint: str
    outcome: str
    attempts: Sequence[Attempt] = ()
    row_count: int | None = None
    declined_reason: str | None = None
    error: str | None = None
    answer_text: str | None = None
    answer_usage: Usage = dataclasses.field(default_factory=Usage)
    translate_ms: float | None = None
    search_ms: float | None = None
    answer_ms: float | None = None
    session_id: str | None = None
    record_id: str = dataclasses.field(default_factory=lambda: str(uuid.uuid4()))
    timestamp: str = dataclasses.field(
        default_factory=lambda: datetime.now(UTC).isoformat()
    )


def event(ask: Ask) -> dict[str, Any]:
    """Returns the JSON-ready ask event.

    Args:
        ask: What happened.

    Returns:
        A dict whose model-authored payloads are strings and whose everything else is
        typed. Nothing derivable is stored: the translation is the last attempt that
        was not rejected, and whether a repair fired is the length of the list.
    """
    return {
        "event": "ask",
        "record_id": ask.record_id,
        "session_id": ask.session_id,
        "timestamp": ask.timestamp,
        "question": ask.question,
        "attempts": [
            {
                "translation": json.dumps(attempt.translation, sort_keys=True),
                "error": attempt.error,
                "usage": dataclasses.asdict(attempt.usage),
            }
            for attempt in ask.attempts
        ],
        "outcome": ask.outcome,
        "row_count": ask.row_count,
        "declined_reason": ask.declined_reason,
        "error": ask.error,
        "answer_text": ask.answer_text,
        "model": ask.model,
        "corpus_fingerprint": ask.corpus_fingerprint,
        "usage": dataclasses.asdict(
            sum((attempt.usage for attempt in ask.attempts), ask.answer_usage)
        ),
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


class NullSink:
    """Drops every record, which is what a caller that asked for no log wants."""

    def write(self, event: Mapping[str, Any]) -> None:
        """Does nothing."""


class JsonlSink:
    """Appends records to a local file, one JSON object per line.

    What a laptop uses: an eval run records to a file with no credentials and no
    network, and the file reads with the same query as the bucket does.
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

    The prefix belongs to a private bucket rather than to one holding publishable
    artifacts: a record is free text somebody typed, against a session identifier, and
    keeping the two apart makes opening the public bucket a policy change rather than a
    standing audit of every prefix.
    """

    def __init__(
        self, store: ObjectStore, *, bucket: str, prefix: str = "nl-log"
    ) -> None:
        """Initializes the sink.

        Args:
            store: Something with boto3's ``put_object``.
            bucket: Bucket to write into.
            prefix: Key prefix; the date partition goes under it.
        """
        self._store = store
        self._bucket = bucket
        self._prefix = prefix.rstrip("/")

    def write(self, event: Mapping[str, Any]) -> None:
        """Stores one event under ``<prefix>/dt=<date>/<record_id>.json``.

        The date comes from the event's own timestamp rather than from the clock, so a
        record lands in the partition it belongs to even when it is written late.
        """
        day = str(event["timestamp"])[:10]
        self._store.put_object(
            Bucket=self._bucket,
            Key=f"{self._prefix}/dt={day}/{event['record_id']}.json",
            Body=json.dumps(event, sort_keys=True).encode("utf-8"),
        )


def emit(sink: Sink | None, ask: Ask) -> None:
    """Records an ask, and never lets the log cost the answer.

    A sink reaches a file or a network, and a reader is waiting on the question rather
    than on the record of it. A failure here is warned about and dropped, which is the
    same contract the search cache holds: an optimization, never a dependency.

    Args:
        sink: Where the record goes; None records nothing.
        ask: What happened.
    """
    if sink is None:
        return
    try:
        sink.write(event(ask))
    except Exception as error:  # noqa: BLE001
        logger.warning("dropping a question log record: %s", error)
