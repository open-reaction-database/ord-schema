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

"""Questions in English, translated into the search grammar and answered.

A model cannot be *constrained* to emit a valid ``Query``. The grammar is recursive, and
both constrained-decoding paths refuse it -- structured outputs and strict tools share
one validator, which rejects circular references and then, once the recursion is
removed, rejects what is left as too large. So translation is ordinary generation,
checked afterwards: the tool schema is documentation the model mostly follows rather
than a rule it cannot break.

Four things follow, and they are the whole design. The predicate tree usually arrives
JSON-encoded inside a string, so it is coerced before validation rather than after a
failure. The compiler's errors name the offending path and suggest a real one, so a
query that does not compile is worth handing back exactly once. A model with no way to
decline invents a query for a question the grammar cannot express, so ``cannot_answer``
is offered beside ``build_query``. And the schema rendering and the grammar are ~15k
tokens of prompt that never change, so they are cached, which is most of what a query
costs.

The measurements behind those choices are in the ord-logbook entry "What constrains a
natural-language layer over the search grammar".
"""

import dataclasses
import functools
import hashlib
import json
import os
import time
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from importlib import resources
from typing import Any

import anthropic
import pyarrow as pa
from anthropic.types import (
    MessageParam,
    TextBlock,
    TextBlockParam,
    ToolParam,
    ToolUseBlock,
)

from ord_schema.artifacts import pivot as pivot_levels
from ord_schema.logging import get_logger
from ord_schema.search import execute, nl_log, placeholders, query, schema

logger = get_logger(__name__)

# Haiku translates most questions correctly and costs a fraction of what Opus does; the
# repair turn is what closes much of the gap. Every entry point takes a model, because
# which one is worth its price is a question for an eval set rather than for this line.
DEFAULT_MODEL = "claude-haiku-5-5"
# Read in preference to ANTHROPIC_API_KEY, which a shell holds for whatever else the
# person at it is doing. A key that answers corpus questions is a deployment's, billed
# to whoever runs the corpus, and naming it apart keeps a general-purpose key from
# quietly paying for a search -- and this one from reaching everything that reads the
# standard name. The SDK's own lookup still applies where this is unset.
API_KEY_VARIABLE = "ORD_ANTHROPIC_API_KEY"
# Room for thinking as well as the reply: models that think by default count it against
# the cap, and a cap sized for the tool call alone ends the turn before the call.
MAX_TOKENS = 16000
ANSWER_MAX_TOKENS = 4096
EFFORTS = ("low", "medium", "high", "xhigh", "max")

SYSTEM_PROMPT = (
    (resources.files("ord_schema.search") / "nl_prompt.md").read_text(encoding="utf-8")
    + "\n"
    + schema.describe()
)

# Built once at import. The tool definition is part of the cached prefix, so a dict
# rebuilt per call -- with any difference in key order -- would silently cost a cache
# miss on every query.
TOOL: ToolParam = {
    "name": "build_query",
    # The input is the query, said explicitly because a model otherwise nests it under
    # the tool's own name: its top-level keys are the Query fields, not a "query" key
    # holding them. _unwrapped recovers from that; saying so avoids it.
    "description": (
        "Build an ORD search query from the user's question. The input is the query "
        "itself, whose top-level keys are where, aggregate, order_by, and limit; do "
        "not nest it under a key of its own."
    ),
    "input_schema": query.Query.model_json_schema(),
}
# Forcing build_query would leave a model with no way to decline, and a model with no
# way to decline invents a query rather than refusing. This is that way. Neither tool is
# forced: newer models refuse a forced tool_choice, so the prompt asks for a call and a
# turn that makes none is handed back the way a query that does not compile is.
REFUSAL_TOOL: ToolParam = {
    "name": "cannot_answer",
    "description": (
        "Say that the question cannot be expressed in this grammar. Use it rather than "
        "building a query that means something else."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "reason": {
                "type": "string",
                "description": "What the question asks for that the grammar lacks.",
            }
        },
        "required": ["reason"],
    },
}
_SYSTEM: list[TextBlockParam] = [
    {
        "type": "text",
        "text": SYSTEM_PROMPT,
        "cache_control": {"type": "ephemeral", "ttl": "1h"},
    }
]
_ANSWER_SYSTEM = (
    "You describe the result of a database query in one or two plain sentences. State "
    "only what the summary shows. Do not invent chemistry, and do not use markdown."
)


@functools.cache
def prompt_fingerprint() -> str:
    """Returns a short digest of the cached prefix.

    Two things at once, which is why it is worth recording on every question. It names
    the translator a record was written by, so a query from months ago is not read as
    evidence about today's prompt. And because the prefix is exactly what the cache
    keys on, the fingerprint changes when and only when a cache read would have missed
    -- so a run where it moves unexpectedly is money already spent.

    Cached: the prompt and the tools are module constants built at import, so a served
    process would otherwise digest the whole prefix again on every question it records.
    Anything that replaces one of them has to call ``cache_clear``.

    Returns:
        Twelve hex characters over the system prompt and both tool definitions.
    """
    digest = hashlib.sha256(SYSTEM_PROMPT.encode())
    for tool in (TOOL, REFUSAL_TOOL):
        digest.update(json.dumps(tool, sort_keys=True).encode())
    return digest.hexdigest()[:12]


class NLQueryError(Exception):
    """A question could not be answered.

    Carries what the failed attempt spent, because a record is most worth having
    exactly when the translation did not work -- and an exception has no return value
    to put it in.

    Attributes:
        attempts: Every ``build_query`` call made before giving up, in order.
        turns: How many times the model was asked to translate before giving up, which
            exceeds the attempts where a turn wrote no query.
        usage: What the whole failed translation cost, including turns that produced no
            attempt.
        elapsed_ms: Wall time spent before the failure.
    """

    def __init__(
        self,
        *args: object,
        attempts: Sequence[nl_log.Attempt] = (),
        turns: int = 0,
        usage: nl_log.Usage | None = None,
        elapsed_ms: float = 0.0,
    ) -> None:
        """Initializes the error.

        Args:
            *args: Passed to ``Exception``.
            attempts: The queries built before giving up.
            turns: How many times the model was asked.
            usage: What those turns cost.
            elapsed_ms: Wall time spent.
        """
        super().__init__(*args)
        self.attempts = tuple(attempts)
        self.turns = turns
        self.usage = usage if usage is not None else nl_log.Usage()
        self.elapsed_ms = elapsed_ms


class ModelUnavailableError(NLQueryError):
    """The model could not be reached."""


class ModelRateLimitedError(NLQueryError):
    """The model refused the request because the caller is over its rate limit."""


class ModelRefusedError(NLQueryError):
    """The model's safety classifiers declined the request.

    Distinct from ``UnanswerableError``, where the model read the question and found the
    grammar lacking: here a classifier declined the request before any answer, and the
    category says which one.
    """


class MalformedQueryError(NLQueryError):
    """The model's query did not compile, and neither did its repair."""


class UnanswerableError(NLQueryError):
    """The question cannot be put to this grammar, and the model said so.

    Distinct from a malformed query: nothing was wrong with the model's reasoning, and
    retrying will not help. Comparing two columns is the standard case -- a value is a
    literal, never another column -- and a layer without a way to say so
    answers with a plausible query that means something else.
    """

    @property
    def attempted(self) -> bool:
        """Returns whether a query came first and failed to compile.

        A caller is served either way, so this changes nothing about the answer; a
        measurement wants it, because recognizing an inexpressible question is not the
        same as backing into it after the compiler refused a guess.
        """
        return bool(self.attempts)


@dataclasses.dataclass(frozen=True)
class Translation:
    """A question turned into a query, and what that took.

    Attributes:
        query: The query, validated and compiled.
        attempts: Every ``build_query`` call, in order, each with the compiler's verdict
            and what the turn cost. The rejected query and the error rejecting it are
            what prompt work runs on, and what pricing a repair turn needs.
        usage: What the whole translation cost.
        elapsed_ms: Wall time spent translating.
        turns: How many times the model was asked: two where the first turn was
            handed back, whether or not it wrote a query to record as an attempt.
    """

    query: query.Query
    attempts: tuple[nl_log.Attempt, ...]
    usage: nl_log.Usage
    elapsed_ms: float
    turns: int = 1


@dataclasses.dataclass(frozen=True)
class Description:
    """A sentence about a result, and what writing it cost.

    Attributes:
        text: The prose.
        usage: What the call cost. The prose is a second model call, so a question's
            price is not the translation's.
        elapsed_ms: Wall time spent.
    """

    text: str
    usage: nl_log.Usage
    elapsed_ms: float


@dataclasses.dataclass(frozen=True)
class Answer:
    """What a question produced, including the query it became.

    Attributes:
        question: The question as asked.
        query: The query that ran, for showing the user what was searched and for
            running it again unchanged.
        table: What the search returned.
        text: A sentence or two describing that result.
        record_id: What this question was logged as. A caller hands it to the reader so
            a thumb has something to reference; without it the log records feedback
            nobody can attach to a question.
        resolved: Each compound name the query gives, mapped to the canonical SMILES
            the search matched it as. The query keeps the name, so this is where a
            reader sees the structure; empty for a query naming no compound.
    """

    question: str
    query: query.Query
    table: pa.Table
    text: str
    record_id: str = ""
    resolved: dict[str, str] = dataclasses.field(default_factory=dict)


def get_client() -> anthropic.Anthropic:
    """Returns a client reading its credentials from the environment.

    Returns:
        A client authenticated with ``ORD_ANTHROPIC_API_KEY`` where the environment
        sets it, and with whatever the SDK finds otherwise.
    """
    key = os.environ.get(API_KEY_VARIABLE)
    return anthropic.Anthropic(api_key=key) if key else anthropic.Anthropic()


def _coerce(value: Any) -> Any:
    """Returns the value with JSON-encoded strings parsed back into objects.

    A model handed a recursive tool schema usually returns the nested predicate as a
    string holding JSON rather than as an object, which the models reject. Parsing it
    back is the normal path rather than a fallback: it is what both measured models do
    most of the time.

    A model encoding the tree that way often closes it with one brace too many, the
    outer query's, so closers after a complete value are dropped rather than failing
    the parse: they add nothing to what the value says.

    Args:
        value: Whatever the tool call carried.

    Returns:
        The same structure, with any JSON-encoded string replaced by what it encodes.
    """
    if isinstance(value, str):
        try:
            decoded, end = json.JSONDecoder().raw_decode(value.strip())
        except json.JSONDecodeError:
            return value
        if value.strip()[end:].strip(" \t\n\r}]"):
            return value
        return _coerce(decoded)
    if isinstance(value, dict):
        return {key: _coerce(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_coerce(item) for item in value]
    return value


def _unwrapped(value: Any) -> Any:
    """Returns a tool input without the envelope a model may nest it in.

    A model that puts its arguments under the tool's parameter name sends
    ``{"query": {...}}`` where the schema asks for the query itself. ``Query`` declares
    no ``query`` field and refuses undeclared ones, so a payload of exactly that shape
    can mean nothing else. Unwrapping it costs no round trip, where handing it to the
    repair turn costs one and does not reliably work: a model that nests once tends to
    nest again when told the result did not validate.

    Args:
        value: A tool input, already coerced.

    Returns:
        What the envelope holds, or ``value`` where it is not one.
    """
    if (
        isinstance(value, dict)
        and set(value) == {"query"}
        and isinstance(value["query"], dict)
    ):
        return value["query"]
    return value


def _validated(coerced: Any) -> query.Query:
    """Returns the query a tool call carries, proven to compile.

    Args:
        coerced: The tool call's input, already run through ``_coerce``.

    Returns:
        The parsed query.

    Raises:
        ValueError: If it does not validate against the grammar, or if it asks nothing.
        QueryError: If it validates but does not compile against the schema.
    """
    parsed = query.Query.model_validate(coerced)
    if parsed.where is None and parsed.aggregate is None and parsed.limit is None:
        # A query with no predicate, no aggregate, and no limit compiles and returns
        # the whole corpus, which is never the answer to a question. It reads as a
        # result rather than as a failure, so it is refused here where the repair turn
        # can still ask again. A limit makes it "show me some reactions", which is a
        # question someone does ask.
        raise ValueError(
            "the query asks nothing: it has no where, no aggregate, and no limit, so "
            "it would return every reaction in the corpus"
        )
    # Every level offered as pivoted, because this asks whether the query is
    # well-formed, not how it would be planned. Which levels a corpus actually holds is
    # a property of that corpus, and refusing here on this process's answer would reject
    # a query the corpus about to run it can serve.
    query.compile_query(
        parsed,
        pivot=lambda path: (
            pivot_levels.table_name(path) if path in pivot_levels.LEVELS else None
        ),
    )
    return parsed


def _spent(response: Any) -> nl_log.Usage:
    """Returns what one response cost.

    Args:
        response: A Messages API response.

    Returns:
        Its usage. The cache fields are read defensively because the API omits them on
        some responses, and a missing one is zero rather than a failed record.
    """
    usage = response.usage
    return nl_log.Usage(
        input=usage.input_tokens,
        output=usage.output_tokens,
        cache_read=getattr(usage, "cache_read_input_tokens", 0) or 0,
        cache_creation=getattr(usage, "cache_creation_input_tokens", 0) or 0,
    )


def _create(client: anthropic.Anthropic, *, effort: str | None, **request: Any) -> Any:
    """Returns one Messages API response, its failures named as this module's.

    Args:
        client: Anthropic client.
        effort: How much the model thinks, or None to leave it to the model's default.
            Omitted rather than sent as a default, because models without effort levels
            reject the parameter outright.
        **request: Everything else the request carries.

    Returns:
        The response.

    Raises:
        ModelRefusedError: If the model's safety classifiers declined the request.
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
    """
    if effort is not None:
        request["output_config"] = {"effort": effort}
    try:
        response = client.messages.create(**request)
    except anthropic.RateLimitError as error:
        raise ModelRateLimitedError(str(error)) from error
    except (anthropic.APIConnectionError, anthropic.APIStatusError) as error:
        raise ModelUnavailableError(str(error)) from error
    if response.stop_reason == "refusal":
        details = getattr(response, "stop_details", None)
        category = getattr(details, "category", None) or "no category given"
        raise ModelRefusedError(f"refused: {category}", usage=_spent(response))
    return response


@dataclasses.dataclass(frozen=True)
class _Turn:
    """One translation turn.

    Attributes:
        content: The response's content blocks, replayed verbatim on a repair: a model
            that thinks needs its thinking blocks back beside the tool call they led to.
        call: The ``build_query`` call, or None where the turn made no tool call.
        usage: What the turn cost.
    """

    content: Any
    call: ToolUseBlock | None
    usage: nl_log.Usage


def _call(
    client: anthropic.Anthropic,
    model: str,
    messages: list[MessageParam],
    effort: str | None,
) -> _Turn:
    """Returns one translation turn.

    Args:
        client: Anthropic client.
        model: Which model to ask.
        messages: The conversation so far.
        effort: How much the model thinks, or None for the model's default.

    Returns:
        The turn, whose ``call`` is None where the model answered without a tool.

    Raises:
        UnanswerableError: If the model calls ``cannot_answer`` instead.
        ModelRefusedError: If the model's safety classifiers declined the request.
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
    """
    response = _create(
        client,
        effort=effort,
        model=model,
        max_tokens=MAX_TOKENS,
        system=_SYSTEM,
        messages=messages,
        tools=[TOOL, REFUSAL_TOOL],
        # At most one call, so a repair has exactly one tool_use to answer.
        tool_choice={"type": "auto", "disable_parallel_tool_use": True},
    )
    spent = _spent(response)
    for block in response.content:
        if isinstance(block, ToolUseBlock):
            if block.name == REFUSAL_TOOL["name"]:
                reason = "no reason given"
                if isinstance(block.input, dict):
                    reason = str(block.input.get("reason", reason))
                raise UnanswerableError(reason, usage=spent)
            return _Turn(response.content, block, spent)
    return _Turn(response.content, None, spent)


def _repair(turn: _Turn, rejection: str) -> MessageParam:
    """Returns the message handing a failed turn back to the model.

    Args:
        turn: The turn that produced no usable query.
        rejection: Why: the compiler's error, or that no tool was called.

    Returns:
        A tool result carrying the error where the turn called ``build_query``, and a
        plain reminder of the two tools where it answered in prose instead.
    """
    if turn.call is None:
        return {
            "role": "user",
            "content": (
                "Answer by calling build_query with the query, or cannot_answer if "
                "the grammar cannot express the question."
            ),
        }
    return {
        "role": "user",
        "content": [
            {
                "type": "tool_result",
                "tool_use_id": turn.call.id,
                "is_error": True,
                "content": (
                    f"That query was rejected: {rejection}. Call build_query again "
                    "with it fixed."
                ),
            }
        ],
    }


def translate(
    question: str,
    *,
    structures: Mapping[str, placeholders.Placeholder] | None = None,
    client: anthropic.Anthropic | None = None,
    model: str = DEFAULT_MODEL,
    effort: str | None = None,
    repair: bool = True,
) -> "Translation":
    """Returns the query a question asks for.

    Args:
        question: The question, in English.
        structures: Each placeholder the question writes, mapped to the structure
            it stands for. The model is told each one's format and never its
            text, which is substituted into the query before it is validated.
        client: Anthropic client; one is built from the environment if omitted.
        model: Which model translates.
        effort: How much the model thinks, one of ``EFFORTS``, or None for the model's
            default. Models without effort levels reject anything but None.
        repair: Hand a failure back once, carrying the compiler's error or, where the
            model called no tool, a reminder to. Turn it off to measure how often a
            model is right without help.

    Returns:
        The query, every ``build_query`` call that led to it, and what they cost.

    Raises:
        MalformedQueryError: If the query does not compile, or the model calls no
            tool, after the repair turn where one was allowed.
        UnanswerableError: If the model says the grammar cannot express the question.
        ModelRefusedError: If the model's safety classifiers declined the request.
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached. Every one of these
            carries the attempts and the usage spent reaching it, because a failed
            translation is a record worth keeping and an exception has no return value.
        PlaceholderError: If the question and its placeholders do not fit each
            other, before the model is asked anything.
    """
    structures = structures or {}
    placeholders.check(question, structures)
    client = client if client is not None else get_client()
    started = time.monotonic()
    attempts: list[nl_log.Attempt] = []
    spent = nl_log.Usage()
    turns = 0

    def elapsed() -> float:
        return (time.monotonic() - started) * 1000

    def failed(error: NLQueryError) -> NLQueryError:
        """Returns the error with what this translation spent attached to it.

        The error's own usage is added rather than replaced: a turn that declined
        carries what it cost and produced no attempt to hold it, so dropping it would
        leave a refusal looking free.
        """
        error.attempts = tuple(attempts)
        error.turns = turns
        error.usage = spent + error.usage
        error.elapsed_ms = elapsed()
        return error

    def call(messages: list[MessageParam]) -> _Turn:
        """Returns one turn, counted, its cost added to what this translation spent."""
        nonlocal spent, turns
        turns += 1
        try:
            turn = _call(client, model, messages, effort)
        except NLQueryError as error:
            raise failed(error) from error
        spent += turn.usage
        return turn

    def judged(turn: _Turn) -> tuple[query.Query | None, str]:
        """Returns the turn's query, or why it has none, recording any attempt."""
        if turn.call is None:
            # Costed in ``spent``, but no attempt: there is no query to record.
            return None, "the model returned no query"
        coerced = _unwrapped(_coerce(turn.call.input))
        try:
            parsed = _validated(placeholders.substitute(coerced, structures))
        except (ValueError, query.QueryError) as error:
            attempts.append(nl_log.Attempt(coerced, str(error), turn.usage))
            return None, str(error)
        attempts.append(nl_log.Attempt(coerced, None, turn.usage))
        return parsed, ""

    described = placeholders.describe(structures)
    content = f"{question}\n\n{described}" if described else question
    messages: list[MessageParam] = [{"role": "user", "content": content}]
    first = call(messages)
    parsed, rejection = judged(first)
    if parsed is not None:
        return Translation(parsed, tuple(attempts), spent, elapsed(), turns)
    if not repair:
        raise failed(MalformedQueryError(rejection))
    logger.info("repairing a translation: %s", rejection)
    messages += [
        {"role": "assistant", "content": first.content},
        _repair(first, rejection),
    ]
    # Declining on the second turn is still the right answer, and the caller gets the
    # same one. An attempt already recorded is what says the model reached it by way
    # of a query the compiler refused rather than by reading the question.
    parsed, rejection = judged(call(messages))
    if parsed is None:
        # One repair, never a loop: a model that missed twice with the compiler's own
        # suggestion in hand is not going to find it on the third paid attempt.
        raise failed(MalformedQueryError(rejection))
    return Translation(parsed, tuple(attempts), spent, elapsed(), turns)


def summarize(table: pa.Table, *, rows: int = 5) -> str:
    """Returns a description of a result small enough to put in a prompt.

    Args:
        table: What the search returned.
        rows: How many sample rows to show.

    Returns:
        The row count, the column names, and up to ``rows`` rows. A result of a hundred
        thousand reactions costs the same prompt as one of three, which is the point:
        the model describes the shape of an answer rather than reading it.
    """
    lines = [f"{table.num_rows} rows, columns: {', '.join(table.column_names)}"]
    lines += [json.dumps(row, default=str) for row in table.slice(0, rows).to_pylist()]
    if table.num_rows > rows:
        lines.append(f"... {table.num_rows - rows} more rows not shown")
    return "\n".join(lines)


def answer(
    question: str,
    table: pa.Table,
    *,
    client: anthropic.Anthropic | None = None,
    model: str = DEFAULT_MODEL,
    effort: str | None = None,
) -> Description:
    """Returns a sentence or two saying what a result shows.

    Args:
        question: The question the result answers.
        table: What the search returned.
        client: Anthropic client; one is built from the environment if omitted.
        model: Which model writes the prose.
        effort: How much the model thinks, or None for the model's default.

    Returns:
        The prose and what it cost. The model sees ``summarize(table)`` rather than the
        rows, so it can describe only what the summary states.

    Raises:
        ModelRefusedError: If the model's safety classifiers declined the request.
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
    """
    client = client if client is not None else get_client()
    started = time.monotonic()
    response = _create(
        client,
        effort=effort,
        model=model,
        max_tokens=ANSWER_MAX_TOKENS,
        system=_ANSWER_SYSTEM,
        messages=[
            {
                "role": "user",
                "content": f"Question: {question}\n\nResult:\n{summarize(table)}",
            }
        ],
    )
    return Description(
        text="".join(
            block.text for block in response.content if isinstance(block, TextBlock)
        ),
        usage=_spent(response),
        elapsed_ms=(time.monotonic() - started) * 1000,
    )


_OUTCOMES: tuple[tuple[type[BaseException], nl_log.Outcome], ...] = (
    (UnanswerableError, nl_log.Outcome.DECLINED),
    (ModelRefusedError, nl_log.Outcome.REFUSED),
    (MalformedQueryError, nl_log.Outcome.MALFORMED),
    (ModelRateLimitedError, nl_log.Outcome.RATE_LIMITED),
    (ModelUnavailableError, nl_log.Outcome.UNAVAILABLE),
    (TimeoutError, nl_log.Outcome.TIMEOUT),
    (query.QueryError, nl_log.Outcome.SEARCH_FAILED),
    # Before the ValueError below, which it inherits from. An unresolvable compound is
    # the question's fault and a corpus whose artifacts disagree is the deployment's,
    # and the whole point of the taxonomy is that a reader can tell them apart.
    (execute.PairingError, nl_log.Outcome.SEARCH_FAILED),
    # Last, because a resolver failure is a ValueError and so are several others; by
    # here the more specific kinds have already matched.
    (ValueError, nl_log.Outcome.UNRESOLVED_COMPOUND),
)


def outcome_of(error: BaseException) -> nl_log.Outcome:
    """Returns the outcome a failure is recorded as.

    Args:
        error: What ended the ask.

    Returns:
        The matching ``nl_log.Outcome``, or ``SEARCH_FAILED`` for anything the search
        raised that has no more specific name.
    """
    for kind, outcome in _OUTCOMES:
        if isinstance(error, kind):
            return outcome
    return nl_log.Outcome.SEARCH_FAILED


def _stamps(corpus: execute.Corpus) -> tuple[str, str]:
    """Returns the fingerprints naming what answered a question.

    Args:
        corpus: What the query will run against.

    Returns:
        The corpus fingerprint and the prompt fingerprint, both empty where they could
        not be read. A record is an optimization and never a dependency, so an artifact
        that has gone unreadable costs the fingerprint rather than the answer somebody
        is waiting on.
    """
    try:
        return corpus.fingerprint, prompt_fingerprint()
    except Exception as error:  # noqa: BLE001
        logger.warning("cannot fingerprint what answered this question: %s", error)
        return "", ""


def ask(
    question: str,
    corpus: execute.Corpus,
    *,
    structures: Mapping[str, placeholders.Placeholder] | None = None,
    client: anthropic.Anthropic | None = None,
    model: str = DEFAULT_MODEL,
    effort: str | None = None,
    timeout_seconds: float = 60.0,
    sink: nl_log.Sink | None = None,
    source: str | None = None,
    session_id: str | None = None,
) -> Answer:
    """Answers a question against a corpus, and records what happened.

    The model's own failures arrive as ``NLQueryError`` subclasses; the search's arrive
    as themselves. Folding both into one taxonomy would cost a caller the distinction
    it answers with -- an unresolvable compound is the question's fault, a corpus whose
    artifacts disagree is the deployment's.

    Every ending is recorded, not only the one where everything worked. A question that
    matched nothing, one the compiler refused twice, and one the model declined are the
    records worth reading, and a log written on the success path holds none of them.

    Args:
        question: The question, in English.
        corpus: The corpus to search.
        structures: Each placeholder the question writes, mapped to the structure
            it stands for.
        client: Anthropic client; one is built from the environment if omitted.
        model: Which model translates and describes.
        effort: How much the model thinks, or None for the model's default.
        timeout_seconds: Passed to the search, so a slow query fails rather than hangs.
        sink: Where the record goes; nothing is recorded if omitted.
        source: Which surface this arrived from -- ``"web"``, ``"eval"``, ``"api"``.
            One prefix holds all of them, and nothing else in the record separates a
            person rephrasing from a harness working through a case file.
        session_id: Groups this question with the rest of a visit, which is what makes
            a rephrasing readable as a correction of the question before it.

    Returns:
        The query that ran, the reactions it found, a description of them, the
        identifier the record was written under, and what each compound the query
        named was searched as.

    Raises:
        MalformedQueryError: If translation produces nothing that compiles.
        UnanswerableError: If the grammar cannot express the question.
        ModelRefusedError: If the model's safety classifiers declined the request.
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
        ValueError: If a compound the query names cannot be resolved.
        PlaceholderError: If the question and its placeholders do not fit each
            other, before the model is asked anything; nothing is recorded.
        TimeoutError: If the search exceeds ``timeout_seconds``.
        QueryError: If the compiled query cannot be run against these artifacts.
        PairingError: If the corpus and its structures do not line up.
    """
    client = client if client is not None else get_client()
    timestamp = datetime.now(UTC).isoformat()
    # Both only where a record is going somewhere: an id nothing wrote would let a
    # caller hang a thumb on an ask that does not exist, and reading the stamps opens a
    # footer per artifact.
    record_id = nl_log.new_identifier() if sink is not None else ""
    stamps = _stamps(corpus) if sink is not None else None

    def record(outcome: nl_log.Outcome, **fields: Any) -> None:
        """Writes one record, filling in what every ending of this ask shares."""
        if stamps is None:
            return
        corpus_fingerprint, translator = stamps
        nl_log.emit(
            sink,
            nl_log.Ask(
                question=question,
                model=model,
                effort=effort,
                corpus_fingerprint=corpus_fingerprint,
                prompt_fingerprint=translator,
                outcome=outcome,
                source=source,
                session_id=session_id,
                record_id=record_id,
                timestamp=timestamp,
                structures=placeholders.dump(structures or {}),
                **fields,
            ),
        )

    try:
        translated = translate(
            question,
            structures=structures,
            client=client,
            model=model,
            effort=effort,
        )
    except NLQueryError as error:
        record(
            outcome_of(error),
            attempts=error.attempts,
            turns=error.turns,
            declined_reason=(
                str(error) if isinstance(error, UnanswerableError) else None
            ),
            error=str(error),
            usage=error.usage,
            translate_ms=error.elapsed_ms,
        )
        raise
    started = time.monotonic()
    try:
        table = corpus.search(translated.query, timeout_seconds=timeout_seconds)
    except Exception as error:
        record(
            outcome_of(error),
            attempts=translated.attempts,
            turns=translated.turns,
            error=str(error),
            usage=translated.usage,
            translate_ms=translated.elapsed_ms,
            search_ms=(time.monotonic() - started) * 1000,
        )
        raise
    search_ms = (time.monotonic() - started) * 1000
    resolved = execute.resolutions(table)
    try:
        described = answer(question, table, client=client, model=model, effort=effort)
    except NLQueryError as error:
        record(
            outcome_of(error),
            attempts=translated.attempts,
            turns=translated.turns,
            row_count=table.num_rows,
            resolved=resolved,
            error=str(error),
            usage=translated.usage + error.usage,
            translate_ms=translated.elapsed_ms,
            search_ms=search_ms,
        )
        raise
    record(
        nl_log.Outcome.ANSWERED if table.num_rows else nl_log.Outcome.EMPTY,
        attempts=translated.attempts,
        turns=translated.turns,
        row_count=table.num_rows,
        resolved=resolved,
        answer_text=described.text,
        usage=translated.usage + described.usage,
        translate_ms=translated.elapsed_ms,
        search_ms=search_ms,
        answer_ms=described.elapsed_ms,
    )
    return Answer(
        question=question,
        query=translated.query,
        table=table,
        text=described.text,
        record_id=record_id,
        resolved=resolved,
    )
