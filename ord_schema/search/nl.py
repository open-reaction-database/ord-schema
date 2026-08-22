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

Three things follow, and they are the whole design. The predicate tree usually arrives
JSON-encoded inside a string, so it is coerced before validation rather than after a
failure. The compiler's errors name the offending path and suggest a real one, so a
query that does not compile is worth handing back exactly once. And the schema rendering
and the grammar are ~15k tokens of prompt that never change, so they are cached, which
is most of what a query costs.

The measurements behind those choices are in the ord-logbook entry "What constrains a
natural-language layer over the search grammar".
"""

import dataclasses
import json
import os
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

from ord_schema.logging import get_logger
from ord_schema.search import execute, query, schema

logger = get_logger(__name__)

# Haiku translates most questions correctly and costs a fraction of what Opus does; the
# repair turn is what closes much of the gap. Every entry point takes a model, because
# which one is worth its price is a question for an eval set rather than for this line.
DEFAULT_MODEL = "claude-haiku-4-5"
# Read in preference to ANTHROPIC_API_KEY, which a shell holds for whatever else the
# person at it is doing. A key that answers corpus questions is a deployment's, billed
# to whoever runs the corpus, and naming it apart keeps a general-purpose key from
# quietly paying for a search -- and this one from reaching everything that reads the
# standard name. The SDK's own lookup still applies where this is unset.
API_KEY_VARIABLE = "ORD_ANTHROPIC_API_KEY"
MAX_TOKENS = 2048
ANSWER_MAX_TOKENS = 512

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
    "description": "Build an ORD search query from the user's question.",
    "input_schema": query.Query.model_json_schema(),
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


class NLQueryError(Exception):
    """A question could not be answered."""


class ModelUnavailableError(NLQueryError):
    """The model could not be reached."""


class ModelRateLimitedError(NLQueryError):
    """The model refused the request because the caller is over its rate limit."""


class MalformedQueryError(NLQueryError):
    """The model's query did not compile, and neither did its repair."""


@dataclasses.dataclass(frozen=True)
class Answer:
    """What a question produced, including the query it became.

    Attributes:
        question: The question as asked.
        query: The query that ran, for showing the user what was searched and for
            running it again unchanged.
        table: What the search returned.
        text: A sentence or two describing that result.
    """

    question: str
    query: query.Query
    table: pa.Table
    text: str


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

    Args:
        value: Whatever the tool call carried.

    Returns:
        The same structure, with any JSON-encoded string replaced by what it encodes.
    """
    if isinstance(value, str):
        try:
            return _coerce(json.loads(value))
        except (json.JSONDecodeError, TypeError):
            return value
    if isinstance(value, dict):
        return {key: _coerce(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_coerce(item) for item in value]
    return value


def _validated(raw: Any) -> query.Query:
    """Returns the query a tool call carries, proven to compile.

    Args:
        raw: The tool call's input.

    Returns:
        The parsed query.

    Raises:
        ValueError: If it does not validate against the grammar.
        QueryError: If it validates but does not compile against the schema.
    """
    parsed = query.Query.model_validate(_coerce(raw))
    query.compile_query(parsed)
    return parsed


def _call(
    client: anthropic.Anthropic, model: str, messages: list[MessageParam]
) -> ToolUseBlock:
    """Returns the tool_use block from one forced call.

    Args:
        client: Anthropic client.
        model: Which model to ask.
        messages: The conversation so far.

    Returns:
        The ``tool_use`` content block.

    Raises:
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
        MalformedQueryError: If the response carries no tool call at all.
    """
    try:
        response = client.messages.create(
            model=model,
            max_tokens=MAX_TOKENS,
            system=_SYSTEM,
            messages=messages,
            tools=[TOOL],
            tool_choice={"type": "tool", "name": "build_query"},
        )
    except anthropic.RateLimitError as error:
        raise ModelRateLimitedError(str(error)) from error
    except (anthropic.APIConnectionError, anthropic.APIStatusError) as error:
        raise ModelUnavailableError(str(error)) from error
    for block in response.content:
        if isinstance(block, ToolUseBlock):
            return block
    raise MalformedQueryError("the model returned no query")


def translate(
    question: str,
    *,
    client: anthropic.Anthropic | None = None,
    model: str = DEFAULT_MODEL,
    repair: bool = True,
) -> query.Query:
    """Returns the query a question asks for.

    Args:
        question: The question, in English.
        client: Anthropic client; one is built from the environment if omitted.
        model: Which model translates.
        repair: Hand a failure back once, carrying the compiler's error. Turn it off to
            measure how often a model is right without help.

    Returns:
        A query that validates and compiles against the projection schema.

    Raises:
        MalformedQueryError: If the query does not compile, after the repair turn where
            one was allowed.
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
    """
    client = client if client is not None else get_client()
    messages: list[MessageParam] = [{"role": "user", "content": question}]
    block = _call(client, model, messages)
    try:
        return _validated(block.input)
    except (ValueError, query.QueryError) as error:
        if not repair:
            raise MalformedQueryError(str(error)) from error
        first = error
    logger.info("repairing a query that did not compile: %s", first)
    messages += [
        {
            "role": "assistant",
            "content": [
                {
                    "type": "tool_use",
                    "id": block.id,
                    "name": block.name,
                    "input": block.input,
                }
            ],
        },
        {
            "role": "user",
            "content": [
                {
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "is_error": True,
                    "content": (
                        f"That query was rejected: {first}. Call build_query again "
                        "with it fixed."
                    ),
                }
            ],
        },
    ]
    retry = _call(client, model, messages)
    try:
        return _validated(retry.input)
    except (ValueError, query.QueryError) as error:
        # One repair, never a loop: a model that missed twice with the compiler's own
        # suggestion in hand is not going to find it on the third paid attempt.
        raise MalformedQueryError(str(error)) from error


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
) -> str:
    """Returns a sentence or two saying what a result shows.

    Args:
        question: The question the result answers.
        table: What the search returned.
        client: Anthropic client; one is built from the environment if omitted.
        model: Which model writes the prose.

    Returns:
        Plain text. The model sees ``summarize(table)`` rather than the rows, so it can
        describe only what the summary states.

    Raises:
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
    """
    client = client if client is not None else get_client()
    try:
        response = client.messages.create(
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
    except anthropic.RateLimitError as error:
        raise ModelRateLimitedError(str(error)) from error
    except (anthropic.APIConnectionError, anthropic.APIStatusError) as error:
        raise ModelUnavailableError(str(error)) from error
    return "".join(
        block.text for block in response.content if isinstance(block, TextBlock)
    )


def ask(
    question: str,
    corpus: execute.Corpus,
    *,
    client: anthropic.Anthropic | None = None,
    model: str = DEFAULT_MODEL,
    timeout_seconds: float = 60.0,
) -> Answer:
    """Answers a question against a corpus.

    The model's own failures arrive as ``NLQueryError`` subclasses; the search's arrive
    as themselves. Folding both into one taxonomy would cost a caller the distinction
    it answers with -- an unresolvable compound is the question's fault, a corpus whose
    artifacts disagree is the deployment's.

    Args:
        question: The question, in English.
        corpus: The corpus to search.
        client: Anthropic client; one is built from the environment if omitted.
        model: Which model translates and describes.
        timeout_seconds: Passed to the search, so a slow query fails rather than hangs.

    Returns:
        The query that ran, the reactions it found, and a description of them.

    Raises:
        MalformedQueryError: If translation produces nothing that compiles.
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached.
        ValueError: If a compound the query names cannot be resolved.
        TimeoutError: If the search exceeds ``timeout_seconds``.
        QueryError: If the compiled query cannot be run against these artifacts.
        PairingError: If the corpus and its structures do not line up.
    """
    client = client if client is not None else get_client()
    translated = translate(question, client=client, model=model)
    table = corpus.search(translated, timeout_seconds=timeout_seconds)
    return Answer(
        question=question,
        query=translated,
        table=table,
        text=answer(question, table, client=client, model=model),
    )
