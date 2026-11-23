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

"""Scoring translations by the reactions they return.

A case states the question and one query that answers it. Scoring compares the
reactions a translation returns against the reactions that reference returns, and they
have to be the same set. That is not the same as pinning the query: several spellings of
a question are right and pinning one would fail a better translation than the one
written the day the case was added, but every right spelling answers with the same
reactions, so the set is what a case can fairly hold a model to.

Both queries run against whatever corpus the command is pointed at, so a case carries
no reaction IDs and nothing in it goes stale when the corpus is rebuilt. The reference
names structures outright rather than by compound, so scoring reaches no resolver, and
it asks for a compound the way the prompt tells a model to -- with ``same_compound``
rather than an equality on a spelling, which would fail a translation for following its
instructions.

Running the cases costs money and reaches the network, so the suite covers the
scoring and the outcomes ``run_case`` maps a translation onto, with the model stubbed.
The measurement itself is a command.
"""

import argparse
import dataclasses
import json
import pathlib
import uuid
from collections.abc import Sequence
from typing import Any

import anthropic
import pandas as pd
import pyarrow as pa
from pandas.testing import assert_frame_equal
from pydantic import BaseModel, model_validator

from ord_schema.logging import get_logger
from ord_schema.search import execute, nl, nl_log, query

logger = get_logger(__name__)

CASES = pathlib.Path(__file__).parent / "nl_cases.json"
# Longer than a served query would wait, because a run is unattended and a case whose
# translation is slow is worth a verdict rather than a stopped harness.
DEFAULT_TIMEOUT_SECONDS = 300.0
# What a record written by a run is marked with. A case file marched through end to end
# is not evidence about what people ask, and the two are indistinguishable in the log
# without this.
SOURCE = "eval"


class EvalCase(BaseModel):
    """One question, and what any correct answer to it must satisfy.

    Attributes:
        question: The question, in English.
        why: What this case is here to catch, for whoever reads a failure.
        reference: One whole query that answers the question, whose result a
            translation's has to match. Not what a translation has to produce --
            several spellings are right -- but the definition of what the right ones
            return: the reactions, or, where the question asks for a number, the
            measured rows.
        compiles: Whether the question should translate at all. False marks one the
            grammar cannot express, which the layer must refuse rather than fudge, and
            which therefore has no reference.
    """

    question: str
    why: str
    reference: dict[str, Any] | None = None
    compiles: bool = True

    @model_validator(mode="after")
    def _expressible_cases_have_a_reference(self) -> "EvalCase":
        if self.compiles == (self.reference is None):
            raise ValueError(
                f"{self.question}: a case the grammar can express needs a reference "
                "query, and one it cannot express has no answer to reference"
            )
        return self


@dataclasses.dataclass(frozen=True)
class CaseResult:
    """How one case came out.

    Attributes:
        case: The case that ran.
        passed: Whether it was satisfied.
        detail: What happened, in a line, whether it passed or not.
    """

    case: EvalCase
    passed: bool
    detail: str


def load_cases(path: pathlib.Path = CASES) -> list[EvalCase]:
    """Returns the cases a file holds.

    Args:
        path: A JSON file holding a list of cases.

    Returns:
        The parsed cases.
    """
    with path.open(encoding="utf-8") as handle:
        return [EvalCase.model_validate(entry) for entry in json.load(handle)]


def _sample(reaction_ids: set[str], limit: int = 3) -> str:
    """Returns a few IDs and a count, for a failure someone has to read."""
    shown = ", ".join(sorted(reaction_ids)[:limit])
    rest = len(reaction_ids) - limit
    return f"{shown} and {rest} more" if rest > 0 else shown


def _frame(table: pa.Table, columns: pd.Index, *, ordered: bool) -> pd.DataFrame:
    """Returns a result as a frame comparable against another with ``columns``.

    Args:
        table: What a query returned.
        columns: The labels to compare under, taken from the reference so that both
            sides carry them. A measure is named by whoever wrote the query, and two
            queries measuring the same thing under different names are the same answer;
            position is what identifies a column, which the compiler fixes as the
            group_by paths followed by the measures.
        ordered: Whether row order is part of the answer. Where it is not, the rows are
            sorted so that two queries returning the same rows compare equal.

    Returns:
        The frame, relabeled and ordered.
    """
    frame = table.to_pandas().set_axis(columns, axis=1)
    if ordered:
        return frame.reset_index(drop=True)
    return frame.sort_values(list(columns)).reset_index(drop=True)


def score_measures(
    case: EvalCase, returned: pa.Table, expected: pa.Table, *, ordered: bool
) -> CaseResult:
    """Returns whether a translation measured what the reference measures.

    Args:
        case: The case.
        returned: What the translated query returned.
        expected: What the reference query returned.
        ordered: Whether the reference sorts, in which case the sequence has to match
            and not only the rows. A question asking for the highest of something is
            answered wrongly by the right rows in the wrong order.

    Returns:
        The result. Values compare within pandas' default relative tolerance, which is
        what lets two spellings of one average agree, and dtypes are not compared: a
        count is the same answer whether it arrives as an integer or a double.
    """
    if returned.num_columns != expected.num_columns:
        return CaseResult(
            case,
            passed=False,
            detail=(
                f"measured {returned.num_columns} columns where the reference "
                f"measures {expected.num_columns}"
            ),
        )
    columns = expected.to_pandas().columns
    try:
        assert_frame_equal(
            _frame(returned, columns, ordered=ordered),
            _frame(expected, columns, ordered=ordered),
            check_dtype=False,
        )
    except AssertionError as difference:
        first = str(difference).strip().splitlines()[0]
        return CaseResult(
            case, passed=False, detail=f"{expected.num_rows} rows expected; {first}"
        )
    return CaseResult(case, passed=True, detail=f"{expected.num_rows} rows")


def score(
    case: EvalCase, returned: Sequence[str], expected: Sequence[str]
) -> CaseResult:
    """Returns whether a translation answered with the reactions the reference does.

    Args:
        case: The case.
        returned: Reaction IDs the translated query returned.
        expected: Reaction IDs the case's reference query returned.

    Returns:
        The result, naming what the translation missed and what it added. Both
        directions matter: a query that is too narrow misses reactions, and one that
        drops a condition adds them.
    """
    found, right = set(returned), set(expected)
    missing, extra = right - found, found - right
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missed {len(missing)} ({_sample(missing)})")
        if extra:
            parts.append(f"wrongly returned {len(extra)} ({_sample(extra)})")
        return CaseResult(
            case, passed=False, detail=f"{len(right)} expected; " + ", ".join(parts)
        )
    return CaseResult(case, passed=True, detail=f"{len(right)} reactions")


def _score_against_reference(
    case: EvalCase,
    corpus: execute.Corpus,
    table: pa.Table,
    *,
    timeout_seconds: float,
) -> CaseResult:
    """Returns how the translated query's rows compare with the reference query's.

    Args:
        case: The case, which carries the reference.
        corpus: The corpus to run the reference against.
        table: What the model's query returned.
        timeout_seconds: Bounds the reference search.

    Returns:
        The result. A question that asks for reactions is scored on the reactions, and
        one that asks for a number is scored on the measured rows -- the reference says
        which by whether it aggregates or sorts. A reference query that will not run
        fails its own case rather than the run: that is the case file being wrong, and
        the other cases still have something to say.
    """
    assert case.reference is not None  # A case that compiles carries one.
    reference = query.Query.model_validate(case.reference)
    try:
        expected = corpus.search(reference, timeout_seconds=timeout_seconds)
    except Exception as error:  # noqa: BLE001
        return CaseResult(
            case, passed=False, detail=f"the reference query failed: {error}"
        )
    if reference.aggregate is not None or reference.order_by:
        return score_measures(case, table, expected, ordered=bool(reference.order_by))
    if "reaction_id" not in table.schema.names:
        # The reference asks which reactions; a translation that grouped them answered
        # a different question, and there is nothing to compare it against.
        return CaseResult(
            case,
            passed=False,
            detail="measured the reactions where the question asks which ones",
        )
    return score(
        case,
        table.column("reaction_id").to_pylist(),
        expected.column("reaction_id").to_pylist(),
    )


def run_case(
    case: EvalCase,
    corpus: execute.Corpus,
    *,
    client: anthropic.Anthropic,
    model: str,
    repair: bool,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    sink: nl_log.Sink | None = None,
    session_id: str | None = None,
) -> CaseResult:
    """Translates one case, runs it, and scores what came back.

    Args:
        case: The case to run.
        corpus: The corpus to search.
        client: Anthropic client.
        model: Which model translates.
        repair: Whether a failure gets the repair turn.
        timeout_seconds: Bounds the search, so a case that translates into something
            slow fails the run rather than stopping it.
        sink: Where each case's record goes; nothing is recorded if omitted. A run is
            worth recording for the same reason a served question is: a real question
            against a real model, priced and timed.
        session_id: Groups this run's cases, the way a visit groups a reader's.

    Returns:
        The result. A case marked ``compiles: false`` passes exactly when the model
        declines outright, since saying so is the right answer to what the grammar
        cannot express. Declining only after a query failed to compile does not pass:
        the case measures whether the model reads the question, not whether it
        eventually stops.

    Raises:
        ModelRateLimitedError: If the caller is over its rate limit.
        ModelUnavailableError: If the model cannot be reached. Neither is scored as a
            failed translation: the measurement did not happen, and calling that a
            wrong answer would understate the model on the next run. Everything else --
            a query that names an unresolvable compound, a search past its timeout, a
            reference query the case file got wrong -- fails its own case and lets the
            rest of the run finish.
    """

    def record(outcome: nl_log.Outcome, **fields: Any) -> None:
        """Writes one record for this case."""
        nl_log.emit(
            sink,
            nl_log.Ask(
                question=case.question,
                model=model,
                corpus_fingerprint=corpus.fingerprint,
                prompt_fingerprint=nl.prompt_fingerprint(),
                outcome=outcome,
                source=SOURCE,
                session_id=session_id,
                **fields,
            ),
        )

    # A rate limit or an unreachable model is deliberately not recorded: the
    # measurement did not happen, and a record of it would read as a translation that
    # went wrong.
    try:
        translated = nl.translate(
            case.question, client=client, model=model, repair=repair
        )
    except nl.UnanswerableError as error:
        record(
            nl_log.Outcome.DECLINED,
            attempts=error.attempts,
            declined_reason=str(error),
            usage=error.usage,
            translate_ms=error.elapsed_ms,
        )
        if error.attempted:
            # The model wrote a query for a question the grammar cannot express and
            # only backed off once the compiler said so. The caller is served, but the
            # case exists to measure whether the model recognizes the question.
            return CaseResult(
                case, passed=False, detail=f"built a query, then declined: {error}"
            )
        return CaseResult(case, passed=not case.compiles, detail=f"declined: {error}")
    except nl.MalformedQueryError as error:
        record(
            nl_log.Outcome.MALFORMED,
            attempts=error.attempts,
            error=str(error),
            usage=error.usage,
            translate_ms=error.elapsed_ms,
        )
        # Never a pass, whatever the case expects. A model that built a query for an
        # inexpressible question reached the right verdict by the wrong road, and
        # counting that as a refusal would hide the day it starts answering instead.
        return CaseResult(case, passed=False, detail=f"did not compile: {error}")
    try:
        table = corpus.search(translated.query, timeout_seconds=timeout_seconds)
    except Exception as error:  # noqa: BLE001
        # One case whose query names an unresolvable compound, or runs past the
        # timeout, is one failed case. Letting it out of here would abort the run and
        # cost every other case's result along with it.
        record(
            nl.outcome_of(error),
            attempts=translated.attempts,
            error=str(error),
            usage=translated.usage,
            translate_ms=translated.elapsed_ms,
        )
        return CaseResult(case, passed=False, detail=f"the search failed: {error}")
    record(
        nl_log.Outcome.ANSWERED if table.num_rows else nl_log.Outcome.EMPTY,
        attempts=translated.attempts,
        row_count=table.num_rows,
        usage=translated.usage,
        translate_ms=translated.elapsed_ms,
    )
    if not case.compiles:
        # Recorded before the verdict, because a model answering a question the
        # grammar supposedly cannot express is the most interesting thing in a run:
        # either the model is wrong, or the case is, and the query is the evidence.
        return CaseResult(
            case, passed=False, detail="compiled, but the grammar cannot express this"
        )
    return _score_against_reference(
        case, corpus, table, timeout_seconds=timeout_seconds
    )


def report(results: Sequence[CaseResult]) -> str:
    """Returns a summary of a run, failures spelled out.

    Args:
        results: What the cases produced.

    Returns:
        A count, then a line per failure naming the question and what went wrong.
    """
    passed = sum(result.passed for result in results)
    lines = [f"{passed}/{len(results)} passed"]
    lines += [
        f"  FAIL {result.case.question}\n"
        f"       {result.detail}\n"
        f"       (this case exists to catch: {result.case.why})"
        for result in results
        if not result.passed
    ]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> None:
    """Runs the cases against a real model and prints the report.

    Args:
        argv: Command-line arguments; ``sys.argv`` when omitted.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projections", required=True, help="Glob for projections")
    parser.add_argument("--structures", required=True, help="Glob for structures")
    parser.add_argument("--model", default=nl.DEFAULT_MODEL, help="Model to translate")
    parser.add_argument("--cases", type=pathlib.Path, default=CASES, help="Cases file")
    parser.add_argument(
        "--no-repair",
        action="store_true",
        help="Measure first-try accuracy, without handing failures back",
    )
    parser.add_argument(
        "--pivots-dir",
        default=None,
        help="Directory of derived pivot artifacts, read rather than built",
    )
    parser.add_argument(
        "--pivot-budget-bytes",
        type=int,
        default=0,
        help=(
            "What pivots built in process may hold. Zero by default: a run is a "
            "handful of queries, and building a pivot over the whole corpus costs "
            "minutes to save milliseconds"
        ),
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=DEFAULT_TIMEOUT_SECONDS,
        help="What one case's search may take before it fails",
    )
    parser.add_argument(
        "--require-current",
        action="store_true",
        help="Refuse artifacts not written by this version of the library",
    )
    parser.add_argument(
        "--log",
        default=None,
        help=(
            "Append a record per case to this JSONL file, so a run's questions, "
            "queries, and token cost outlive the printed report"
        ),
    )
    args = parser.parse_args(argv)
    cases = load_cases(args.cases)
    client = nl.get_client()
    sink = None if args.log is None else nl_log.JsonlSink(args.log)
    # One run is one session, so its cases group the way a reader's questions do.
    session_id = str(uuid.uuid4())
    with execute.Corpus(
        args.projections,
        args.structures,
        require_current=args.require_current,
        pivots_dir=args.pivots_dir,
        pivot_budget_bytes=args.pivot_budget_bytes,
    ) as corpus:
        results = [
            run_case(
                case,
                corpus,
                client=client,
                model=args.model,
                repair=not args.no_repair,
                timeout_seconds=args.timeout_seconds,
                sink=sink,
                session_id=session_id,
            )
            for case in cases
        ]
    logger.info(
        "%s on %s", "with repair" if not args.no_repair else "first try", args.model
    )
    print(report(results))


if __name__ == "__main__":
    main()
