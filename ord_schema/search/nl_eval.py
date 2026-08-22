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

A case does not state the query it expects, because several spellings of a question are
right and pinning one would fail a better translation than the one written the day the
case was added. It states reactions any correct query returns, and reactions a plausible
wrong one returns -- usually the near-miss where two conditions on the same element
become two quantifiers, which matches a solvent here and a reactant there.

Running the cases costs money and reaches the network, so the suite covers the
scoring and the outcomes ``run_case`` maps a translation onto, with the model stubbed.
The measurement itself is a command.
"""

import argparse
import dataclasses
import json
import pathlib
from collections.abc import Sequence

import anthropic
from pydantic import BaseModel, Field

from ord_schema.logging import get_logger
from ord_schema.search import execute, nl

logger = get_logger(__name__)

CASES = pathlib.Path(__file__).parent / "nl_cases.json"
# Longer than a served query would wait, because a run is unattended and a case whose
# translation is slow is worth a verdict rather than a stopped harness.
DEFAULT_TIMEOUT_SECONDS = 300.0


class EvalCase(BaseModel):
    """One question, and what any correct answer to it must satisfy.

    Attributes:
        question: The question, in English.
        why: What this case is here to catch, for whoever reads a failure.
        compiles: Whether the question should translate at all. False marks one the
            grammar cannot express, which the layer must refuse rather than fudge.
        must_return: Reaction IDs any correct query returns.
        must_not_return: Reaction IDs a plausible wrong query returns and a correct one
            does not.
    """

    question: str
    why: str
    compiles: bool = True
    must_return: list[str] = Field(default_factory=list)
    must_not_return: list[str] = Field(default_factory=list)


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


def score(case: EvalCase, returned: Sequence[str]) -> CaseResult:
    """Returns whether the reactions a query returned satisfy a case.

    Args:
        case: The case.
        returned: Reaction IDs the translated query returned.

    Returns:
        The result, naming what was missing or wrongly present.
    """
    found = set(returned)
    missing = [value for value in case.must_return if value not in found]
    forbidden = [value for value in case.must_not_return if value in found]
    if missing:
        return CaseResult(case, passed=False, detail=f"must_return absent: {missing}")
    if forbidden:
        return CaseResult(
            case, passed=False, detail=f"must_not_return present: {forbidden}"
        )
    return CaseResult(case, passed=True, detail=f"{len(found)} reactions")


def run_case(
    case: EvalCase,
    corpus: execute.Corpus,
    *,
    client: anthropic.Anthropic,
    model: str,
    repair: bool,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
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
            wrong answer would understate the model on the next run.
    """
    try:
        translated = nl.translate(
            case.question, client=client, model=model, repair=repair
        )
    except nl.UnanswerableError as error:
        if error.attempted:
            # The model wrote a query for a question the grammar cannot express and
            # only backed off once the compiler said so. The caller is served, but the
            # case exists to measure whether the model recognizes the question.
            return CaseResult(
                case, passed=False, detail=f"built a query, then declined: {error}"
            )
        return CaseResult(case, passed=not case.compiles, detail=f"declined: {error}")
    except nl.MalformedQueryError as error:
        # Never a pass, whatever the case expects. A model that built a query for an
        # inexpressible question reached the right verdict by the wrong road, and
        # counting that as a refusal would hide the day it starts answering instead.
        return CaseResult(case, passed=False, detail=f"did not compile: {error}")
    if not case.compiles:
        return CaseResult(
            case, passed=False, detail="compiled, but the grammar cannot express this"
        )
    table = corpus.search(translated, timeout_seconds=timeout_seconds)
    return score(case, table.column("reaction_id").to_pylist())


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
    args = parser.parse_args(argv)
    cases = load_cases(args.cases)
    client = nl.get_client()
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
            )
            for case in cases
        ]
    logger.info(
        "%s on %s", "with repair" if not args.no_repair else "first try", args.model
    )
    print(report(results))


if __name__ == "__main__":
    main()
