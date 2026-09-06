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

"""Checking a built corpus before anyone publishes it.

The unit tests say the deriving code is right. They say nothing about a corpus someone
actually built, which is what gets published and what a wrong answer comes out of. Three
questions, and a command that answers them:

**Is it faithful to its sources?** Every dataset has a projection, every projection
holds exactly the reactions its source does, and a sample of rows re-projects to what
is stored. A projection that quietly dropped a field would pass every other check here.

**Would a rebuild answer the same?** A baseline records what a set of queries returned,
one digest per query over the reactions it found. Rerun against a new build and any
answer that moved is named. The baseline also records which corpus it was taken from, so
running it against a *different* corpus reports that rather than reporting every query
as broken -- an ID from one corpus is not evidence about another, and a check that
forgets so becomes a check nobody believes.

**What is in it?** Counts that decide what a release can promise: how many reactions
carry a temperature, a yield, a solvent; how many structures RDKit cannot read; what the
yields actually range over. These never fail. They are here because the answer to
"should we publish this" is not a boolean, and the numbers belong beside the checks
rather than in someone's terminal history.
"""

import argparse
import dataclasses
import glob
import hashlib
import json
import pathlib
import random
from collections.abc import Iterable, Sequence
from typing import Any

import pyarrow.parquet as pq

from ord_schema import parquet
from ord_schema.artifacts import base, projection
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.search import execute, query

logger = get_logger(__name__)

BASELINE = pathlib.Path(__file__).parent / "corpus_baseline.json"

# One query per thing the corpus can be asked, so a build that broke any of them shows
# up as an answer that moved rather than as a passing check. Each says what it covers,
# because a digest that changes is only useful if someone can tell what changed.
#
# Python rather than a data file, which is what the eval cases beside this module are.
# The rule is who changes it and when: an eval case is a question somebody thought to
# ask a model, and arrives without touching the code, so it is data. This is a checklist
# of code paths, and an entry is owed the day an operator is added to the grammar -- in
# that same change, by that same person. ``check_test`` holds it to that: every operator
# the grammar can express has to appear here.
QUERIES: list[dict[str, Any]] = [
    {
        "name": "scalar_comparison",
        "covers": "a leaf predicate needing no quantifier",
        "query": {
            "where": {
                "op": "gt",
                "path": "conditions.temperature.setpoint_kelvin",
                "value": {"literal": 350},
            }
        },
    },
    {
        "name": "existential_quantifier",
        "covers": "an exists over a repeated level",
        "query": {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "eq",
                    "path": "reaction_role",
                    "value": {"literal": "SOLVENT"},
                },
            }
        },
    },
    {
        "name": "universal_quantifier",
        "covers": "a forall, which the occurrence index cannot answer",
        "query": {
            "where": {
                "op": "forall",
                "path": "inputs.components",
                "where": {
                    "op": "ne",
                    "path": "reaction_role",
                    "value": {"literal": "SOLVENT"},
                },
            }
        },
    },
    {
        "name": "correlated_conditions",
        "covers": "two conditions bound to one element, which the pivots serve",
        "query": {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "and",
                    "clauses": [
                        {
                            "op": "eq",
                            "path": "reaction_role",
                            "value": {"literal": "SOLVENT"},
                        },
                        {
                            "op": "same_compound",
                            "path": "smiles",
                            "smiles": "c1ccncc1",
                        },
                    ],
                },
            }
        },
    },
    {
        "name": "nested_quantifiers",
        "covers": "a quantifier inside a quantifier",
        "query": {
            "where": {
                "op": "exists",
                "path": "outcomes.products",
                "where": {
                    "op": "and",
                    "clauses": [
                        {
                            "op": "eq",
                            "path": "is_desired_product",
                            "value": {"literal": True},
                        },
                        {
                            "op": "exists",
                            "path": "measurements",
                            "where": {
                                "op": "and",
                                "clauses": [
                                    {
                                        "op": "eq",
                                        "path": "type",
                                        "value": {"literal": "YIELD"},
                                    },
                                    {
                                        "op": "gt",
                                        "path": "percentage.value",
                                        "value": {"literal": 50},
                                    },
                                ],
                            },
                        },
                    ],
                },
            }
        },
    },
    {
        "name": "reduction_ordering",
        "covers": "a list aggregate over a repeated path, and ORDER BY over it",
        "query": {
            "order_by": [
                {
                    "key": {
                        "reduce": "max",
                        "path": "outcomes.products.measurements.percentage.value",
                    },
                    "descending": True,
                }
            ],
            "limit": 25,
        },
    },
    {
        # Averaged rather than ordered, because every ordering over this path answers
        # the same 25 reactions whatever the reducer and whatever the filter: the top
        # is held by a handful of percentages recorded as 9.0e19, which are yields. A
        # query whose digest cannot move is coverage in name only.
        "name": "filtered_reduction",
        "covers": "a reduction narrowed to the elements its where keeps",
        "query": {
            "aggregate": {
                "measures": [
                    {
                        "fn": "avg",
                        "path": {
                            "reduce": "count",
                            "path": "outcomes.products.measurements.percentage.value",
                            "where": {
                                "op": "eq",
                                "path": "type",
                                "value": {"literal": "YIELD"},
                            },
                        },
                        "name": "yields_per_reaction",
                    }
                ]
            }
        },
    },
    {
        "name": "temporal_range",
        "covers": "an ordered comparison against a typed timestamp, at both ends",
        "query": {
            "where": {
                "op": "and",
                "clauses": [
                    {
                        "op": "ge",
                        "path": "provenance.record_created.time.timestamp",
                        "value": {"literal": "2020-01-01"},
                    },
                    {
                        "op": "lt",
                        "path": "provenance.record_created.time.timestamp",
                        "value": {"literal": "2021-01-01"},
                    },
                ],
            }
        },
    },
    {
        "name": "aggregate_grouping",
        "covers": "GROUP BY with a measure",
        "query": {
            "aggregate": {
                "group_by": ["conditions.stirring.type"],
                "measures": [{"fn": "count", "name": "reactions"}],
            }
        },
    },
    {
        "name": "substructure_screen",
        "covers": "the pattern-fingerprint screen and its verification",
        "query": {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "substructure",
                    "path": "smiles",
                    "smarts": "c1ccncc1",
                },
            }
        },
    },
    {
        "name": "similarity_screen",
        "covers": "the Morgan fingerprint and the popcount bound",
        "query": {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "similarity",
                    "path": "smiles",
                    "smiles": "CCO",
                    "threshold": 0.6,
                },
            }
        },
    },
    {
        "name": "compound_identity",
        "covers": "the mol_hash column",
        "query": {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "same_compound",
                    "path": "smiles",
                    "smiles": "CC(=O)O",
                },
            }
        },
    },
    {
        "name": "parent_identity",
        "covers": "the parent_hash column, and the salts mol_hash does not reach",
        "query": {
            "where": {
                "op": "exists",
                "path": "inputs.components",
                "where": {
                    "op": "same_parent",
                    "path": "smiles",
                    "smiles": "CC(=O)O",
                },
            }
        },
    },
    {
        "name": "bounded_range",
        "covers": "ge and lt, the comparisons a range is written from",
        "query": {
            "where": {
                "op": "and",
                "clauses": [
                    {
                        "op": "ge",
                        "path": "conditions.temperature.setpoint_kelvin",
                        "value": {"literal": 300},
                    },
                    {
                        "op": "lt",
                        "path": "conditions.temperature.setpoint_kelvin",
                        "value": {"literal": 350},
                    },
                ],
            }
        },
    },
    {
        "name": "at_most",
        "covers": "le, which a bound reaches through and a range does not",
        "query": {
            "where": {
                "op": "le",
                "path": "conditions.temperature.setpoint_kelvin",
                "value": {"literal": 273},
            }
        },
    },
    {
        "name": "disjunction",
        "covers": "or, whose clauses compile the same way and combine differently",
        "query": {
            "where": {
                "op": "or",
                "clauses": [
                    {
                        "op": "eq",
                        "path": "conditions.stirring.type",
                        "value": {"literal": "STIR_BAR"},
                    },
                    {
                        "op": "eq",
                        "path": "conditions.stirring.type",
                        "value": {"literal": "OVERHEAD_MIXER"},
                    },
                ],
            }
        },
    },
    {
        "name": "negation",
        "covers": "not, which the quantifier answers a different way",
        "query": {
            "where": {
                "op": "not",
                "clause": {
                    "op": "exists",
                    "path": "inputs.components",
                    "where": {
                        "op": "eq",
                        "path": "reaction_role",
                        "value": {"literal": "SOLVENT"},
                    },
                },
            }
        },
    },
    {
        "name": "presence",
        "covers": "not_null, the other half of the is_null pair",
        "query": {"where": {"op": "not_null", "path": "provenance.doi"}},
    },
    {
        "name": "text_edges",
        "covers": "starts_with and ends_with, which contains does not exercise",
        "query": {
            "where": {
                "op": "or",
                "clauses": [
                    {
                        "op": "starts_with",
                        "path": "provenance.doi",
                        "value": {"literal": "10."},
                    },
                    {
                        "op": "ends_with",
                        "path": "provenance.doi",
                        "value": {"literal": "x"},
                    },
                ],
            }
        },
    },
    {
        "name": "text_predicate",
        "covers": "a substring match on a text column",
        "query": {
            "where": {
                "op": "contains",
                "path": "provenance.city",
                "value": {"literal": "Cambridge"},
            }
        },
    },
    {
        "name": "null_check",
        "covers": "the difference between an unset field and a zero",
        "query": {
            "where": {
                "op": "is_null",
                "path": "conditions.temperature.setpoint_kelvin",
            }
        },
    },
]


@dataclasses.dataclass(frozen=True)
class Finding:
    """One thing checked, and how it came out.

    Attributes:
        check: What was checked, in a few words.
        passed: Whether it holds. A coverage count is always True; it reports rather
            than judges.
        detail: What was found, whether it passed or not.
    """

    check: str
    passed: bool
    detail: str


def digest(reaction_ids: Iterable[str]) -> str:
    """Returns a digest of a result, independent of the order it came back in.

    Args:
        reaction_ids: The reactions a query returned.

    Returns:
        A hex digest over the sorted IDs. Sorted because two runs of one query need not
        agree on order and a difference in order is not a difference in the answer.
    """
    hasher = hashlib.sha256()
    for value in sorted(reaction_ids):
        hasher.update(value.encode())
        hasher.update(b"\x00")
    return hasher.hexdigest()


def corpus_digest(projection_pattern: str) -> str:
    """Returns a digest of which datasets a corpus is built from.

    Recorded beside a baseline so that running it against another corpus is reported as
    another corpus rather than as every answer having changed.

    Args:
        projection_pattern: Glob matching the projections.

    Returns:
        A hex digest over the sorted source hashes the projections stamp.
    """
    hashes = sorted(
        base.load_stamps(path).source_md5
        for path in glob.glob(projection_pattern, recursive=True)
    )
    hasher = hashlib.sha256()
    for value in hashes:
        hasher.update(value.encode())
        hasher.update(b"\x00")
    return hasher.hexdigest()


def _reaction_ids(path: str) -> set[str]:
    """Returns the reaction IDs a projection holds."""
    return set(pq.read_table(path, columns=["reaction_id"]).column(0).to_pylist())


def check_sources(
    projection_pattern: str,
    source_pattern: str,
    *,
    datasets: int,
    rows: int,
    seed: int,
) -> list[Finding]:
    """Returns whether the projections restate the datasets they were derived from.

    Args:
        projection_pattern: Glob matching the projections.
        source_pattern: Glob matching the source datasets.
        datasets: How many datasets to sample rows from. Re-projecting a row costs what
            projecting it did, so this is the knob between a check and a rebuild.
        rows: How many rows to re-project per sampled dataset.
        seed: Fixes which datasets and rows are sampled, so a failure can be rerun.

    Returns:
        One finding for the pairing, one for the reaction IDs, and one for the sampled
        rows.
    """
    stamped = {
        path: base.load_stamps(path).source_dataset_id
        for path in sorted(glob.glob(projection_pattern, recursive=True))
    }
    named = {
        path: parquet.DatasetView(path).dataset_id
        for path in sorted(glob.glob(source_pattern, recursive=True))
    }
    projections = {value: path for path, value in stamped.items() if value is not None}
    sources = {value: path for path, value in named.items() if value is not None}
    findings = []
    anonymous = sorted(
        [path for path, value in stamped.items() if value is None]
        + [path for path, value in named.items() if value is None]
    )
    findings.append(
        Finding(
            "every file names the dataset it holds",
            passed=not anonymous,
            detail="all named"
            if not anonymous
            # Nothing can be paired with a file that does not say which dataset it is,
            # so it is not merely unchecked here -- it is unreachable by any check.
            else f"no dataset ID, so nothing can be paired with them: {anonymous}",
        )
    )
    unpaired = sorted(set(sources) ^ set(projections))
    findings.append(
        Finding(
            "every dataset has a projection",
            not unpaired,
            f"{len(sources)} datasets, {len(projections)} projections"
            + (f"; unpaired: {unpaired}" if unpaired else ""),
        )
    )
    paired = sorted(set(sources) & set(projections))
    mismatched = []
    total = 0
    for dataset_id in paired:
        view = parquet.DatasetView(sources[dataset_id])
        expected = set(view.iter_reaction_ids())
        found = _reaction_ids(projections[dataset_id])
        total += len(expected)
        if expected != found:
            mismatched.append(
                f"{dataset_id}: {len(expected - found)} missing, "
                f"{len(found - expected)} unexpected"
            )
    findings.append(
        Finding(
            "every reaction is projected exactly once",
            not mismatched,
            f"{total:,} reactions across {len(paired)} datasets"
            + (f"; {mismatched}" if mismatched else ""),
        )
    )
    findings.append(
        _check_round_trip(paired, sources, projections, datasets, rows, seed)
    )
    return findings


def _check_round_trip(
    paired: Sequence[str],
    sources: dict[str, str],
    projections: dict[str, str],
    datasets: int,
    rows: int,
    seed: int,
) -> Finding:
    """Returns whether sampled reactions re-project to what the projection stores.

    Args:
        paired: Dataset IDs present on both sides.
        sources: Dataset ID to source path.
        projections: Dataset ID to projection path.
        datasets: How many datasets to sample.
        rows: How many rows per dataset.
        seed: Fixes the sample.

    Returns:
        The finding, naming the first few rows that came out differently.
    """
    # S311: choosing which rows to spot-check, which wants reproducibility from a seed
    # rather than unpredictability.
    chooser = random.Random(seed)  # noqa: S311
    chosen = chooser.sample(list(paired), min(datasets, len(paired)))
    differing: list[str] = []
    checked = 0
    for dataset_id in sorted(chosen):
        view = parquet.DatasetView(sources[dataset_id])
        stored = pq.read_table(projections[dataset_id]).to_pylist()
        positions = chooser.sample(range(len(stored)), min(rows, len(stored)))
        for position in positions:
            row = stored[position]
            reaction = view.reactions[position]
            checked += 1
            try:
                # structure_id is assigned in first-seen order across a whole dataset,
                # so one row cannot reproduce it; the projection's own tests cover it.
                expected = projection.reaction_row(reaction, row["reaction_id"])
            except ValueError as error:
                # Projecting refuses a row whose ID disagrees with the reaction it
                # stores. That is the check finding something, not the check breaking.
                differing.append(f"{dataset_id}[{position}] does not project: {error}")
                continue
            if _without_structure_ids(expected) != _without_structure_ids(row):
                differing.append(f"{dataset_id}[{position}] {row['reaction_id']}")
    return Finding(
        "sampled reactions re-project to what is stored",
        not differing,
        f"{checked} rows across {len(chosen)} datasets"
        + (f"; differing: {differing[:5]}" if differing else ""),
    )


def _without_structure_ids(value: Any) -> Any:
    """Returns the value with every ``structure_id`` dropped, at any depth.

    Args:
        value: A projected row, or any part of one.

    Returns:
        The same structure without the key. A structure ID is a position in a sibling
        artifact rather than a fact about the reaction, and it is assigned across a
        whole dataset, so a single row re-projected on its own cannot carry one.
    """
    if isinstance(value, dict):
        return {
            key: _without_structure_ids(item)
            for key, item in value.items()
            if key != "structure_id"
        }
    if isinstance(value, list):
        return [_without_structure_ids(item) for item in value]
    if isinstance(value, tuple):
        # A map column reads back as (key, value) pairs, and the values under it are
        # where a component's structure_id lives.
        return tuple(_without_structure_ids(item) for item in value)
    return value


def measure(corpus: execute.Corpus, *, timeout_seconds: float) -> dict[str, Any]:
    """Returns what each canonical query answers on this corpus.

    Args:
        corpus: The corpus to ask.
        timeout_seconds: Bounds each query, so a slow one fails the run rather than
            stopping it.

    Returns:
        A mapping from query name to the reactions found and their digest.
    """
    results = {}
    for entry in QUERIES:
        table = corpus.search(
            query.Query.model_validate(entry["query"]), timeout_seconds=timeout_seconds
        )
        column = (
            table.column("reaction_id").to_pylist()
            if "reaction_id" in table.column_names
            else [json.dumps(row, default=str) for row in table.to_pylist()]
        )
        results[entry["name"]] = {"rows": len(column), "digest": digest(column)}
        logger.info("%s answered %d rows", entry["name"], len(column))
    return results


def check_answers(measured: dict[str, Any], baseline: dict[str, Any]) -> list[Finding]:
    """Returns whether a corpus answers as the baseline recorded.

    Args:
        measured: What ``measure`` found.
        baseline: A baseline file's contents.

    Returns:
        One finding per query, plus one naming any query the baseline does not carry.
        A baseline taken from another corpus produces a single finding saying so: its
        digests are not evidence about this corpus, and comparing them anyway would
        report every query as broken.
    """
    recorded = baseline.get("queries", {})
    unrecorded = sorted(set(measured) - set(recorded))
    findings = []
    for name in sorted(set(measured) & set(recorded)):
        was, now = recorded[name], measured[name]
        same = was["digest"] == now["digest"]
        findings.append(
            Finding(
                f"answers unchanged: {name}",
                passed=same,
                detail=f"{now['rows']:,} rows"
                if same
                else f"{was['rows']:,} rows before, {now['rows']:,} now",
            )
        )
    if unrecorded:
        findings.append(
            Finding(
                "every query is in the baseline",
                passed=False,
                detail=f"not recorded, so nothing to compare: {unrecorded}",
            )
        )
    return findings


def check_coverage(corpus: execute.Corpus, *, timeout_seconds: float) -> list[Finding]:
    """Returns counts describing what the corpus holds.

    These never fail. What a release can promise is a judgment about the numbers, and
    the numbers belong beside the checks rather than in someone's terminal history.

    Args:
        corpus: The corpus to describe.
        timeout_seconds: Bounds each query.

    Returns:
        One finding per count, all passing.
    """
    counts: list[Finding] = []

    def count(label: str, where: dict[str, Any] | None) -> int:
        payload: dict[str, Any] = {
            "aggregate": {"measures": [{"fn": "count", "name": "n"}]}
        }
        if where is not None:
            payload["where"] = where
        table = corpus.search(
            query.Query.model_validate(payload), timeout_seconds=timeout_seconds
        )
        found = int(table.column("n").to_pylist()[0])
        counts.append(Finding(label, passed=True, detail=f"{found:,}"))
        return found

    total = count("reactions", None)
    if not total:
        # Every check below this reports a share of nothing, and every check above it
        # passed by having nothing to disagree about. A corpus nobody can query is not
        # a corpus that passed: say so rather than divide by it.
        counts.append(
            Finding(
                "the corpus holds reactions",
                passed=False,
                detail="no reactions, so every other check passed vacuously",
            )
        )
        return counts

    def share(label: str, where: dict[str, Any]) -> None:
        found = count(label, where)
        counts[-1] = Finding(
            label, passed=True, detail=f"{found:,} ({100 * found / total:.1f}%)"
        )

    share(
        "with a temperature setpoint",
        {"op": "not_null", "path": "conditions.temperature.setpoint_kelvin"},
    )
    share(
        "with a solvent",
        {
            "op": "exists",
            "path": "inputs.components",
            "where": {
                "op": "eq",
                "path": "reaction_role",
                "value": {"literal": "SOLVENT"},
            },
        },
    )
    share(
        "with a yield",
        {
            "op": "exists",
            "path": "outcomes.products.measurements",
            "where": {"op": "eq", "path": "type", "value": {"literal": "YIELD"}},
        },
    )
    extremes = corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "measures": [
                        {
                            "fn": "min",
                            "path": {
                                "reduce": "min",
                                "path": (
                                    "outcomes.products.measurements.percentage.value"
                                ),
                            },
                            "name": "lowest",
                        },
                        {
                            "fn": "max",
                            "path": {
                                "reduce": "max",
                                "path": (
                                    "outcomes.products.measurements.percentage.value"
                                ),
                            },
                            "name": "highest",
                        },
                    ]
                }
            }
        ),
        timeout_seconds=timeout_seconds,
    ).to_pylist()[0]
    # Reported rather than asserted: the corpus holds percentages that are not
    # percentages, and anything that ranks by yield will surface them. Whether that is
    # a release blocker is a judgment, and this is the number it is made from.
    counts.append(
        Finding(
            "yield range, as recorded",
            passed=True,
            detail=f"{extremes['lowest']} to {extremes['highest']} percent",
        )
    )
    return counts


def check_structures(structures_pattern: str) -> list[Finding]:
    """Returns counts of what RDKit could and could not read.

    Args:
        structures_pattern: Glob matching the structures artifacts.

    Returns:
        One finding for the readable share, one for the hashed share. Both pass; a
        corpus with unreadable structures is normal and the number is what matters.
    """
    total = readable = hashed = 0
    for path in sorted(glob.glob(structures_pattern, recursive=True)):
        table = pq.read_table(path, columns=["mol_binary", "mol_hash"])
        total += table.num_rows
        readable += table.num_rows - table.column("mol_binary").null_count
        hashed += table.num_rows - table.column("mol_hash").null_count
    if not total:
        return [
            Finding(
                "structures readable",
                passed=False,
                detail="no structures artifacts matched",
            )
        ]
    return [
        Finding(
            "structures RDKit can read",
            passed=True,
            detail=f"{readable:,} of {total:,} ({100 * readable / total:.2f}%)",
        ),
        Finding(
            "structures with a compound hash",
            passed=True,
            detail=f"{hashed:,} of {total:,} ({100 * hashed / total:.2f}%)",
        ),
    ]


def report(findings: Sequence[Finding]) -> str:
    """Returns the findings as lines, failures first.

    Args:
        findings: What the checks produced.

    Returns:
        A summary line, then every finding.
    """
    failures = [finding for finding in findings if not finding.passed]
    lines = [f"{len(findings) - len(failures)}/{len(findings)} checks passed"]
    lines += [f"  FAIL {value.check}\n        {value.detail}" for value in failures]
    lines += [
        f"  ok   {value.check}: {value.detail}" for value in findings if value.passed
    ]
    return "\n".join(lines)


def load_baseline(path: pathlib.Path) -> dict[str, Any]:
    """Returns a baseline file's contents.

    Args:
        path: The file to read.

    Returns:
        The parsed baseline.

    Raises:
        FileNotFoundError: If no baseline has been recorded yet, which is a different
            thing from a corpus that answers differently.
    """
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def main(argv: Sequence[str] | None = None) -> int:
    """Checks a built corpus and prints what it found.

    Args:
        argv: Command-line arguments; ``sys.argv`` when omitted.

    Returns:
        Zero if every check passed, one otherwise, so a pipeline can gate on it.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projections", required=True, help="Glob for projections")
    parser.add_argument("--structures", required=True, help="Glob for structures")
    parser.add_argument(
        "--sources",
        default=None,
        help="Glob for the source datasets. Without it the fidelity checks are skipped",
    )
    parser.add_argument(
        "--baseline", type=pathlib.Path, default=BASELINE, help="Recorded answers"
    )
    parser.add_argument(
        "--write-baseline",
        action="store_true",
        help="Record what this corpus answers instead of comparing against a record",
    )
    parser.add_argument(
        "--sample-datasets",
        type=int,
        default=5,
        help="How many datasets to re-project rows from",
    )
    parser.add_argument(
        "--sample-rows",
        type=int,
        default=20,
        help="How many rows to re-project per sampled dataset",
    )
    parser.add_argument("--seed", type=int, default=0, help="Fixes which rows are read")
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=600.0,
        help="What one query may take before it fails",
    )
    parser.add_argument(
        "--pivot-budget-bytes",
        type=int,
        default=0,
        help="What pivots built in process may hold; zero builds none",
    )
    parser.add_argument(
        "--pivots",
        default=None,
        help=(
            "Directory holding pivot artifacts. Without it a reduction reads the "
            "projection's lists, which answers the same but leaves the pivoted route "
            "unmeasured"
        ),
    )
    args = parser.parse_args(argv)
    # Re-projecting a sampled reaction runs the same RDKit the derivation did, and it
    # says the same things about atom maps it said then.
    silence_rdkit_logs()

    findings: list[Finding] = []
    if args.sources:
        findings += check_sources(
            args.projections,
            args.sources,
            datasets=args.sample_datasets,
            rows=args.sample_rows,
            seed=args.seed,
        )
    findings += check_structures(args.structures)
    with execute.Corpus(
        args.projections,
        args.structures,
        require_current=True,
        pivot_budget_bytes=args.pivot_budget_bytes,
        pivots_dir=args.pivots,
        # Unbounded: a baseline is a count and a digest over every matching reaction,
        # and under the corpus's own default it would be a digest over an arbitrary
        # page of them -- which compares unequal against the recorded one and reports a
        # regression in the corpus rather than in the bound that produced it.
        max_rows=None,
    ) as corpus:
        measured = measure(corpus, timeout_seconds=args.timeout_seconds)
        if args.write_baseline:
            args.baseline.write_text(
                json.dumps(
                    {
                        "corpus": corpus_digest(args.projections),
                        "queries": measured,
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            logger.info("recorded %d answers to %s", len(measured), args.baseline)
        else:
            baseline = load_baseline(args.baseline)
            recorded_for = baseline.get("corpus")
            here = corpus_digest(args.projections)
            if recorded_for != here:
                findings.append(
                    Finding(
                        "the baseline is about this corpus",
                        passed=False,
                        detail=(
                            "recorded against another set of datasets, so its answers "
                            "say nothing about this one; rerun with --write-baseline "
                            "if this corpus is the one to record"
                        ),
                    )
                )
            else:
                findings += check_answers(measured, baseline)
        findings += check_coverage(corpus, timeout_seconds=args.timeout_seconds)
    print(report(findings))
    return 0 if all(finding.passed for finding in findings) else 1


if __name__ == "__main__":
    raise SystemExit(main())
