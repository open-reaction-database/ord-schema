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

"""Tests for ord_schema.search.check.

Built on a two-dataset corpus small enough to derive in the suite, so the checks are
exercised against artifacts rather than against mocks -- which is the whole point of
them: they exist to catch what a built corpus gets wrong.
"""

import json
import pathlib
import typing

import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import projection, structures
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.search import check, execute, query

_ROLE = reaction_pb2.ReactionRole


def _reaction(reaction_id: str, smiles: str, role) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    component = reaction.inputs["in"].components.add()
    component.identifiers.add(type="SMILES", value=smiles)
    component.reaction_role = role
    reaction.conditions.temperature.setpoint.value = 400
    reaction.conditions.temperature.setpoint.units = 3  # KELVIN
    return reaction


@pytest.fixture(scope="module")
def corpus_root(tmp_path_factory) -> pathlib.Path:
    """A two-dataset corpus, derived the way a real one is."""
    root = tmp_path_factory.mktemp("check")
    shards = {
        "aa": [
            _reaction("ord-aa01", "c1ccncc1", _ROLE.SOLVENT),
            _reaction("ord-aa02", "CCO", _ROLE.REACTANT),
        ],
        "bb": [_reaction("ord-bb01", "CC(=O)[O-].[Na+]", _ROLE.REAGENT)],
    }
    for shard, reactions in shards.items():
        source = root / "data" / f"ord_dataset-{shard}.parquet"
        source.parent.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            dataset_pb2.Dataset(
                dataset_id=f"ord_dataset-{shard}",
                name="test",
                description="test",
                reactions=reactions,
            ),
            str(source),
        )
        projected = root / "projections" / source.name
        projected.parent.mkdir(parents=True, exist_ok=True)
        projection.write_projection(source, projected)
        structured = root / "structures" / source.name
        structured.parent.mkdir(parents=True, exist_ok=True)
        structures.write_structures(projected, structured)
    return root


def _globs(root: pathlib.Path) -> tuple[str, str, str]:
    return (
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        str(root / "data" / "*.parquet"),
    )


def _failures(findings) -> list[str]:
    return [finding.check for finding in findings if not finding.passed]


def test_a_faithful_corpus_passes_every_fidelity_check(corpus_root):
    projections, _, sources = _globs(corpus_root)
    findings = check.check_sources(projections, sources, datasets=2, rows=5, seed=0)
    assert _failures(findings) == []


def test_a_projection_missing_a_reaction_is_caught(corpus_root, tmp_path):
    # The failure a row count alone would miss: the right number of reactions, one of
    # them the wrong one.
    _, _, sources = _globs(corpus_root)
    altered = tmp_path / "projections"
    altered.mkdir()
    for path in sorted(pathlib.Path(corpus_root / "projections").glob("*.parquet")):
        table = pq.read_table(path)
        if "aa" in path.name:
            ids = table.column("reaction_id").to_pylist()
            ids[0] = "ord-nowhere"
            table = table.set_column(
                table.schema.get_field_index("reaction_id"),
                "reaction_id",
                [ids],
            )
        pq.write_table(
            table.replace_schema_metadata(pq.read_schema(path).metadata),
            altered / path.name,
        )
    findings = check.check_sources(
        str(altered / "*.parquet"), sources, datasets=2, rows=5, seed=0
    )
    assert "every reaction is projected exactly once" in _failures(findings)


def test_a_source_with_no_projection_is_caught(corpus_root, tmp_path):
    _, _, sources = _globs(corpus_root)
    only_one = tmp_path / "projections"
    only_one.mkdir()
    first = sorted(pathlib.Path(corpus_root / "projections").glob("*.parquet"))[0]
    (only_one / first.name).write_bytes(first.read_bytes())
    findings = check.check_sources(
        str(only_one / "*.parquet"), sources, datasets=1, rows=5, seed=0
    )
    assert "every dataset has a projection" in _failures(findings)


def test_the_structure_counts_are_reported(corpus_root):
    _, structures_glob, _ = _globs(corpus_root)
    findings = check.check_structures(structures_glob)
    assert _failures(findings) == []
    assert all("100.00%" in finding.detail for finding in findings)


def test_a_digest_ignores_order():
    assert check.digest(["b", "a"]) == check.digest(["a", "b"])


def test_a_digest_separates_its_values():
    # Without a separator "ab" + "c" and "a" + "bc" would hash alike, and two different
    # answers would compare equal.
    assert check.digest(["ab", "c"]) != check.digest(["a", "bc"])


def test_an_unchanged_corpus_answers_as_recorded(corpus_root):
    projections, structures_glob, _ = _globs(corpus_root)
    with execute.Corpus(projections, structures_glob, pivot_budget_bytes=0) as corpus:
        measured = check.measure(corpus, timeout_seconds=60)
    baseline = {"corpus": check.corpus_digest(projections), "queries": measured}
    assert _failures(check.check_answers(measured, baseline)) == []


def test_an_answer_that_moved_is_named(corpus_root):
    projections, structures_glob, _ = _globs(corpus_root)
    with execute.Corpus(projections, structures_glob, pivot_budget_bytes=0) as corpus:
        measured = check.measure(corpus, timeout_seconds=60)
    baseline = {"corpus": "x", "queries": json.loads(json.dumps(measured))}
    baseline["queries"]["scalar_comparison"]["digest"] = "something else"
    failures = _failures(check.check_answers(measured, baseline))
    assert failures == ["answers unchanged: scalar_comparison"]


def test_a_query_the_baseline_does_not_carry_is_named():
    measured = {"new_query": {"rows": 1, "digest": "a"}}
    findings = check.check_answers(measured, {"queries": {}})
    assert _failures(findings) == ["every query is in the baseline"]


def test_every_canonical_query_compiles_and_runs(corpus_root):
    # A query naming a path the schema lacks would otherwise fail only on the day
    # someone ran the check against a real corpus.
    projections, structures_glob, _ = _globs(corpus_root)
    with execute.Corpus(projections, structures_glob, pivot_budget_bytes=0) as corpus:
        measured = check.measure(corpus, timeout_seconds=60)
    assert sorted(measured) == sorted(entry["name"] for entry in check.QUERIES)


def test_every_canonical_query_says_what_it_covers():
    for entry in check.QUERIES:
        assert entry["covers"], entry["name"]


def test_the_coverage_counts_never_fail(corpus_root):
    projections, structures_glob, _ = _globs(corpus_root)
    with execute.Corpus(projections, structures_glob, pivot_budget_bytes=0) as corpus:
        findings = check.check_coverage(corpus, timeout_seconds=60)
    assert _failures(findings) == []
    assert any(finding.check == "reactions" for finding in findings)


def test_the_report_puts_failures_first():
    findings = [
        check.Finding("fine", passed=True, detail="ok"),
        check.Finding("broken", passed=False, detail="not ok"),
    ]
    text = check.report(findings)
    assert text.splitlines()[0] == "1/2 checks passed"
    assert text.index("FAIL broken") < text.index("ok   fine")


def _grammar_operators() -> set[str]:
    """Returns every ``op`` the predicate grammar can express."""
    members = typing.get_args(typing.get_args(query.Predicate)[0])
    return {
        value
        for member in members
        for value in typing.get_args(member.model_fields["op"].annotation)
    }


def _operators_used(node) -> set[str]:
    """Returns the ``op`` values appearing anywhere in a query payload."""
    if isinstance(node, dict):
        used = {node["op"]} if "op" in node else set()
        for value in node.values():
            used |= _operators_used(value)
        return used
    if isinstance(node, list):
        return {value for item in node for value in _operators_used(item)}
    return set()


def _reductions_used(node) -> list[dict]:
    """Returns every reduction appearing anywhere in a query payload."""
    if isinstance(node, dict):
        found = [node] if "reduce" in node else []
        for value in node.values():
            found += _reductions_used(value)
        return found
    if isinstance(node, list):
        return [found for item in node for found in _reductions_used(item)]
    return []


def test_every_grammar_operator_is_covered_by_a_canonical_query():
    # The baseline is only as good as what it asks. An operator nothing asks about can
    # break without moving a single digest, and the day that matters is the day someone
    # adds one and forgets this list.
    covered = set()
    for entry in check.QUERIES:
        covered |= _operators_used(entry["query"])
    assert _grammar_operators() - covered == set()


def test_every_reducer_is_covered_by_a_canonical_query():
    reducers = set(query._REDUCERS)
    covered = set()
    for entry in check.QUERIES:
        payload = json.dumps(entry["query"])
        covered |= {name for name in reducers if f'"reduce": "{name}"' in payload}
    # One is enough to prove the list-aggregate path compiles and runs; the rest differ
    # only in which DuckDB function the same expression wraps, which query_test pins.
    assert covered


def test_every_field_of_a_reduction_is_covered_by_a_canonical_query():
    # A reduction compiles differently for each field it carries -- a where narrows the
    # elements on both routes -- so a field nothing asks about can break without moving
    # a digest. The same argument as the operator list above, for the same reason.
    covered = set()
    for entry in check.QUERIES:
        for reduction in _reductions_used(entry["query"]):
            covered |= set(reduction)
    assert set(query.Reduction.model_fields) - covered == set()


@pytest.fixture(scope="module")
def empty_corpus(corpus_root, tmp_path_factory):
    """A corpus whose projections hold no rows.

    Built by emptying a real one rather than by deriving an empty dataset, which
    ``save_dataset`` refuses. Unreachable through the deriving path today, which is the
    point: the coverage counts are a division, and the day a corpus can be empty is not
    the day to discover that.
    """
    root = tmp_path_factory.mktemp("empty")
    for kind in ("projections", "structures"):
        (root / kind).mkdir()
        for path in sorted((corpus_root / kind).glob("*.parquet")):
            table = pq.read_table(path)
            pq.write_table(
                table.slice(0, 0).replace_schema_metadata(table.schema.metadata),
                root / kind / path.name,
            )
    with execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        pivot_budget_bytes=0,
    ) as corpus:
        yield corpus


def test_an_empty_corpus_fails_rather_than_dividing_by_nothing(empty_corpus):
    findings = check.check_coverage(empty_corpus, timeout_seconds=60)
    assert "the corpus holds reactions" in _failures(findings)


def test_an_empty_corpus_says_why_its_other_checks_are_worthless(empty_corpus):
    # A run that reports "everything passed" over an empty corpus is the worst outcome
    # available: it looks like the answer someone wanted.
    findings = check.check_coverage(empty_corpus, timeout_seconds=60)
    failure = next(f for f in findings if not f.passed)
    assert "vacuously" in failure.detail


def test_the_baseline_is_measured_over_every_match(monkeypatch):
    # A baseline is a count and a digest over every matching reaction. Under the
    # corpus's own default bound it would be a digest over an arbitrary page, which
    # compares unequal against the recorded one and reports a regression in the corpus
    # rather than in the bound that produced it. The fixtures here hold three
    # reactions, far under any bound, so only the argument itself can say this.
    seen = {}

    class _RecordedError(Exception):
        pass

    def _corpus(*args, **kwargs):
        seen.update(kwargs)
        raise _RecordedError

    monkeypatch.setattr(check.execute, "Corpus", _corpus)
    with pytest.raises(_RecordedError):
        check.main(["--projections=x", "--structures=y"])
    assert seen["max_rows"] is None


@pytest.mark.parametrize(
    ("argv", "expected"),
    [([], None), (["--pivots=/artifacts/pivots"], "/artifacts/pivots")],
)
def test_the_pivots_directory_reaches_the_corpus(monkeypatch, argv, expected):
    # A reduction routes to a pivot only where the corpus holds one, so a baseline
    # recorded without this argument measures the list spelling and says nothing about
    # the pivoted route -- for either the answers or the coverage counts.
    seen = {}

    class _RecordedError(Exception):
        pass

    def _corpus(*args, **kwargs):
        seen.update(kwargs)
        raise _RecordedError

    monkeypatch.setattr(check.execute, "Corpus", _corpus)
    with pytest.raises(_RecordedError):
        check.main(["--projections=x", "--structures=y", *argv])
    assert seen["pivots_dir"] == expected
