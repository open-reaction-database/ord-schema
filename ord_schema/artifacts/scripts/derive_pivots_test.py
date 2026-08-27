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

"""Tests for ord_schema.scripts.derive_pivots.

Artifact behavior is covered in ord_schema.artifacts.pivot_test; these cover the CLI:
the per-level layout a Corpus reads, which matches count as sources, the skip-if-current
shortcut, how an unknown level is refused, that a projection rewritten under a running
process is counted again, that each level is written the count taken for it rather than
another level's, and that a level empty in every artifact says so and keeps saying so.
"""

import logging
import pathlib

import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import base, pivot
from ord_schema.artifacts.scripts import (
    derive_pivots,
    derive_projection,
    derive_structures,
)
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.search import execute, query


def _dataset(dataset_id: str, *, observations: bool = False):
    """Returns a one-reaction dataset.

    Args:
        dataset_id: The dataset ID, which also names the reaction.
        observations: Whether to record an observation. Off in both shards of the
            default tree, which is what makes ``observations`` a level the whole tree
            is empty at; on in one shard where a test needs the levels to differ, since
            two identical shards make "empty in every artifact" and "empty in any of
            them" the same question.

    Returns:
        The dataset.
    """
    reaction = reaction_pb2.Reaction(reaction_id=f"ord-{dataset_id}-01")
    # Three components against one workup, so a run deriving both levels pairs two
    # counts that are not interchangeable: a query built in one order and read back in
    # another would answer each level with the other's count, and levels holding the
    # same number of elements make every such pairing the identity.
    for smiles in ("c1ccccc1", "CCO", "CC(=O)O"):
        reaction.inputs["a"].components.add().identifiers.add(
            type="SMILES", value=smiles
        )
    reaction.workups.add(type="EXTRACTION")
    reaction.outcomes.add().products.add(isolated_color="white")
    if observations:
        reaction.observations.add(comment="recorded")
    return dataset_pb2.Dataset(
        dataset_id=dataset_id, name="test", description="desc", reactions=[reaction]
    )


def _reproject(root: pathlib.Path) -> None:
    """Derives the projections under ``root`` from the datasets under it."""
    derive_projection.main(
        derive_projection.parse_args(
            [
                f"--input_pattern={root / 'data' / '*' / '*.parquet'}",
                f"--output_dir={root / 'projections'}",
                "--force",
            ]
        )
    )


@pytest.fixture
def projected(tmp_path) -> pathlib.Path:
    """Writes two shards of projections, and returns the root holding them."""
    for shard in ("aa", "bb"):
        directory = tmp_path / "data" / shard
        directory.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            _dataset(f"ord_dataset-{shard}"),
            str(directory / f"ord_dataset-{shard}.parquet"),
        )
    _reproject(tmp_path)
    return tmp_path


def _run(root: pathlib.Path, *extra: str) -> None:
    derive_pivots.main(
        derive_pivots.parse_args(
            [
                f"--input_pattern={root / 'projections' / '*' / '*.parquet'}",
                f"--output_dir={root / 'pivots'}",
                *extra,
            ]
        )
    )


def test_each_level_gets_its_own_tree(projected):
    _run(projected, "--levels", "workups", "outcomes.products")
    for level in ("workups", "outcomes.products"):
        for shard in ("aa", "bb"):
            written = (
                projected / "pivots" / level / shard / f"ord_dataset-{shard}.parquet"
            )
            assert written.exists()
            assert pivot.pivot_path(written) == level
            assert base.load_stamps(written).artifact == pivot.ARTIFACT
    # A level not asked for is not derived, which is what makes naming them cheaper
    # than deriving all 39.
    assert not (projected / "pivots" / "inputs.components").exists()


def test_each_level_is_written_its_own_count(projected):
    # One query carries an aggregate per level and the row is paired back to the levels
    # by position, so the whole batching optimization rests on that pairing. The two
    # levels hold different numbers of elements, or any mispairing would be invisible.
    _run(projected, "--levels", "workups", "inputs.components")
    for level, rows in (("workups", 1), ("inputs.components", 3)):
        for shard in ("aa", "bb"):
            written = (
                projected / "pivots" / level / shard / f"ord_dataset-{shard}.parquet"
            )
            assert pq.read_table(written).num_rows == rows, level


def _warnings(caplog) -> list[str]:
    """Returns the WARNING messages recorded so far."""
    return [
        record.message for record in caplog.records if record.levelno == logging.WARNING
    ]


def test_a_level_empty_in_every_artifact_is_reported(projected, caplog):
    # Ordinary for a level this corpus never records, and identical on disk to a level
    # whose count came back wrong. Nothing else in the chain aggregates rows across a
    # tree, so without this an empty artifact is published, stamped current, and never
    # looked at again.
    with caplog.at_level(logging.WARNING):
        # The empty level named second, and derived after a level whose artifacts are
        # already on disk: named first it is also levels[0] and the only level in the
        # tree, so a report that ignored level_path, or that scanned the whole output
        # directory rather than this level's, would read the same.
        _run(projected, "--levels", "workups", "observations")
    warnings = _warnings(caplog)
    assert any("observations: every artifact is empty" in m for m in warnings)
    assert not any("workups" in m for m in warnings)


def test_a_level_empty_in_only_some_artifacts_is_not_reported(tmp_path, caplog):
    # The two shards differ at this level, which is the whole point: with identical
    # shards, "every artifact is empty" and "any artifact is empty" are the same
    # question and a warning that fired on an ordinary mixed corpus would pass.
    for shard, observations in (("aa", True), ("bb", False)):
        directory = tmp_path / "data" / shard
        directory.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            _dataset(f"ord_dataset-{shard}", observations=observations),
            str(directory / f"ord_dataset-{shard}.parquet"),
        )
    _reproject(tmp_path)
    with caplog.at_level(logging.INFO):
        _run(tmp_path, "--levels", "observations")
    rows = [
        pq.read_metadata(path).num_rows
        for path in sorted((tmp_path / "pivots" / "observations").rglob("*.parquet"))
    ]
    assert sorted(rows) == [0, 1], "the shards must differ for this to prove anything"
    assert not _warnings(caplog)
    # The tally an operator reads, with the two numbers distinct so their order shows.
    assert any(
        "observations: 2 written, 0 already current, 1 of 2 artifacts empty" in record
        for record in [r.getMessage() for r in caplog.records]
    )


def test_a_level_empty_in_every_artifact_is_reported_again_on_a_re_run(
    projected, caplog
):
    # Counted from the artifacts rather than from what the run wrote, so a re-run that
    # skips every artifact says the same thing. The signal would otherwise appear only
    # on the run that created them and never again -- vanishing exactly once the empty
    # artifacts are entrenched.
    _run(projected, "--levels", "observations")
    with caplog.at_level(logging.WARNING):
        _run(projected, "--levels", "observations")
    assert any("observations: every artifact is empty" in m for m in _warnings(caplog))


def test_a_second_run_skips_what_is_already_current(projected):
    _run(projected, "--levels", "workups")
    written = projected / "pivots" / "workups" / "aa" / "ord_dataset-aa.parquet"
    before = written.stat().st_mtime_ns
    _run(projected, "--levels", "workups")
    assert written.stat().st_mtime_ns == before


def test_force_rewrites_a_current_artifact(projected):
    _run(projected, "--levels", "workups")
    written = projected / "pivots" / "workups" / "aa" / "ord_dataset-aa.parquet"
    before = written.stat().st_mtime_ns
    _run(projected, "--levels", "workups", "--force")
    assert written.stat().st_mtime_ns != before


def test_a_rewritten_projection_is_counted_again(projected):
    # The count decides whether a level is unnested at all, and it is cached across the
    # levels of one run. A process that outlives a projection must not pivot the new
    # content against the old count: counting the workups this dataset does not yet
    # have would answer zero, skip the unnest, and publish an empty artifact stamped
    # current -- which every later run skips and every quantifier reads as no matches.
    without = _dataset("ord_dataset-aa")
    del without.reactions[0].workups[:]
    parquet.save_dataset(
        without, str(projected / "data" / "aa" / "ord_dataset-aa.parquet")
    )
    _reproject(projected)
    _run(projected, "--levels", "workups")
    written = projected / "pivots" / "workups" / "aa" / "ord_dataset-aa.parquet"
    assert pq.read_table(written).num_rows == 0

    parquet.save_dataset(
        _dataset("ord_dataset-aa"),
        str(projected / "data" / "aa" / "ord_dataset-aa.parquet"),
    )
    _reproject(projected)
    _run(projected, "--levels", "workups")
    assert pq.read_table(written).num_rows == 1


def test_an_unknown_level_is_refused_before_anything_is_written(projected):
    with pytest.raises(ValueError, match="not repeated levels"):
        _run(projected, "--levels", "workups", "conditions.temperature")
    assert not (projected / "pivots").exists()


def test_a_pattern_matching_only_non_projections_is_an_error(tmp_path):
    # The matches are readable Parquet that derive_tree ignores, so the run finishes
    # having written nothing -- which a pipeline step downstream would otherwise take
    # for artifacts that were built.
    directory = tmp_path / "projections" / "aa"
    directory.mkdir(parents=True)
    parquet.save_dataset(_dataset("ord_dataset-aa"), str(directory / "source.parquet"))
    with pytest.raises(ValueError, match="no pivots derived for workups"):
        _run(tmp_path, "--levels", "workups")


def test_a_corpus_reads_the_tree_this_script_writes(projected, caplog):
    # The layout is a contract between this script and Corpus(pivots_dir=...), and each
    # side tested against its own idea of it would agree with nobody. Derived here and
    # read there, end to end.
    derive_structures.main(
        derive_structures.parse_args(
            [
                f"--input_pattern={projected / 'projections' / '*' / '*.parquet'}",
                f"--output_dir={projected / 'structures'}",
            ]
        )
    )
    _run(projected, "--levels", "workups")
    where = {
        "op": "exists",
        "path": "workups",
        "where": {"op": "eq", "path": "type", "value": {"literal": "EXTRACTION"}},
    }
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(projected / "projections" / "*" / "*.parquet"),
            str(projected / "structures" / "*" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(projected / "pivots"),
        ) as corpus,
    ):
        found = corpus.search(query.Query.model_validate({"where": where}))
    assert set(found.column("reaction_id").to_pylist()) == {
        "ord-ord_dataset-aa-01",
        "ord-ord_dataset-bb-01",
    }
    messages = [record.message for record in caplog.records]
    assert any("read 2 pivot artifacts for workups" in message for message in messages)
    assert not any("building the pivot" in message for message in messages)
