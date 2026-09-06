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
another level's, that a level empty in every artifact says so and keeps saying so, that
what is counted is what a reader reads -- and nothing else in the tree -- and that an
artifact short a column is derived again rather than left current, from here and from a
Corpus.
"""

import logging
import pathlib
import shutil

import pyarrow as pa
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


def _messages(caplog) -> list[str]:
    """Returns every message recorded so far, formatted."""
    return [record.getMessage() for record in caplog.records]


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
    with caplog.at_level(logging.INFO):
        # The empty level named second, and derived after a level whose artifacts are
        # already on disk: named first it is also levels[0] and the only level in the
        # tree, so a report that ignored level_path, or that scanned the whole output
        # directory rather than this level's, would read the same.
        _run(projected, "--levels", "workups", "observations")
    # The premise above, asserted rather than assumed: workups really was derived
    # first, so its artifacts were on disk while observations was being measured.
    reports = [m for m in _messages(caplog) if "artifacts empty" in m]
    assert [m.split(":")[0] for m in reports] == ["workups", "observations"]
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
    messages = [record.getMessage() for record in caplog.records]
    # The tally an operator reads, with the two numbers distinct so their order shows.
    assert any(
        "observations: 2 written, 0 already current, 1 of 2 artifacts empty" in message
        for message in messages
    )
    # Named, not just counted: the projection that came out empty is what an operator
    # would go and look at, and a count alone sends them through the whole tree.
    named = [
        message for message in messages if "no elements at observations" in message
    ]
    assert len(named) == 1
    assert "bb" in named[0]
    assert "aa" not in named[0]


def test_files_that_are_not_this_level_s_pivots_are_not_counted(projected, caplog):
    # The tally decides whether a level is reported as empty, so anything it miscounts
    # either invents that report or suppresses it. All of these turn up in a real tree:
    # atomic_io writes its temp as a sibling inside the level's own directory, so a run
    # killed mid-write leaves one behind, and an output directory is a place people put
    # things.
    _run(projected, "--levels", "workups", "outcomes.products")
    beside = projected / "pivots" / "workups" / "aa"
    # A leaked temp is not named like an artifact, so no reader sees it and neither
    # does the tally -- which matters because this one holds rows, and counting it
    # would say the level has something no query can reach.
    (beside / "ord_dataset-aa.parquet.q7x1.tmp").write_bytes(b"")
    # These three a reader would pick up, and each is not this level's pivot.
    (beside / "corrupt.parquet").write_bytes(b"")
    pq.write_table(pa.table({"x": [1]}), beside / "unstamped.parquet")
    shutil.copy(
        projected / "pivots" / "outcomes.products" / "aa" / "ord_dataset-aa.parquet",
        beside / "another-level.parquet",
    )
    with caplog.at_level(logging.INFO):
        _run(projected, "--levels", "workups")
    messages = [record.getMessage() for record in caplog.records]
    assert any(
        "workups: 0 written, 2 already current, 0 of 2 artifacts empty" in message
        for message in messages
    )
    warnings = _warnings(caplog)
    assert sorted(warning.split("/")[-1] for warning in warnings) == [
        "another-level.parquet holds outcomes.products: not counted at workups",
        "corrupt.parquet cannot be read as Parquet: not counted at workups",
        "unstamped.parquet carries no pivot level: not counted at workups",
    ]


def test_a_shard_a_reader_cannot_descend_is_an_error(projected, tmp_path):
    # A shard directory symlinked onto other storage is written straight through and
    # then has to be listed back through the link. Half a level going missing while the
    # run reports success is the shape that has no other alarm: the denominator shrinks
    # to match, so the level reads as complete.
    _run(projected, "--levels", "workups")
    shard = projected / "pivots" / "workups" / "bb"
    elsewhere = tmp_path / "elsewhere"
    shard.rename(elsewhere)
    shard.symlink_to(elsewhere, target_is_directory=True)
    assert len(pivot.artifact_paths(projected / "pivots", "workups")) == 2
    # Nothing is rewritten, so a run that could not descend would find one of two.
    _run(projected, "--levels", "workups")


def test_a_level_whose_artifacts_a_reader_cannot_find_is_an_error(tmp_path):
    # A tree derived from projections a reader's glob does not match is written where
    # nothing looks: derive_tree reports every artifact current, and every quantifier
    # over the level falls back to unnesting the projection, forever.
    for shard in ("aa", "bb"):
        directory = tmp_path / "data" / shard
        directory.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            _dataset(f"ord_dataset-{shard}"),
            str(directory / f"ord_dataset-{shard}.parquet"),
        )
    _reproject(tmp_path)
    for projection_path in (tmp_path / "projections").rglob("*.parquet"):
        projection_path.replace(projection_path.with_suffix(".pq"))
    arguments = derive_pivots.parse_args(
        [
            f"--input_pattern={tmp_path / 'projections' / '*' / '*.pq'}",
            f"--output_dir={tmp_path / 'pivots'}",
            "--levels",
            "workups",
        ]
    )
    with pytest.raises(ValueError, match="holds 0 pivot artifacts for workups"):
        derive_pivots.main(arguments)
    # Again, where the artifacts are now current and the run writes nothing: that is
    # the steady state of every build after the first, and a guard that counted only
    # what this run wrote would find nothing missing and report the level empty.
    with pytest.raises(ValueError, match="holds 0 pivot artifacts for workups"):
        derive_pivots.main(arguments)


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
    # current -- which every later run skips, and which answers no exists and
    # every forall.
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


def _drop_column(path: pathlib.Path, column: str) -> None:
    """Rewrites the artifact at ``path`` without ``column``, keeping its footer.

    Stands in for an artifact written before the schema carried that column, by a
    library version that would not have bumped for it: the stamps are what a run reads
    to decide the file is current, and they say nothing about its columns.

    Args:
        path: The artifact to rewrite.
        column: The column to drop.
    """
    table = pq.read_table(path)
    metadata = table.schema.metadata
    table = table.drop_columns([column])
    pq.write_table(table.replace_schema_metadata(metadata), path)


def test_an_artifact_missing_an_ordinal_is_derived_again(projected):
    # A repeated level added above this one gives it another ordinal, and the artifacts
    # written before it carry stamps a later run cannot tell apart, absent a version
    # bump -- so the columns are the only difference, and without declaring them the
    # short artifact is current forever. An uncorrelated query still answers off it,
    # and only a correlation joining on the missing ordinal fails, so this asserts on
    # the file rather than on a result.
    _run(projected, "--levels", "outcomes.products")
    written = (
        projected / "pivots" / "outcomes.products" / "aa" / "ord_dataset-aa.parquet"
    )
    _drop_column(written, "product_index")
    assert "product_index" not in pq.read_schema(written).names
    _run(projected, "--levels", "outcomes.products")
    assert "product_index" in pq.read_schema(written).names


def test_a_corpus_deriving_pivots_also_derives_an_artifact_missing_an_ordinal(
    projected,
):
    # The script and the Corpus write into one tree, so they have to agree on which
    # artifacts are current; one of them accepting a column-short artifact the other
    # rebuilds would leave the answer to whichever ran last.
    derive_structures.main(
        derive_structures.parse_args(
            [
                f"--input_pattern={projected / 'projections' / '*' / '*.parquet'}",
                f"--output_dir={projected / 'structures'}",
            ]
        )
    )
    _run(projected, "--levels", "workups")
    written = projected / "pivots" / "workups" / "aa" / "ord_dataset-aa.parquet"
    _drop_column(written, "workup_index")
    with execute.Corpus(
        str(projected / "projections" / "*" / "*.parquet"),
        str(projected / "structures" / "*" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(projected / "pivots"),
        derive_pivots=True,
        warm=False,
    ) as corpus:
        corpus.search(
            query.Query.model_validate(
                {
                    "where": {
                        "op": "exists",
                        "path": "workups",
                        "where": {
                            "op": "eq",
                            "path": "type",
                            "value": {"literal": "EXTRACTION"},
                        },
                    }
                }
            )
        )
    assert "workup_index" in pq.read_schema(written).names


def test_a_level_whose_artifacts_cannot_be_found_is_an_error(projected, monkeypatch):
    # Without this, no artifacts found means none of them hold rows, and the run
    # announces that the level is empty everywhere -- the report this exists to make,
    # made about a scan that failed rather than about the corpus.
    def unreadable(path):
        raise pa.ArrowInvalid(f"{path} is not Parquet")

    monkeypatch.setattr(derive_pivots.pq, "read_metadata", unreadable)
    # Naming what it could not read is the difference between an actionable error and
    # "2 were written and a reader found none of them".
    with pytest.raises(
        ValueError,
        match=r"holds 0 pivot artifacts for workups.*"
        r"2 not counted.*ord_dataset-aa\.parquet cannot be read",
    ):
        _run(projected, "--levels", "workups")


def test_artifacts_a_shrunken_corpus_left_behind_are_not_an_error(projected, caplog):
    # A corpus that loses a dataset -- removed, or consolidated into another shard --
    # leaves its pivots behind with no projection to match. The guard against artifacts
    # a reader cannot find has to stay one-directional through that: surplus is not the
    # same failure as absence, and a build that stopped on it would stop on every run
    # from then on, sending an operator to look for files that are not missing.
    _run(projected, "--levels", "workups")
    orphan = projected / "pivots" / "workups" / "cc"
    orphan.mkdir(parents=True)
    shutil.copy(
        projected / "pivots" / "workups" / "aa" / "ord_dataset-aa.parquet",
        orphan / "ord_dataset-cc.parquet",
    )
    with caplog.at_level(logging.INFO):
        _run(projected, "--levels", "workups")
    assert any(
        "workups: 0 written, 2 already current, 0 of 3 artifacts empty" in message
        for message in _messages(caplog)
    )


def test_an_empty_level_is_not_also_named_artifact_by_artifact(projected, caplog):
    # The warning already says every artifact is empty. Naming each one as well would
    # be a line per shard across the levels a corpus records nothing at, which on ORD
    # is most of them -- burying the case the lines exist for, a level empty in only
    # some of its artifacts.
    with caplog.at_level(logging.INFO):
        _run(projected, "--levels", "observations")
    assert any("observations: every artifact is empty" in m for m in _warnings(caplog))
    assert not [m for m in _messages(caplog) if "no elements at" in m]


def test_each_level_is_derived_against_its_own_columns(projected):
    # The columns a level declares are its own, and a level's are a strict subset of
    # every level beneath it: "outcomes" has no product_index. A run deriving both
    # against one of their schemas would find a product artifact short an ordinal
    # current and leave it, which only a correlation joining on that ordinal notices,
    # by over-returning. Two levels, the shallower named first.
    _run(projected, "--levels", "outcomes", "outcomes.products")
    written = (
        projected / "pivots" / "outcomes.products" / "aa" / "ord_dataset-aa.parquet"
    )
    _drop_column(written, "product_index")
    _run(projected, "--levels", "outcomes", "outcomes.products")
    assert "product_index" in pq.read_schema(written).names


def test_a_level_this_run_did_not_name_is_checked_against_the_schema(projected):
    # Naming a subset is the documented cheaper run, and the tree it writes into keeps
    # whatever an earlier, larger run put there. Those artifacts stamp current forever,
    # so a projection that grew a column leaves every level nobody named short of it --
    # which is the whole tree a corpus reads, not just the part this run touched.
    _run(projected, "--levels", "workups", "outcomes.products")
    left_behind = projected / "pivots" / "workups" / "aa" / "ord_dataset-aa.parquet"
    _drop_column(left_behind, "workup_index")
    with pytest.raises(ValueError, match="did not derive"):
        _run(projected, "--levels", "outcomes.products")


def test_a_level_this_run_did_not_name_and_did_not_outgrow_is_not_an_error(projected):
    # The subset run has to stay usable. A level left alone whose artifacts still match
    # the schema is exactly what --levels is for, and stopping on it would make every
    # deployment derive all 39 levels to answer questions over four.
    _run(projected, "--levels", "workups", "outcomes.products")
    _run(projected, "--levels", "outcomes.products")


def test_a_level_the_tree_does_not_hold_is_not_left_behind(projected):
    # A level nobody has ever derived has no artifacts to be stale, and a run that
    # stopped on one would refuse every first build of a subset.
    _run(projected, "--levels", "outcomes.products")
