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
shortcut, and how an unknown level is refused.
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


def _dataset(dataset_id: str):
    reaction = reaction_pb2.Reaction(reaction_id=f"ord-{dataset_id}-01")
    reaction.inputs["a"].components.add().identifiers.add(
        type="SMILES", value="c1ccccc1"
    )
    reaction.workups.add(type="EXTRACTION")
    reaction.outcomes.add().products.add(isolated_color="white")
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
