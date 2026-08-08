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

"""Tests for ord_schema.scripts.derive_views.

View behavior is covered in ord_schema.views_test; these cover the CLI: path mapping,
which matches count as sources, the skip-if-current shortcut, and --force.
"""

import pathlib

import pyarrow.parquet as pq
import pytest

from ord_schema import artifacts, parquet, views
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.scripts import derive_projection, derive_views


def _dataset(dataset_id: str, *, reaction_ids=("ord-0001",)):
    reactions = []
    for reaction_id in reaction_ids:
        reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
        reaction.identifiers.add(type="REACTION_SMILES", value="c1ccccc1>>Cc1ccccc1")
        reactions.append(reaction)
    return dataset_pb2.Dataset(
        dataset_id=dataset_id, name="test", description="desc", reactions=reactions
    )


def _corpus(root: pathlib.Path, shards=("aa", "bb")) -> None:
    """Writes data/<shard>/ord_dataset-<shard>.parquet under root."""
    for shard in shards:
        directory = root / "data" / shard
        directory.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            _dataset(f"ord_dataset-{shard}"),
            str(directory / f"ord_dataset-{shard}.parquet"),
        )


def _project(root: pathlib.Path) -> None:
    """Projects the corpus, which is what views are derived from."""
    derive_projection.main(
        derive_projection.parse_args(
            [
                f"--input_pattern={root / 'data' / '*' / '*.parquet'}",
                f"--output_dir={root / 'projections'}",
            ]
        )
    )


def _run(root: pathlib.Path, *extra: str) -> None:
    derive_views.main(
        derive_views.parse_args(
            [
                f"--input_pattern={root / 'projections' / '*' / '*.parquet'}",
                f"--output_dir={root / 'views'}",
                *extra,
            ]
        )
    )


def test_main_writes_one_view_per_dataset(tmp_path):
    _corpus(tmp_path)
    _project(tmp_path)
    _run(tmp_path)
    for shard in ("aa", "bb"):
        view = tmp_path / "views" / shard / f"ord_dataset-{shard}.parquet"
        assert view.exists()
        table = pq.read_table(view)
        assert table.schema.names == views.SCHEMA.names
        assert table.column("reaction_id").to_pylist() == ["ord-0001"]


def test_main_stamps_views_with_their_source(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    source = tmp_path / "data" / "aa" / "ord_dataset-aa.parquet"
    stamps = artifacts.load_stamps(tmp_path / "views" / "aa" / "ord_dataset-aa.parquet")
    assert stamps.source_dataset_id == "ord_dataset-aa"
    assert stamps.source_md5 == parquet.DatasetView(source).md5()


# write_view publishes by renaming a fresh temp file into place, so the inode changes
# if and only if the view was rewritten. That is exact, unlike comparing mtimes.
def test_main_skips_views_that_are_current(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    view = tmp_path / "views" / "aa" / "ord_dataset-aa.parquet"
    before = view.stat().st_ino
    _run(tmp_path)
    assert view.stat().st_ino == before


def test_force_rewrites_a_current_view(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    view = tmp_path / "views" / "aa" / "ord_dataset-aa.parquet"
    before = view.stat().st_ino
    _run(tmp_path, "--force")
    assert view.stat().st_ino != before


def test_main_rederives_when_the_source_changes(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    source = tmp_path / "data" / "aa" / "ord_dataset-aa.parquet"
    parquet.save_dataset(
        _dataset("ord_dataset-aa", reaction_ids=("ord-0001", "ord-0002")), str(source)
    )
    _project(tmp_path)
    _run(tmp_path)
    view = tmp_path / "views" / "aa" / "ord_dataset-aa.parquet"
    assert pq.read_table(view).column("reaction_id").to_pylist() == [
        "ord-0001",
        "ord-0002",
    ]


def test_main_raises_when_every_match_is_already_a_view(tmp_path):
    # Aiming --input_pattern at the output tree would otherwise exit 0 having written
    # nothing, and a downstream step would proceed as though the views were built.
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    _run(tmp_path)
    with pytest.raises(ValueError, match="are projections"):
        derive_views.main(
            derive_views.parse_args(
                [
                    f"--input_pattern={tmp_path}/views/*/*.parquet",
                    f"--output_dir={tmp_path}/again",
                ]
            )
        )


def test_main_ignores_source_datasets(tmp_path):
    # A view is a narrowing of the projection, so a pattern aimed at the sources has
    # nothing to narrow -- rather than silently deriving one the slow way.
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    with pytest.raises(ValueError, match="are projections"):
        derive_views.main(
            derive_views.parse_args(
                [
                    f"--input_pattern={tmp_path}/data/*/*.parquet",
                    f"--output_dir={tmp_path}/views",
                ]
            )
        )


def test_main_raises_when_nothing_matches(tmp_path):
    with pytest.raises(ValueError, match="no datasets matched"):
        derive_views.main(
            derive_views.parse_args(
                [
                    f"--input_pattern={tmp_path / 'absent' / '*.parquet'}",
                    f"--output_dir={tmp_path / 'views'}",
                ]
            )
        )


def test_main_accepts_an_exact_filename(tmp_path):
    _corpus(tmp_path, shards=("aa",))
    _project(tmp_path)
    source = tmp_path / "projections" / "aa" / "ord_dataset-aa.parquet"
    derive_views.main(
        derive_views.parse_args(
            [f"--input_pattern={source}", f"--output_dir={tmp_path / 'views'}"]
        )
    )
    view = tmp_path / "views" / "ord_dataset-aa.parquet"
    assert view.is_file()
    assert pq.read_table(view).column("reaction_id").to_pylist() == ["ord-0001"]
