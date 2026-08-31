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

"""Tests for ord_schema.artifacts.occurrences."""

import pathlib
import shutil

import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import base, occurrences, pivot, projection, structures
from ord_schema.proto import dataset_pb2, reaction_pb2

_ROLE = reaction_pb2.ReactionRole.ReactionRoleType


def _component(reaction_input, smiles: str, role) -> None:
    component = reaction_input.components.add(reaction_role=role)
    component.identifiers.add(type="SMILES", value=smiles)


@pytest.fixture(scope="module")
def projected(tmp_path_factory) -> pathlib.Path:
    """A projection carrying a structure at every indexed path, and one at none.

    The second reaction has inputs and nothing else, so the deeper paths are what the
    projection writes for a level the source never recorded. The third component of the
    first reaction carries no identifier at all, which is the element a pivot holds a
    row for and this artifact must not.
    """
    reaction = reaction_pb2.Reaction(reaction_id="ord-oc01")
    _component(reaction.inputs["in"], "CCO", _ROLE.REACTANT)
    _component(reaction.inputs["in"], "c1ccncc1", _ROLE.SOLVENT)
    reaction.inputs["in"].components.add(reaction_role=_ROLE.CATALYST)
    _component(reaction.workups.add(type="EXTRACTION").input, "CCOCC", _ROLE.SOLVENT)
    product = reaction.outcomes.add().products.add()
    product.identifiers.add(type="SMILES", value="Cc1ccccc1")
    standard = product.measurements.add(analysis_key="nmr").authentic_standard
    standard.identifiers.add(type="SMILES", value="c1ccccc1")
    shallow = reaction_pb2.Reaction(reaction_id="ord-oc02")
    _component(shallow.inputs["in"], "CCO", _ROLE.REACTANT)
    root = tmp_path_factory.mktemp("occurrences")
    source = root / "ord_dataset-oc.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-oc",
            name="test",
            description="test",
            reactions=[reaction, shallow],
        ),
        str(source),
    )
    output = root / "projection.parquet"
    projection.write_projection(source, output)
    return output


@pytest.fixture(scope="module")
def pivots(projected, tmp_path_factory) -> pathlib.Path:
    """A pivot per level the indexed paths are read from."""
    root = tmp_path_factory.mktemp("pivots")
    for level in {level.path for level, _ in occurrences.PATHS.values()}:
        target = root / level / projected.name
        target.parent.mkdir(parents=True, exist_ok=True)
        pivot.write_pivot(projected, target, level_path=level)
    return root


def _rows(path: pathlib.Path) -> list[tuple]:
    table = pq.read_table(path)
    return sorted(
        zip(
            *(table.column(name).to_pylist() for name in table.schema.names),
            strict=True,
        )
    )


def _derive(pivots: pathlib.Path, into: pathlib.Path, path: str) -> pathlib.Path:
    level = occurrences.PATHS[path][0].path
    source = next((pivots / level).glob("*.parquet"))
    target = into / path / source.name
    target.parent.mkdir(parents=True, exist_ok=True)
    occurrences.write_occurrences(source, target, path=path)
    return target


# What the artifact covers


def test_the_covered_paths_are_the_structure_bearing_ones():
    # Walked rather than listed, so a compound field added upstream is covered wherever
    # it lands. A path dropped here is one whose structures no query can find, which
    # reads as a corpus holding none of them.
    walked = {path for _, path, _ in structures.structure_levels()}
    assert set(occurrences.PATHS) == walked


def test_a_path_is_read_from_the_level_it_ranges_within():
    # Three of the four name a repeated level of their own; the authentic standard is
    # one compound per measurement, so it is reached through the remainder.
    level, remainder = occurrences.PATHS["inputs.components"]
    assert (level.path, remainder) == ("inputs.components", ())
    level, remainder = occurrences.PATHS[
        "outcomes.products.measurements.authentic_standard"
    ]
    assert (level.path, remainder) == (
        "outcomes.products.measurements",
        ("authentic_standard",),
    )


# Derivation


def test_an_occurrence_carries_the_reaction_the_structure_and_the_role(
    pivots, tmp_path
):
    written = _derive(pivots, tmp_path, "inputs.components")
    assert _rows(written) == [
        ("ord-oc01", 0, "REACTANT"),
        ("ord-oc01", 1, "SOLVENT"),
        ("ord-oc02", 0, "REACTANT"),
    ]


def test_an_element_carrying_no_structure_is_no_occurrence(pivots, tmp_path):
    # The pivot is complete -- every element gets a row, including one whose fields are
    # all NULL -- and this is not: a row with no structure_id would join to no molecule
    # and count against the structures the index has to reach.
    level = occurrences.PATHS["inputs.components"][0].path
    source = next((pivots / level).glob("*.parquet"))
    assert pq.read_metadata(source).num_rows == 4
    written = _derive(pivots, tmp_path, "inputs.components")
    assert pq.read_metadata(written).num_rows == 3


def test_a_path_below_its_level_is_read_through_the_remainder(pivots, tmp_path):
    # The authentic standard sits inside the measurement's element rather than being an
    # element of its own, so this is the one path whose read is not the level itself.
    # Its role is NULL, which is what the source records: a standard is named to
    # identify a peak rather than to play a part in the reaction, so a role-bound query
    # over this path matches nothing, and that is the right answer.
    written = _derive(
        pivots, tmp_path, "outcomes.products.measurements.authentic_standard"
    )
    assert _rows(written) == [("ord-oc01", 4, None)]


def test_a_level_the_source_never_recorded_writes_an_empty_artifact(pivots, tmp_path):
    # Empty rather than missing: a reader globs a path and finds a file holding nothing,
    # which is what a corpus recording no workups looks like from here.
    for level in ("workups.input.components", "inputs.components"):
        (pivots / level).exists()
    written = _derive(pivots, tmp_path, "workups.input.components")
    assert _rows(written) == [("ord-oc01", 2, "SOLVENT")]


def test_the_artifact_carries_the_schema_it_promises(pivots, tmp_path):
    written = _derive(pivots, tmp_path, "inputs.components")
    assert pq.read_schema(written).remove_metadata() == occurrences.SCHEMA


def test_the_artifact_names_the_path_it_holds(pivots, tmp_path):
    written = _derive(pivots, tmp_path, "outcomes.products")
    assert occurrences.occurrence_path(written) == "outcomes.products"


def test_a_file_that_is_not_this_artifact_names_no_path(projected):
    assert occurrences.occurrence_path(projected) is None


def test_the_source_dataset_is_passed_through_rather_than_rehashed(pivots, tmp_path):
    # Every artifact derived from a dataset names that dataset, however many
    # derivations away it sits, so one comparison answers "is this current for it?"
    level = occurrences.PATHS["inputs.components"][0].path
    source = next((pivots / level).glob("*.parquet"))
    written = _derive(pivots, tmp_path, "inputs.components")
    assert base.load_stamps(written).source_md5 == base.load_stamps(source).source_md5
    assert base.load_stamps(written).artifact == occurrences.ARTIFACT


# What it refuses


def test_a_path_the_schema_carries_no_structure_at_is_refused(pivots, tmp_path):
    level = occurrences.PATHS["inputs.components"][0].path
    source = next((pivots / level).glob("*.parquet"))
    with pytest.raises(ValueError, match="not a path this artifact covers"):
        occurrences.write_occurrences(source, tmp_path / "x.parquet", path="workups")


def test_a_parent_that_is_not_a_pivot_is_refused(projected, tmp_path):
    # A projection holds the same elements nested, so deriving from one would produce
    # something -- against a schema walk this module does not do.
    with pytest.raises(ValueError, match="is a projection, not a pivot"):
        occurrences.write_occurrences(
            projected, tmp_path / "x.parquet", path="inputs.components"
        )


def test_a_pivot_over_another_level_is_refused(pivots, tmp_path):
    # The rows would be the wrong elements, and every column this reads exists on both.
    source = next((pivots / "outcomes.products").glob("*.parquet"))
    with pytest.raises(ValueError, match=r"holds the pivot over outcomes\.products"):
        occurrences.write_occurrences(
            source, tmp_path / "x.parquet", path="inputs.components"
        )


def test_a_stale_pivot_is_refused(pivots, tmp_path):
    # The stamps are the only thing saying an artifact was written by the library
    # reading it, and occurrences derived from a stale one would claim its provenance.
    level = occurrences.PATHS["inputs.components"][0].path
    source = next((pivots / level).glob("*.parquet"))
    staled = tmp_path / "stale.parquet"
    table = pq.read_table(source)
    metadata = dict(table.schema.metadata)
    metadata[b"ord.rdkit_version"] = b"0000.00.0"
    pq.write_table(table.replace_schema_metadata(metadata), staled)
    with pytest.raises(ValueError, match="is stale"):
        occurrences.write_occurrences(
            staled, tmp_path / "x.parquet", path="inputs.components"
        )


# Finding them


def test_artifacts_are_found_beneath_a_shard_directory(pivots, tmp_path):
    # The tree mirrors the projections it descends from, so a sharded corpus puts them
    # a directory deeper and a reader that does not descend finds nothing.
    written = _derive(pivots, tmp_path, "inputs.components")
    sharded = tmp_path / "sharded" / "inputs.components" / "aa"
    sharded.mkdir(parents=True)
    shutil.copy(written, sharded / written.name)
    found = occurrences.artifact_paths(tmp_path / "sharded", "inputs.components")
    assert [path.name for path in found] == [written.name]


def test_a_path_with_no_directory_has_no_artifacts(tmp_path):
    assert occurrences.artifact_paths(tmp_path, "inputs.components") == []
