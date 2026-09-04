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

"""Tests for ord_schema.artifacts.pivot."""

import os
import pathlib
import re
import shutil

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ord_schema import parquet
from ord_schema.artifacts import base, pivot, projection
from ord_schema.proto import dataset_pb2, reaction_pb2


def test_covers_the_levels_a_query_quantifies_over():
    levels = pivot.repeated_levels()
    for path in (
        "workups",
        "inputs",
        "inputs.components",
        "outcomes",
        "outcomes.products",
        "outcomes.products.measurements",
        "conditions.temperature.measurements",
    ):
        assert path in levels


def test_ordinals_accumulate_down_the_path():
    levels = pivot.repeated_levels()
    assert levels["workups"].ordinals == ("workup_index",)
    assert levels["outcomes"].ordinals == ("outcome_index",)
    assert levels["outcomes.products"].ordinals == ("outcome_index", "product_index")
    assert levels["outcomes.products.measurements"].ordinals == (
        "outcome_index",
        "product_index",
        "measurement_index",
    )
    # A map is a repeated level too, so its components sit one ordinal deeper.
    assert levels["inputs.components"].ordinals == ("input_index", "component_index")


def test_element_type_drops_repeated_fields_and_keeps_structs():
    element = pivot.repeated_levels()["outcomes.products"].element_type
    names = [field.name for field in element]
    assert "isolated_color" in names
    assert "texture" in names
    assert "measurements" not in names
    assert "identifiers" not in names
    assert "analyses" not in names
    assert [field.name for field in element.field("texture").type] == [
        "type",
        "details",
    ]


def test_a_repeated_field_inside_a_kept_struct_is_dropped_too():
    element = pivot.repeated_levels()["outcomes.products.measurements"].element_type
    standard = element.field("authentic_standard").type
    names = [field.name for field in standard]
    assert "smiles" in names
    assert "identifiers" not in names
    assert "preparations" not in names


def test_the_structure_id_is_kept():
    # A pivot can answer a structure predicate, which needs the corpus-wide ID the
    # schema marks internal for a model's benefit rather than this walk's.
    element = pivot.repeated_levels()["inputs.components"].element_type
    assert "structure_id" in [field.name for field in element]


def test_a_level_whose_elements_are_all_repeated_is_not_covered():
    schema = pa.schema(
        [
            pa.field(
                "things",
                pa.list_(pa.struct([pa.field("bits", pa.list_(pa.int32()))])),
            )
        ]
    )
    assert pivot.repeated_levels(schema) == {}


def test_a_list_of_scalars_is_not_a_level():
    schema = pa.schema([pa.field("names", pa.list_(pa.string()))])
    assert pivot.repeated_levels(schema) == {}


def test_enum_metadata_survives_pruning():
    # The compiler validates an enum comparison against the members on the field, so a
    # pruned type that dropped them would refuse a query the projection accepts.
    element = pivot.repeated_levels()["inputs.components"].element_type
    field = element.field("reaction_role")
    assert field.metadata is not None


def test_table_name_is_a_bare_identifier():
    assert pivot.table_name("outcomes.products") == "pivot_outcomes_products"
    assert pivot.table_name("workups") == "pivot_workups"


def test_ordinals_that_would_collide_are_refused():
    # Two segments singularizing to one name would give a row two columns of that name,
    # and the second would silently win.
    schema = pa.schema(
        [
            pa.field(
                "analyses",
                pa.list_(
                    pa.struct(
                        [
                            pa.field(
                                "analysis",
                                pa.list_(pa.struct([pa.field("value", pa.int32())])),
                            )
                        ]
                    )
                ),
            )
        ]
    )
    with pytest.raises(ValueError, match="analysis_index"):
        pivot.repeated_levels(schema)


@pytest.fixture(scope="module")
def projected(tmp_path_factory) -> pathlib.Path:
    """A written projection with two outcomes on one reaction and none on another."""
    root = tmp_path_factory.mktemp("pivot")
    split = reaction_pb2.Reaction(reaction_id="ord-pv01")
    first = split.outcomes.add()
    first.products.add(isolated_color="white", is_desired_product=False)
    second = split.outcomes.add()
    second.products.add(isolated_color="yellow", is_desired_product=True)
    second.products.add(isolated_color="red", is_desired_product=False)
    source = root / "ord_dataset-pv.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-pv",
            name="test",
            description="test",
            reactions=[split, reaction_pb2.Reaction(reaction_id="ord-pv02")],
        ),
        str(source),
    )
    output = root / "projection.parquet"
    projection.write_projection(source, output)
    return output


def _fill_compound(compound) -> None:
    """Populates every repeated level a compound carries, several elements per list.

    Args:
        compound: The compound to fill.
    """
    for name in ("ethanol", "ethyl alcohol"):
        compound.identifiers.add(type=reaction_pb2.CompoundIdentifier.NAME, value=name)
    compound.preparations.add(type=reaction_pb2.CompoundPreparation.DRIED)
    compound.preparations.add(type=reaction_pb2.CompoundPreparation.SYNTHESIZED)
    compound.features["first"].string_value = "recorded"
    compound.features["second"].string_value = "recorded"
    for key in ("first", "second", "third"):
        analysis = compound.analyses[key]
        analysis.type = reaction_pb2.Analysis.LC
        analysis.data["datum"].string_value = "recorded"
        analysis.data["another"].string_value = "recorded"


def _fill_input(reaction_input, components: int) -> None:
    """Populates a reaction input, its components, and their levels.

    Args:
        reaction_input: The input to fill.
        components: How many components to add, which is what makes a count of
            elements differ from a count of the lists holding them.
    """
    for _ in range(components):
        _fill_compound(reaction_input.components.add())
    reaction_input.crude_components.add(reaction_id="ord-pvfull", includes_workup=True)


def _second_reaction() -> reaction_pb2.Reaction:
    """A second reaction recording something at every repeated level.

    What separates a total from a per-row maximum, or from the first row's count, is a
    second row recording something at the same level -- not two rows recording
    different numbers, since (1, 1) already sums to more than either. Every level this
    fills is one where taking a maximum falls short.

    Returns:
        The reaction.
    """
    reaction = reaction_pb2.Reaction(reaction_id="ord-pvsecond")
    reaction.identifiers.add(
        type=reaction_pb2.ReactionIdentifier.REACTION_TYPE, value="reduction"
    )
    _fill_input(reaction.inputs["only"], components=1)
    reaction.setup.vessel.attachments.add(type=reaction_pb2.VesselAttachment.CAP)
    reaction.setup.vessel.preparations.add(type=reaction_pb2.VesselPreparation.PURGED)
    reaction.setup.automation_code["other"].string_value = "recorded"
    reaction.conditions.temperature.measurements.add(
        type=reaction_pb2.TemperatureConditions.TemperatureMeasurement.THERMOCOUPLE_EXTERNAL
    )
    reaction.conditions.pressure.measurements.add(
        type=reaction_pb2.PressureConditions.PressureMeasurement.CUSTOM,
        details="recorded",
    )
    reaction.conditions.electrochemistry.measurements.add(
        current={"value": 2.0, "units": "AMPERE"}
    )
    reaction.observations.add(comment="also recorded")
    workup = reaction.workups.add(type=reaction_pb2.ReactionWorkup.FILTRATION)
    _fill_input(workup.input, components=1)
    workup.temperature.measurements.add(
        type=reaction_pb2.TemperatureConditions.TemperatureMeasurement.THERMOCOUPLE_EXTERNAL
    )
    outcome = reaction.outcomes.add()
    outcome_analysis = outcome.analyses["other"]
    outcome_analysis.type = reaction_pb2.Analysis.GC
    outcome_analysis.data["datum"].string_value = "recorded"
    product = outcome.products.add(is_desired_product=True)
    product.identifiers.add(type=reaction_pb2.CompoundIdentifier.NAME, value="ethanol")
    product.features["only"].string_value = "recorded"
    measurement = product.measurements.add(type=reaction_pb2.ProductMeasurement.YIELD)
    _fill_compound(measurement.authentic_standard)
    reaction.provenance.record_modified.add(details="also recorded")
    reaction.provenance.reaction_metadata["other"].string_value = "recorded"
    return reaction


@pytest.fixture(scope="module")
def populated(tmp_path_factory) -> pathlib.Path:
    """A projection recording elements at every repeated level, over three reactions.

    The count decides whether a level is unnested at all, so a wrong count is only
    catchable where the level holds something: against a projection recording nothing
    below the second level, a count that answered zero for everything deeper would
    agree with the unnest on every level there is.

    Three shapes, because three different wrong counts each agree with the unnest on a
    projection missing one of them. Lists holding more than one element, at enough
    levels that counting the lists rather than summing their lengths is caught -- not
    at every level, since a level whose lists happen to hold one element apiece ties
    either way. Two reactions recording something at every level, so a maximum or a
    first-row count falls short of the total everywhere. And a third recording nothing,
    so every level is also counted over a chain that is NULL rather than empty.
    """
    root = tmp_path_factory.mktemp("pivot_full")
    reaction = reaction_pb2.Reaction(reaction_id="ord-pvfull")
    reaction.identifiers.add(
        type=reaction_pb2.ReactionIdentifier.REACTION_TYPE, value="oxidation"
    )
    _fill_input(reaction.inputs["first"], components=3)
    _fill_input(reaction.inputs["second"], components=2)
    reaction.setup.vessel.attachments.add(type=reaction_pb2.VesselAttachment.SEPTUM)
    reaction.setup.vessel.preparations.add(
        type=reaction_pb2.VesselPreparation.OVEN_DRIED
    )
    reaction.setup.automation_code["code"].string_value = "recorded"
    reaction.conditions.temperature.measurements.add(
        type=reaction_pb2.TemperatureConditions.TemperatureMeasurement.THERMOCOUPLE_INTERNAL
    )
    reaction.conditions.pressure.measurements.add(
        type=reaction_pb2.PressureConditions.PressureMeasurement.PRESSURE_TRANSDUCER
    )
    reaction.conditions.electrochemistry.measurements.add(
        current={"value": 1.0, "units": "AMPERE"}
    )
    reaction.observations.add(comment="recorded")
    reaction.workups.add(type=reaction_pb2.ReactionWorkup.FILTRATION)
    workup = reaction.workups.add(type=reaction_pb2.ReactionWorkup.EXTRACTION)
    _fill_input(workup.input, components=2)
    workup.temperature.measurements.add(
        type=reaction_pb2.TemperatureConditions.TemperatureMeasurement.THERMOCOUPLE_INTERNAL
    )
    reaction.outcomes.add()
    outcome = reaction.outcomes.add()
    outcome_analysis = outcome.analyses["analysis"]
    outcome_analysis.type = reaction_pb2.Analysis.LC
    outcome_analysis.data["datum"].string_value = "recorded"
    outcome.products.add(is_desired_product=False)
    product = outcome.products.add(is_desired_product=True)
    for name in ("acetaldehyde", "ethanal"):
        product.identifiers.add(type=reaction_pb2.CompoundIdentifier.NAME, value=name)
    product.features["first"].string_value = "recorded"
    product.features["second"].string_value = "recorded"
    product.measurements.add(type=reaction_pb2.ProductMeasurement.YIELD)
    measurement = product.measurements.add(
        type=reaction_pb2.ProductMeasurement.YIELD, analysis_key="first"
    )
    _fill_compound(measurement.authentic_standard)
    reaction.provenance.record_modified.add(details="recorded")
    reaction.provenance.reaction_metadata["note"].string_value = "recorded"
    source = root / "ord_dataset-pvfull.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-pvfull",
            name="test",
            description="test",
            reactions=[
                reaction,
                _second_reaction(),
                reaction_pb2.Reaction(reaction_id="ord-pvempty"),
            ],
        ),
        str(source),
    )
    output = root / "projection.parquet"
    projection.write_projection(source, output)
    return output


def test_a_written_pivot_carries_its_rows_and_its_ordinals(projected, tmp_path):
    output = tmp_path / "products.parquet"
    written = pivot.write_pivot(projected, output, level_path="outcomes.products")
    assert written == 3
    table = pq.read_table(output)
    assert table.column_names == [
        "reaction_id",
        "outcome_index",
        "product_index",
        "element",
    ]
    rows = [
        (
            row["reaction_id"],
            row["outcome_index"],
            row["product_index"],
            row["element"]["isolated_color"],
        )
        for row in table.to_pylist()
    ]
    assert sorted(rows) == [
        ("ord-pv01", 1, 1, "white"),
        ("ord-pv01", 2, 1, "yellow"),
        ("ord-pv01", 2, 2, "red"),
    ]


def test_a_written_pivot_is_stamped_with_its_level(projected, tmp_path):
    output = tmp_path / "products.parquet"
    pivot.write_pivot(projected, output, level_path="outcomes.products")
    stamps = base.load_stamps(output)
    assert stamps.artifact == pivot.ARTIFACT
    assert stamps.source_dataset_id == "ord_dataset-pv"
    assert stamps.source_md5 == base.load_stamps(projected).source_md5
    assert base.stamps_are_current(stamps, pivot.ARTIFACT)
    assert pivot.pivot_path(output) == "outcomes.products"


def test_a_level_with_no_elements_still_writes_a_readable_file(projected, tmp_path):
    # A reader globs the level it wants; a file that is missing and a file that is
    # empty are different things, and only the second says the build ran.
    output = tmp_path / "workups.parquet"
    assert pivot.write_pivot(projected, output, level_path="workups") == 0
    assert pq.read_table(output).num_rows == 0
    assert pivot.pivot_path(output) == "workups"


def test_a_level_with_no_elements_is_never_unnested(projected, tmp_path, monkeypatch):
    # The saving: discovering that a level holds nothing by unnesting it costs a full
    # pass over every ancestor, and this projection records nothing at most levels.
    # Patched by name, which is the only way to assert that work did *not* happen --
    # and which stops holding the moment the query is built any way but through one of
    # these module globals. Both are patched: write_pivot builds through stage, and
    # select is the unstaged form it stands in for.
    def refuse(*args: object, **kwargs: object) -> str:
        raise AssertionError("the level was unnested")

    monkeypatch.setattr(pivot, "select", refuse)
    monkeypatch.setattr(pivot, "stage", refuse)
    output = tmp_path / "workups.parquet"
    assert pivot.write_pivot(projected, output, level_path="workups") == 0
    assert pq.read_table(output).num_rows == 0


def _unnested(source: pathlib.Path, level_path: str) -> int:
    """Returns the rows unnesting ``level_path`` out of ``source`` produces."""
    connection = duckdb.connect()
    try:
        connection.read_parquet(str(source)).create_view("reactions")
        return (
            connection.execute(pivot.select(pivot.LEVELS[level_path], "reactions"))
            .to_arrow_table()
            .num_rows
        )
    finally:
        connection.close()


@pytest.mark.parametrize("level_path", sorted(pivot.LEVELS))
def test_staging_a_level_produces_what_unnesting_it_in_one_query_does(
    populated, level_path
):
    # What makes staging safe is that it answers identically, and the part most easily
    # got wrong is invisible in a row count: each temp table carries the ordinals
    # accumulated above it, and a level whose prefix is dropped or misordered produces
    # exactly as many rows, with the wrong ones. Only a correlation joining on that
    # prefix would notice, in production rather than here, so compare the rows.
    #
    # Sorted reprs, because neither form promises an order and the artifact does not
    # need one -- and a row is a dict, which sorted() cannot order directly.
    level = pivot.LEVELS[level_path]
    connection = duckdb.connect()
    try:
        connection.read_parquet(str(populated)).create_view("reactions")
        direct = connection.execute(pivot.select(level, "reactions")).to_arrow_table()
        staged = connection.execute(
            pivot.stage(connection, level, "reactions")
        ).to_arrow_table()
    finally:
        connection.close()
    assert staged.schema == direct.schema
    assert staged.num_rows > 0, f"{level_path} is unpopulated; this proves nothing"
    assert sorted(map(repr, staged.to_pylist())) == sorted(
        map(repr, direct.to_pylist())
    )


def _staged_tables(connection) -> list[str]:
    """Returns the names of the temp tables staging left in ``connection``."""
    rows = connection.execute(
        "SELECT table_name FROM duckdb_tables() ORDER BY table_name"
    ).fetchall()
    return [name for (name,) in rows if name.startswith("stage_")]


def test_staging_leaves_a_table_behind_for_every_step_but_the_last(populated):
    # A table, specifically: a view would read as staged and plan as though it were
    # not, since the engine inlines it and the pipeline collapses back into one. The
    # saving would be gone with nothing in the SQL to show it.
    connection = duckdb.connect()
    try:
        connection.read_parquet(str(populated)).create_view("reactions")
        level = pivot.LEVELS["inputs.components.analyses"]
        connection.execute(pivot.stage(connection, level, "reactions")).fetchall()
        assert _staged_tables(connection) == ["stage_inputs", "stage_inputs_components"]
    finally:
        connection.close()


def test_a_level_with_one_step_is_staged_rather_than_fused(populated):
    # Staging is what lets an unnest run on its own, and a fused step runs with
    # whatever follows it. A single-step level fused is therefore the whole saving
    # given back, in the case that looks like it has nothing to stage.
    connection = duckdb.connect()
    try:
        connection.read_parquet(str(populated)).create_view("reactions")
        one_step = pivot.LEVELS["workups"]
        assert len(one_step.steps) == 1
        connection.execute(pivot.stage(connection, one_step, "reactions")).fetchall()
        assert _staged_tables(connection) == ["stage_workups"]
    finally:
        connection.close()


def test_staging_a_second_level_on_one_connection_is_refused(populated):
    # The tables are named for the level's ancestors, so two levels of one chain ask
    # for one name. Replacing a table a reader is still streaming from ends its rows
    # early and raises nothing, so the second call has to fail instead.
    connection = duckdb.connect()
    try:
        connection.read_parquet(str(populated)).create_view("reactions")
        pivot.stage(connection, pivot.LEVELS["inputs"], "reactions")
        with pytest.raises(duckdb.Error, match="stage_inputs"):
            pivot.stage(connection, pivot.LEVELS["inputs.components"], "reactions")
    finally:
        connection.close()


def test_write_pivot_stages_rather_than_unnesting_in_one_query(
    populated, tmp_path, monkeypatch
):
    # The inverse of the zero-element test: there select must not be called at all,
    # here it must not be what write_pivot builds with. Nothing else can tell the two
    # apart -- they answer identically, which is what the equivalence test proves -- so
    # every other assertion holds whichever one runs, and only this fails when a merge
    # or an edit quietly puts the level back in a single query.
    def refuse(*args: object, **kwargs: object) -> str:
        raise AssertionError("write_pivot built its query with select")

    monkeypatch.setattr(pivot, "select", refuse)
    output = tmp_path / "components.parquet"
    written = pivot.write_pivot(populated, output, level_path="inputs.components")
    assert written > 0
    assert pq.read_table(output).num_rows == written


def test_write_pivot_does_not_hold_the_order_of_the_unnest(
    populated, tmp_path, monkeypatch
):
    # Asserted on the connection rather than on a duration, which here would be a
    # flake; the setting is the mechanism, and what it buys is in pivot.write_pivot.
    #
    # Read as the connection closes, since write_pivot closes it before returning and
    # the setting is per-instance: any other connection answers for another database.
    seen = []

    class Watched:
        """Delegates to a connection, reading the setting as it closes."""

        def __init__(self, connection) -> None:
            self._connection = connection

        def __getattr__(self, name: str):
            return getattr(self._connection, name)

        def close(self) -> None:
            seen.append(
                self._connection.execute(
                    "SELECT current_setting('preserve_insertion_order')"
                ).fetchone()[0]
            )
            self._connection.close()

    original = duckdb.connect
    monkeypatch.setattr(
        pivot.duckdb, "connect", lambda *a, **k: Watched(original(*a, **k))
    )
    pivot.write_pivot(populated, tmp_path / "workups.parquet", level_path="workups")
    assert seen == [False]


@pytest.mark.parametrize("level_path", sorted(pivot.LEVELS))
def test_counting_a_level_agrees_with_unnesting_it(populated, level_path):
    # Compared against the unnest itself, never against write_pivot: its row count is
    # zero *because* the count said zero, so those two agree however wrong the count
    # is. A count wrongly zero publishes an empty artifact over a level that holds
    # rows, and every forall over it is then vacuously true of the whole corpus.
    counted = pivot.count_levels(populated, [level_path])[level_path]
    # LEVELS is walked from the schema, so a repeated field added upstream joins this
    # parametrize on its own -- and would be compared against an unnest of the nothing
    # the fixture records for it, which proves nothing. Whoever adds the field is the
    # one who has to populate it.
    assert counted > 0, f"{level_path} is unpopulated; the comparison proves nothing"
    assert counted == _unnested(populated, level_path)


def test_counting_many_levels_answers_each_one_for_itself(populated):
    # The aggregates come back as one row and are paired with the levels by position,
    # so a query built in one order and read in another answers every level with some
    # other level's count. Levels whose counts are pairwise distinct, or a permutation
    # would be the identity and there would be nothing here to fail.
    wanted = [
        "workups",
        "inputs.components",
        "inputs.components.analyses",
        "inputs.components.analyses.data",
    ]
    counted = pivot.count_levels(populated, wanted)
    assert len(counted) == len(wanted)
    assert len(set(counted.values())) == len(wanted), counted
    assert dict(counted) == {path: _unnested(populated, path) for path in wanted}


def test_each_part_of_a_file_identity_catches_a_replacement_the_others_miss(tmp_path):
    # Three components because each covers what the others let through, and the ways a
    # file is replaced without appearing to change are ordinary rather than exotic:
    # restoring timestamps is what rsync --times, tar -p, and a snapshot rollback all
    # do, and a filesystem with coarse mtime granularity does it by accident.
    path = tmp_path / "file"
    path.write_bytes(b"first")
    identity = pivot.file_identity(path)

    # Replaced by rename, with the timestamp put back: same size, same mtime, and only
    # the inode says a different file is there.
    replacement = tmp_path / "replacement"
    replacement.write_bytes(b"secnd")
    os.utime(replacement, ns=(identity.modified, identity.modified))
    replacement.replace(path)
    renamed = pivot.file_identity(path)
    assert (renamed.size, renamed.modified) == (identity.size, identity.modified)
    assert renamed != identity

    # Rewritten in place with the timestamp put back: the inode is unchanged and the
    # mtime restored, so the size is the only difference left.
    path.write_bytes(b"longer than before")
    os.utime(path, ns=(renamed.modified, renamed.modified))
    resized = pivot.file_identity(path)
    assert (resized.inode, resized.modified) == (renamed.inode, renamed.modified)
    assert resized != renamed

    # Rewritten in place at the same length: inode and size both hold, and the
    # timestamp is what is left to notice.
    path.write_bytes(b"longer than beforE")
    retimed = pivot.file_identity(path)
    assert (retimed.inode, retimed.size) == (resized.inode, resized.size)
    assert retimed != resized


def test_two_files_alike_but_for_their_filesystem_are_told_apart(tmp_path, monkeypatch):
    # Inode numbers are unique per filesystem, not across them, so without the device
    # two projections on different mounts can agree on all three of the rest -- and a
    # memo keyed on the identity would answer one with the other's counts. Staged
    # through stat, since a second filesystem is not a thing a test can mount.
    path = tmp_path / "file"
    path.write_bytes(b"first")
    here = pivot.file_identity(path)

    class _Elsewhere:
        """The same file as stat would report it from another filesystem."""

        st_dev = here.device + 1
        st_ino = here.inode
        st_size = here.size
        st_mtime_ns = here.modified

    monkeypatch.setattr(pathlib.Path, "stat", lambda self, **kwargs: _Elsewhere())
    there = pivot.file_identity(path)
    assert (there.inode, there.size, there.modified) == (
        here.inode,
        here.size,
        here.modified,
    )
    assert there != here


def test_artifact_paths_are_sorted(tmp_path):
    # The order artifacts are read in is the order a view holds their rows, and the
    # order a report names them. A glob answers in directory order, which is whatever
    # the filesystem last did.
    directory = tmp_path / "workups"
    for shard in ("cc", "aa", "bb"):
        (directory / shard).mkdir(parents=True)
        (directory / shard / f"ord_dataset-{shard}.parquet").write_bytes(b"")
    found = pivot.artifact_paths(tmp_path, "workups")
    assert [path.parent.name for path in found] == ["aa", "bb", "cc"]


def test_artifact_paths_reads_a_directory_whose_name_looks_like_a_pattern(tmp_path):
    # Only the last two segments are a pattern; everything above them is a name someone
    # chose. A bracket in it would be read as a character class and match nothing, and
    # the build would then reject a level whose artifacts are all present.
    directory = tmp_path / "pivots[v2]"
    (directory / "workups" / "aa").mkdir(parents=True)
    written = directory / "workups" / "aa" / "ord_dataset-aa.parquet"
    written.write_bytes(b"")
    assert pivot.artifact_paths(directory, "workups") == [written]


def test_artifact_paths_finds_nothing_where_there_is_no_level(tmp_path):
    assert pivot.artifact_paths(tmp_path, "workups") == []


def _miscounted(counts: pivot.Counts, level_path: str, value: int) -> pivot.Counts:
    """Returns ``counts`` with ``level_path`` answered wrongly, everything else kept."""
    return pivot.Counts(counts.identity, dict(counts) | {level_path: value})


@pytest.mark.parametrize("wrong", [99, 1])
def test_a_count_that_disagrees_with_the_unnest_is_refused(populated, tmp_path, wrong):
    # Wrong high and wrong low, since a check comparing one way would let the other
    # publish a short artifact. A count of zero skips the unnest and has nothing to
    # disagree with; what holds that case is that the count came from Counts, which
    # names the file it was read from and the level it answers for --
    # test_a_count_of_the_bytes_that_are_gone_is_refused and
    # test_a_count_of_another_level_is_not_an_answer_for_this_one.
    output = tmp_path / "products.parquet"
    counts = pivot.count_levels(populated, ["outcomes.products"])
    assert counts["outcomes.products"] not in (1, 99), "the miscount must be a miscount"
    with pytest.raises(ValueError, match="disagree"):
        pivot.write_pivot(
            populated,
            output,
            level_path="outcomes.products",
            counts=_miscounted(counts, "outcomes.products", wrong),
        )
    # Nothing published: an artifact that reached the destination would be stamped
    # current, skipped by every later run, and read as the truth about the level.
    assert not output.exists()


def test_a_count_of_the_bytes_that_are_gone_is_refused(populated, tmp_path):
    # The case the agreement check cannot reach, at the same path rather than a
    # different one -- an identity that merely told two paths apart would satisfy a
    # test using two files, and this is the accident the check exists for. Zero is the
    # count that skips the unnest, so nothing downstream contradicts it and the empty
    # artifact is published stamped current.
    source = tmp_path / "projection.parquet"
    shutil.copy(populated, source)
    counts = pivot.count_levels(source, ["outcomes.products"])
    assert counts["outcomes.products"] > 0
    stale = _miscounted(counts, "outcomes.products", 0)
    # Republished the way atomic_io does it, a sibling renamed over the destination,
    # with the timestamp put back the way restoring a backup would: a new inode is the
    # only thing separating the two, which is the component a same-path rewrite would
    # otherwise never exercise.
    replacement = tmp_path / "replacement.parquet"
    shutil.copy(populated, replacement)
    os.utime(replacement, ns=(counts.identity.modified, counts.identity.modified))
    replacement.replace(source)
    output = tmp_path / "products.parquet"
    with pytest.raises(ValueError, match="has changed since it was counted"):
        pivot.write_pivot(source, output, level_path="outcomes.products", counts=stale)
    assert not output.exists()


def test_a_count_of_another_level_is_not_an_answer_for_this_one(populated, tmp_path):
    # The count is taken from Counts by the level being written rather than passed
    # beside it, so an honest count of one level cannot be spent on another. Zero is
    # the one that would go unnoticed: the unnest is skipped and the empty artifact
    # published.
    counts = pivot.count_levels(populated, ["observations"])
    output = tmp_path / "components.parquet"
    with pytest.raises(KeyError, match=re.escape("inputs.components")):
        pivot.write_pivot(
            populated, output, level_path="inputs.components", counts=counts
        )
    assert not output.exists()


def test_a_projection_replaced_while_it_is_counted_is_refused(
    populated, tmp_path, monkeypatch
):
    # The counts would describe neither version. Caught by the identity being taken on
    # both sides of the read rather than only before it. Republished as atomic_io does
    # it, timestamp restored, so a new inode is the only difference and no one
    # component of the identity can carry this on its own.
    source = tmp_path / "projection.parquet"
    shutil.copy(populated, source)
    before = pivot.file_identity(source)
    original = pivot._count_view

    def replace_then_count(*args, **kwargs):
        result = original(*args, **kwargs)
        replacement = tmp_path / "replacement.parquet"
        shutil.copy(populated, replacement)
        os.utime(replacement, ns=(before.modified, before.modified))
        replacement.replace(source)
        return result

    monkeypatch.setattr(pivot, "_count_view", replace_then_count)
    with pytest.raises(ValueError, match="replaced while it was being counted"):
        pivot.count_levels(source, ["outcomes.products"])


def test_a_pivot_larger_than_one_batch_counts_every_batch(populated, tmp_path):
    # The rows are accumulated across batches and then checked against the count, so a
    # writer that streams more than one batch is the only shape where the accumulation
    # can be wrong. One write_table per batch means one row group per batch, so the
    # file says how many batches there were -- where the row count would only say that
    # row_group_size was honored, which is DuckDB's business and not this test's.
    output = tmp_path / "data.parquet"
    level_path = "inputs.components.analyses.data"
    counts = pivot.count_levels(populated, [level_path])
    written = pivot.write_pivot(
        populated, output, level_path=level_path, row_group_size=7, counts=counts
    )
    assert pq.ParquetFile(output).num_row_groups > 1, "one batch proves nothing here"
    assert pq.read_table(output).num_rows == written == counts[level_path]


def test_a_rejected_pivot_leaves_an_existing_artifact_alone(populated, tmp_path):
    # The rejected write is over a *different* level than the one already there. Aimed
    # at the same level it would unnest the same rows, so a destination overwritten by
    # the artifact being rejected would hold a file equal to the one it replaced, and
    # comparing them would pass whether or not anything was published.
    output = tmp_path / "pivot.parquet"
    pivot.write_pivot(populated, output, level_path="workups")
    good = pq.read_table(output)
    with pytest.raises(ValueError, match="disagree"):
        pivot.write_pivot(
            populated,
            output,
            level_path="outcomes.products",
            counts=_miscounted(
                pivot.count_levels(populated, ["outcomes.products"]),
                "outcomes.products",
                99,
            ),
        )
    assert pq.read_table(output).equals(good)
    assert pivot.pivot_path(output) == "workups"


def test_counting_a_level_twice_is_refused(populated):
    # zip would pair both, and the dict would collapse them, leaving the caller one
    # fewer answer than it asked for and no sign of it.
    with pytest.raises(ValueError, match="named more than once"):
        pivot.count_levels(populated, ["workups", "workups"])


def test_counting_a_stale_projection_is_refused(populated, tmp_path, monkeypatch):
    # An artifact derived from a stale projection inherits the dataset hash and so
    # claims a provenance it does not have; nothing would mark it stale again.
    monkeypatch.setattr(
        base, "ARTIFACT_VERSIONS", base.ARTIFACT_VERSIONS | {"projection": "later"}
    )
    with pytest.raises(ValueError, match="stale projection"):
        pivot.count_levels(populated, ["workups"])
    with pytest.raises(ValueError, match="stale projection"):
        pivot.write_pivot(populated, tmp_path / "workups.parquet", level_path="workups")


def test_counting_another_kind_of_artifact_is_refused(populated, tmp_path):
    # Distinct from a file carrying no stamps at all, which load_stamps refuses one
    # line earlier; this is the branch that reads the stamps and finds the wrong kind.
    other = tmp_path / "a-pivot.parquet"
    pivot.write_pivot(populated, other, level_path="workups")
    with pytest.raises(ValueError, match="is a pivot, not a projection"):
        pivot.count_levels(other, ["workups"])


def test_counting_something_that_is_not_a_projection_is_refused(tmp_path):
    # Otherwise the failure is a DuckDB binder error naming a view the caller never
    # created, for a file it should have been told is the wrong kind.
    other = tmp_path / "not-a-projection.parquet"
    pq.write_table(pa.table({"x": [1]}), other)
    with pytest.raises(ValueError, match="not a derived artifact"):
        pivot.count_levels(other, ["workups"])


def test_a_path_that_is_not_a_level_is_refused(projected, tmp_path):
    with pytest.raises(ValueError, match="not a repeated level"):
        pivot.write_pivot(
            projected, tmp_path / "x.parquet", level_path="conditions.temperature"
        )


def test_pivoting_something_that_is_not_a_projection_is_refused(projected, tmp_path):
    output = tmp_path / "products.parquet"
    pivot.write_pivot(projected, output, level_path="outcomes.products")
    with pytest.raises(ValueError, match="not a projection"):
        pivot.write_pivot(
            output, tmp_path / "again.parquet", level_path="outcomes.products"
        )


def test_reach_finds_a_level_that_is_the_path_itself():
    reached = pivot.reach("outcomes.products")
    assert reached is not None
    level, remainder, dtype = reached
    assert level.path == "outcomes.products"
    assert remainder == ()
    assert dtype is level.element_type


def test_reach_descends_to_a_singular_struct_under_a_level():
    # One authentic standard per measurement rather than a list of its own, so the
    # level it ranges over is the measurements, and the pivot over those carries it.
    reached = pivot.reach("outcomes.products.measurements.authentic_standard")
    assert reached is not None
    level, remainder, dtype = reached
    assert level.path == "outcomes.products.measurements"
    assert remainder == ("authentic_standard",)
    assert "smiles" in [field.name for field in dtype]


def test_reach_declines_a_path_leaving_the_pruned_element():
    # measurements is a repeated field of a product, so it is not on the pruned type;
    # it is a level of its own, and reach finds that level rather than descending.
    reached = pivot.reach("outcomes.products.measurements")
    assert reached is not None
    level, remainder, _ = reached
    assert level.path == "outcomes.products.measurements"
    assert remainder == ()
    # A repeated field whose elements are scalars is no level of its own -- there is no
    # element struct to pivot -- and pruning removed it from the measurement, so
    # descending to it finds nothing either way.
    assert (
        pivot.reach("outcomes.products.measurements.mass_spec_details.eic_masses")
        is None
    )
    # Whereas a repeated field of structs is a level, reached as itself rather than by
    # descending into the measurement that holds it.
    nested = pivot.reach(
        "outcomes.products.measurements.authentic_standard.identifiers"
    )
    assert nested is not None
    standard, remainder, _ = nested
    assert standard.path == (
        "outcomes.products.measurements.authentic_standard.identifiers"
    )
    assert remainder == ()


def test_reach_declines_a_path_with_no_repeated_ancestor():
    assert pivot.reach("conditions.temperature") is None
    assert pivot.reach("reaction_id") is None


def _write_through(mapping, key, value) -> None:
    """Assigns into ``mapping``, for a mapping that is not supposed to allow it.

    Args:
        mapping: The mapping to write into.
        key: The key to assign.
        value: The value to assign.
    """
    mapping[key] = value


def test_counts_of_the_same_levels_hash_alike_whatever_order_they_were_asked_in(
    populated,
):
    # A frozen dataclass advertises hashability, and the generated hash would reach the
    # mapping and raise. Equality does not depend on the order the levels were asked
    # in, so a hash that did would put equal values in different buckets.
    wanted = ["workups", "outcomes.products"]
    counts = pivot.count_levels(populated, wanted)
    reversed_counts = pivot.count_levels(populated, list(reversed(wanted)))
    assert counts == reversed_counts
    assert hash(counts) == hash(reversed_counts)
    assert len({counts, reversed_counts}) == 1


def test_counts_own_their_mapping(populated):
    # One build shares a projection's counts across every level it derives, so a write
    # that landed would answer for the rest of the run -- and move the hash under
    # anything holding one. Refused through the handle, and not reachable through the
    # mapping the Counts was built from, which a proxy alone would leave open.
    counts = pivot.count_levels(populated, ["workups"])
    with pytest.raises(TypeError):
        # Through an unannotated helper: a subscript here is refused statically too,
        # and what this asserts is that it is refused at runtime.
        _write_through(counts.at, "workups", 99)
    assert counts["workups"] == _unnested(populated, "workups")

    mutable = {"workups": 7}
    built = pivot.Counts(pivot.file_identity(populated), mutable)
    before = hash(built)
    mutable["workups"] = 99
    assert built["workups"] == 7
    assert hash(built) == before
