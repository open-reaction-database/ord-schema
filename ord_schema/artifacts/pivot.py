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

"""The repeated levels of the projection, and what a pivot over one holds.

A quantifier asks whether some element of a repeated level satisfies a body. Answering
that from the projection means decoding the whole nested column and filtering its
elements; answering it from a *pivot* -- one row per element -- means a semi-join
against a flat table. Over ORD the second is milliseconds where the first is seconds,
because the cost was never the predicate but the reconstruction of lists of structs.

A pivot row carries three things:

* ``reaction_id``, which is what the semi-join joins on.
* The ordinal of every repeated level from the root down to and including this one.
  Binding an element in the nested form gives co-membership structurally -- a
  measurement reached through a product *belongs* to it. A flat row has to say so, and
  a correlation joining on anything short of the full prefix over-returns.
* ``element``, a struct holding the element's own fields with every repeated field
  removed, recursively. Struct nesting is kept, so a path reaches the same leaf by the
  same spelling on either route, and a body reaching a dropped field simply fails to
  resolve -- which is how a quantifier the pivot cannot answer declines to the
  projection without a flag anyone has to keep honest.

Keeping the element's fields inside one struct rather than hoisting them is what makes
the key columns safe: ``inputs.components.preparations`` carries a ``reaction_id`` of
its own, and hoisted it would collide with the row's.

A pivot is unfiltered -- every element gets a row, including one whose fields are all
NULL. That completeness is what lets it answer ``forall``, which an index holding only
the elements that match something never can.
"""

import dataclasses
import os
import pathlib

import duckdb
import inflection
import pyarrow as pa
import pyarrow.parquet as pq

from ord_schema import atomic_io
from ord_schema.artifacts import base, projection
from ord_schema.logging import get_logger

logger = get_logger(__name__)

# The artifact name a pivot is stamped with, and the footer key naming which level it
# holds. One artifact name covers every level, because they are one kind of thing
# derived one way; the path is what tells them apart.
ARTIFACT = "pivot"
_META_PIVOT_PATH = "ord.pivot_path"

# The column the element's own fields sit under, and the row's join key.
ELEMENT = "element"
REACTION_ID = "reaction_id"

# What an ordinal is stored as. A reaction cannot hold four billion of anything at one
# level, and the narrower column is the one a pivot carries once per row.
ORDINAL_TYPE = pa.uint32()
ORDINAL_SQL = "UINTEGER"


@dataclasses.dataclass(frozen=True)
class Step:
    """One repeated level on the way to another, and how to reach it.

    Attributes:
        segments: The dotted hops from the previous element -- or from the row, for the
            first step -- down to the repeated field itself.
        is_map: Whether the repeated field is a map, whose values are the elements.
        ordinal: The column naming this element's position among its siblings.
    """

    segments: tuple[str, ...]
    is_map: bool
    ordinal: str

    def expression(self, source: str | None) -> str:
        """Returns the list to unnest at this step.

        Args:
            source: The variable bound to the previous element, or None at the row.

        Returns:
            A DuckDB expression evaluating to a list. A map contributes its values,
            since those are the elements a query quantifies over.
        """
        parts = self.segments if source is None else (source, *self.segments)
        reached = ".".join(parts)
        return f"map_values({reached})" if self.is_map else reached


@dataclasses.dataclass(frozen=True)
class RepeatedLevel:
    """A repeated level of the projection, and the shape of a pivot over it.

    Attributes:
        path: The level as the query grammar names it, with no wrapper segments.
        steps: One per repeated level from the root down to and including this one,
            outermost first. Unnesting them in order is what a build walks, and their
            ordinals are what a row carries to say which element it was.
        element_type: The element's struct type with repeated fields removed
            recursively, which is what a body must resolve against to be answerable
            from the pivot.
    """

    path: str
    steps: tuple[Step, ...]
    element_type: pa.DataType

    @property
    def ordinals(self) -> tuple[str, ...]:
        """Returns the ordinal column names, outermost first."""
        return tuple(step.ordinal for step in self.steps)


def _prune(dtype: pa.DataType) -> pa.DataType | None:
    """Returns ``dtype`` without its repeated fields, or None if none are left.

    Field metadata is carried through, since a comparison against an enum is checked
    against the members recorded there and a pruned type that dropped them would refuse
    a query the projection accepts.

    Args:
        dtype: A struct type to prune.

    Returns:
        The struct holding only the scalar and struct fields, recursively pruned, or
        None if pruning empties it -- a level whose elements are entirely repeated has
        nothing a pivot could hold.
    """
    fields = []
    for field in dtype:
        if pa.types.is_list(field.type) or pa.types.is_map(field.type):
            continue
        if pa.types.is_struct(field.type):
            pruned = _prune(field.type)
            if pruned is not None:
                fields.append(
                    pa.field(
                        field.name,
                        pruned,
                        nullable=field.nullable,
                        metadata=field.metadata,
                    )
                )
        else:
            fields.append(field)
    return pa.struct(fields) if fields else None


def repeated_levels(
    schema: pa.Schema = projection.SCHEMA,
) -> dict[str, RepeatedLevel]:
    """Returns every repeated level of ``schema``, by the path a query names it.

    Walked rather than listed by hand, in the same spirit as
    ``structures.structure_levels``: a repeated field added upstream becomes a level
    nobody had to remember to add, and one walk cannot disagree with itself the way two
    lists can.

    A list or a map is a repeated level; both are levels a query quantifies over, and
    the grammar spells neither with a wrapper segment, so ``inputs.components`` names
    the components under a map of inputs.

    Args:
        schema: Schema to walk; the projection schema by default.

    Returns:
        Each level's path mapped to what a pivot over it holds. A list of scalars is
        absent -- there is no element struct to pivot -- and so is a level whose
        elements carry nothing but repeated fields.

    Raises:
        ValueError: If two levels on one path yield the same ordinal name, which would
            give a row two columns of that name and let the inner one silently win.
            Raised at import, naming the ordinal, because the answer is to name it
            differently here rather than to retry.
    """
    levels: dict[str, RepeatedLevel] = {}

    def walk(
        dtype: pa.DataType,
        path: str,
        steps: tuple[Step, ...],
        pending: tuple[str, ...],
    ) -> None:
        if pa.types.is_list(dtype) or pa.types.is_map(dtype):
            is_map = pa.types.is_map(dtype)
            value = dtype.item_type if is_map else dtype.value_type
            ordinal = f"{inflection.singularize(path.rsplit('.', 1)[-1])}_index"
            if ordinal in [step.ordinal for step in steps]:
                raise ValueError(
                    f"{path} would carry two ordinals named {ordinal}, so a row's "
                    "outer position would be overwritten by its inner one"
                )
            steps = (*steps, Step(pending, is_map, ordinal))
            if pa.types.is_struct(value):
                element = _prune(value)
                if element is not None:
                    levels[path] = RepeatedLevel(path, steps, element)
            walk(value, path, steps, ())
        elif pa.types.is_struct(dtype):
            for field in dtype:
                walk(
                    field.type,
                    f"{path}.{field.name}" if path else field.name,
                    steps,
                    (*pending, field.name),
                )

    for field in schema:
        walk(field.type, field.name, (), (field.name,))
    return levels


def select(level: RepeatedLevel, table: str) -> str:
    """Returns the query building the pivot over ``level``.

    Unnests one repeated level at a time rather than taking the flattened traversal the
    compiler uses, because the ordinals are what each unnest contributes and a list
    already flattened has lost them.

    Args:
        level: The level to pivot.
        table: The relation holding the reactions.

    Returns:
        A SELECT over ``table``, one row per element of the level.
    """
    froms, source = [table], None
    for depth, step in enumerate(level.steps):
        froms.append(
            f"unnest({step.expression(source)}) WITH ORDINALITY "
            f"AS u{depth}(e{depth}, {step.ordinal})"
        )
        source = f"e{depth}"
    # Cast, because WITH ORDINALITY counts in BIGINT and an artifact's schema is a
    # promise rather than whatever the build happened to produce.
    ordinals = [f"{ordinal}::{ORDINAL_SQL} AS {ordinal}" for ordinal in level.ordinals]
    selected = [
        REACTION_ID,
        *ordinals,
        f"{element_expression(level.element_type, str(source))} AS {ELEMENT}",
    ]
    # S608: every fragment is this module's own walk of the projection schema.
    return f"SELECT {', '.join(selected)} FROM {', '.join(froms)}"  # noqa: S608


def schema(level: RepeatedLevel) -> pa.Schema:
    """Returns the schema of a pivot artifact over ``level``.

    Args:
        level: The level the artifact holds.

    Returns:
        ``reaction_id``, one ordinal per repeated level above and including this one,
        and the pruned element.
    """
    fields = [pa.field(REACTION_ID, pa.string())]
    fields += [pa.field(ordinal, ORDINAL_TYPE) for ordinal in level.ordinals]
    fields.append(pa.field(ELEMENT, level.element_type))
    return pa.schema(fields)


def element_expression(dtype: pa.DataType, source: str) -> str:
    """Returns a struct literal rebuilding ``source`` without its repeated fields.

    The pivot stores what a predicate can read rather than the whole element, which is
    where the size comes down: the lists a build already unnested into their own levels
    would otherwise be carried again on every row above them.

    Args:
        dtype: The pruned element type, as ``RepeatedLevel.element_type`` gives it.
        source: The DuckDB expression bound to the element.

    Returns:
        A struct literal holding each kept field, recursing into kept structs.
    """
    parts = []
    for field in dtype:
        reached = f"{source}.{field.name}"
        value = (
            element_expression(field.type, reached)
            if pa.types.is_struct(field.type)
            else reached
        )
        parts.append(f"'{field.name}': {value}")
    return "{" + ", ".join(parts) + "}"


def reach(
    path: str, levels: dict[str, RepeatedLevel] | None = None
) -> tuple[RepeatedLevel, tuple[str, ...], pa.DataType] | None:
    """Returns the pivot a quantifier over ``path`` ranges within, if one does.

    A quantifier's path need not name a repeated level itself. Descending from one
    through singular struct fields reaches one value per element rather than a list of
    its own -- ``outcomes.products.measurements.authentic_standard`` is one compound per
    measurement -- so the level it ranges over is the nearest repeated ancestor, and the
    pivot over that level already carries the struct.

    Args:
        path: The dotted path the quantifier ranges over.
        levels: Levels to search; the projection's by default.

    Returns:
        The level whose pivot holds it, the remaining segments from that level's
        element down to the value, and the type they reach -- or None if no repeated
        ancestor covers the path, or the descent leaves the element's pruned type,
        which is where every repeated field it dropped would have been.
    """
    found = LEVELS if levels is None else levels
    segments = path.split(".")
    for cut in range(len(segments), 0, -1):
        level = found.get(".".join(segments[:cut]))
        if level is None:
            continue
        remainder = tuple(segments[cut:])
        dtype = level.element_type
        for segment in remainder:
            if not pa.types.is_struct(dtype) or segment not in [
                field.name for field in dtype
            ]:
                return None
            dtype = dtype.field(segment).type
        return level, remainder, dtype
    return None


def table_name(path: str) -> str:
    """Returns the relation name holding the pivot over ``path``.

    Args:
        path: The level as the query grammar names it.

    Returns:
        A bare identifier, so it reaches SQL without quoting.
    """
    return f"pivot_{path.replace('.', '_')}"


LEVELS = repeated_levels()


def pivot_path(path: str | os.PathLike[str]) -> str | None:
    """Returns the level a pivot artifact holds, or None if it holds none.

    Args:
        path: Path to a Parquet file.

    Returns:
        The dotted level path stamped in the footer, or None if the file is not a
        pivot artifact.
    """
    metadata = pq.read_schema(path).metadata or {}
    value = metadata.get(_META_PIVOT_PATH.encode())
    return value.decode() if value is not None else None


def count_expression(level: RepeatedLevel) -> str:
    """Returns SQL counting the elements a projection records at ``level``.

    Reads list offsets rather than element bodies: ``len`` needs to know how long a
    list is, not what is in it, so this touches a fraction of what unnesting the level
    would. Most levels of this schema are never recorded, and discovering that by
    unnesting costs a full pass over every ancestor.

    Args:
        level: The level to count.

    Returns:
        A scalar expression over the reactions relation, zero where nothing is
        recorded. Flattened at each step, so the count is of elements rather than of
        the lists holding them.
    """
    expression = level.steps[0].expression(None)
    for depth, step in enumerate(level.steps[1:], start=1):
        element = f"e{depth}"
        expression = (
            f"flatten(list_transform({expression}, "
            f"{element} -> {step.expression(element)}))"
        )
    return f"coalesce(sum(len({expression})), 0)"


def write_pivot(
    source: str | os.PathLike[str],
    output: str | os.PathLike[str],
    /,
    *,
    level_path: str,
    source_md5: str | None = None,
    source_dataset_id: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 100_000,
) -> int:
    """Derives a pivot artifact over one repeated level and writes it.

    One artifact holds one level of one projection, because the levels have different
    schemas and a reader globs the one it wants. The rows stream out of DuckDB a batch
    at a time, so a level that explodes -- a corpus of long measurement lists -- costs
    a batch of memory rather than the whole level.

    The output is published atomically, so a failure partway leaves any existing
    artifact untouched.

    Args:
        source: Path to the projection to pivot.
        output: Path to write the artifact to.
        level_path: The repeated level to pivot, as the query grammar names it.
        source_md5: Hash of the *source dataset* to stamp, if the caller already read
            one. Taken from the projection's own stamps when omitted, so the artifact
            names the dataset it reflects rather than the file it read.
        source_dataset_id: Source dataset ID to stamp, if the caller already read one.
            Taken from the projection's stamps when omitted.
        compression: Parquet codec, any name ``pq.ParquetWriter`` accepts.
        row_group_size: Rows per output row group, which is also how many elements are
            held in memory at a time.

    Returns:
        Number of rows written: the number of elements at that level.

    Raises:
        ValueError: If ``source`` is not a current projection, or the schema has no
            such repeated level.
    """
    level = LEVELS.get(level_path)
    if level is None:
        raise ValueError(
            f"{level_path} is not a repeated level of the projection schema; "
            f"a pivot has nothing to hold. Known levels: {sorted(LEVELS)}"
        )
    parent = base.load_stamps(source)
    if parent.artifact != projection.ARTIFACT:
        raise ValueError(
            f"{source} is a {parent.artifact}, not a {projection.ARTIFACT}; a pivot "
            "unnests a projection"
        )
    # derive_tree refuses stale parents, but this writer is public and its output
    # inherits the dataset hash: an artifact derived from a stale projection would
    # claim a provenance it does not have and nothing would ever mark it stale again.
    if not base.stamps_are_current(parent, projection.ARTIFACT):
        raise ValueError(
            f"{source} is a stale {projection.ARTIFACT}; derive it again first"
        )
    if source_md5 is None:
        source_md5 = parent.source_md5
    if source_dataset_id is None:
        source_dataset_id = parent.source_dataset_id
    stamps = base.current_stamps(ARTIFACT, source_dataset_id, source_md5)
    metadata = base.to_metadata(stamps) | {_META_PIVOT_PATH: level_path}
    target = schema(level).with_metadata(metadata)
    connection = duckdb.connect()
    written = 0
    try:
        # Through the relational API rather than a path interpolated into SQL, which a
        # filename holding a quote would otherwise close.
        connection.read_parquet(str(pathlib.Path(source))).create_view("reactions")
        # S608: the expression is this module's own walk of the projection schema.
        counted = connection.execute(
            f"SELECT {count_expression(level)} FROM reactions"  # noqa: S608
        ).fetchone()
        if counted is not None and not counted[0]:
            # Counting read the offsets; unnesting would have walked every ancestor of
            # a level this projection never records, to arrive at the same nothing.
            with (
                atomic_io.atomic_path(output) as temp_path,
                pq.ParquetWriter(temp_path, target, compression=compression),
            ):
                pass
            logger.info("%s: no elements at %s", source, level_path)
            return 0
        result = connection.execute(select(level, "reactions"))
        with (
            atomic_io.atomic_path(output) as temp_path,
            pq.ParquetWriter(temp_path, target, compression=compression) as writer,
        ):
            # A projection whose reactions record nothing at this level leaves the loop
            # empty; the writer's close still produces a valid zero-row file carrying
            # the schema and stamps, which is what a reader globs and finds nothing in
            # rather than a file that is missing.
            for batch in result.to_arrow_reader(row_group_size):
                # Cast rather than trust: DuckDB's own widths are its business, and the
                # artifact's schema is a promise to whoever reads it later.
                writer.write_table(
                    pa.Table.from_batches([batch])
                    .cast(schema(level))
                    .replace_schema_metadata(metadata)
                )
                written += batch.num_rows
    finally:
        connection.close()
    logger.info("%s: %d elements at %s", source, written, level_path)
    return written
