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

"""Exemplar reactions: the handful a reviewer should read, and why each.

Sampling a large dataset at random mostly finds reactions just like the last one.
Instead this picks:

- one reaction of each **shape**, the set of fields a reaction fills, the most
  common shape first, so every way the dataset records a reaction is read once;
- the reactions with the **lowest and highest** value of each numeric field of the
  experiment, such as the yield, where transcription errors and edge cases show;
- under the ``uniform`` profile, a reaction holding an **input value no other
  reaction uses**, which in a screen is either a deliberate one-off or a mistake.

Each exemplar says why it was picked, and a reaction picked for several reasons is
listed once. :func:`write_txtpb` writes them as a readable subset of the dataset.
"""

import dataclasses
import os
import pathlib
from typing import Protocol

import numpy as np
from google.protobuf import text_format

from ord_schema.analysis import columns, flat
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.review import structure

# The most exemplars picked: a handful to read, not a sample to analyze.
MAX_EXEMPLARS = 20

# Reactions picked per field for holding a value no other reaction uses.
MAX_SINGLETONS_PER_FIELD = 2

# Fields named in a shape's description before the rest are counted.
_FIELDS_NAMED = 4


class _ByPosition(Protocol):
    """Reactions fetched by position: a ``Dataset``'s or a ``DatasetView``'s."""

    def __getitem__(self, position: int, /) -> reaction_pb2.Reaction: ...


@dataclasses.dataclass
class Exemplar:
    """A reaction picked for the reviewer to read.

    Attributes:
        position: Where the reaction is in the dataset, from 0.
        reaction_id: Its ID, or empty if it has none.
        reasons: Why it was picked, one line each.
        reaction: The reaction itself.
    """

    position: int
    reaction_id: str
    reasons: list[str]
    reaction: reaction_pb2.Reaction

    @property
    def label(self) -> str:
        """The reaction's ID, or its position if it has none."""
        return self.reaction_id or f"reaction #{self.position + 1}"


# A repeated field's index in a column name; map keys are kept.
def _keyed_path(column: str) -> str:
    """A column's field path, keeping which input or map entry it is in.

    The position in a list does not change what a reaction records, but which input
    it records it in does.
    """
    return flat.field_path(column, keep_keys=True)


def _shapes(table: flat.FlatTable) -> list[tuple[frozenset[str], list[int]]]:
    """Each shape the table's rows take, with the rows that take it, most common first.

    A shape is the set of keyed field paths a row fills, so a control that leaves
    out the ligand's amount has a different shape from the reactions that give it.
    This is finer than :attr:`DatasetSummary.shapes`, whose paths ignore which input
    a field is in. A path that some constant column holds is filled in every
    reaction, so only paths held in varying columns alone can tell shapes apart.
    """
    frame = table.varying
    if not table.reaction_count:
        return []
    always = {_keyed_path(c) for c in table.constants}
    by_path_columns: dict[str, list[str]] = {}
    for column in frame.columns:
        path = _keyed_path(column)
        if path not in always:
            by_path_columns.setdefault(path, []).append(column)
    paths = sorted(by_path_columns)
    if not paths:
        return [(frozenset(), list(range(table.reaction_count)))]
    filled = frame.notna()
    by_path = np.column_stack(
        [filled[by_path_columns[path]].any(axis=1).to_numpy() for path in paths]
    )
    rows: dict[bytes, list[int]] = {}
    for position, row in enumerate(by_path):
        rows.setdefault(np.packbits(row).tobytes(), []).append(position)
    shapes = [
        (
            frozenset(
                path
                for path, set_ in zip(paths, by_path[members[0]], strict=True)
                if set_
            ),
            members,
        )
        for members in rows.values()
    ]
    return sorted(shapes, key=lambda shape: (-len(shape[1]), shape[1][0]))


def _describe(shape: frozenset[str], common: frozenset[str]) -> str:
    parts = []
    for verb, fields in (("adds", shape - common), ("lacks", common - shape)):
        if fields:
            named = ", ".join(sorted(fields)[:_FIELDS_NAMED])
            more = len(fields) - _FIELDS_NAMED
            parts.append(f"{verb} {named}" + (f" and {more} more" if more > 0 else ""))
    return "; ".join(parts)


def _short(column: str) -> str:
    """A column, named well enough for a reason.

    Says which field, in which input or product, without the parts every column of
    its kind shares.
    """
    name = column.removeprefix("outcomes[0].").removesuffix(".value")
    return name.removesuffix(".float_value")


def pick(
    table: flat.FlatTable,
    stats: list[columns.ColumnStats],
    reactions: _ByPosition,
    profile: structure.ValidationProfile = structure.ValidationProfile.UNIFORM,
    *,
    max_exemplars: int = MAX_EXEMPLARS,
) -> list[Exemplar]:
    """Picks the reactions a reviewer should read.

    Args:
        table: The dataset's flattened table.
        stats: Its column statistics.
        reactions: The dataset's reactions, indexable by position: the
            ``reactions`` of a ``Dataset`` or of a Parquet ``DatasetView``.
        profile: The dataset's profile; one-off values are only picked under
            ``uniform``.
        max_exemplars: The most to pick.

    Returns:
        The exemplars, in the order picked: shapes, then extremes, then one-offs.
    """
    reasons: dict[int, list[str]] = {}

    def add(position: int, reason: str) -> None:
        if position in reasons:
            reasons[position].append(reason)
        elif len(reasons) < max_exemplars:
            reasons[position] = [reason]

    shapes = _shapes(table)
    for index, (shape, members) in enumerate(shapes):
        if index == 0 and len(members) > 1:
            add(members[0], f"the most common shape ({len(members):,} reactions)")
        elif index == 0:
            add(members[0], "the first reaction; no two reactions share a shape")
        elif len(members) == 1:
            add(members[0], f"a shape of its own: {_describe(shape, shapes[0][0])}")
        else:
            add(
                members[0],
                f"a shape {len(members):,} reactions share: "
                f"{_describe(shape, shapes[0][0])}",
            )

    frame = table.varying
    for stat in stats:
        if stat.kind != "numeric" or stat.section != "experiment":
            continue
        values = frame[stat.column]
        low, high = int(values.idxmin()), int(values.idxmax())
        add(low, f"lowest {_short(stat.column)} ({columns.format_value(values[low])})")
        add(
            high,
            f"highest {_short(stat.column)} ({columns.format_value(values[high])})",
        )

    if profile is structure.ValidationProfile.UNIFORM:
        for stat in stats:
            if (
                stat.kind != "categorical"
                or stat.section != "experiment"
                or not stat.column.startswith("inputs")
            ):
                continue
            column = frame[stat.column]
            counts = column.value_counts()
            for value in list(counts[counts == 1].index)[:MAX_SINGLETONS_PER_FIELD]:
                position = int(column.index[column == value][0])
                add(
                    position,
                    f"the only reaction with {_short(stat.column)} = "
                    f"{columns.format_value(value)}",
                )

    exemplars = []
    for position, why in reasons.items():
        reaction = reactions[position]
        exemplars.append(
            Exemplar(
                position=position,
                reaction_id=reaction.reaction_id,
                reasons=why,
                reaction=reaction,
            )
        )
    return exemplars


def to_dataset(exemplars: list[Exemplar], *, name: str) -> dataset_pb2.Dataset:
    """The exemplars as a dataset, each reaction exactly as the original holds it.

    The description lists why each was picked, in order, since a reaction has no
    field for a reviewer's note.
    """
    lines = [f"Exemplar reactions picked for review from {name or 'a dataset'}:"]
    lines.extend(
        f"{i}. {exemplar.label}: {'; '.join(exemplar.reasons)}"
        for i, exemplar in enumerate(exemplars, start=1)
    )
    return dataset_pb2.Dataset(
        name=f"{name} exemplars" if name else "Exemplars",
        description="\n".join(lines),
        reactions=[exemplar.reaction for exemplar in exemplars],
    )


def write_txtpb(
    exemplars: list[Exemplar], path: str | os.PathLike[str], *, name: str
) -> None:
    """Writes the exemplars as a text-format Dataset, to read or load.

    Named ``.txtpb`` by the report, the suffix the ORD reaction editor's dataset
    import accepts, so a reviewer can import the file as it is.
    """
    pathlib.Path(path).write_text(
        text_format.MessageToString(to_dataset(exemplars, name=name)), encoding="utf-8"
    )
