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

"""A dataset as a table: one row per reaction, one column per field that varies.

Columns are the ones :func:`ord_schema.message_helpers.message_to_row` writes, e.g.
``inputs["base"].components[0].identifiers[0].value``: a field at one position in
one input. A column is *constant* when every reaction holds it with the same value;
those are left out of the table and listed separately, with their values, so a reader
sees what the dataset varies and, beside it, what it held fixed.

The columns kept and dropped are exactly those of
``message_helpers.messages_to_dataframe(reactions, drop_constant_columns=True)``,
with ``reaction_id`` moved to the front.
That helper builds every column of every reaction before dropping the constant ones,
and does not return them. This reads the reactions once, keeping values only for
columns that have varied so far: while a column is constant, all it needs to remember
is which reactions hold it, since they all hold the same value. A large dataset
therefore never has its constant columns in memory::

    table = flatten(datasets.load_dataset("dataset.parquet").reactions)
    table.varying  # pandas DataFrame, one row per reaction
    table.constants  # {column: value}

Columns are positional, so a dataset that lists the same components in a different
order from one reaction to the next shows them as varying columns; see
:meth:`ord_schema.analysis.summary.DatasetSummary.reordered_fields`.
"""

import array
import dataclasses
import os
import pathlib
from collections.abc import Iterable

import pandas as pd

import ord_schema
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2


@dataclasses.dataclass
class _Column:
    """One column, as read so far.

    Attributes:
        first: The first value read.
        rows: While every value read has equalled ``first``, the positions of the
            reactions that hold the column; empty once it has varied.
        varied: Once a value has differed from ``first``, every value read, by
            reaction position; None until then.
    """

    first: ord_schema.ScalarType
    rows: array.array = dataclasses.field(default_factory=lambda: array.array("L"))
    varied: dict[int, ord_schema.ScalarType] | None = None

    def add(self, position: int, value: ord_schema.ScalarType) -> None:
        """Records that the reaction at ``position`` holds ``value``."""
        if self.varied is None:
            if value == self.first:
                self.rows.append(position)
                return
            # The first difference: every reaction read so far held ``first``.
            self.varied = dict.fromkeys(self.rows, self.first)
            self.rows = array.array("L")
        self.varied[position] = value

    def all_values(self) -> dict[int, ord_schema.ScalarType]:
        """Every value read, by reaction position."""
        if self.varied is not None:
            return self.varied
        return dict.fromkeys(self.rows, self.first)


@dataclasses.dataclass
class FlatTable:
    """A dataset's varying fields as a table, and its constant fields beside it.

    Attributes:
        varying: One row per reaction, in dataset order, and one column per
            ``message_to_row`` column that varies: ``reaction_id`` first, then the
            rest in the order columns first appear.
            A reaction that does not populate a column holds NaN there.
        constants: Each column every reaction holds with one and the same value,
            mapped to that value, in the order columns first appear.
        reaction_count: How many reactions were read.
    """

    varying: pd.DataFrame
    constants: dict[str, ord_schema.ScalarType]
    reaction_count: int

    def constants_frame(self) -> pd.DataFrame:
        """The constant fields as a two-column table: ``field`` and ``value``."""
        return pd.DataFrame(
            {"field": list(self.constants), "value": list(self.constants.values())}
        )

    def write_csv(self, path: str | os.PathLike[str]) -> None:
        """Writes the varying fields, one row per reaction."""
        with pathlib.Path(path).open("w", newline="", encoding="utf-8") as f:
            self.varying.to_csv(f, index=False)

    def write_constants_csv(self, path: str | os.PathLike[str]) -> None:
        """Writes the constant fields, one row each: ``field`` and ``value``."""
        with pathlib.Path(path).open("w", newline="", encoding="utf-8") as f:
            self.constants_frame().to_csv(f, index=False)


def flatten(reactions: Iterable[reaction_pb2.Reaction]) -> FlatTable:
    """Returns the varying and constant fields of ``reactions``.

    Args:
        reactions: The dataset's reactions, read once, in order.

    Returns:
        The table of varying fields and the constant fields beside it.
    """
    columns: dict[str, _Column] = {}
    reaction_count = 0
    for position, reaction in enumerate(reactions):
        reaction_count += 1
        for name, value in message_helpers.message_to_row(reaction).items():
            column = columns.get(name)
            if column is None:
                columns[name] = column = _Column(first=value)
            column.add(position, value)

    # A column some reaction leaves empty varies, even if every value it has agrees.
    varying = [
        name
        for name, column in columns.items()
        if column.varied is not None or len(column.rows) != reaction_count
    ]
    keep = set(varying)
    constants = {
        name: column.first for name, column in columns.items() if name not in keep
    }
    # Built as messages_to_dataframe builds its frame, from one dict per reaction, so
    # pandas infers the same dtypes.
    rows: list[dict[str, ord_schema.ScalarType]] = [{} for _ in range(reaction_count)]
    for name in varying:
        for position, value in columns[name].all_values().items():
            rows[position][name] = value
    # The reaction ID first, so each row says which reaction it is; the rest in the
    # order message_to_row writes them.
    if "reaction_id" in keep:
        varying.remove("reaction_id")
        varying.insert(0, "reaction_id")
    frame = pd.DataFrame(rows, columns=pd.Index(varying))
    return FlatTable(varying=frame, constants=constants, reaction_count=reaction_count)
