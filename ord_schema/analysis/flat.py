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

With ``date_resolution=None``, the columns kept and dropped are exactly those of
``message_helpers.messages_to_dataframe(reactions, drop_constant_columns=True)``,
with ``reaction_id`` moved to the front. At the default ``"day"``, two times on one
day are one value here, so a column that varies only within a day is constant.
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

Three things make the table easier to read without hiding anything:

- **Grouping.** Columns that always change together, such as a compound's SMILES
  and its name, form one group, so each design variable is one column of
  :meth:`FlatTable.grouped`. Listing a group's values with
  :meth:`FlatTable.unique_values` shows every column of the group side by side, so
  a name and the structure it stands for stay paired.
- **Dates.** ``DateTime`` fields are compared at a resolution, a day by default, so
  records created minutes apart do not count as different values.
- **Provenance.** Who recorded the data and when stays in the table, but
  :func:`is_provenance` tells it apart from the experiment, so a report can show it
  in its own section.
"""

import array
import dataclasses
import datetime
import functools
import os
import pathlib
import re
from collections.abc import Iterable
from typing import TYPE_CHECKING, Literal, cast

import pandas as pd

import ord_schema
from ord_schema import message_helpers, validations
from ord_schema.proto import reaction_pb2

if TYPE_CHECKING:
    from google.protobuf.descriptor import Descriptor, FieldDescriptor

# How finely DateTime values are compared, as strftime formats.
DateResolution = Literal["month", "day", "hour", "minute", "second"]
_DATE_FORMATS: dict[str, str] = {
    "month": "%Y-%m",
    "day": "%Y-%m-%d",
    "hour": "%Y-%m-%d %H:00",
    "minute": "%Y-%m-%d %H:%M",
    "second": "%Y-%m-%d %H:%M:%S",
}

# A map key or repeated-field index in a message_to_row column name.
_SUBSCRIPT = re.compile(r'\["[^"]*"\]|\[\d+\]')


def is_provenance(column: str) -> bool:
    """Whether ``column`` records who made the data and when, not the experiment."""
    return column.startswith("provenance.")


@functools.lru_cache(maxsize=4096)
def field_for(column: str) -> "FieldDescriptor | None":
    """Returns the proto field a column or field path names, or None if none.

    Accepts ``message_to_row`` columns and summary paths alike: map keys and
    repeated-field indices are ignored, and a map's entries are looked up through
    its value type. Read from the proto descriptors, so it follows the schema.
    """
    descriptor = cast("Descriptor", reaction_pb2.Reaction.DESCRIPTOR)
    names = _SUBSCRIPT.sub("", column).split(".")
    field = None
    for index, name in enumerate(names):
        field = descriptor.fields_by_name.get(name)
        if field is None:
            return None
        if index == len(names) - 1:
            break
        if field.message_type is None:
            return None
        descriptor = field.message_type
        if descriptor.GetOptions().map_entry:
            value = descriptor.fields_by_name["value"]
            if value.message_type is None:
                # A map of scalars: the next name is the entry's own "value".
                return value if names[index + 1 :] == ["value"] else None
            descriptor = value.message_type
    return field


def is_datetime(column: str) -> bool:
    """Whether ``column`` is the value of a ``DateTime`` message.

    Read from the proto descriptors, so a ``DateTime`` field added to the schema is
    recognized without a change here.
    """
    field = field_for(column)
    return (
        field is not None
        and field.containing_type is reaction_pb2.DateTime.DESCRIPTOR
        and field.name == "value"
    )


def _bin_date(
    value: ord_schema.ScalarType, resolution: DateResolution
) -> ord_schema.ScalarType:
    """Returns ``value`` cut to ``resolution``, or unchanged if it is not a DateTime.

    Read with :func:`ord_schema.validations.parse_date_time`, the schema's definition
    of a DateTime, so a value validation rejects is left as written: it is for a
    reviewer to see, not to be hidden by binning. A time is cut as written, in the
    zone it was recorded in. A date alone has no time to cut, so it stays a date at
    any resolution finer than a day.
    """
    if not isinstance(value, str):
        return value
    parsed = validations.parse_date_time(value)
    if parsed is None:
        return value
    if not isinstance(parsed, datetime.datetime) and resolution != "month":
        return parsed.isoformat()
    return parsed.strftime(_DATE_FORMATS[resolution])


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
        date_resolution: The resolution ``DateTime`` columns were compared at, and
            are shown at, or None if they were left as written.
    """

    varying: pd.DataFrame
    constants: dict[str, ord_schema.ScalarType]
    reaction_count: int
    date_resolution: DateResolution | None = "day"

    def constants_frame(self) -> pd.DataFrame:
        """The constant fields as a table: ``field``, ``value`` and ``section``.

        ``section`` is ``provenance`` or ``experiment``; see :func:`is_provenance`.
        """
        return pd.DataFrame(
            {
                "field": list(self.constants),
                "value": list(self.constants.values()),
                "section": [_section(name) for name in self.constants],
            }
        )

    def experiment_columns(self) -> list[str]:
        """The varying columns that describe the experiment, in table order."""
        return [c for c in self.varying.columns if not is_provenance(c)]

    def provenance_columns(self) -> list[str]:
        """The varying columns that record who made the data and when."""
        return [c for c in self.varying.columns if is_provenance(c)]

    @functools.cached_property
    def column_groups(self) -> list[list[str]]:
        """Varying columns that always change together, in table order.

        Two columns are grouped when each value of one always goes with the same
        value of the other, and the other way round; a column's absence counts as a
        value. A compound's SMILES and name group this way, and so do all the fields
        of a component that only some reactions have. Two kinds of column are never
        grouped: ``reaction_id``, and a column whose every value is different, since
        pairing them shows nothing but that both are unique.

        Returns:
            One list per group, the first column standing for the group.
        """
        groups: dict[bytes, list[str]] = {}
        order: list[list[str]] = []
        for name in self.varying.columns:
            codes, uniques = pd.factorize(self.varying[name], use_na_sentinel=False)
            if name == "reaction_id" or len(uniques) == len(codes):
                order.append([name])
                continue
            # factorize numbers values in order of first appearance, so two columns
            # get the same codes exactly when they pair one to one.
            key = codes.tobytes()
            if key in groups:
                groups[key].append(name)
            else:
                groups[key] = [name]
                order.append(groups[key])
        return order

    def grouped(self) -> pd.DataFrame:
        """The table with one column per group of :attr:`column_groups`."""
        return self.varying[[group[0] for group in self.column_groups]]

    def unique_values(self, column: str) -> pd.DataFrame:
        """The distinct values of ``column``'s group, every column side by side.

        Args:
            column: Any column of the group.

        Returns:
            One row per distinct combination, with a ``reactions`` count, most
            common first. Absent values show as NaN.

        Raises:
            KeyError: If ``column`` is not a varying column.
        """
        for group in self.column_groups:
            if column in group:
                counts = self.varying[group].value_counts(dropna=False, sort=True)
                return counts.rename("reactions").reset_index()
        raise KeyError(column)

    def write_csv(self, path: str | os.PathLike[str]) -> None:
        """Writes the varying fields, one row per reaction."""
        with pathlib.Path(path).open("w", newline="", encoding="utf-8") as f:
            self.varying.to_csv(f, index=False)

    def write_constants_csv(self, path: str | os.PathLike[str]) -> None:
        """Writes the constant fields, one row each: field, value and section."""
        with pathlib.Path(path).open("w", newline="", encoding="utf-8") as f:
            self.constants_frame().to_csv(f, index=False)


def _section(column: str) -> str:
    return "provenance" if is_provenance(column) else "experiment"


def flatten(
    reactions: Iterable[reaction_pb2.Reaction],
    *,
    date_resolution: DateResolution | None = "day",
) -> FlatTable:
    """Returns the varying and constant fields of ``reactions``.

    Args:
        reactions: The dataset's reactions, read once, in order.
        date_resolution: How finely to compare ``DateTime`` values, which the table
            then shows at that resolution: by default the day, so records made
            minutes apart agree. None compares them as written.

    Returns:
        The table of varying fields and the constant fields beside it.
    """
    columns: dict[str, _Column] = {}
    reaction_count = 0
    for position, reaction in enumerate(reactions):
        reaction_count += 1
        for name, value in message_helpers.message_to_row(reaction).items():
            if date_resolution is not None and is_datetime(name):
                value = _bin_date(value, date_resolution)  # noqa: PLW2901
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
    return FlatTable(
        varying=frame,
        constants=constants,
        reaction_count=reaction_count,
        date_resolution=date_resolution,
    )
