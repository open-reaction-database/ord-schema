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

"""Statistics and charts for each column a dataset varies.

One entry per group of :attr:`~ord_schema.analysis.flat.FlatTable.column_groups`, so
a compound's SMILES and name make one entry rather than two. A group whose values
are numbers, with more than a handful of distinct ones, is *numeric*: it gets
quartiles and a histogram. Everything else is *categorical*: it gets each distinct
value with its count, every column of the group side by side, so a name stays
paired with its structure.

Charts are inline SVG, drawn here rather than with a plotting library so the
``analysis`` extra needs nothing beyond the base install and a report embedding them
stays one self-contained file. Colors are CSS variables with light-theme fallbacks
(``--series-1``, ``--text-secondary``, ``--grid``), so a page can restyle them, for a
dark theme say, without redrawing::

    table = flat.flatten(dataset.reactions)
    stats = column_statistics(table)
    write_csv(stats, "column_stats.csv")
    svg = stats[0].to_svg()
"""

import csv
import dataclasses
import html
import math
import os
import pathlib
from typing import cast

import numpy as np
import pandas as pd

from ord_schema.analysis import flat

# A numeric column with at most this many distinct values is shown as categories:
# a temperature set at 25 and 60 is two levels, not a distribution.
MAX_DISCRETE_NUMERIC = 12

# Categories drawn before the rest fold into one "other" bar; the CSV keeps more.
MAX_CHART_CATEGORIES = 12
MAX_LISTED_CATEGORIES = 50

HISTOGRAM_BINS = 20

# Chart geometry, in pixels.
_WIDTH = 360
_BAR = 18  # Categorical bar thickness; marks stay under 24px.
_ROW = 26
_LABEL_WIDTH = 150
_HIST_HEIGHT = 120
_GAP = 2  # The surface gap between touching marks.
_RADIUS = 4  # Rounded data-end; the baseline end stays square.
_FONT = 11


@dataclasses.dataclass
class ColumnStats:
    """What one group of varying columns holds.

    Attributes:
        column: The column that stands for the group.
        columns: Every column of the group, ``column`` first.
        section: ``experiment`` or ``provenance``; see :func:`flat.is_provenance`.
        kind: ``numeric`` or ``categorical``.
        filled: How many reactions hold a value.
        empty: How many do not.
        distinct: How many distinct values (combinations, for a group) there are,
            not counting empty.
        summary: For a numeric group: ``min``, ``q1``, ``median``, ``q3``, ``max``,
            ``mean`` and ``std``. Empty for a categorical one.
        histogram: For a numeric group: the bin edges and the count in each bin.
        categories: For a categorical group: each distinct combination of the
            group's values, with its count, most common first, up to
            ``MAX_LISTED_CATEGORIES``; an empty value shows as None.
        other: How many reactions hold values beyond those in ``categories``.
    """

    column: str
    columns: list[str]
    section: str
    kind: str
    filled: int
    empty: int
    distinct: int
    summary: dict[str, float] = dataclasses.field(default_factory=dict)
    histogram: tuple[list[float], list[int]] | None = None
    categories: list[tuple[tuple[object, ...], int]] = dataclasses.field(
        default_factory=list
    )
    other: int = 0

    def values_not_listed(self, listed: int) -> int:
        """How many distinct values the first ``listed`` categories leave out.

        ``distinct`` does not count the empty value, but the categories can include
        it, so only the non-empty categories listed are taken off.
        """
        shown = sum(
            1
            for values, _ in self.categories[:listed]
            if any(v is not None for v in values)
        )
        return self.distinct - shown

    def to_svg(self, *, title: str | None = None) -> str:
        """Returns the group's chart as an inline SVG element.

        Args:
            title: The chart's accessible name; by default the column.
        """
        label = title or self.column
        if self.kind == "numeric":
            return _histogram_svg(self, label)
        return _categories_svg(self, label)


def _is_number(value: object) -> bool:
    return isinstance(value, int | float | np.integer | np.floating) and not isinstance(
        value, bool | np.bool_
    )


def _format(value: object) -> str:
    """A value as it reads in a chart or CSV: numbers to 4 significant figures."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "(empty)"
    if _is_number(value):
        number = float(cast("float", value))
        if number.is_integer() and abs(number) < 1e15:
            return f"{int(number):,}"
        return f"{number:,.4g}"
    return str(value)


def _stats_for(table: flat.FlatTable, group: list[str]) -> ColumnStats:
    frame = table.varying
    column = frame[group[0]]
    filled = int(column.notna().sum())
    values = column.dropna()
    # Judged on the group's first column; a group of numbers pairs one to one, as
    # a yield stored twice does, so the others follow it.
    numeric = (
        filled > 0
        and all(_is_number(v) for v in values)
        and values.nunique() > MAX_DISCRETE_NUMERIC
    )
    stats = ColumnStats(
        column=group[0],
        columns=list(group),
        section="provenance" if flat.is_provenance(group[0]) else "experiment",
        kind="numeric" if numeric else "categorical",
        filled=filled,
        empty=table.reaction_count - filled,
        distinct=0,
    )
    if numeric:
        data = values.astype(float).to_numpy()
        counts, edges = np.histogram(data, bins=HISTOGRAM_BINS)
        q1, median, q3 = np.percentile(data, [25, 50, 75])
        stats.distinct = int(values.nunique())
        stats.summary = {
            "min": float(data.min()),
            "q1": float(q1),
            "median": float(median),
            "q3": float(q3),
            "max": float(data.max()),
            "mean": float(data.mean()),
            "std": float(data.std(ddof=1)) if len(data) > 1 else 0.0,
        }
        stats.histogram = ([float(e) for e in edges], [int(c) for c in counts])
        return stats
    unique = table.unique_values(group[0])
    rows = [
        (
            tuple(None if pd.isna(v) else v for v in row[: len(group)]),
            int(row[len(group)]),
        )
        for row in unique.itertuples(index=False, name=None)
    ]
    stats.distinct = sum(1 for kept, _ in rows if any(v is not None for v in kept))
    stats.categories = rows[:MAX_LISTED_CATEGORIES]
    stats.other = sum(count for _, count in rows[MAX_LISTED_CATEGORIES:])
    return stats


def column_statistics(table: flat.FlatTable) -> list[ColumnStats]:
    """Returns statistics for each group of the table's varying columns.

    ``reaction_id`` is left out: it identifies rows rather than describing them.

    Args:
        table: The dataset's flattened table.

    Returns:
        One entry per group, in table order.
    """
    return [
        _stats_for(table, group)
        for group in table.column_groups
        if group[0] != "reaction_id"
    ]


_CSV_FIELDS = [
    "column",
    "columns",
    "section",
    "kind",
    "filled",
    "empty",
    "distinct",
    "min",
    "q1",
    "median",
    "q3",
    "max",
    "mean",
    "std",
    "values",
]


def write_csv(stats: list[ColumnStats], path: str | os.PathLike[str]) -> None:
    """Writes one row per group: its counts, quartiles, and most common values.

    A categorical group's ``values`` lists each combination as its columns' values
    joined by `` | ``, with its count, separated by ``; ``.
    """
    with pathlib.Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_CSV_FIELDS)
        writer.writeheader()
        for stat in stats:
            row: dict[str, object] = {
                "column": stat.column,
                "columns": " | ".join(stat.columns),
                "section": stat.section,
                "kind": stat.kind,
                "filled": stat.filled,
                "empty": stat.empty,
                "distinct": stat.distinct,
                "values": "; ".join(
                    f"{' | '.join(_format(v) for v in values)} ({count})"
                    for values, count in stat.categories
                ),
            }
            row.update(stat.summary)
            writer.writerow(row)


def _style() -> str:
    return (
        f'font-family="system-ui, -apple-system, Segoe UI, sans-serif" '
        f'font-size="{_FONT}"'
    )


def _bar_path(x: float, y: float, width: float, height: float, *, up: bool) -> str:
    """A bar with a rounded data-end and a square baseline end.

    ``up`` draws a column growing upward from its bottom edge; otherwise a bar
    growing rightward from its left edge.
    """
    if up:
        r = min(_RADIUS, width / 2, height)
        return (
            f"M{x:.1f},{y + height:.1f}V{y + r:.1f}Q{x:.1f},{y:.1f} {x + r:.1f},{y:.1f}"
            f"H{x + width - r:.1f}Q{x + width:.1f},{y:.1f} {x + width:.1f},{y + r:.1f}"
            f"V{y + height:.1f}Z"
        )
    r = min(_RADIUS, height / 2, width)
    return (
        f"M{x:.1f},{y:.1f}H{x + width - r:.1f}Q{x + width:.1f},{y:.1f} "
        f"{x + width:.1f},{y + r:.1f}V{y + height - r:.1f}Q{x + width:.1f},"
        f"{y + height:.1f} {x + width - r:.1f},{y + height:.1f}H{x:.1f}Z"
    )


def _text(
    x: float, y: float, value: str, *, anchor: str = "start", muted: bool = True
) -> str:
    color = (
        "var(--text-secondary, #52514e)" if muted else "var(--text-primary, #0b0b0b)"
    )
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" fill="{color}">'
        f"{html.escape(value)}</text>"
    )


def _text_width(value: str) -> float:
    """Roughly how wide ``value`` renders at the chart's font size, in pixels.

    Errs wide, so a label cut to fit does not reach the chart's edge: capitals,
    digits and the brackets, ``@`` and ``=`` that SMILES are full of count wider
    than lowercase letters.
    """
    narrow = set("ijlrtf.,:;|!()[]'` ")
    wide = set("@%#=MWmwQOGDHNUCB")
    return sum(
        0.42 * _FONT if c in narrow else 0.82 * _FONT if c in wide else 0.62 * _FONT
        for c in value
    )


def _fit(value: str, room: float) -> str:
    """``value`` cut, with an ellipsis, to fit in ``room`` pixels."""
    if _text_width(value) <= room:
        return value
    while value and _text_width(value + "…") > room:
        value = value[:-1]
    return value + "…"


# Room for a category label, in pixels; a longer label is cut, and the whole value
# stays in the bar's tooltip.
_LABEL_ROOM = _LABEL_WIDTH - 12


def _categories_svg(stat: ColumnStats, label: str) -> str:
    shown = stat.categories[:MAX_CHART_CATEGORIES]
    rest = (
        sum(count for _, count in stat.categories[MAX_CHART_CATEGORIES:]) + stat.other
    )
    bars = [(" | ".join(_format(v) for v in values), count) for values, count in shown]
    # Scaled to the categories drawn: the folded rest is a line of text, not a bar,
    # since its total would dwarf them.
    top = max((count for _, count in bars), default=1)
    plot = _WIDTH - _LABEL_WIDTH - 48
    rows = len(bars) + (1 if rest else 0)
    height = rows * _ROW + 8
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {_WIDTH} {height}" '
        f'width="{_WIDTH}" height="{height}" role="img" {_style()}>',
        f"<title>{html.escape(label)}</title>",
        f'<line x1="{_LABEL_WIDTH}" y1="2" x2="{_LABEL_WIDTH}" '
        f'y2="{len(bars) * _ROW + 4}" stroke="var(--grid, #e2e1dc)" stroke-width="1"/>',
    ]
    for index, (name, count) in enumerate(bars):
        y = 4 + index * _ROW + (_ROW - _BAR) / 2
        width = max(1.0, plot * count / top)
        tooltip = f"{name}: {count:,} reactions"
        parts.append(
            f"<g><title>{html.escape(tooltip)}</title>"
            f'<path d="{_bar_path(_LABEL_WIDTH, y, width, _BAR, up=False)}" '
            'fill="var(--series-1, #2a78d6)"/>'
            + _text(
                _LABEL_WIDTH - 6,
                y + _BAR / 2 + 4,
                _fit(name, _LABEL_ROOM),
                anchor="end",
            )
            + _text(
                _LABEL_WIDTH + width + 4, y + _BAR / 2 + 4, f"{count:,}", muted=False
            )
            + "</g>"
        )
    if rest:
        folded = stat.values_not_listed(len(shown))
        y = 4 + len(bars) * _ROW + _ROW / 2 + 4
        parts.append(
            _text(
                _LABEL_WIDTH + 2,
                y,
                f"+ {folded:,} more values, {rest:,} reactions",
            )
        )
    parts.append("</svg>")
    return "".join(parts)


def _histogram_svg(stat: ColumnStats, label: str) -> str:
    assert stat.histogram is not None
    edges, counts = stat.histogram
    left, bottom, top_pad = 36, 22, 14
    plot_width = _WIDTH - left - 8
    height = _HIST_HEIGHT + bottom + top_pad
    top = max(counts) or 1
    slot = plot_width / len(counts)
    baseline = top_pad + _HIST_HEIGHT
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {_WIDTH} {height}" '
        f'width="{_WIDTH}" height="{height}" role="img" {_style()}>',
        f"<title>{html.escape(label)}</title>",
        f'<line x1="{left}" y1="{baseline}" x2="{_WIDTH - 8}" y2="{baseline}" '
        'stroke="var(--grid, #e2e1dc)" stroke-width="1"/>',
        _text(left - 4, top_pad + 4, f"{top:,}", anchor="end"),
        _text(left - 4, baseline, "0", anchor="end"),
    ]
    for index, count in enumerate(counts):
        if not count:
            continue
        bar_height = _HIST_HEIGHT * count / top
        x = left + index * slot + _GAP / 2
        tooltip = (
            f"{_format(edges[index])} to {_format(edges[index + 1])}: "
            f"{count:,} reactions"
        )
        path = _bar_path(x, baseline - bar_height, slot - _GAP, bar_height, up=True)
        parts.append(
            f"<g><title>{html.escape(tooltip)}</title>"
            f'<path d="{path}" fill="var(--series-1, #2a78d6)"/></g>'
        )
    parts.append(_text(left, baseline + 15, _format(edges[0])))
    parts.append(_text(_WIDTH - 8, baseline + 15, _format(edges[-1]), anchor="end"))
    parts.append(
        _text(
            left + plot_width / 2,
            baseline + 15,
            f"median {_format(stat.summary['median'])}",
            anchor="middle",
        )
    )
    parts.append("</svg>")
    return "".join(parts)
