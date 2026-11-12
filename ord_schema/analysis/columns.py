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
from collections.abc import Callable, Collection
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

    def to_html(
        self,
        *,
        title: str | None = None,
        drawing: Callable[[object], str | None] | None = None,
    ) -> str:
        """Returns the group's chart, with each value in full where its bar can't.

        A bar whose label is cut short, whose value runs over several lines (a
        MOLBLOCK, say) or that has a drawing shows a panel while hovered or
        focused: the value in full, its count and the drawing. Under such a chart,
        a folded list gives every bar's value in full, to read and copy.

        Args:
            title: The chart's accessible name; by default the column.
            drawing: Returns the HTML to show for a value, e.g. a structure drawing
                for a SMILES, or None for none. A bar of a group's values shows the
                first of them that has one.

        Returns:
            The bare SVG when every bar shows its value whole and none has a
            drawing. Otherwise the SVG and its panels in a ``div.chart``, and the
            values list; a page styles them with :data:`TOOLTIP_CSS`, and must let
            a panel overflow the chart's box.
        """
        if self.kind == "numeric":
            return self.to_svg(title=title)
        panels = _panels(self, drawing)
        svg = _categories_svg(self, title or self.column, set(panels))
        if not panels:
            return svg
        tips = "".join(
            f'<div class="moltip t{index}" aria-hidden="true" '
            f'style="left: 8px; top: {4 + (index + 1) * _ROW}px">{panel}</div>'
            for index, panel in panels.items()
        )
        return f'<div class="chart">{svg}{tips}</div>{_values_list(self)}'


def _is_number(value: object) -> bool:
    return isinstance(value, int | float | np.integer | np.floating) and not isinstance(
        value, bool | np.bool_
    )


def format_value(value: object) -> str:
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
                    f"{' | '.join(format_value(v) for v in values)} ({count})"
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


# How many lines of a multi-line value, such as a MOLBLOCK, a hover panel shows;
# the values list under the chart shows them all.
MAX_TIP_LINES = 15

# The hover panel and the values list. In a chart the panel opens below its bar;
# the page must let it overflow the chart's box. It ignores the pointer, so it
# never takes the hover from the bar it covers. Each bar i shows the panel t{i}.
TOOLTIP_CSS = (
    ".moltip { display: none; padding: 6px; background: #fff; color: #0b0b0b; "
    "border: 1px solid var(--grid, #e2e1dc); border-radius: 6px; font-size: 11px; "
    "box-shadow: 0 2px 8px rgba(0, 0, 0, 0.18); width: max-content; "
    "max-width: 344px; }\n"
    ".moltip svg { display: block; }\n"
    ".moltip .count { color: #52514e; margin: 2px 0 4px; }\n"
    ".moltip pre { overflow: hidden; }\n"
    ".value { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; "
    "overflow-wrap: anywhere; }\n"
    "pre.value { margin: 0; font-size: 10px; line-height: 1.25; "
    "white-space: pre; overflow-wrap: normal; }\n"
    ".chart { position: relative; }\n"
    ".chart .moltip { position: absolute; z-index: 2; pointer-events: none; }\n"
    + "".join(
        f".chart:has(.b{i}:hover, .b{i}:focus) .t{i} {{ display: block; }}\n"
        for i in range(MAX_CHART_CATEGORIES)
    )
    + "details.values { margin-top: 6px; font-size: 12px; }\n"
    "details.values table { table-layout: fixed; width: 100%; }\n"
    "details.values th:first-child { width: 6em; }\n"
    "details.values pre.value { overflow-x: auto; font-size: 11px; }\n"
)


def _value_html(values: tuple[object, ...], max_lines: int | None = None) -> str:
    """A bar's values in full, one per line; a multi-line one as preformatted text.

    Args:
        values: The bar's value in each column of the group.
        max_lines: The most lines of a multi-line value to show; the rest are
            counted. None shows them all.
    """
    parts = []
    for value in values:
        text = _format(value)
        if "\n" not in text:
            parts.append(f'<div class="value">{html.escape(text)}</div>')
            continue
        lines = text.split("\n")
        if max_lines is not None and len(lines) > max_lines:
            hidden = len(lines) - max_lines
            lines = [*lines[:max_lines], f"... {hidden} more lines"]
        parts.append(f'<pre class="value">{html.escape(chr(10).join(lines))}</pre>')
    return "".join(parts)


def _panels(
    stat: ColumnStats, drawing: Callable[[object], str | None] | None
) -> dict[int, str]:
    """The hover panel for each bar whose label doesn't show it all, by position.

    A bar gets a panel when its label is cut short, its value runs over more than
    one line, or it has a drawing: the panel shows the value in full, with the
    drawing if there is one.
    """
    panels = {}
    for index, (values, count) in enumerate(stat.categories[:MAX_CHART_CATEGORIES]):
        name = " | ".join(_format(v) for v in values)
        tip = None
        if drawing is not None:
            tip = next((t for v in values if v is not None and (t := drawing(v))), None)
        if tip is None and "\n" not in name and _fit(name, _LABEL_ROOM) == name:
            continue
        panels[index] = (
            _value_html(values, MAX_TIP_LINES)
            + f'<div class="count">{count:,} reactions</div>'
            + (tip or "")
        )
    return panels


def _values_list(stat: ColumnStats) -> str:
    """Each bar's values in full, to read and copy, folded under the chart."""
    shown = stat.categories[:MAX_CHART_CATEGORIES]
    rows = "".join(
        f"<tr><td>{count:,}</td><td>{_value_html(values)}</td></tr>"
        for values, count in shown
    )
    return (
        f'<details class="values"><summary>Values shown, in full ({len(shown)})'
        "</summary><table><thead><tr><th>Reactions</th><th>Value</th></tr></thead>"
        f"<tbody>{rows}</tbody></table></details>"
    )


def _categories_svg(
    stat: ColumnStats, label: str, with_panel: Collection[int] = ()
) -> str:
    shown = stat.categories[:MAX_CHART_CATEGORIES]
    rest = (
        sum(count for _, count in stat.categories[MAX_CHART_CATEGORIES:]) + stat.other
    )
    bars = [
        (" | ".join(format_value(v) for v in values), count) for values, count in shown
    ]
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
        # A bar with a panel can be focused, so the keyboard can show it too.
        marked = f' class="b{index}" tabindex="0"' if index in with_panel else ""
        parts.append(
            f"<g{marked}><title>{html.escape(tooltip)}</title>"
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
            f"{format_value(edges[index])} to {format_value(edges[index + 1])}: "
            f"{count:,} reactions"
        )
        path = _bar_path(x, baseline - bar_height, slot - _GAP, bar_height, up=True)
        parts.append(
            f"<g><title>{html.escape(tooltip)}</title>"
            f'<path d="{path}" fill="var(--series-1, #2a78d6)"/></g>'
        )
    parts.append(_text(left, baseline + 15, format_value(edges[0])))
    parts.append(
        _text(_WIDTH - 8, baseline + 15, format_value(edges[-1]), anchor="end")
    )
    parts.append(
        _text(
            left + plot_width / 2,
            baseline + 15,
            f"median {format_value(stat.summary['median'])}",
            anchor="middle",
        )
    )
    parts.append("</svg>")
    return "".join(parts)
