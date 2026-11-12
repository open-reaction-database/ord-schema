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

"""The review report: everything a reviewer reads about one dataset, in one bundle.

:func:`build_report` reads a dataset once and gathers:

- validation findings for every reaction, grouped so that one message repeated
  across thousands of reactions is one row with a count;
- the dataset-level structural findings of :mod:`ord_schema.review.structure`,
  under the dataset's profile;
- the compound table, the flattened table, and the column statistics from
  :mod:`ord_schema.analysis`.

:func:`write_bundle` writes it as::

    <output_dir>/
        <dataset>_review.html   one self-contained page; long lists are capped
        tables/
            findings.csv        every finding, validation and structural
            compounds.csv       every compound
            flat.csv            the varying fields, one row per reaction
            constants.csv       the fields every reaction holds the same way
            column_stats.csv    statistics for each group of varying columns

The CSVs are always complete; only the page caps what it shows. Nothing here blocks
a submission: the report is for the reviewer, who decides.
"""

import csv
import dataclasses
import datetime
import html
import importlib.metadata
import os
import pathlib
import re
from collections.abc import Iterable, Iterator

from ord_schema import parquet, validations
from ord_schema.analysis import columns
from ord_schema.analysis.compounds import StructureDrawings
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.review import structure
from ord_schema.validations import Severity

# What the page shows before pointing to the CSVs.
MAX_FINDINGS_SHOWN = 200
MAX_COMPOUNDS_SHOWN = 200
MAX_CONSTANTS_SHOWN = 300
MAX_EXAMPLES = 5

# Digits in a validation message, which differ from reaction to reaction (amounts,
# times, indices) while the problem they describe is the same.
_DIGITS = re.compile(r"\d+(?:\.\d+)?")


@dataclasses.dataclass
class ValidationGroup:
    """One validation message, as repeated across a dataset's reactions.

    Attributes:
        severity: How serious the message is.
        message: The message as first seen.
        reactions: How many reactions raised it.
        examples: A few of those reactions.
    """

    severity: Severity
    message: str
    reactions: int = 0
    examples: list[str] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class ReviewReport:
    """Everything the review found about one dataset.

    Attributes:
        name: The dataset's name.
        dataset_id: Its ID, if it has one.
        description: Its description.
        profile: The profile the structural checks ran under.
        validation: Validation messages, grouped, most serious and most common first.
        structure: The structural checks' findings, and which checks ran.
        evidence: The summary, compound table and flattened table they read.
        column_stats: Statistics for each group of varying columns.
        generated: When the report was made, in UTC.
        ord_schema_version: The ord-schema that made it.
    """

    name: str
    dataset_id: str
    description: str
    profile: structure.ValidationProfile
    validation: list[ValidationGroup]
    structure: structure.StructureReport
    evidence: structure.DatasetEvidence
    column_stats: list[columns.ColumnStats]
    generated: str
    ord_schema_version: str

    @property
    def reaction_count(self) -> int:
        """How many reactions the dataset holds."""
        return self.evidence.table.reaction_count

    def severity_counts(self) -> dict[Severity, int]:
        """Findings at each severity: structural findings, and validation messages."""
        counts = dict.fromkeys(Severity, 0)
        for finding in self.structure.findings:
            counts[finding.severity] += 1
        for group in self.validation:
            counts[group.severity] += 1
        return counts


class _Validator:
    """Validates reactions one at a time, grouping their messages."""

    def __init__(self, options: validations.ValidationOptions) -> None:
        self.options = options
        self.groups: dict[tuple[Severity, str], ValidationGroup] = {}

    def watch(
        self, reactions: Iterable[reaction_pb2.Reaction]
    ) -> Iterator[reaction_pb2.Reaction]:
        """Yields ``reactions`` unchanged, validating each on the way past."""
        for position, reaction in enumerate(reactions):
            # validate_message may edit what it validates, so it gets a copy.
            copy = reaction_pb2.Reaction()
            copy.CopyFrom(reaction)
            output = validations.validate_message(
                copy, raise_on_error=False, options=self.options
            )
            label = reaction.reaction_id or f"reaction #{position + 1}"
            for severity, messages in (
                (Severity.ERROR, output.errors),
                (Severity.WARNING, output.warnings),
                (Severity.INFO, output.info),
            ):
                for message in dict.fromkeys(messages):
                    key = (severity, _DIGITS.sub("#", message))
                    group = self.groups.get(key)
                    if group is None:
                        group = self.groups[key] = ValidationGroup(severity, message)
                    group.reactions += 1
                    if len(group.examples) < MAX_EXAMPLES:
                        group.examples.append(label)
            yield reaction

    def sorted_groups(self) -> list[ValidationGroup]:
        return sorted(
            self.groups.values(), key=lambda g: (-g.severity, -g.reactions, g.message)
        )


def build_report(
    dataset: dataset_pb2.Dataset | parquet.DatasetView,
    profile: structure.ValidationProfile = structure.ValidationProfile.UNIFORM,
    *,
    options: validations.ValidationOptions | None = None,
) -> ReviewReport:
    """Reviews a dataset, reading its reactions once.

    Args:
        dataset: The dataset: a ``Dataset`` message, or a Parquet ``DatasetView``,
            which is streamed rather than loaded whole.
        profile: The dataset's profile; ``uniform`` unless it declares otherwise.
        options: Options for validating each reaction; the defaults if omitted.

    Returns:
        The report.
    """
    validator = _Validator(options or validations.ValidationOptions())
    evidence = structure.gather(validator.watch(dataset.reactions))
    return ReviewReport(
        name=dataset.name,
        dataset_id=dataset.dataset_id,
        description=dataset.description,
        profile=profile,
        validation=validator.sorted_groups(),
        structure=structure.run_checks(evidence, profile),
        evidence=evidence,
        column_stats=columns.column_statistics(evidence.table),
        generated=datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M UTC"),
        ord_schema_version=importlib.metadata.version("ord-schema"),
    )


def _write_findings_csv(report: ReviewReport, path: pathlib.Path) -> None:
    fields = [
        "source",
        "check",
        "severity",
        "message",
        "fields",
        "reactions",
        "examples",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for finding in report.structure.findings:
            writer.writerow(
                {
                    "source": "structure",
                    "check": finding.check,
                    "severity": finding.severity.name,
                    "message": finding.message,
                    "fields": " | ".join(finding.fields),
                    "reactions": "" if finding.reactions is None else finding.reactions,
                    "examples": " | ".join(finding.examples),
                }
            )
        for group in report.validation:
            writer.writerow(
                {
                    "source": "validation",
                    "check": "",
                    "severity": group.severity.name,
                    "message": group.message,
                    "fields": "",
                    "reactions": group.reactions,
                    "examples": " | ".join(group.examples),
                }
            )


def _file_stem(report: ReviewReport) -> str:
    """A filename for the page, from the dataset's ID or name."""
    stem = report.dataset_id or report.name or "dataset"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_") or "dataset"


def write_bundle(
    report: ReviewReport, output_dir: str | os.PathLike[str]
) -> pathlib.Path:
    """Writes the report's page and tables into ``output_dir``.

    Args:
        report: The report.
        output_dir: Where to write; created if it does not exist.

    Returns:
        The path of the HTML page.
    """
    root = pathlib.Path(output_dir)
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    evidence = report.evidence
    _write_findings_csv(report, tables / "findings.csv")
    evidence.summary.compounds.write_csv(tables / "compounds.csv")
    evidence.table.write_csv(tables / "flat.csv")
    evidence.table.write_constants_csv(tables / "constants.csv")
    columns.write_csv(report.column_stats, tables / "column_stats.csv")
    page = root / f"{_file_stem(report)}_review.html"
    page.write_text(render_html(report), encoding="utf-8")
    return page


_STYLE = """
:root {
  --surface: #fcfcfb; --card: #ffffff; --text-primary: #0b0b0b;
  --text-secondary: #52514e; --grid: #e2e1dc; --series-1: #2a78d6;
  --error: #d03b3b; --error-bg: #fbeaea; --warning: #8a5a00; --warning-bg: #fdf3dc;
  --info: #2a5f9e; --info-bg: #e8f0fb; --flag-bg: #fbeaea;
  color-scheme: light;
}
@media (prefers-color-scheme: dark) {
  :root {
    --surface: #1a1a19; --card: #222221; --text-primary: #f0efec;
    --text-secondary: #c3c2b7; --grid: #383835; --series-1: #3987e5;
    --error: #f08080; --error-bg: #3a1f1f; --warning: #f0c060; --warning-bg: #3a3020;
    --info: #8fb8ef; --info-bg: #1f2a3a; --flag-bg: #3a1f1f;
    color-scheme: dark;
  }
}
* { box-sizing: border-box; }
body { margin: 0; padding: 24px; background: var(--surface);
  color: var(--text-primary); font: 14px/1.5 system-ui, -apple-system, sans-serif; }
main { max-width: 1280px; margin: 0 auto; display: grid; gap: 32px; }
h1 { font-size: 26px; margin: 0; text-wrap: balance; }
h2 { font-size: 20px; margin: 0 0 8px; padding-bottom: 4px;
  border-bottom: 2px solid var(--grid); }
h3 { font-size: 16px; margin: 16px 0 6px; }
p { margin: 0 0 8px; max-width: 78ch; }
.meta, .muted { color: var(--text-secondary); }
code { font-size: 12px; overflow-wrap: break-word; }
td code { overflow-wrap: anywhere; }
td.check code { white-space: nowrap; overflow-wrap: normal; }
.examples { display: flex; flex-wrap: wrap; gap: 2px 8px; min-width: 150px; }
.tiles { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 12px; }
.tile { background: var(--card); border: 1px solid var(--grid); border-radius: 8px;
  padding: 8px 14px; min-width: 120px; }
.tile b { display: block; font-size: 22px; }
.tile span { color: var(--text-secondary); font-size: 12px; }
.notice { background: var(--card); border: 1px solid var(--grid); border-radius: 8px;
  padding: 10px 14px; color: var(--text-secondary); }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; background: var(--card); }
th, td { border: 1px solid var(--grid); padding: 6px 8px; vertical-align: top;
  text-align: left; }
th { font-size: 12px; color: var(--text-secondary); }
td ul { margin: 0; padding-left: 16px; }
.pill { display: inline-block; font-size: 11px; font-weight: 600; border-radius: 999px;
  padding: 1px 8px; white-space: nowrap; }
.pill.ERROR { color: var(--error); background: var(--error-bg); }
.pill.WARNING { color: var(--warning); background: var(--warning-bg); }
.pill.INFO { color: var(--info); background: var(--info-bg); }
.grid { display: grid; gap: 14px;
  grid-template-columns: repeat(auto-fill, minmax(390px, 1fr)); }
figure { margin: 0; padding: 12px; background: var(--card);
  border: 1px solid var(--grid); border-radius: 8px; overflow-x: auto; }
figcaption { font-size: 12px; margin-bottom: 8px; display: grid; gap: 2px; }
td.structure { width: 190px; background: #ffffff; }
td.structure svg { display: block; }
tr.flagged td { background: var(--flag-bg); }
tr.flagged td.structure { background: #ffffff; }
.conflict { color: var(--error); font-weight: 600; }
.note { color: var(--warning); }
details summary { cursor: pointer; color: var(--text-secondary); }
@media print {
  body { padding: 0; font-size: 10pt; }
  section { break-before: page; }
  tr, figure { break-inside: avoid; }
  thead { display: table-header-group; }
}
"""


def _e(value: object) -> str:
    return html.escape(str(value))


def _pill(severity: Severity) -> str:
    return f'<span class="pill {severity.name}">{severity.name}</span>'


def _list(values: Iterable[str], drawings: StructureDrawings | None = None) -> str:
    """Example reaction IDs or values, kept whole and wrapped between items.

    With ``drawings``, a value that is a recorded compound's SMILES shows its
    structure on hover.
    """
    if drawings is None:
        items = "".join(f"<code>{_e(v)}</code>" for v in values)
    else:
        items = "".join(drawings.code(str(v)) for v in values)
    return f'<div class="examples">{items}</div>' if items else ""


def _header(report: ReviewReport) -> str:
    counts = report.severity_counts()
    compounds = report.evidence.summary.compounds
    tiles = [
        (f"{report.reaction_count:,}", "reactions"),
        (f"{counts[Severity.ERROR]:,}", "errors"),
        (f"{counts[Severity.WARNING]:,}", "warnings"),
        (f"{counts[Severity.INFO]:,}", "notes"),
        (f"{len(compounds):,}", "compounds"),
        (f"{len(report.column_stats):,}", "varying fields"),
    ]
    meta = [
        f"profile <b>{report.profile.value}</b>",
        f"generated {_e(report.generated)}",
        f"ord-schema {_e(report.ord_schema_version)}",
    ]
    if report.dataset_id:
        meta.insert(0, f"<code>{_e(report.dataset_id)}</code>")
    return (
        "<header>"
        f"<h1>Review: {_e(report.name or 'unnamed dataset')}</h1>"
        f'<p class="meta">{" · ".join(meta)}</p>'
        + (f"<p>{_e(report.description)}</p>" if report.description else "")
        + '<div class="tiles">'
        + "".join(
            f'<div class="tile"><b>{v}</b><span>{k}</span></div>' for v, k in tiles
        )
        + "</div></header>"
        '<p class="notice">This report is for the reviewer. Nothing in it blocks a '
        "submission; the reviewer decides what needs changing. Every table is complete "
        "in <code>tables/</code>; this page caps long lists. Not yet included: "
        "completeness scores (#1054) and exemplar reactions (#1052).</p>"
    )


def _findings_section(report: ReviewReport, drawings: StructureDrawings) -> str:
    findings = report.structure.findings
    parts = ["<section><h2>Dataset-level findings</h2>"]
    if findings:
        parts.append(
            '<div class="scroll"><table><thead><tr><th>Severity</th><th>Check</th>'
            "<th>Finding</th><th>Reactions</th><th>Examples</th></tr></thead><tbody>"
        )
        parts.extend(
            f"<tr><td>{_pill(f.severity)}</td>"
            f'<td class="check"><code>{_e(f.check)}</code></td>'
            f"<td>{_e(f.message)}"
            + (
                f'<div class="muted">{_e(", ".join(f.fields[:3]))}'
                + (f" and {len(f.fields) - 3} more" if len(f.fields) > 3 else "")
                + "</div>"
                if f.fields
                else ""
            )
            + f"</td><td>{'' if f.reactions is None else f'{f.reactions:,}'}</td>"
            f"<td>{_list(f.examples, drawings)}</td></tr>"
            for f in findings[:MAX_FINDINGS_SHOWN]
        )
        parts.append("</tbody></table></div>")
        if len(findings) > MAX_FINDINGS_SHOWN:
            parts.append(
                f'<p class="muted">{len(findings) - MAX_FINDINGS_SHOWN:,} more in '
                "tables/findings.csv.</p>"
            )
    else:
        parts.append("<p>None of the checks that ran found anything.</p>")
    applied = ", ".join(f"<code>{_e(c)}</code>" for c in report.structure.applied)
    parts.append(f'<p class="muted">Checks run: {applied}.</p>')
    if report.structure.not_applied:
        skipped = ", ".join(
            f"<code>{_e(c)}</code>" for c in report.structure.not_applied
        )
        parts.append(
            f'<p class="muted">Not applied under the <b>{report.profile.value}</b> '
            f"profile: {skipped}.</p>"
        )
    parts.append("</section>")
    return "".join(parts)


def _validation_section(report: ReviewReport) -> str:
    groups = report.validation
    parts = [
        "<section><h2>Reaction validation</h2>",
        '<p class="muted">Each message from validating the reactions one at a time, '
        "with how many reactions raised it. Messages that differ only in numbers are "
        "counted together; the first is shown.</p>",
    ]
    if not groups:
        parts.append("<p>Every reaction validates with nothing to report.</p>")
    else:
        parts.append(
            '<div class="scroll"><table><thead><tr><th>Severity</th><th>Message</th>'
            "<th>Reactions</th><th>Examples</th></tr></thead><tbody>"
        )
        parts.extend(
            f"<tr><td>{_pill(g.severity)}</td><td>{_e(g.message)}</td>"
            f"<td>{g.reactions:,}</td><td>{_list(g.examples)}</td></tr>"
            for g in groups[:MAX_FINDINGS_SHOWN]
        )
        parts.append("</tbody></table></div>")
        if len(groups) > MAX_FINDINGS_SHOWN:
            parts.append(
                f'<p class="muted">{len(groups) - MAX_FINDINGS_SHOWN:,} more in '
                "tables/findings.csv.</p>"
            )
    parts.append("</section>")
    return "".join(parts)


def _shapes_section(report: ReviewReport) -> str:
    summary = report.evidence.summary
    shapes = summary.ordered_shapes()
    parts = [
        "<section><h2>How reactions are recorded</h2>",
        f"<p>{len(shapes):,} distinct shape{'s' if len(shapes) != 1 else ''}: sets of "
        "fields a reaction fills. In a uniform dataset one shape is expected; more "
        "can be controls, or fields some reactions miss.</p>",
    ]
    if len(shapes) > 1:
        common = shapes[0][0]
        rows = []
        for shape, count in shapes[:20]:
            extra = sorted(shape - common)
            missing = sorted(common - shape)
            examples = summary.shape_examples.get(shape, [])
            diff = []
            if extra:
                diff.append(f"adds {len(extra)}: " + ", ".join(extra[:4]))
            if missing:
                diff.append(f"lacks {len(missing)}: " + ", ".join(missing[:4]))
            rows.append(
                f"<tr><td>{count:,}</td><td>{_e('; '.join(diff) or 'most common')}</td>"
                f"<td>{_list(examples)}</td></tr>"
            )
        parts.append(
            '<div class="scroll"><table><thead><tr><th>Reactions</th>'
            "<th>Compared with the most common shape</th><th>Examples</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table></div>"
        )
    parts.append("</section>")
    return "".join(parts)


def _figure(stat: columns.ColumnStats, drawings: StructureDrawings) -> str:
    caption = [f"<code>{_e(stat.column)}</code>"]
    if len(stat.columns) > 1:
        caption.append(
            f'<span class="muted">with {_e(", ".join(stat.columns[1:]))}</span>'
        )
    detail = (
        f"{stat.kind}, {stat.filled:,} filled, {stat.empty:,} empty, "
        f"{stat.distinct:,} distinct"
    )
    caption.append(f'<span class="muted">{detail}</span>')
    body = "".join(f"<div>{c}</div>" for c in caption)
    chart = stat.to_html(drawing=drawings.tip)
    return f"<figure><figcaption>{body}</figcaption>{chart}</figure>"


def _figures(stats: list[columns.ColumnStats], drawings: StructureDrawings) -> str:
    return "".join(_figure(stat, drawings) for stat in stats)


def _columns_section(report: ReviewReport, drawings: StructureDrawings) -> str:
    table = report.evidence.table
    parts = ["<section><h2>What the dataset varies</h2>"]
    if report.profile is structure.ValidationProfile.VARIED:
        parts.append(
            "<p>Under the <b>varied</b> profile the reactions are not expected to "
            "share a template, so a table of what varies between them says little. The "
            "full table is in <code>tables/flat.csv</code>.</p></section>"
        )
        return "".join(parts)
    experiment = [s for s in report.column_stats if s.section == "experiment"]
    provenance = [s for s in report.column_stats if s.section == "provenance"]
    groups = len(report.column_stats)
    parts.append(
        f"<p>{len(table.varying.columns):,} varying columns in {groups:,} "
        "groups; columns that always change together, such as a compound's name and "
        "SMILES, are one group. Hover over a bar for its full value.</p>"
    )
    if experiment:
        parts.append(f'<div class="grid">{_figures(experiment, drawings)}</div>')
    else:
        parts.append("<p>No experimental field varies.</p>")
    parts.append("<h3>Provenance</h3>")
    if provenance:
        parts.append(
            '<p class="muted">Who recorded the data and when. Dates are compared by '
            f"{_e(table.date_resolution or 'exact value')}.</p>"
            f'<div class="grid">{_figures(provenance, drawings)}</div>'
        )
    else:
        parts.append('<p class="muted">The same for every reaction; see below.</p>')
    constants = table.constants_frame()
    shown = constants.head(MAX_CONSTANTS_SHOWN)
    rows = "".join(
        f"<tr><td><code>{_e(field)}</code></td><td>{_e(value)}</td>"
        f"<td>{_e(section)}</td></tr>"
        for field, value, section in zip(
            shown["field"], shown["value"], shown["section"], strict=True
        )
    )
    parts.append(
        f"<details><summary>{len(constants):,} fields every reaction holds the same "
        'way</summary><div class="scroll"><table><thead><tr><th>Field</th>'
        f"<th>Value</th><th>Section</th></tr></thead><tbody>{rows}</tbody></table></div>"
        + (
            f'<p class="muted">{len(constants) - MAX_CONSTANTS_SHOWN:,} more in '
            "tables/constants.csv.</p>"
            if len(constants) > MAX_CONSTANTS_SHOWN
            else ""
        )
        + "</details></section>"
    )
    return "".join(parts)


def _compounds_section(report: ReviewReport, drawings: StructureDrawings) -> str:
    compounds = report.evidence.summary.compounds
    return (
        "<section><h2>Compounds</h2>"
        f"<p>{len(compounds):,} distinct compounds, grouped by role. Values that "
        "disagree are listed, not resolved; rows with something to check are "
        "shaded.</p>"
        + compounds.to_html_sections(
            heading="h3", max_rows=MAX_COMPOUNDS_SHOWN, drawings=drawings
        )
        + "</section>"
    )


def render_html(report: ReviewReport) -> str:
    """Returns the report as one self-contained HTML page."""
    title = f"Review: {report.name or 'unnamed dataset'}"
    # Each structure shown on hover is drawn once, in defs, after the sections that
    # use it have been written.
    drawings = StructureDrawings(report.evidence.summary.compounds)
    body = (
        _header(report)
        + _findings_section(report, drawings)
        + _validation_section(report)
        + _shapes_section(report)
        + _columns_section(report, drawings)
        + _compounds_section(report, drawings)
    )
    return (
        # The marker in Go's form, which tools such as addlicense recognize, so a
        # report kept in a repository is not mistaken for source missing a license.
        "<!DOCTYPE html>\n<!--\n// Code generated by ord_schema.review.report. "
        "DO NOT EDIT.\n-->\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{_e(title)}</title><style>{_STYLE}"
        f"{columns.TOOLTIP_CSS}{StructureDrawings.CSS}</style></head><body>"
        + drawings.defs()
        + "<main>"
        + body
        + "</main></body></html>\n"
    )


__all__ = [
    "ReviewReport",
    "ValidationGroup",
    "build_report",
    "render_html",
    "write_bundle",
]
