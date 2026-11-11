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
            findings.csv        every finding: validation, structural, completeness
            compounds.csv       every compound
            flat.csv            the varying fields, one row per reaction
            constants.csv       the fields every reaction holds the same way
            column_stats.csv    statistics for each group of varying columns
            exemplars.txtpb     the reactions picked for the reviewer to read
            completeness.csv    each reaction's completeness score, and what it lacks
        packet/
            review_packet.md    the review for an AI reviewer, sized to fit its context
            review.json         the same, structured

The CSVs are always complete; only the page caps what it shows. Nothing here blocks
a submission: the report is for the reviewer, who decides.
"""

import base64
import csv
import dataclasses
import datetime
import html
import importlib.metadata
import os
import pathlib
import re
from collections.abc import Iterable, Iterator

from google.protobuf import text_format

from ord_schema import parquet, validations
from ord_schema.analysis import columns, flat
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.review import chemistry, completeness, exemplars, packet, structure
from ord_schema.validations import Severity

# What the page shows before pointing to the CSVs.
MAX_FINDINGS_SHOWN = 200
MAX_COMPOUNDS_SHOWN = 200
MAX_CONSTANTS_SHOWN = 300
MAX_EXAMPLES = 5
# Characters of an exemplar's text shown on the page; the pbtxt holds the rest.
MAX_EXEMPLAR_TEXT = 6000

# Where a reviewer can import an exemplar to see it rendered, until the report can
# link to a viewer page instead (#1183).
EDITOR_URL = "https://editor.open-reaction-database.org"

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
        exemplars: The reactions picked for the reviewer to read, and why.
        completeness: The reactions' completeness scores, rolled up.
        chemistry: What the chemistry sanity checks found, rolled up.
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
    exemplars: list[exemplars.Exemplar]
    completeness: completeness.CompletenessRollup
    chemistry: list[structure.Finding]
    generated: str
    ord_schema_version: str

    @property
    def reaction_count(self) -> int:
        """How many reactions the dataset holds."""
        return self.evidence.table.reaction_count

    def severity_counts(self) -> dict[Severity, int]:
        """Findings at each severity: structural findings, and validation messages."""
        counts = dict.fromkeys(Severity, 0)
        for finding in [*self.structure.findings, *self.chemistry]:
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
    scores = completeness.CompletenessRollup(profile.value)
    sanity = chemistry.ChemistryRollup()
    evidence = structure.gather(
        sanity.watch(scores.watch(validator.watch(dataset.reactions)))
    )
    stats = columns.column_statistics(evidence.table)
    return ReviewReport(
        name=dataset.name,
        dataset_id=dataset.dataset_id,
        description=dataset.description,
        profile=profile,
        validation=validator.sorted_groups(),
        structure=structure.run_checks(evidence, profile),
        evidence=evidence,
        column_stats=stats,
        # Fetched back by position, which a Dataset and a DatasetView both allow.
        exemplars=exemplars.pick(evidence.table, stats, dataset.reactions, profile),
        completeness=scores,
        chemistry=sanity.findings(),
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
        for source, finding in [
            *(("structure", f) for f in report.structure.findings),
            *(("chemistry", f) for f in report.chemistry),
        ]:
            writer.writerow(
                {
                    "source": source,
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
        for rollup in report.completeness.lacking():
            writer.writerow(
                {
                    "source": "completeness",
                    "check": rollup.item.id,
                    "severity": Severity.INFO.name,
                    "message": rollup.suggestion(),
                    "fields": " | ".join(rollup.item.fields),
                    "reactions": rollup.lacking,
                    "examples": " | ".join(rollup.examples),
                }
            )


def _file_stem(report: ReviewReport) -> str:
    """A filename for the page, from the dataset's ID or name."""
    stem = report.dataset_id or report.name or "dataset"
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_") or "dataset"


def write_bundle(
    report: ReviewReport,
    output_dir: str | os.PathLike[str],
    *,
    packet_budget_tokens: int = packet.DEFAULT_BUDGET_TOKENS,
) -> pathlib.Path:
    """Writes the report's page and tables into ``output_dir``.

    Args:
        report: The report.
        output_dir: Where to write; created if it does not exist.
        packet_budget_tokens: The most estimated tokens each file of the AI review
            packet may hold; see :mod:`ord_schema.review.packet`.

    Returns:
        The path of the HTML page.
    """
    root = pathlib.Path(output_dir)
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    evidence = report.evidence
    _write_findings_csv(report, tables / "findings.csv")
    report.completeness.write_csv(tables / "completeness.csv")
    evidence.summary.compounds.write_csv(tables / "compounds.csv")
    evidence.table.write_csv(tables / "flat.csv")
    evidence.table.write_constants_csv(tables / "constants.csv")
    columns.write_csv(report.column_stats, tables / "column_stats.csv")
    exemplars.write_txtpb(
        report.exemplars, tables / "exemplars.txtpb", name=report.name
    )
    packet.write_packet(report, root, packet_budget_tokens)
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
body { margin: 0; padding: 24px 16px; background: var(--surface);
  color: var(--text-primary); font: 14px/1.5 system-ui, -apple-system, sans-serif; }
main { max-width: 1280px; margin: 0 auto; display: grid; gap: 32px; }
/* A grid item is as wide as its widest content unless told otherwise; this lets
   long tables and reaction text scroll in place instead of widening the page. */
main > * { min-width: 0; }
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
  grid-template-columns: repeat(auto-fill, minmax(min(390px, 100%), 1fr)); }
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
a { color: var(--series-1); }
a.download { font-weight: 600; }
details pre { background: var(--card); border: 1px solid var(--grid);
  border-radius: 6px; padding: 10px; font-size: 12px; overflow-x: auto;
  max-height: 480px; }
section > details { margin: 6px 0; }
details summary, details li { overflow-wrap: anywhere; }
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


def _list(values: Iterable[str]) -> str:
    """Example reaction IDs or values, kept whole and wrapped between items."""
    items = "".join(f"<code>{_e(v)}</code>" for v in values)
    return f'<div class="examples">{items}</div>' if items else ""


def _header(report: ReviewReport) -> str:
    counts = report.severity_counts()
    compounds = report.evidence.summary.compounds
    tiles = [
        (f"{report.reaction_count:,}", "reactions"),
        (f"{counts[Severity.ERROR]:,}", "errors"),
        (f"{counts[Severity.WARNING]:,}", "warnings"),
        (f"{counts[Severity.INFO]:,}", "notes"),
        (_percent(report.completeness.mean_overall()), "completeness"),
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
        "in <code>tables/</code>; this page caps long lists.</p>"
    )


def _findings_table(findings: list[structure.Finding]) -> str:
    """Findings as a table, capped, pointing to findings.csv for the rest."""
    parts = [
        '<div class="scroll"><table><thead><tr><th>Severity</th><th>Check</th>'
        "<th>Finding</th><th>Reactions</th><th>Examples</th></tr></thead><tbody>"
    ]
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
        f"<td>{_list(f.examples)}</td></tr>"
        for f in findings[:MAX_FINDINGS_SHOWN]
    )
    parts.append("</tbody></table></div>")
    if len(findings) > MAX_FINDINGS_SHOWN:
        parts.append(
            f'<p class="muted">{len(findings) - MAX_FINDINGS_SHOWN:,} more in '
            "tables/findings.csv.</p>"
        )
    return "".join(parts)


def _findings_section(report: ReviewReport) -> str:
    findings = report.structure.findings
    parts = ["<section><h2>Dataset-level findings</h2>"]
    if findings:
        parts.append(_findings_table(findings))
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


def _chemistry_section(report: ReviewReport) -> str:
    checks = ", ".join(f"<code>{_e(c)}</code>" for c in chemistry.CHECKS)
    return (
        "<section><h2>Chemistry checks</h2>"
        "<p>Whether each reaction makes chemical sense: elements conserved, a stated "
        "reaction SMILES matching the recorded compounds, no stray radicals or "
        "unbalanced charges, stereochemistry where a name implies it, and plausible "
        "yields, temperatures and times. Messages that differ only in numbers are "
        "counted together.</p>"
        + (
            _findings_table(report.chemistry)
            if report.chemistry
            else "<p>No reaction raised a chemistry check.</p>"
        )
        + f'<p class="muted">Checks run: {checks}.</p></section>'
    )


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


def _percent(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.0%}"


def _completeness_section(report: ReviewReport) -> str:
    rollup = report.completeness
    levels = rollup.mean_levels()
    parts = [
        "<section><h2>Completeness</h2>",
        "<p>How completely the reactions record details that can almost always be "
        "recorded, scored against rubric version "
        f"{_e(rollup.rubric.version)} under the <b>{_e(rollup.profile)}</b> profile. "
        "Advisory: a low score never blocks a submission. Each reaction's score is in "
        "<code>tables/completeness.csv</code>.</p>",
        '<div class="tiles">',
        f'<div class="tile"><b>{_percent(rollup.mean_overall())}</b>'
        "<span>overall, mean of reactions</span></div>",
    ]
    parts.extend(
        f'<div class="tile"><b>{_percent(value)}</b><span>{_e(level)}</span></div>'
        for level, value in levels.items()
    )
    parts.append("</div>")
    lacking = rollup.lacking()
    if lacking:
        parts.append(
            '<div class="scroll"><table><thead><tr><th>Suggestion</th><th>Level</th>'
            "<th>Reactions lacking it</th><th>Examples</th></tr></thead><tbody>"
        )
        parts.extend(
            f"<tr><td>{_e(r.item.suggestion[0].upper() + r.item.suggestion[1:])}"
            f'<div class="muted">{_e(", ".join(r.item.fields))}</div></td>'
            f"<td>{_e(r.item.level)}</td>"
            f"<td>{r.lacking:,} of {r.applicable:,}</td>"
            f"<td>{_list(r.examples)}</td></tr>"
            for r in lacking
        )
        parts.append("</tbody></table></div>")
    else:
        parts.append("<p>Every reaction is complete on every item that applies.</p>")
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


def _download(text: str, filename: str, label: str) -> str:
    """A link that saves ``text`` as ``filename``, carried in the page itself."""
    data = base64.b64encode(text.encode()).decode()
    return (
        f'<a class="download" download="{_e(filename)}" '
        f'href="data:text/plain;charset=utf-8;base64,{data}">{_e(label)}</a>'
    )


def _filename(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("_") or "reaction"


def _exemplars_section(report: ReviewReport) -> str:
    picked = report.exemplars
    parts = [
        "<section><h2>Exemplar reactions</h2>",
        f"<p>{len(picked):,} reaction{'s' if len(picked) != 1 else ''} to read in "
        "full: one of each shape, the extremes of each measured value and, in a "
        "uniform dataset, inputs no other reaction uses. All are in "
        "<code>tables/exemplars.txtpb</code>.</p>",
    ]
    # The page carries each exemplar without its binary payloads: an image in every
    # reaction would make the page tens of megabytes. tables/exemplars.txtpb has them.
    stripped = [exemplars.without_binary(e.reaction) for e in picked]
    picked = [
        dataclasses.replace(e, reaction=reaction)
        for e, (reaction, _) in zip(picked, stripped, strict=True)
    ]
    omitted = sum(size for _, size in stripped)
    if picked:
        everything = text_format.MessageToString(
            exemplars.to_dataset(picked, name=report.name)
        )
        parts.append(
            '<p class="notice">To see an exemplar rendered, download it and import it '
            f'into the <a href="{EDITOR_URL}">ORD reaction editor</a>: one reaction '
            "with <b>Create Reaction from File</b>, or all of them as a dataset with "
            "<b>Create Dataset from File</b>. Imported reactions are saved to your "
            "editor account, so import only data you may share there."
            + (
                f" Images and other binary attachments ({omitted:,} bytes in all) are "
                "left out of these downloads and the text below; "
                "<code>tables/exemplars.txtpb</code> holds them."
                if omitted
                else ""
            )
            + "<br>"
            + _download(
                everything,
                f"{_file_stem(report)}_exemplars.txtpb",
                f"Download all {len(picked):,} exemplars (.txtpb)",
            )
            + "</p>"
        )
    for exemplar in picked:
        full = text_format.MessageToString(exemplar.reaction)
        text = full
        if len(text) > MAX_EXEMPLAR_TEXT:
            text = (
                text[:MAX_EXEMPLAR_TEXT]
                + "\n# ... continued in the download and tables/exemplars.txtpb\n"
            )
        reasons = "".join(f"<li>{_e(reason)}</li>" for reason in exemplar.reasons)
        parts.append(
            f"<details><summary><b>{_e(exemplar.label)}</b>: "
            f"{_e(exemplar.reasons[0])}"
            + (
                f" (+{len(exemplar.reasons) - 1} more "
                f"reason{'s' if len(exemplar.reasons) > 2 else ''})"
                if len(exemplar.reasons) > 1
                else ""
            )
            + f"</summary><ul>{reasons}</ul><p>"
            + _download(
                full,
                f"{_filename(exemplar.label)}.txtpb",
                "Download this reaction (.txtpb)",
            )
            + f"</p><pre>{_e(text)}</pre></details>"
        )
    parts.append("</section>")
    return "".join(parts)


def _figure(stat: columns.ColumnStats) -> str:
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
    return f"<figure><figcaption>{body}</figcaption>{stat.to_svg()}</figure>"


def _columns_section(report: ReviewReport) -> str:
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
        parts.append(f'<div class="grid">{"".join(map(_figure, experiment))}</div>')
    else:
        parts.append("<p>No experimental field varies.</p>")
    parts.append("<h3>Provenance</h3>")
    if provenance:
        parts.append(
            '<p class="muted">Who recorded the data and when. Dates are compared by '
            f"{_e(table.date_resolution or 'exact value')}.</p>"
            f'<div class="grid">{"".join(map(_figure, provenance))}</div>'
        )
    else:
        parts.append('<p class="muted">The same for every reaction; see below.</p>')
    constants = table.constants_frame()
    shown = constants.head(MAX_CONSTANTS_SHOWN)
    rows = "".join(
        f"<tr><td><code>{_e(field)}</code></td><td>{_e(flat.display_value(value))}</td>"
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


def _compounds_section(report: ReviewReport) -> str:
    compounds = report.evidence.summary.compounds
    return (
        "<section><h2>Compounds</h2>"
        f"<p>{len(compounds):,} distinct compounds, grouped by role. Values that "
        "disagree are listed, not resolved; rows with something to check are "
        "shaded.</p>"
        '<div class="scroll">'
        + compounds.to_html_sections(heading="h3", max_rows=MAX_COMPOUNDS_SHOWN)
        + "</div></section>"
    )


def render_html(report: ReviewReport) -> str:
    """Returns the report as one self-contained HTML page."""
    title = f"Review: {report.name or 'unnamed dataset'}"
    return (
        # The marker in Go's form, which tools such as addlicense recognize, so a
        # report kept in a repository is not mistaken for source missing a license.
        "<!DOCTYPE html>\n<!--\n// Code generated by ord_schema.review.report. "
        "DO NOT EDIT.\n-->\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{_e(title)}</title><style>{_STYLE}</style></head><body><main>"
        + _header(report)
        + _findings_section(report)
        + _chemistry_section(report)
        + _validation_section(report)
        + _completeness_section(report)
        + _shapes_section(report)
        + _exemplars_section(report)
        + _columns_section(report)
        + _compounds_section(report)
        + "</main></body></html>\n"
    )


__all__ = [
    "ReviewReport",
    "ValidationGroup",
    "build_report",
    "render_html",
    "write_bundle",
]
