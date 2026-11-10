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

"""A table of every compound in a dataset, for a reviewer to read.

Each row is one compound, identified by its canonical SMILES, or by its name where it
has no readable structure. A row lists the roles the compound is used in, every
identifier recorded for it, how many reactions use it, and anything that disagrees:
two CAS numbers for one structure, one name used for two structures, identifiers in
one compound that describe different molecules, or a structure that will not parse.
Conflicts are listed, never resolved; deciding which value is right is the
reviewer's job.

Compounds are collected from reaction inputs, workup inputs and outcomes. A crude
component refers to another reaction's product rather than naming a compound, so it
appears as its own row, keyed by that reaction's ID.

The table is built one reaction at a time and two tables merge, so a large dataset
can be read as a stream, or in parallel by row group::

    table = CompoundTable()
    for reaction in datasets.load_dataset("dataset.parquet").reactions:
        table.add_reaction(reaction)
    table.write_csv("compounds.csv")
    table.write_html("compounds.html", title="My dataset")
"""

import copy
import csv
import dataclasses
import functools
import html
import os
import pathlib
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

_IdentifierType = reaction_pb2.CompoundIdentifier.CompoundIdentifierType
_RoleType = reaction_pb2.ReactionRole.ReactionRoleType

# Where in a reaction a compound was recorded.
INPUT = "input"
WORKUP = "workup"
OUTCOME = "outcome"

# The role group a crude component is listed under. Not a ReactionRoleType: a crude
# component carries no role, since it stands for another reaction's whole product.
CRUDE = "CRUDE"
# Key prefix for a crude component's row; the ID after it is the upstream reaction's.
_CRUDE_PREFIX = "crude:"

# Reaction IDs kept per row, so a reviewer can find an example without the table
# growing with the dataset.
_EXAMPLES_KEPT = 3

# Identifier types whose values are too long to read in a table cell.
_BULKY_TYPES = frozenset({_IdentifierType.MOLBLOCK, _IdentifierType.XYZ})

# Role groups in the order a reader expects them: what goes in, then what comes out.
_ROLE_ORDER = (
    "REACTANT",
    "REAGENT",
    "CATALYST",
    "SOLVENT",
    "INTERNAL_STANDARD",
    "AUTHENTIC_STANDARD",
    "WORKUP",
    CRUDE,
    "PRODUCT",
    "BYPRODUCT",
    "SIDE_PRODUCT",
    "UNSPECIFIED",
)


def _role_name(role: int) -> str:
    """Returns the name of a ReactionRoleType value, from the proto descriptor."""
    return _RoleType.Name(role)


def _identifier_name(identifier_type: int) -> str:
    """Returns the name of a CompoundIdentifierType value, from the proto descriptor."""
    return _IdentifierType.Name(identifier_type)


def _normalize_name(name: str) -> str:
    """Returns ``name`` in the form two spellings of the same name share."""
    return re.sub(r"\s+", " ", name).strip().casefold()


@dataclasses.dataclass
class CompoundRecord:
    """Everything recorded about one compound across a dataset.

    Attributes:
        key: What identifies the compound: ``smiles:<canonical SMILES>``,
            ``name:<normalized name>``, ``crude:<reaction ID>``, or for a compound
            with neither structure nor name, its first identifier.
        smiles: Canonical SMILES, or empty where no structure could be read.
        roles: How often the compound appears in each role, by role name.
        locations: How often it appears in inputs, workups and outcomes.
        identifiers: Every value recorded for each identifier type, with counts.
        smiles_variants: The distinct structural identifier values that canonicalize
            to ``smiles``. More than one means the same structure was written
            differently, which is worth a reviewer's glance but is not an error.
        details: Free-text details recorded on identifiers and preparations.
        problems: What disagrees within the compound's own messages, one line each.
        reaction_count: How many reactions use the compound.
        example_reaction_ids: A few of those reactions, for finding an example.
    """

    key: str
    smiles: str = ""
    roles: Counter[str] = dataclasses.field(default_factory=Counter)
    locations: Counter[str] = dataclasses.field(default_factory=Counter)
    identifiers: dict[str, Counter[str]] = dataclasses.field(
        default_factory=lambda: defaultdict(Counter)
    )
    smiles_variants: set[str] = dataclasses.field(default_factory=set)
    details: set[str] = dataclasses.field(default_factory=set)
    problems: set[str] = dataclasses.field(default_factory=set)
    reaction_count: int = 0
    example_reaction_ids: list[str] = dataclasses.field(default_factory=list)

    @property
    def names(self) -> list[str]:
        """The names recorded for the compound, most used first."""
        return [
            name for name, _ in self.identifiers.get("NAME", Counter()).most_common()
        ]

    @property
    def label(self) -> str:
        """A short, readable label: the most used name, else the SMILES or key."""
        if self.key.startswith(_CRUDE_PREFIX):
            # Names the upstream reaction that made the crude, not the one using it.
            return f"crude product of {self.key.removeprefix(_CRUDE_PREFIX)}"
        names = self.names
        return names[0] if names else (self.smiles or self.key)

    def merge(self, other: "CompoundRecord") -> None:
        """Adds ``other``'s observations of the same compound to this record."""
        if other.key != self.key:
            raise ValueError(f"cannot merge {other.key!r} into {self.key!r}")
        self.smiles = self.smiles or other.smiles
        self.roles.update(other.roles)
        self.locations.update(other.locations)
        for identifier_type, values in other.identifiers.items():
            self.identifiers[identifier_type].update(values)
        self.smiles_variants |= other.smiles_variants
        self.details |= other.details
        self.problems |= other.problems
        self.reaction_count += other.reaction_count
        for reaction_id in other.example_reaction_ids:
            if len(self.example_reaction_ids) >= _EXAMPLES_KEPT:
                break
            if reaction_id not in self.example_reaction_ids:
                self.example_reaction_ids.append(reaction_id)

    def conflicts(self) -> list[str]:
        """What disagrees within this record, one line each, in a stable order."""
        found = sorted(self.problems)
        cas = sorted(self.identifiers.get("CAS_NUMBER", Counter()))
        if len(cas) > 1:
            found.append(f"more than one CAS number: {', '.join(cas)}")
        return found


class CompoundTable:
    """The compounds of a dataset, built one reaction at a time.

    Add reactions with :meth:`add_reaction`, or combine tables built over separate
    parts of a dataset with :meth:`merge`. Nothing is resolved while adding: the
    table only records what each reaction says.
    """

    def __init__(self) -> None:
        """Starts an empty table."""
        self._records: dict[str, CompoundRecord] = {}
        self._name_conflicts: dict[str, list[str]] | None = None
        self.reaction_count = 0

    def __len__(self) -> int:
        """The number of distinct compounds."""
        return len(self._records)

    def add_reaction(self, reaction: reaction_pb2.Reaction) -> None:
        """Records every compound ``reaction`` uses or produces.

        Args:
            reaction: The reaction to read.
        """
        self.reaction_count += 1
        seen: set[str] = set()
        for location, compound, crude in _compounds_in(reaction):
            record = (
                _crude_record(crude) if crude is not None else _record_for(compound)
            )
            record.locations[location] += 1
            if record.key not in seen:
                seen.add(record.key)
                record.reaction_count = 1
                if reaction.reaction_id:
                    record.example_reaction_ids.append(reaction.reaction_id)
            self._absorb(record)

    def add_reactions(self, reactions: Iterable[reaction_pb2.Reaction]) -> None:
        """Records every compound in each of ``reactions``."""
        for reaction in reactions:
            self.add_reaction(reaction)

    def merge(self, other: "CompoundTable") -> None:
        """Adds the compounds of a table built over another part of the dataset."""
        self.reaction_count += other.reaction_count
        for record in other.records():
            # A copy, so later additions to either table do not reach into the other.
            self._absorb(copy.deepcopy(record))

    def _absorb(self, record: CompoundRecord) -> None:
        self._name_conflicts = None
        existing = self._records.get(record.key)
        if existing is None:
            self._records[record.key] = record
        else:
            existing.merge(record)

    def records(self) -> list[CompoundRecord]:
        """Every compound, most used first, then by key so the order is stable."""
        return sorted(self._records.values(), key=lambda r: (-r.reaction_count, r.key))

    def by_role(self) -> dict[str, list[CompoundRecord]]:
        """Compounds grouped under each role they appear in, in reading order.

        A compound used in several roles appears under each of them.
        """
        groups: dict[str, list[CompoundRecord]] = defaultdict(list)
        for record in self.records():
            for role in record.roles:
                groups[role].append(record)
        ordered = {role: groups[role] for role in _ROLE_ORDER if role in groups}
        ordered.update(
            {role: groups[role] for role in sorted(groups) if role not in ordered}
        )
        return ordered

    def cross_role(self) -> list[CompoundRecord]:
        """Compounds used in more than one role, e.g. REACTANT here, SOLVENT there."""
        return [record for record in self.records() if len(record.roles) > 1]

    def name_conflicts(self) -> dict[str, list[str]]:
        """Names used for more than one structure, with the structures they name.

        Returns:
            Each such name, normalized, mapped to the sorted SMILES it is used for.
        """
        if self._name_conflicts is None:
            self._name_conflicts = self._find_name_conflicts()
        return self._name_conflicts

    def _find_name_conflicts(self) -> dict[str, list[str]]:
        structures: dict[str, set[str]] = defaultdict(set)
        for record in self._records.values():
            if not record.smiles:
                continue
            for name in record.identifiers.get("NAME", Counter()):
                structures[_normalize_name(name)].add(record.smiles)
        return {
            name: sorted(smiles)
            for name, smiles in sorted(structures.items())
            if len(smiles) > 1
        }

    def conflicts(self, record: CompoundRecord) -> list[str]:
        """Everything that disagrees about ``record``, its own and across the table."""
        found = record.conflicts()
        clashes = self.name_conflicts()
        reported = set()
        for name in record.names:
            normalized = _normalize_name(name)
            if normalized in clashes and normalized not in reported:
                reported.add(normalized)
                others = [s for s in clashes[normalized] if s != record.smiles]
                found.append(f"name {name!r} is also used for: {', '.join(others)}")
        return found

    def identifier_types(self) -> list[str]:
        """The identifier types in use anywhere in the table, in schema order."""
        in_use = {name for r in self._records.values() for name in r.identifiers}
        # Schema order, so columns line up the same way in every report.
        ordered = [_identifier_name(number) for number in _IdentifierType.values()]
        return [name for name in ordered if name in in_use]

    def summary(self) -> dict[str, int]:
        """Headline counts for the top of a report."""
        records = self.records()
        return {
            "reactions": self.reaction_count,
            "compounds": len(records),
            "with_structure": sum(1 for r in records if r.smiles),
            "without_cas": sum(
                1 for r in records if r.smiles and "CAS_NUMBER" not in r.identifiers
            ),
            "in_more_than_one_role": len(self.cross_role()),
            "with_conflicts": sum(1 for r in records if self.conflicts(r)),
            "written_more_than_one_way": sum(
                1 for r in records if len(r.smiles_variants) > 1
            ),
        }

    def rows(self) -> Iterator[dict[str, str]]:
        """One flat row per compound, as written to the CSV.

        Every value of every identifier type in use gets a column; values are joined
        with `` | `` in order of use, so nothing recorded is dropped.
        """
        types = self.identifier_types()
        for record in self.records():
            row = {
                "key": record.key,
                "label": record.label,
                "smiles": record.smiles,
                "roles": _counts(record.roles),
                "locations": _counts(record.locations),
                "reaction_count": str(record.reaction_count),
            }
            for identifier_type in types:
                values = record.identifiers.get(identifier_type, Counter())
                row[identifier_type] = " | ".join(v for v, _ in values.most_common())
            row["smiles_variants"] = (
                " | ".join(sorted(record.smiles_variants))
                if len(record.smiles_variants) > 1
                else ""
            )
            row["details"] = " | ".join(sorted(record.details))
            row["conflicts"] = " | ".join(self.conflicts(record))
            row["example_reaction_ids"] = " | ".join(record.example_reaction_ids)
            yield row

    def write_csv(self, path: str | os.PathLike[str]) -> None:
        """Writes the table as CSV, one row per compound."""
        rows = list(self.rows())
        fields = (
            list(rows[0])
            if rows
            else ["key", "label", "smiles", "roles", "locations", "reaction_count"]
        )
        with pathlib.Path(path).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    def write_html(self, path: str | os.PathLike[str], *, title: str = "") -> None:
        """Writes a self-contained HTML page of the table, one section per role."""
        with pathlib.Path(path).open("w", encoding="utf-8") as f:
            f.write(self.to_html(title=title))

    def to_html(self, *, title: str = "") -> str:
        """Returns a self-contained HTML page of the table, one section per role."""
        return _render_html(self, title or "Compounds")


def _counts(counter: Counter[str]) -> str:
    """Returns ``counter`` as ``A (3); B (1)``, most common first."""
    return "; ".join(f"{name} ({count})" for name, count in counter.most_common())


def _compounds_in(
    reaction: reaction_pb2.Reaction,
) -> Iterator[
    tuple[
        str,
        reaction_pb2.Compound | reaction_pb2.ProductCompound | None,
        reaction_pb2.CrudeComponent | None,
    ]
]:
    """Yields each compound in ``reaction`` with where it was recorded.

    Yields:
        ``(location, compound, None)`` for a compound, or ``(location, None, crude)``
        for a crude component.
    """
    inputs = [(INPUT, value) for value in reaction.inputs.values()]
    inputs += [(WORKUP, workup.input) for workup in reaction.workups]
    for location, reaction_input in inputs:
        for compound in reaction_input.components:
            yield location, compound, None
        for crude in reaction_input.crude_components:
            yield location, None, crude
    for outcome in reaction.outcomes:
        for product in outcome.products:
            yield OUTCOME, product, None


def _crude_record(crude: reaction_pb2.CrudeComponent) -> CompoundRecord:
    """Returns a one-observation record for a crude component."""
    record = CompoundRecord(key=f"{_CRUDE_PREFIX}{crude.reaction_id}")
    record.roles[CRUDE] += 1
    return record


def _record_for(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound | None,
) -> CompoundRecord:
    """Returns a one-observation record for a compound, keyed by its identity."""
    assert compound is not None
    structures: dict[str, str] = {}
    problems: set[str] = set()
    for identifier in message_helpers.structural_identifiers(compound):
        canonical = message_helpers.canonical_smiles_for_identifier(
            identifier.type, identifier.value
        )
        type_name = _identifier_name(identifier.type)
        if canonical is None:
            problems.add(f"unreadable {type_name}: {identifier.value}")
        else:
            structures[identifier.value] = canonical
    distinct = sorted(set(structures.values()))
    if len(distinct) > 1:
        problems.add(
            "identifiers in one compound describe different structures: "
            + ", ".join(distinct)
        )
    smiles = message_helpers.smiles_from_compound(compound) or ""

    names = [
        identifier.value
        for identifier in compound.identifiers
        if identifier.type == _IdentifierType.NAME and identifier.value
    ]
    if smiles:
        key = f"smiles:{smiles}"
    elif names:
        key = f"name:{_normalize_name(names[0])}"
    elif compound.identifiers:
        first = compound.identifiers[0]
        key = f"{_identifier_name(first.type).lower()}:{first.value}"
    else:
        key = "unidentified"

    record = CompoundRecord(key=key, smiles=smiles, problems=problems)
    record.smiles_variants = {
        value for value, canonical in structures.items() if canonical == smiles
    }
    role = compound.reaction_role
    if (
        isinstance(compound, reaction_pb2.ProductCompound)
        and role == _RoleType.UNSPECIFIED
    ):
        # An outcome's product with no role recorded is still a product; listing it
        # under UNSPECIFIED would put it beside inputs nobody gave a role.
        role = _RoleType.PRODUCT
    record.roles[_role_name(role)] += 1
    for identifier in compound.identifiers:
        value = (
            f"({_identifier_name(identifier.type).lower()} present)"
            if identifier.type in _BULKY_TYPES
            else identifier.value
        )
        record.identifiers[_identifier_name(identifier.type)][value] += 1
        if identifier.details:
            record.details.add(identifier.details)
    if isinstance(compound, reaction_pb2.Compound):
        for preparation in compound.preparations:
            if preparation.details:
                record.details.add(preparation.details)
    return record


@functools.lru_cache(maxsize=4096)
def _structure_svg(smiles: str) -> str:
    """Returns an inline SVG drawing of ``smiles``, or "" if it cannot be drawn."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    drawer = rdMolDraw2D.MolDraw2DSVG(180, 130)
    drawer.drawOptions().clearBackground = False
    try:
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    except (RuntimeError, ValueError):
        return ""
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # Drop the XML prolog: the drawing is embedded in an HTML page, not a file.
    return svg[svg.index("<svg") :] if "<svg" in svg else ""


_STYLE = """
:root { --ink: #1b2430; --muted: #5d6b7a; --line: #d9dee3; --bg: #ffffff;
  --head: #f2f4f6; --flag: #a33a2b; --flag-bg: #fbecea; --note: #7a5a00; }
* { box-sizing: border-box; }
body { margin: 0; padding: 24px; color: var(--ink); background: var(--bg);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
h1 { font-size: 24px; margin: 0 0 4px; }
h2 { font-size: 18px; margin: 32px 0 8px; padding-bottom: 4px;
  border-bottom: 2px solid var(--ink); }
p.meta { color: var(--muted); margin: 0 0 16px; }
dl.summary { display: grid; gap: 8px 16px; margin: 0 0 8px;
  grid-template-columns: repeat(auto-fill, minmax(170px, 1fr)); }
dl.summary div { border: 1px solid var(--line); border-radius: 6px; padding: 8px 10px; }
dl.summary dt { color: var(--muted); font-size: 12px; }
dl.summary dd { margin: 0; font-size: 20px; font-variant-numeric: tabular-nums; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid var(--line); padding: 6px 8px; vertical-align: top;
  text-align: left; }
th { background: var(--head); font-size: 12px; }
td.structure { width: 190px; }
td.structure svg { display: block; }
td code, li code { font-size: 12px; word-break: break-all; }
td ul { margin: 0; padding-left: 16px; }
tr.flagged td { background: var(--flag-bg); }
.conflict { color: var(--flag); font-weight: 600; }
.note { color: var(--note); }
.muted { color: var(--muted); }
@media print {
  body { padding: 0; font-size: 10pt; }
  h2 { break-after: avoid; }
  section { break-before: page; }
  section:first-of-type { break-before: auto; }
  tr { break-inside: avoid; }
  thead { display: table-header-group; }
  tr.flagged td { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
}
"""


def _cell_list(values: Iterable[str], css: str = "") -> str:
    items = [f"<li><code>{html.escape(v)}</code></li>" for v in values]
    if not items:
        return ""
    cls = f' class="{css}"' if css else ""
    return f"<ul{cls}>{''.join(items)}</ul>"


def _render_html(table: CompoundTable, title: str) -> str:
    summary = table.summary()
    types = table.identifier_types()
    labels = {
        "reactions": "reactions",
        "compounds": "distinct compounds",
        "with_structure": "with a readable structure",
        "without_cas": "structures without a CAS number",
        "in_more_than_one_role": "used in more than one role",
        "with_conflicts": "with something to check",
        "written_more_than_one_way": "structures written more than one way",
    }
    parts = [
        "<!DOCTYPE html>",
        '<html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{html.escape(title)}</title><style>{_STYLE}</style></head><body>",
        f"<h1>{html.escape(title)}</h1>",
        '<p class="meta">Every compound in the dataset, grouped by role. Values that '
        "disagree are listed, not resolved.</p>",
        '<dl class="summary">',
    ]
    parts.extend(
        f"<div><dt>{label}</dt><dd>{summary[key]}</dd></div>"
        for key, label in labels.items()
    )
    parts.append("</dl>")

    cross = table.cross_role()
    if cross:
        parts.append("<section><h2>Compounds used in more than one role</h2>")
        parts.append("<table><thead><tr><th>Compound</th><th>Roles</th><th>SMILES</th>")
        parts.append("</tr></thead><tbody>")
        parts.extend(
            f"<tr><td>{html.escape(record.label)}</td>"
            f"<td>{html.escape(_counts(record.roles))}</td>"
            f"<td><code>{html.escape(record.smiles)}</code></td></tr>"
            for record in cross
        )
        parts.append("</tbody></table></section>")

    for role, records in table.by_role().items():
        parts.append(f"<section><h2>{html.escape(role)} ({len(records)})</h2>")
        parts.append("<table><thead><tr><th>Structure</th><th>Compound</th>")
        parts.extend(f"<th>{html.escape(name)}</th>" for name in types)
        parts.append("<th>Reactions</th><th>Details</th><th>To check</th>")
        parts.append("</tr></thead><tbody>")
        for record in records:
            conflicts = table.conflicts(record)
            row_class = ' class="flagged"' if conflicts else ""
            svg = _structure_svg(record.smiles) if record.smiles else ""
            structure = svg or '<span class="muted">no structure</span>'
            compound = [f"<strong>{html.escape(record.label)}</strong>"]
            if record.smiles and record.label != record.smiles:
                compound.append(f"<div><code>{html.escape(record.smiles)}</code></div>")
            if len(record.roles) > 1:
                roles = html.escape(_counts(record.roles))
                compound.append(f'<div class="muted">roles: {roles}</div>')
            parts.append(f'<tr{row_class}><td class="structure">{structure}</td>')
            parts.append(f"<td>{''.join(compound)}</td>")
            for identifier_type in types:
                values = record.identifiers.get(identifier_type, Counter())
                parts.append(
                    f"<td>{_cell_list(v for v, _ in values.most_common())}</td>"
                )
            examples = ", ".join(record.example_reaction_ids)
            parts.append(
                f"<td>{record.reaction_count}"
                + (
                    f'<div class="muted"><code>{html.escape(examples)}</code></div>'
                    if examples
                    else ""
                )
                + "</td>"
            )
            parts.append(f"<td>{_cell_list(sorted(record.details))}</td>")
            to_check = _cell_list(conflicts, "conflict")
            if len(record.smiles_variants) > 1:
                to_check += '<div class="note">written more than one way:</div>'
                to_check += _cell_list(sorted(record.smiles_variants), "note")
            parts.append(f"<td>{to_check}</td></tr>")
        parts.append("</tbody></table></section>")
    parts.append("</body></html>")
    return "\n".join(parts)
