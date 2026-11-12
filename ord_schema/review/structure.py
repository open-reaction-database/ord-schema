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

"""Dataset-level structural checks, and the dataset profile that selects them.

Single-reaction validation cannot see a dataset whose reactions are each fine but
disagree with one another: a field most reactions fill and a few leave empty, a unit
that changes halfway through, a well used twice. These checks look across the whole
dataset for that kind of problem and report it to a reviewer.

A dataset declares a **profile**. ``uniform`` (high-throughput screens, DoE and BO
tables, automated synthesis) is the default and runs every check. ``varied``
(literature extraction, ELN exports, collections of single syntheses) skips the
checks that assume every reaction follows one template. A profile never changes how
serious a check is, only whether it applies; checks it skips are reported as not
applied, so what a profile relaxed stays visible.

Every check has a stable ID, e.g. ``coverage-gap``, so findings can be counted by
check rather than by message, and so a check can later move into
:mod:`ord_schema.validations` unchanged. Severities are
:class:`ord_schema.validations.Severity`. Nothing here blocks a submission: the
findings are for the reviewer, who decides.

Prototyped here rather than in :mod:`ord_schema.validations` until where
dataset-level checks live is agreed (see #1043)::

    report = check_dataset(datasets.load_dataset("dataset.pb.gz").reactions)
    for finding in report.findings:
        print(finding.severity.name, finding.check, finding.message)
"""

import dataclasses
import enum
import hashlib
import itertools
import math
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable

import pandas as pd

from ord_schema.analysis import flat, summary
from ord_schema.proto import reaction_pb2
from ord_schema.validations import Severity


class ValidationProfile(enum.Enum):
    """How a dataset's reactions relate to one another, which selects the checks."""

    # Reactions share one template: screens, DoE/BO tables, automated synthesis.
    UNIFORM = "uniform"
    # Reactions do not: literature extraction, ELN exports, single syntheses.
    VARIED = "varied"


# A field this share of reactions fills, but not all, is a coverage gap.
COVERAGE_GAP_THRESHOLD = 0.8

# A grid with more cells than this is not checked for gaps; listing them would not
# help a reviewer, and a design that large is unlikely to be a full grid.
MAX_GRID_CELLS = 1_000_000

# Below this share of its cells filled, a set of factors is not treated as a grid.
MIN_GRID_FILL = 0.5

# With more candidate factors than this, a dataset is not a designed screen but an
# aggregate, such as a year of USPTO grants with hundreds of distinct input keys;
# comparing every pair of candidates there would take hours and find no grid.
MAX_FACTOR_CANDIDATES = 40

# Examples kept per finding.
MAX_EXAMPLES = 5

# A check that reports compound by compound lists at most this many, the most
# frequent first, and sums up the rest in one finding: an aggregate such as a
# reaction-mining dataset can otherwise yield one finding for each of 100,000
# compounds.
MAX_COMPOUND_FINDINGS = 20

# Values recorded where a real value was not known.
PLACEHOLDER_STRINGS = frozenset(
    {"", "-", "--", "?", "n/a", "na", "nan", "none", "null", "tbd", "unknown"}
)

_AMOUNT_KINDS = ("mass", "moles", "volume")


@dataclasses.dataclass
class Finding:
    """One thing a check found.

    Attributes:
        check: The stable ID of the check that found it.
        severity: How serious it is.
        message: What was found and, where it helps, what to change.
        fields: The field paths or columns concerned.
        reactions: How many reactions it affects, where that is countable.
        examples: A few reaction IDs or values, for finding it in the data.
    """

    check: str
    severity: Severity
    message: str
    fields: list[str] = dataclasses.field(default_factory=list)
    reactions: int | None = None
    examples: list[str] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class _Observations:
    """What the checks need from each reaction beyond the summary and table."""

    fingerprints: dict[str, list[str]] = dataclasses.field(
        default_factory=lambda: defaultdict(list)
    )
    wells: dict[tuple[str, str], list[str]] = dataclasses.field(
        default_factory=lambda: defaultdict(list)
    )
    with_source: int = 0
    without_source: list[str] = dataclasses.field(default_factory=list)
    with_source_examples: list[str] = dataclasses.field(default_factory=list)

    def observe(self, reaction: reaction_pb2.Reaction, label: str) -> None:
        self.fingerprints[_fingerprint(reaction)].append(label)
        vessel = reaction.setup.vessel
        position = vessel.position or (
            f"{vessel.row}{vessel.col}" if vessel.row or vessel.col else ""
        )
        if position:
            self.wells[(vessel.vessel_id, position)].append(label)
        provenance = reaction.provenance
        if provenance.doi or provenance.patent or provenance.publication_url:
            self.with_source += 1
            if len(self.with_source_examples) < MAX_EXAMPLES:
                self.with_source_examples.append(label)
        else:
            self.without_source.append(label)


def _fingerprint(reaction: reaction_pb2.Reaction) -> str:
    """Returns a digest of ``reaction`` with its IDs and record history removed.

    Two reactions with the same fingerprint record the same experiment; they differ
    at most in their reaction ID and in who recorded them and when.
    """
    copy = reaction_pb2.Reaction()
    copy.CopyFrom(reaction)
    copy.ClearField("reaction_id")
    copy.provenance.ClearField("record_created")
    copy.provenance.ClearField("record_modified")
    return hashlib.sha256(copy.SerializeToString(deterministic=True)).hexdigest()


@dataclasses.dataclass
class DatasetEvidence:
    """Everything the checks read: the summary, the flat table, and observations."""

    summary: summary.DatasetSummary
    table: flat.FlatTable
    observations: _Observations

    def label(self, row: int) -> str:
        """The reaction ID of the table's ``row``, or its position if it has none."""
        if "reaction_id" in self.table.varying.columns:
            value = self.table.varying["reaction_id"].iloc[row]
            if isinstance(value, str) and value:
                return value
        return f"reaction #{row + 1}"


def gather(reactions: Iterable[reaction_pb2.Reaction]) -> DatasetEvidence:
    """Reads ``reactions`` once and returns what the checks need."""
    dataset_summary = summary.DatasetSummary()
    observations = _Observations()

    def each() -> Iterable[reaction_pb2.Reaction]:
        for position, reaction in enumerate(reactions):
            dataset_summary.observe(reaction)
            observations.observe(
                reaction, reaction.reaction_id or f"reaction #{position + 1}"
            )
            yield reaction

    table = flat.flatten(each())
    return DatasetEvidence(
        summary=dataset_summary, table=table, observations=observations
    )


@dataclasses.dataclass(frozen=True)
class Check:
    """A dataset-level check.

    Attributes:
        id: Stable ID, used in findings and to decide which profiles skip it.
        description: What it looks for, for documentation and the report.
        profiles: The profiles it applies to.
        run: Returns its findings from the evidence.
    """

    id: str
    description: str
    profiles: frozenset[ValidationProfile]
    run: Callable[[DatasetEvidence], list[Finding]]


_UNIFORM_ONLY = frozenset({ValidationProfile.UNIFORM})
_BOTH = frozenset(ValidationProfile)


def _examples(labels: Iterable[str]) -> list[str]:
    return list(itertools.islice(labels, MAX_EXAMPLES))


def _check_coverage_gaps(evidence: DatasetEvidence) -> list[Finding]:
    """Fields most reactions fill and a few leave empty.

    Fields empty in exactly the same reactions make one finding, since they are
    usually one omission: a measurement, or a component, those reactions lack.
    """
    table = evidence.table
    total = table.reaction_count
    gaps: dict[tuple[int, ...], list[str]] = {}
    for group in table.column_groups:
        column = table.varying[group[0]]
        filled = int(column.notna().sum())
        if total and COVERAGE_GAP_THRESHOLD <= filled / total < 1:
            empty = tuple(int(i) for i in column.index[column.isna()])
            gaps.setdefault(empty, []).extend(group)
    findings = []
    for empty, fields in gaps.items():
        named = (
            fields[0] if len(fields) == 1 else f"{fields[0]} and {len(fields) - 1} more"
        )
        findings.append(
            Finding(
                check="coverage-gap",
                severity=Severity.WARNING,
                message=(
                    f"{named} filled in {total - len(empty)} of {total} reactions "
                    f"and empty in the same {len(empty)}; check whether those are "
                    "missing data or deliberate, such as controls"
                ),
                fields=fields,
                reactions=len(empty),
                examples=_examples(evidence.label(row) for row in empty),
            )
        )
    return findings


def _is_placeholder(path: str, value: object) -> bool:
    if isinstance(value, str):
        # An enum's NONE or UNSPECIFIED is a real value (stirring.type NONE means not
        # stirred), so only free text is checked.
        field = flat.field_for(path)
        if field is None or field.type != field.TYPE_STRING:
            return False
        # An identifier's value is a name or structure, and "Na" or "None" can be a
        # real one. A placeholder name is a compound with no structure, which the
        # compound table shows.
        if path.endswith("identifiers.value"):
            return False
        return value.strip().casefold() in PLACEHOLDER_STRINGS
    if isinstance(value, int | float) and not isinstance(value, bool):
        # A zero amount stands in for an unknown one; a zero elsewhere (a
        # temperature, an addition order) can be real.
        return value == 0 and ".amount." in f".{path}."
    return False


def _check_placeholder_values(evidence: DatasetEvidence) -> list[Finding]:
    """Values that stand in for unknown ones, repeated across reactions."""
    dataset_summary = evidence.summary
    findings = []
    # DatasetSummary.values, a dict of counters, not a pandas attribute.
    for path, counter in sorted(dataset_summary.values.items()):  # noqa: PD011
        if summary.is_keyed(path):
            continue
        for value, count in counter.items():
            if count > 1 and _is_placeholder(path, value):
                findings.append(
                    Finding(
                        check="placeholder-value",
                        severity=Severity.WARNING,
                        message=(
                            f"{path} is {value!r} {count} times; if that stands for "
                            "an unknown value, leave the field empty instead"
                        ),
                        fields=[path],
                        reactions=None,
                        examples=[repr(value)],
                    )
                )
    return findings


# Two factors whose combinations fill less than this share of their joint grid are
# one factor at paired levels: reactant pairs, say, rather than every reactant with
# every other.
PAIRED_FILL = 0.75


def _factor_groups(table: flat.FlatTable) -> list[list[str]]:
    """Column groups that look like a design factor: an input set at a few levels.

    A factor is an input group whose values each recur, and which no other factor
    determines: a field that only says whether the ligand is present is decided by
    which ligand it is, so it adds no cells to the grid.
    """
    frame = table.varying
    candidates = []
    for group in table.column_groups:
        name = group[0]
        if not name.startswith("inputs"):
            continue
        counts = frame[name].value_counts(dropna=False)
        if len(counts) < 2 or counts.min() < 2:
            continue
        candidates.append(group)
        if len(candidates) > MAX_FACTOR_CANDIDATES:
            return []
    factors: list[list[str]] = []
    for group in sorted(candidates, key=lambda g: -frame[g[0]].nunique(dropna=False)):
        determined = any(
            frame[[kept[0], group[0]]].drop_duplicates().shape[0]
            == frame[kept[0]].nunique(dropna=False)
            for kept in factors
        )
        if not determined:
            factors.append(group)
    return sorted(factors, key=lambda g: list(frame.columns).index(g[0]))


def _paired_factors(codes: pd.DataFrame) -> list[list[str]]:
    """Merges factors used in fixed pairings into one, until none remain.

    Args:
        codes: One column per factor, holding each reaction's level as an integer.

    Returns:
        The factors, each a list of the columns it combines.
    """
    factors = [[name] for name in codes.columns]

    def levels(factor: list[str]) -> int:
        return len(codes[factor].drop_duplicates())

    merged = True
    while merged and len(factors) > 1:
        merged = False
        for a, b in itertools.combinations(range(len(factors)), 2):
            joint = levels(factors[a] + factors[b])
            if joint < PAIRED_FILL * levels(factors[a]) * levels(factors[b]):
                factors[a] = factors[a] + factors[b]
                del factors[b]
                merged = True
                break
    return factors


def _check_grid_gaps(evidence: DatasetEvidence) -> list[Finding]:
    """Combinations of a screen's factors that no reaction covers."""
    table = evidence.table
    names = [group[0] for group in _factor_groups(table)]
    if len(names) < 2:
        return []
    columns = table.varying[names]
    codes = pd.DataFrame(
        {name: pd.factorize(columns[name], use_na_sentinel=False)[0] for name in names}
    )
    factors = _paired_factors(codes)
    if len(factors) < 2:
        return []
    # Each factor's levels, as the tuples of codes that occur, and the cells filled.
    levels = [
        list(dict.fromkeys(map(tuple, codes[f].to_numpy().tolist()))) for f in factors
    ]
    cells = math.prod(len(level) for level in levels)
    positions = [[codes.columns.get_loc(name) for name in f] for f in factors]
    filled = {
        tuple(tuple(row[i] for i in where) for where in positions)
        for row in codes.drop_duplicates().to_numpy().tolist()
    }
    if cells > MAX_GRID_CELLS or len(filled) / cells < MIN_GRID_FILL:
        return []
    missing = [cell for cell in itertools.product(*levels) if cell not in filled]
    if not missing:
        return []
    labels = {
        name: pd.factorize(columns[name], use_na_sentinel=False)[1] for name in names
    }

    def describe(cell: tuple[tuple[int, ...], ...]) -> str:
        parts = []
        for factor, level in zip(factors, cell, strict=True):
            parts.append(
                " + ".join(
                    str(labels[name][code])
                    for name, code in zip(factor, level, strict=True)
                )
            )
        return " | ".join(parts)

    shape = " × ".join(str(len(level)) for level in levels)
    return [
        Finding(
            check="grid-gap",
            severity=Severity.INFO,
            message=(
                f"{len(missing)} of the {cells} combinations of the {len(factors)} "
                f"varied inputs ({shape}) have no reaction; check whether they were "
                "left out deliberately"
            ),
            fields=[name for factor in factors for name in factor],
            examples=[describe(cell) for cell in missing[:MAX_EXAMPLES]],
        )
    ]


def _check_duplicate_wells(evidence: DatasetEvidence) -> list[Finding]:
    """A plate position more than one reaction claims."""
    findings = []
    for (plate, position), labels in sorted(evidence.observations.wells.items()):
        if len(labels) > 1:
            where = f"plate {plate!r}, " if plate else ""
            findings.append(
                Finding(
                    check="duplicate-well",
                    severity=Severity.WARNING,
                    message=f"{len(labels)} reactions share {where}position {position}",
                    fields=["setup.vessel.position"],
                    reactions=len(labels),
                    examples=_examples(labels),
                )
            )
    return findings


def _check_mixed_units(evidence: DatasetEvidence) -> list[Finding]:
    """One position recorded in different units, or as different kinds of amount."""
    table = evidence.table
    frame = table.varying
    findings = []
    for name in frame.columns:
        if name.endswith(".units"):
            counts = frame[name].value_counts()
            if len(counts) > 1:
                used = ", ".join(f"{unit} ({count})" for unit, count in counts.items())
                findings.append(
                    Finding(
                        check="mixed-units",
                        severity=Severity.WARNING,
                        message=f"{name} is recorded in {used}",
                        fields=[name],
                        reactions=int(counts.sum() - counts.iloc[0]),
                    )
                )
    # The same position given as a mass in some reactions and moles in others.
    kinds: dict[str, dict[str, pd.Series]] = defaultdict(dict)
    for name in list(frame.columns) + list(table.constants):
        head, sep, tail = name.partition(".amount.")
        kind = tail.split(".")[0]
        if sep and kind in _AMOUNT_KINDS and tail.endswith(".value"):
            column = (
                frame[name].notna()
                if name in frame.columns
                else pd.Series(data=True, index=range(table.reaction_count))
            )
            kinds[head][kind] = column
    for head, by_kind in kinds.items():
        if len(by_kind) < 2:
            continue
        present = pd.DataFrame(by_kind)
        patterns = present.apply(lambda row: "+".join(row.index[row]), axis=1)
        counts = patterns[patterns != ""].value_counts()
        if len(counts) > 1:
            used = ", ".join(f"{p} ({c})" for p, c in counts.items())
            findings.append(
                Finding(
                    check="mixed-units",
                    severity=Severity.WARNING,
                    message=f"{head}.amount is given as {used}",
                    fields=[f"{head}.amount"],
                    reactions=int(counts.sum() - counts.iloc[0]),
                )
            )
    return findings


def _capped(findings: list[Finding], noun: str) -> list[Finding]:
    """The MAX_COMPOUND_FINDINGS most frequent findings, and one for the rest.

    Args:
        findings: One check's findings, one per compound.
        noun: What each finding is about, in the plural, for the summary.

    Returns:
        The findings, capped; the summary's examples are the first messages it covers.
        The summary leaves ``reactions`` unset: one reaction can hold several of the
        compounds it covers, so their counts don't add up to a count of reactions.
        Its message gives the most any one of them is in instead.
    """
    if len(findings) <= MAX_COMPOUND_FINDINGS:
        return findings
    ranked = sorted(findings, key=lambda f: -(f.reactions or 0))
    kept, rest = ranked[:MAX_COMPOUND_FINDINGS], ranked[MAX_COMPOUND_FINDINGS:]
    most = max((f.reactions for f in rest if f.reactions is not None), default=None)
    message = f"and {len(rest):,} more {noun}"
    if most is not None:
        message += f", each in at most {most:,} reactions"
    summary = Finding(
        check=rest[0].check,
        severity=rest[0].severity,
        message=message,
        fields=rest[0].fields,
        reactions=None,
        examples=[f.message for f in rest[:MAX_EXAMPLES]],
    )
    return [*kept, summary]


def _check_name_structure(evidence: DatasetEvidence) -> list[Finding]:
    """One name for several structures, or one structure under several names."""
    compounds = evidence.summary.compounds
    conflicts = [
        Finding(
            check="name-structure-conflict",
            severity=Severity.WARNING,
            message=f"the name {name!r} is used for {len(smiles)} structures",
            fields=["NAME"],
            examples=smiles[:MAX_EXAMPLES],
        )
        for name, smiles in compounds.name_conflicts().items()
    ]
    variants = []
    for record in compounds.records():
        names = record.names
        if record.smiles and len(names) > 1:
            variants.append(
                Finding(
                    check="structure-name-variants",
                    severity=Severity.INFO,
                    message=(
                        f"{record.smiles} is named {len(names)} ways: "
                        + ", ".join(repr(n) for n in names[:MAX_EXAMPLES])
                    ),
                    fields=["NAME"],
                    reactions=record.reaction_count,
                    examples=record.example_reaction_ids,
                )
            )
    return _capped(conflicts, "names used for several structures") + _capped(
        variants, "structures under several names"
    )


def _check_conflicting_roles(evidence: DatasetEvidence) -> list[Finding]:
    """The same compound in different roles in different reactions."""
    findings = [
        Finding(
            check="conflicting-roles",
            severity=Severity.INFO,
            message=(
                f"{record.label} appears as "
                + ", ".join(f"{role} ({n})" for role, n in record.roles.most_common())
            ),
            fields=["reaction_role"],
            reactions=record.reaction_count,
            examples=record.example_reaction_ids,
        )
        for record in evidence.summary.compounds.cross_role()
    ]
    return _capped(findings, "compounds in more than one role")


def _check_duplicate_reactions(evidence: DatasetEvidence) -> list[Finding]:
    """Reactions identical once their IDs and record history are removed."""
    duplicates = [
        labels
        for labels in evidence.observations.fingerprints.values()
        if len(labels) > 1
    ]
    if not duplicates:
        return []
    extra = sum(len(labels) - 1 for labels in duplicates)
    return [
        Finding(
            check="duplicate-reaction",
            severity=Severity.WARNING,
            message=(
                f"{len(duplicates)} sets of identical reactions ({extra} more than "
                "needed); check for repeated entries, or record replicates so they "
                "can be told apart"
            ),
            reactions=extra + len(duplicates),
            examples=[" = ".join(labels[:3]) for labels in duplicates[:MAX_EXAMPLES]],
        )
    ]


def _check_mixed_sourcing(evidence: DatasetEvidence) -> list[Finding]:
    """Some reactions cite a source and others do not."""
    observations = evidence.observations
    if not observations.with_source or not observations.without_source:
        return []
    return [
        Finding(
            check="mixed-sourcing",
            severity=Severity.WARNING,
            message=(
                f"{observations.with_source} reactions cite a doi, patent or "
                f"publication_url and {len(observations.without_source)} do not"
            ),
            fields=[
                "provenance.doi",
                "provenance.patent",
                "provenance.publication_url",
            ],
            reactions=len(observations.without_source),
            examples=_examples(observations.without_source),
        )
    ]


CHECKS: tuple[Check, ...] = (
    Check(
        "coverage-gap",
        "A field most reactions fill and a few leave empty.",
        _UNIFORM_ONLY,
        _check_coverage_gaps,
    ),
    Check(
        "placeholder-value",
        "A value that stands for an unknown one, such as a zero amount, repeated.",
        _UNIFORM_ONLY,
        _check_placeholder_values,
    ),
    Check(
        "grid-gap",
        "Combinations of a screen's varied inputs that no reaction covers.",
        _UNIFORM_ONLY,
        _check_grid_gaps,
    ),
    Check(
        "duplicate-well",
        "A plate position more than one reaction claims.",
        _UNIFORM_ONLY,
        _check_duplicate_wells,
    ),
    Check(
        "mixed-units",
        "One field recorded in different units, or as different kinds of amount.",
        _BOTH,
        _check_mixed_units,
    ),
    Check(
        "name-structure-conflict",
        "One name used for several structures, or one structure under several names.",
        _BOTH,
        _check_name_structure,
    ),
    Check(
        "conflicting-roles",
        "The same compound in different roles in different reactions.",
        _BOTH,
        _check_conflicting_roles,
    ),
    Check(
        "duplicate-reaction",
        "Reactions identical once their IDs and record history are removed.",
        _BOTH,
        _check_duplicate_reactions,
    ),
    Check(
        "mixed-sourcing",
        "Some reactions cite a source and others do not.",
        _BOTH,
        _check_mixed_sourcing,
    ),
)


@dataclasses.dataclass
class StructureReport:
    """The findings of the structural checks for one dataset.

    Attributes:
        profile: The profile the checks ran under.
        findings: What the applied checks found, most serious first.
        applied: The IDs of the checks that ran.
        not_applied: The IDs of the checks the profile skipped.
    """

    profile: ValidationProfile
    findings: list[Finding]
    applied: list[str]
    not_applied: list[str]

    def by_check(self) -> Counter[str]:
        """How many findings each check made."""
        return Counter(finding.check for finding in self.findings)


def run_checks(
    evidence: DatasetEvidence,
    profile: ValidationProfile = ValidationProfile.UNIFORM,
) -> StructureReport:
    """Runs the checks ``profile`` applies over already-gathered evidence."""
    findings: list[Finding] = []
    applied, not_applied = [], []
    for check in CHECKS:
        if profile in check.profiles:
            applied.append(check.id)
            findings.extend(check.run(evidence))
        else:
            not_applied.append(check.id)
    findings.sort(key=lambda f: -f.severity)
    return StructureReport(
        profile=profile, findings=findings, applied=applied, not_applied=not_applied
    )


def check_dataset(
    reactions: Iterable[reaction_pb2.Reaction],
    profile: ValidationProfile = ValidationProfile.UNIFORM,
) -> StructureReport:
    """Runs the structural checks ``profile`` applies over ``reactions``.

    Args:
        reactions: The dataset's reactions, read once.
        profile: The dataset's profile; ``uniform`` unless it declares otherwise.

    Returns:
        The findings, and which checks ran and which the profile skipped.
    """
    return run_checks(gather(reactions), profile)
