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

"""Chemistry sanity checks: whether a reaction makes chemical sense.

Validation checks that a reaction is well formed. These check what a chemist reading
it would notice: a product with an element no input brings, a stated reaction SMILES
that names different molecules from the ones recorded, a radical where none belongs,
an ion without its counterion, a product named as one stereoisomer but drawn without
stereochemistry, yields that add up past 100%, a temperature or time that cannot be
right, a product nothing measured.

Each check looks at one reaction and has a stable ID. :class:`ChemistryRollup` runs
them over a dataset as it streams past and groups what they find, so one problem in
3,000 reactions is one finding with a count and examples. Findings use
:class:`ord_schema.review.structure.Finding` and
:class:`ord_schema.validations.Severity`, and are advisory like the rest of the review.

Not repeated here, because validation already reports them: a single yield outside
0-100%, and a measurement with no ``analysis_key`` or one naming no recorded analysis.
Not checked, because the schema has no field for it: a stated concentration that
disagrees with a volume, since ``Reaction`` records no concentrations.
"""

import dataclasses
import re
from collections.abc import Callable, Iterable, Iterator

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.review.structure import Finding
from ord_schema.validations import Severity

# A temperature setpoint outside this range, in Celsius, is more likely a typo or a
# wrong unit than a real condition. Liquid nitrogen (-196 C) is inside it.
PLAUSIBLE_CELSIUS = (-200.0, 400.0)

# Product yields in one outcome may add up to a little over 100% from rounding.
YIELD_SUM_TOLERANCE = 0.5

# Examples kept per finding.
MAX_EXAMPLES = 5

# Stereodescriptors in a name: (R), (S), (1R,2S), (E), (Z), (+), (-), cis-, trans-,
# and amino-acid style L- and D-. Racemic markers ("rac-", "(±)", "RS") say the
# opposite, so a name with one is not taken to imply a single stereoisomer.
_STEREO_NAME = re.compile(
    r"\((?:\d+[a-z']*)?[RSEZ](?:,\s*(?:\d+[a-z']*)?[RSEZ])*\)"
    r"|\([+-]\)"
    r"|\b(?:cis|trans)-"
    r"|(?:^|[\s(])[LD]-"
)
_RACEMIC_NAME = re.compile(r"\brac-|\(±\)|\(\+/-\)|\(RS\)|\(SR\)")

# Digits in a message, which differ from reaction to reaction while the problem is
# the same, so messages are grouped with them masked.
_DIGITS = re.compile(r"\d+(?:\.\d+)?")

_Compound = reaction_pb2.Compound | reaction_pb2.ProductCompound
# What one check found in one reaction: (check ID, severity, message).
Problem = tuple[str, Severity, str]


def _mol(compound: _Compound) -> Chem.Mol | None:
    """The compound's structure, or None if it has none RDKit can read."""
    try:
        mol = message_helpers.mol_from_compound(compound)
    except ValueError:
        return None
    assert isinstance(mol, Chem.Mol)
    return mol


def _name(compound: _Compound) -> str:
    # message_helpers.get_compound_name is annotated for Compound alone, though a
    # ProductCompound carries identifiers the same way.
    return next(
        (
            i.value
            for i in compound.identifiers
            if i.type == reaction_pb2.CompoundIdentifier.NAME and i.value
        ),
        "",
    )


def _label(compound: _Compound) -> str:
    """A compound named well enough for a message: its name, else its SMILES."""
    return _name(compound) or message_helpers.smiles_from_compound(compound) or "?"


def _inputs(reaction: reaction_pb2.Reaction) -> list[reaction_pb2.Compound]:
    """Every compound added: reaction inputs and workup inputs."""
    compounds = [c for value in reaction.inputs.values() for c in value.components]
    compounds += [c for workup in reaction.workups for c in workup.input.components]
    return compounds


def _products(reaction: reaction_pb2.Reaction) -> list[reaction_pb2.ProductCompound]:
    return [p for outcome in reaction.outcomes for p in outcome.products]


def _elements(mol: Chem.Mol) -> set[str]:
    # Hydrogen is left out: it is implicit in most structures, and any input brings it.
    return {atom.GetSymbol() for atom in mol.GetAtoms()} - {"H"}


def _check_element_conservation(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """A product with an element that no input or workup compound brings."""
    inputs = _inputs(reaction)
    mols = [_mol(c) for c in inputs]
    # An input without a structure could be what brings the element.
    if not inputs or any(m is None for m in mols):
        return []
    available = set().union(*(_elements(m) for m in mols if m is not None))
    problems = []
    for product in _products(reaction):
        mol = _mol(product)
        if mol is None:
            continue
        missing = sorted(_elements(mol) - available)
        if missing:
            problems.append(
                (
                    "element-conservation",
                    Severity.WARNING,
                    f"product {_label(product)} contains {', '.join(missing)}, which "
                    "no input or workup compound brings",
                )
            )
    return problems


def _canonical_molecules(mols: Iterable[Chem.Mol]) -> set[str]:
    """Each molecule and salt fragment as canonical SMILES, without atom maps."""
    found = set()
    for mol in mols:
        copy = Chem.Mol(mol)
        for atom in copy.GetAtoms():
            atom.SetAtomMapNum(0)
        try:
            Chem.SanitizeMol(copy)
        except (ValueError, RuntimeError):
            continue
        found.update(Chem.MolToSmiles(copy).split("."))
    return found


def _check_reaction_smiles(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """A stated reaction SMILES that names other molecules than the ones recorded."""
    stated = next(
        (
            i.value
            for i in reaction.identifiers
            if i.type
            in (
                reaction_pb2.ReactionIdentifier.REACTION_SMILES,
                reaction_pb2.ReactionIdentifier.REACTION_CXSMILES,
            )
            and i.value
        ),
        None,
    )
    if stated is None:
        return []
    try:
        parsed = rdChemReactions.ReactionFromSmarts(stated, useSmiles=True)
    except ValueError:
        return []  # Validation reports an unreadable reaction SMILES.
    if parsed is None:
        return []
    stated_inputs = _canonical_molecules([*parsed.GetReactants(), *parsed.GetAgents()])
    stated_products = _canonical_molecules(parsed.GetProducts())
    recorded_inputs = _canonical_molecules(
        m for m in (_mol(c) for c in _inputs(reaction)) if m is not None
    )
    recorded_products = _canonical_molecules(
        m for m in (_mol(p) for p in _products(reaction)) if m is not None
    )
    problems = []
    unknown_inputs = sorted(stated_inputs - recorded_inputs)
    if recorded_inputs and unknown_inputs:
        problems.append(
            (
                "reaction-smiles-mismatch",
                Severity.WARNING,
                "the stated reaction SMILES names reactants or agents that no input "
                f"records: {', '.join(unknown_inputs)}",
            )
        )
    if recorded_products:
        unstated = sorted(recorded_products - stated_products)
        unrecorded = sorted(stated_products - recorded_products)
        if unstated or unrecorded:
            parts = []
            if unrecorded:
                parts.append(
                    f"states {', '.join(unrecorded)}, which no outcome records"
                )
            if unstated:
                parts.append(f"leaves out recorded {', '.join(unstated)}")
            problems.append(
                (
                    "reaction-smiles-mismatch",
                    Severity.WARNING,
                    f"the stated reaction SMILES {'; '.join(parts)}",
                )
            )
    return problems


def _check_radicals_and_charges(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """A structure with unpaired electrons, or a net charge without a counterion."""
    problems = []
    for compound in [*_inputs(reaction), *_products(reaction)]:
        mol = _mol(compound)
        if mol is None:
            continue
        # Open-shell metals are expected: copper(II) has an unpaired electron.
        radicals = sum(
            atom.GetNumRadicalElectrons()
            for atom in mol.GetAtoms()
            if not message_helpers.is_transition_metal(atom)
        )
        if radicals:
            problems.append(
                (
                    "radical",
                    Severity.WARNING,
                    f"{_label(compound)} has {radicals} unpaired electron(s); check "
                    "for missing hydrogens or a mistyped atom",
                )
            )
        charge = Chem.GetFormalCharge(mol)
        if charge:
            problems.append(
                (
                    "net-charge",
                    Severity.WARNING,
                    f"{_label(compound)} has a net charge of {charge:+d}; is a "
                    "counterion missing?",
                )
            )
    return problems


def _check_stereo_in_name(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """A product named as one stereoisomer whose structure leaves stereo unassigned."""
    problems = []
    for product in _products(reaction):
        name = _name(product)
        if not name or not _STEREO_NAME.search(name) or _RACEMIC_NAME.search(name):
            continue
        mol = _mol(product)
        if mol is None:
            continue
        # Stereocentres and stereogenic double bonds alike, with whether each is
        # given; a double bond written without / or \\ counts as not given.
        unassigned = sum(
            1
            for info in Chem.FindPotentialStereo(mol)
            if info.specified == Chem.StereoSpecified.Unspecified
        )
        if unassigned:
            problems.append(
                (
                    "stereo-unspecified",
                    Severity.INFO,
                    f"product {name!r} is named as one stereoisomer, but its structure "
                    f"leaves {unassigned} stereo element(s) unassigned",
                )
            )
    return problems


def _check_yield_sum(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """Product yields in one outcome that add up past 100%."""
    problems = []
    for index, outcome in enumerate(reaction.outcomes):
        total = sum(
            m.percentage.value
            for p in outcome.products
            for m in p.measurements
            if m.type == reaction_pb2.ProductMeasurement.YIELD
            and m.HasField("percentage")
        )
        if total > 100 + YIELD_SUM_TOLERANCE:
            problems.append(
                (
                    "yield-sum",
                    Severity.WARNING,
                    f"product yields in outcome {index + 1} add up to {total:g}%",
                )
            )
    return problems


def _celsius(temperature: reaction_pb2.Temperature) -> float | None:
    if temperature.units == reaction_pb2.Temperature.CELSIUS:
        return temperature.value
    if temperature.units == reaction_pb2.Temperature.KELVIN:
        return temperature.value - 273.15
    if temperature.units == reaction_pb2.Temperature.FAHRENHEIT:
        return (temperature.value - 32) * 5 / 9
    return None


def _check_temperature(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """A temperature setpoint outside the range lab reactions run at."""
    setpoint = reaction.conditions.temperature.setpoint
    if not reaction.conditions.temperature.HasField("setpoint"):
        return []
    celsius = _celsius(setpoint)
    low, high = PLAUSIBLE_CELSIUS
    if celsius is None or low <= celsius <= high:
        return []
    unit = reaction_pb2.Temperature.TemperatureUnit.Name(setpoint.units)
    return [
        (
            "implausible-temperature",
            Severity.WARNING,
            f"the temperature setpoint is {setpoint.value:g} {unit}; check the value "
            "and the unit",
        )
    ]


def _check_reaction_time(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """A reaction time of zero."""
    return [
        (
            "zero-reaction-time",
            Severity.WARNING,
            "a reaction time is zero; leave it unset if it is not known",
        )
        for outcome in reaction.outcomes
        if outcome.HasField("reaction_time") and outcome.reaction_time.value == 0
    ]


def _check_unmeasured_products(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """A product recorded with no measurement."""
    return [
        (
            "product-without-measurement",
            Severity.INFO,
            f"product {_label(product)} has no measurement",
        )
        for product in _products(reaction)
        if not product.measurements
    ]


CHECKS: dict[str, Callable[[reaction_pb2.Reaction], list[Problem]]] = {
    "element-conservation": _check_element_conservation,
    "reaction-smiles-mismatch": _check_reaction_smiles,
    "radicals-and-charges": _check_radicals_and_charges,
    "stereo-unspecified": _check_stereo_in_name,
    "yield-sum": _check_yield_sum,
    "implausible-temperature": _check_temperature,
    "zero-reaction-time": _check_reaction_time,
    "product-without-measurement": _check_unmeasured_products,
}


def check_reaction(reaction: reaction_pb2.Reaction) -> list[Problem]:
    """Runs every chemistry check on one reaction.

    Args:
        reaction: The reaction; it is only read.

    Returns:
        What each check found, one entry per problem.
    """
    found: list[Problem] = []
    for check in CHECKS.values():
        found.extend(check(reaction))
    return found


@dataclasses.dataclass
class _Group:
    check: str
    severity: Severity
    message: str
    reactions: int = 0
    examples: list[str] = dataclasses.field(default_factory=list)


class ChemistryRollup:
    """Runs the chemistry checks over a dataset as it streams past, and groups them."""

    def __init__(self) -> None:
        """Starts an empty rollup."""
        self._groups: dict[tuple[str, str], _Group] = {}
        self.reaction_count = 0

    def observe(self, reaction: reaction_pb2.Reaction, label: str) -> None:
        """Checks one reaction and adds what it finds."""
        self.reaction_count += 1
        seen = set()
        for check, severity, message in check_reaction(reaction):
            key = (check, _DIGITS.sub("#", message))
            if key in seen:
                continue
            seen.add(key)
            group = self._groups.get(key)
            if group is None:
                group = self._groups[key] = _Group(check, severity, message)
            group.reactions += 1
            if len(group.examples) < MAX_EXAMPLES:
                group.examples.append(label)

    def watch(
        self, reactions: Iterable[reaction_pb2.Reaction]
    ) -> Iterator[reaction_pb2.Reaction]:
        """Yields ``reactions`` unchanged, checking each on the way past."""
        for position, reaction in enumerate(reactions):
            self.observe(reaction, reaction.reaction_id or f"reaction #{position + 1}")
            yield reaction

    def findings(self) -> list[Finding]:
        """What the checks found, most serious and most widespread first.

        Messages that differ only in numbers are one finding; its message is the
        first seen, and ``reactions`` counts every reaction that raised it.
        """
        groups = sorted(
            self._groups.values(),
            key=lambda g: (-g.severity, -g.reactions, g.check, g.message),
        )
        return [
            Finding(
                check=g.check,
                severity=g.severity,
                message=g.message
                + (
                    f" (and {g.reactions - 1:,} more reactions like it)"
                    if g.reactions > 1
                    else ""
                ),
                reactions=g.reactions,
                examples=g.examples,
            )
            for g in groups
        ]
