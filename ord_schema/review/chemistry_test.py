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

"""Tests for ord_schema.review.chemistry."""

import csv
import json
import pathlib

import pytest

from ord_schema import datasets, message_helpers
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.review import chemistry, report
from ord_schema.validations import Severity

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)


def _reaction(
    *inputs: str, product: str = "CC(=O)Nc1ccccc1", name: str | None = None
) -> reaction_pb2.Reaction:
    """A reaction with each input SMILES in its own input, and one product."""
    reaction = reaction_pb2.Reaction(reaction_id="r1")
    for index, smiles in enumerate(inputs or ("Nc1ccccc1", "CC(=O)OC(C)=O")):
        reaction.inputs[f"input {index}"].components.append(
            message_helpers.build_compound(smiles=smiles, role="reactant")
        )
    outcome = reaction.outcomes.add()
    outcome.reaction_time.value = 1
    outcome.reaction_time.units = reaction_pb2.Time.HOUR
    compound = outcome.products.add()
    compound.identifiers.add(type="SMILES", value=product)
    if name:
        compound.identifiers.add(type="NAME", value=name)
    compound.measurements.add(type="YIELD").percentage.value = 80
    return reaction


def _checks(reaction: reaction_pb2.Reaction) -> list[str]:
    return [check for check, _, _ in chemistry.check_reaction(reaction)]


def test_a_sound_reaction_raises_nothing():
    assert chemistry.check_reaction(_reaction()) == []


def test_a_product_element_no_input_brings_is_flagged():
    (problem,) = chemistry.check_reaction(_reaction(product="Brc1ccccc1NC(C)=O"))
    check, severity, message = problem
    assert (check, severity) == ("element-conservation", Severity.WARNING)
    assert "contains Br" in message


def test_an_element_brought_in_a_workup_is_conserved():
    reaction = _reaction(product="Brc1ccccc1NC(C)=O")
    reaction.workups.add().input.components.append(
        message_helpers.build_compound(smiles="BrBr", role="workup")
    )
    assert "element-conservation" not in _checks(reaction)


def test_conservation_is_not_judged_when_an_input_has_no_structure():
    reaction = _reaction(product="Brc1ccccc1NC(C)=O")
    reaction.inputs["mystery"].components.append(
        message_helpers.build_compound(name="brominating agent")
    )
    assert "element-conservation" not in _checks(reaction)


def _with_reaction_smiles(reaction: reaction_pb2.Reaction, smiles: str):
    reaction.identifiers.add(type="REACTION_SMILES", value=smiles)
    return reaction


def test_a_stated_reaction_smiles_that_matches_is_fine_even_when_mapped():
    reaction = _with_reaction_smiles(
        _reaction(), "[NH2:1]c1ccccc1.CC(=O)OC(C)=O>>CC(=O)[NH:1]c1ccccc1"
    )
    assert "reaction-smiles-mismatch" not in _checks(reaction)


def test_a_stated_reaction_smiles_with_another_product_is_flagged():
    reaction = _with_reaction_smiles(
        _reaction(), "Nc1ccccc1.CC(=O)OC(C)=O>>CC(=O)Nc1ccccn1"
    )
    problems = [
        p
        for p in chemistry.check_reaction(reaction)
        if p[0] == "reaction-smiles-mismatch"
    ]
    (problem,) = problems
    assert "CC(=O)Nc1ccccn1" in problem[2]
    assert "which no outcome records" in problem[2]


def test_a_stated_reactant_no_input_records_is_flagged():
    reaction = _with_reaction_smiles(_reaction(), "Nc1ccccc1.CC(=O)Cl>>CC(=O)Nc1ccccc1")
    (problem,) = [
        p
        for p in chemistry.check_reaction(reaction)
        if p[0] == "reaction-smiles-mismatch"
    ]
    assert "CC(=O)Cl" in problem[2]


def test_an_unreadable_reaction_smiles_is_left_to_validation():
    reaction = _with_reaction_smiles(_reaction(), "not a reaction")
    assert "reaction-smiles-mismatch" not in _checks(reaction)


def test_a_radical_is_flagged_but_an_open_shell_metal_is_not():
    # Dimethoxy(methyl)silane without its Si-H, as the Liu example records it.
    assert "radical" in _checks(_reaction("CO[Si](C)OC", product="CO[SiH](C)OC"))
    copper_acetate = "CC(=O)[O-].CC(=O)[O-].[Cu+2]"
    assert "radical" not in _checks(_reaction(copper_acetate, product="CC(=O)O"))


def test_an_ion_without_its_counterion_is_flagged():
    problems = chemistry.check_reaction(_reaction("[NH4+]", product="N"))
    assert ("net-charge", Severity.WARNING) in [(c, s) for c, s, _ in problems]
    assert "net-charge" not in _checks(_reaction("[NH4+].[Cl-]", product="N"))


@pytest.mark.parametrize(
    ("name", "smiles", "flagged"),
    [
        ("(R)-1-phenylethanol", "CC(O)c1ccccc1", True),
        ("(R)-1-phenylethanol", "C[C@@H](O)c1ccccc1", False),
        ("rac-1-phenylethanol", "CC(O)c1ccccc1", False),
        ("(E)-stilbene", "C(=Cc1ccccc1)c1ccccc1", True),
        ("(E)-stilbene", "C(/c1ccccc1)=C/c1ccccc1", False),
        ("L-alanine", "CC(N)C(=O)O", True),
        ("1-phenylethanol", "CC(O)c1ccccc1", False),
    ],
)
def test_a_stereo_name_needs_a_stereo_structure(name, smiles, flagged):
    reaction = _reaction("c1ccccc1", "CC(O)C(N)C(=O)O", product=smiles, name=name)
    assert ("stereo-unspecified" in _checks(reaction)) is flagged


def test_yields_adding_up_past_100_percent_are_flagged():
    reaction = _reaction()
    second = reaction.outcomes[0].products.add()
    second.identifiers.add(type="SMILES", value="CC(=O)O")
    second.measurements.add(type="YIELD").percentage.value = 30
    assert "yield-sum" in _checks(reaction)
    reaction.outcomes[0].products[0].measurements[0].percentage.value = 70
    assert "yield-sum" not in _checks(reaction)


@pytest.mark.parametrize(
    ("value", "units", "flagged"),
    [
        (1000, "CELSIUS", True),
        (-196, "CELSIUS", False),
        (298, "KELVIN", False),
        (1500, "FAHRENHEIT", True),
    ],
)
def test_implausible_temperatures_are_flagged(value, units, flagged):
    reaction = _reaction()
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = value
    setpoint.units = reaction_pb2.Temperature.TemperatureUnit.Value(units)
    assert ("implausible-temperature" in _checks(reaction)) is flagged


def test_a_zero_reaction_time_is_flagged():
    reaction = _reaction()
    reaction.outcomes[0].reaction_time.value = 0
    assert "zero-reaction-time" in _checks(reaction)


def test_a_product_with_no_measurement_is_noted():
    reaction = _reaction()
    del reaction.outcomes[0].products[0].measurements[:]
    (problem,) = chemistry.check_reaction(reaction)
    assert problem[:2] == ("product-without-measurement", Severity.INFO)


def test_checking_does_not_change_the_reaction():
    reaction = _with_reaction_smiles(
        _reaction(product="Brc1ccccc1NC(C)=O"),
        "[NH2:1]c1ccccc1>>Brc1ccccc1[NH:1]C(C)=O",
    )
    before = reaction.SerializeToString(deterministic=True)
    chemistry.check_reaction(reaction)
    assert reaction.SerializeToString(deterministic=True) == before


def test_the_rollup_counts_each_problem_once_per_reaction():
    rollup = chemistry.ChemistryRollup()
    reactions = []
    for i in range(3):
        reaction = _reaction(product="Brc1ccccc1NC(C)=O")
        reaction.reaction_id = f"r{i}"
        reaction.conditions.temperature.setpoint.value = 1000 + i
        reaction.conditions.temperature.setpoint.units = (
            reaction_pb2.Temperature.CELSIUS
        )
        reactions.append(reaction)
    for _ in rollup.watch(reactions):
        pass
    findings = {f.check: f for f in rollup.findings()}
    assert findings["element-conservation"].reactions == 3
    assert findings["element-conservation"].examples == ["r0", "r1", "r2"]
    # The temperatures differ only in their numbers, so they are one finding.
    assert findings["implausible-temperature"].reactions == 3
    assert "and 2 more reactions like it" in findings["implausible-temperature"].message


def test_each_check_keeps_its_most_widespread_findings_and_sums_up_the_rest(
    monkeypatch,
):
    monkeypatch.setattr(chemistry, "MAX_GROUPS_PER_CHECK", 2)
    rollup = chemistry.ChemistryRollup()
    products = ["CCO", "CCCO", "CCCCO", "CCCCCO"]
    for i, smiles in enumerate(products):
        for j in range(4 - i):  # CCO in 4 reactions, CCCO in 3, ...
            reaction = _reaction(product=smiles)
            del reaction.outcomes[0].products[0].measurements[:]
            rollup.observe(reaction, f"r{i}-{j}")
    findings = [
        f for f in rollup.findings() if f.check == "product-without-measurement"
    ]
    assert [f.reactions for f in findings] == [4, 3, 10]
    summary = findings[-1]
    assert summary.message.startswith("and 2 more findings like these")
    assert "10 reactions raise this check in all" in summary.message
    assert summary.examples == [
        "product CCCCO has no measurement",
        "product CCCCCO has no measurement",
    ]


def test_a_consistent_real_screen_raises_nothing():
    rollup = chemistry.ChemistryRollup()
    for _ in rollup.watch(datasets.load_dataset(_NIELSEN).reactions):
        pass
    assert rollup.findings() == []


def test_the_report_and_packet_carry_the_chemistry_findings(tmp_path):
    reaction = _reaction(product="Brc1ccccc1NC(C)=O")
    review = report.build_report(dataset_pb2.Dataset(name="x", reactions=[reaction]))
    page = report.write_bundle(review, tmp_path).read_text()
    assert "<h2>Chemistry checks</h2>" in page
    with (tmp_path / "tables" / "findings.csv").open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["source"] == "chemistry"]
    assert [r["check"] for r in rows] == ["element-conservation"]
    data = json.loads((tmp_path / "packet" / "review.json").read_text())
    assert [f["check"] for f in data["chemistry_findings"]] == ["element-conservation"]
    assert review.severity_counts()[Severity.WARNING] >= 1


def test_a_crude_input_could_bring_the_element():
    # Step 2 uses step 1's crude product, which carries the bromine.
    reaction = _reaction("OB(O)c1ccccc1", product="Brc1ccc(-c2ccccc2)cc1")
    reaction.inputs["crude"].crude_components.add(reaction_id="step-1")
    assert "element-conservation" not in _checks(reaction)


def test_a_yield_measured_two_ways_counts_once():
    reaction = _reaction()
    product = reaction.outcomes[0].products[0]
    del product.measurements[:]
    product.measurements.add(type="YIELD", analysis_key="nmr").percentage.value = 82
    product.measurements.add(
        type="YIELD", analysis_key="isolated"
    ).percentage.value = 74
    assert "yield-sum" not in _checks(reaction)


def test_every_check_is_keyed_by_the_id_its_findings_carry():
    assert {"radical", "net-charge"} <= set(chemistry.CHECKS)
    assert "radicals-and-charges" not in chemistry.CHECKS
