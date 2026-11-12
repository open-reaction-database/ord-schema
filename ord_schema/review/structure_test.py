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

"""Tests for ord_schema.review.structure."""

import itertools
import pathlib

import pytest

from ord_schema import datasets, message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.review import structure
from ord_schema.validations import Severity

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)
_UNIFORM = structure.ValidationProfile.UNIFORM
_VARIED = structure.ValidationProfile.VARIED


def _reaction(
    reaction_id: str,
    *,
    base: str = "CCN(CC)CC",
    halide: str = "Brc1ccccc1",
    base_amount: str = "1 mmol",
) -> reaction_pb2.Reaction:
    """A small screen reaction: a solvent, a base and an aryl halide."""
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.inputs["solvent"].components.append(
        message_helpers.build_compound(smiles="CS(C)=O", amount="1 mL", role="solvent")
    )
    reaction.inputs["base"].components.append(
        message_helpers.build_compound(smiles=base, amount=base_amount, role="reagent")
    )
    reaction.inputs["halide"].components.append(
        message_helpers.build_compound(smiles=halide, amount="1 mmol", role="reactant")
    )
    return reaction


def _findings(
    reactions: list[reaction_pb2.Reaction],
    check: str,
    profile: structure.ValidationProfile = _UNIFORM,
) -> list[structure.Finding]:
    report = structure.check_dataset(reactions, profile)
    return [finding for finding in report.findings if finding.check == check]


def test_a_consistently_recorded_screen_has_nothing_to_report():
    report = structure.check_dataset(datasets.load_dataset(_NIELSEN).reactions)
    assert report.findings == []
    assert report.not_applied == []


def test_varied_skips_the_template_checks_and_says_so():
    report = structure.check_dataset([_reaction("r1")], _VARIED)
    assert set(report.not_applied) == {
        "coverage-gap",
        "placeholder-value",
        "grid-gap",
        "duplicate-well",
    }
    assert set(report.applied) | set(report.not_applied) == {
        check.id for check in structure.CHECKS
    }


def test_check_ids_are_unique():
    ids = [check.id for check in structure.CHECKS]
    assert len(ids) == len(set(ids))


def _with_measurement(reaction: reaction_pb2.Reaction, value: float | None):
    product = reaction.outcomes.add().products.add()
    product.identifiers.add(type="SMILES", value="c1ccccc1N(CC)CC")
    if value is not None:
        measurement = product.measurements.add(type="YIELD")
        measurement.percentage.value = value
    return reaction


def test_fields_empty_in_the_same_reactions_are_one_coverage_gap():
    reactions = [
        _with_measurement(_reaction(f"r{i}"), None if i in (3, 7) else 10 + i)
        for i in range(10)
    ]
    (finding,) = _findings(reactions, "coverage-gap")
    assert finding.severity is Severity.WARNING
    assert finding.reactions == 2
    assert finding.examples == ["r3", "r7"]
    assert len(finding.fields) > 1


def test_a_field_half_the_reactions_fill_is_not_a_coverage_gap():
    reactions = [
        _with_measurement(_reaction(f"r{i}"), None if i % 2 else 10 + i)
        for i in range(10)
    ]
    assert _findings(reactions, "coverage-gap") == []


def test_coverage_gaps_are_not_checked_for_varied_datasets():
    reactions = [
        _with_measurement(_reaction(f"r{i}"), None if i == 3 else 10 + i)
        for i in range(10)
    ]
    assert _findings(reactions, "coverage-gap", _VARIED) == []


def test_a_repeated_zero_amount_is_a_placeholder():
    reactions = [_reaction(f"r{i}", base_amount="0 mmol") for i in range(3)]
    (finding,) = _findings(reactions, "placeholder-value")
    assert finding.fields == ["inputs.components.amount.moles.value"]


def test_placeholder_text_is_found_in_free_text_fields():
    reactions = [_reaction(f"r{i}") for i in range(3)]
    for reaction in reactions:
        reaction.notes.procedure_details = "N/A"
    (finding,) = _findings(reactions, "placeholder-value")
    assert finding.fields == ["notes.procedure_details"]


def test_an_enum_none_is_a_real_value_not_a_placeholder():
    # stirring.type NONE means the reaction was not stirred.
    reactions = [_reaction(f"r{i}") for i in range(3)]
    for reaction in reactions:
        reaction.conditions.stirring.type = reaction_pb2.StirringConditions.NONE
    assert _findings(reactions, "placeholder-value") == []


_BASES = ["CCN(CC)CC", "CN(C)C", "C1CCNCC1"]
_HALIDES = ["Brc1ccccc1", "Clc1ccccc1", "Ic1ccccc1", "Brc1ccncc1"]


def test_a_missing_combination_in_a_screen_is_a_grid_gap():
    cells = list(itertools.product(_BASES, _HALIDES))
    reactions = [
        _reaction(f"r{i}", base=base, halide=halide)
        for i, (base, halide) in enumerate(cells)
        if (base, halide) != ("CN(C)C", "Ic1ccccc1")
    ]
    (finding,) = _findings(reactions, "grid-gap")
    assert finding.severity is Severity.INFO
    assert "1 of the 12 combinations" in finding.message
    (example,) = finding.examples
    assert set(example.split(" | ")) == {"CN(C)C", "Ic1ccccc1"}


def test_a_full_grid_has_no_gap():
    reactions = [
        _reaction(f"r{i}", base=base, halide=halide)
        for i, (base, halide) in enumerate(itertools.product(_BASES, _HALIDES))
    ]
    assert _findings(reactions, "grid-gap") == []


def test_inputs_used_in_fixed_pairs_are_one_factor_not_a_sparse_grid():
    # As in the Perera screen: bases and halides are used in seven fixed pairs, not
    # every base with every halide, and neither decides the other. The pairs are
    # crossed with the solvent. Read as separate factors they would fill 7 of 12
    # cells, enough to look like a grid with 10 gaps; they are one factor, and the
    # grid is complete.
    pairs = [
        (_BASES[0], _HALIDES[0]),
        (_BASES[0], _HALIDES[1]),
        (_BASES[0], _HALIDES[2]),
        (_BASES[1], _HALIDES[1]),
        (_BASES[1], _HALIDES[2]),
        (_BASES[2], _HALIDES[2]),
        (_BASES[2], _HALIDES[3]),
    ]
    reactions = []
    for i, ((base, halide), solvent) in enumerate(
        itertools.product(pairs, ["CS(C)=O", "CN(C)C=O"])
    ):
        reaction = _reaction(f"r{i}", base=base, halide=halide)
        reaction.inputs["solvent"].components[0].identifiers[0].value = solvent
        reactions.append(reaction)
    assert _findings(reactions, "grid-gap") == []


def test_an_aggregate_with_too_many_candidate_factors_is_not_read_as_a_grid(
    monkeypatch,
):
    monkeypatch.setattr(structure, "MAX_FACTOR_CANDIDATES", 1)
    cells = list(itertools.product(_BASES, _HALIDES))
    reactions = [
        _reaction(f"r{i}", base=base, halide=halide)
        for i, (base, halide) in enumerate(cells)
        if (base, halide) != ("CN(C)C", "Ic1ccccc1")
    ]
    assert _findings(reactions, "grid-gap") == []


def test_two_reactions_in_one_well_are_flagged():
    reactions = [_reaction(f"r{i}") for i in range(3)]
    for reaction, position in zip(reactions, ["A1", "A2", "A1"], strict=True):
        reaction.setup.vessel.vessel_id = "plate 1"
        reaction.setup.vessel.position = position
    (finding,) = _findings(reactions, "duplicate-well")
    assert finding.examples == ["r0", "r2"]
    assert "plate 1" in finding.message


def test_the_same_position_on_different_plates_is_fine():
    reactions = [_reaction(f"r{i}") for i in range(2)]
    for reaction, plate in zip(reactions, ["plate 1", "plate 2"], strict=True):
        reaction.setup.vessel.vessel_id = plate
        reaction.setup.vessel.row = "B"
        reaction.setup.vessel.col = "3"
    assert _findings(reactions, "duplicate-well") == []


def test_one_position_in_two_units_is_mixed():
    reactions = [
        _reaction("r1", base_amount="1 mmol"),
        _reaction("r2", base_amount="2 umol"),
    ]
    (finding,) = _findings(reactions, "mixed-units")
    assert "MILLIMOLE (1)" in finding.message
    assert "MICROMOLE (1)" in finding.message


def test_one_position_as_mass_in_some_reactions_and_moles_in_others_is_mixed():
    reactions = [
        _reaction("r1", base_amount="1 mmol"),
        _reaction("r2", base_amount="5 mg"),
    ]
    findings = _findings(reactions, "mixed-units")
    assert any("given as" in finding.message for finding in findings)


def test_one_name_for_two_structures_is_a_conflict():
    reactions = [_reaction("r1"), _reaction("r2", base="CN(C)C")]
    for reaction in reactions:
        reaction.inputs["base"].components[0].identifiers.add(
            type="NAME", value="amine"
        )
    (finding,) = _findings(reactions, "name-structure-conflict")
    assert finding.severity is Severity.WARNING
    assert "'amine'" in finding.message


def test_one_structure_under_two_names_is_noted():
    reactions = [_reaction("r1"), _reaction("r2")]
    for reaction, name in zip(reactions, ["Et3N", "triethylamine"], strict=True):
        reaction.inputs["base"].components[0].identifiers.add(type="NAME", value=name)
    (finding,) = _findings(reactions, "structure-name-variants")
    assert finding.severity is Severity.INFO


def test_a_compound_in_two_roles_is_noted():
    reactions = [_reaction("r1"), _reaction("r2")]
    reactions[1].inputs["base"].components[
        0
    ].reaction_role = reaction_pb2.ReactionRole.REACTANT
    (finding,) = _findings(reactions, "conflicting-roles")
    assert "REAGENT (1)" in finding.message
    assert "REACTANT (1)" in finding.message


def test_compound_by_compound_findings_are_capped_with_a_summary(monkeypatch):
    monkeypatch.setattr(structure, "MAX_COMPOUND_FINDINGS", 2)
    findings = [
        structure.Finding(
            check="conflicting-roles",
            severity=Severity.INFO,
            message=f"compound {n} appears as ...",
            reactions=n,
        )
        for n in (3, 9, 1, 5)
    ]
    capped = structure._capped(findings, "compounds in more than one role")
    assert [f.reactions for f in capped] == [9, 5, 4]
    summary = capped[-1]
    assert (
        summary.message
        == "and 2 more compounds in more than one role, each in fewer reactions"
    )
    assert summary.examples == [
        "compound 3 appears as ...",
        "compound 1 appears as ...",
    ]
    assert structure._capped(findings[:2], "x") == findings[:2]


def test_reactions_identical_but_for_ids_and_record_times_are_duplicates():
    reactions = [_reaction("r1"), _reaction("r2"), _reaction("r3", base="CN(C)C")]
    reactions[0].provenance.record_created.time.value = "2026-01-05"
    reactions[1].provenance.record_created.time.value = "2026-01-06"
    (finding,) = _findings(reactions, "duplicate-reaction")
    assert finding.reactions == 2
    assert finding.examples == ["r1 = r2"]


def test_a_dataset_where_only_some_reactions_cite_a_source_is_flagged():
    reactions = [_reaction(f"r{i}") for i in range(3)]
    reactions[0].provenance.doi = "10.1000/example"
    reactions[1].provenance.patent = "US1234567"
    (finding,) = _findings(reactions, "mixed-sourcing")
    assert finding.reactions == 1
    assert finding.examples == ["r2"]


@pytest.mark.parametrize("cited", [0, 3])
def test_all_or_none_citing_a_source_is_consistent(cited):
    reactions = [_reaction(f"r{i}") for i in range(3)]
    for reaction in reactions[:cited]:
        reaction.provenance.doi = "10.1000/example"
    assert _findings(reactions, "mixed-sourcing") == []


def test_findings_are_ordered_most_serious_first():
    reactions = [_reaction("r1"), _reaction("r2")]
    reactions[1].inputs["base"].components[
        0
    ].reaction_role = reaction_pb2.ReactionRole.REACTANT
    reactions[0].provenance.doi = "10.1000/example"
    report = structure.check_dataset(reactions)
    severities = [finding.severity for finding in report.findings]
    assert severities == sorted(severities, reverse=True)
    assert report.by_check()["mixed-sourcing"] == 1
