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

"""Tests for ord_schema.completeness."""

import dataclasses

import pytest

from ord_schema import completeness, message_helpers
from ord_schema.proto import reaction_pb2


def _complete_reaction() -> reaction_pb2.Reaction:
    """A reaction that records every detail the rubric asks for."""
    reaction = reaction_pb2.Reaction(reaction_id="r1")
    amine = message_helpers.build_compound(
        smiles="CCN", name="ethylamine", amount="1 mmol", role="reactant"
    )
    amine.is_limiting = True
    reaction.inputs["amine"].components.append(amine)
    reaction.inputs["amine"].addition_order = 1
    reaction.inputs["solvent"].components.append(
        message_helpers.build_compound(
            smiles="CS(C)=O", name="DMSO", amount="1 mL", role="solvent"
        )
    )
    reaction.inputs["solvent"].addition_order = 2
    conditions = reaction.conditions
    conditions.temperature.control.type = (
        reaction_pb2.TemperatureConditions.TemperatureControl.AMBIENT
    )
    conditions.pressure.atmosphere.type = (
        reaction_pb2.PressureConditions.Atmosphere.NITROGEN
    )
    conditions.stirring.type = reaction_pb2.StirringConditions.STIR_BAR
    reaction.workups.add(type=reaction_pb2.ReactionWorkup.CONCENTRATION)
    outcome = reaction.outcomes.add()
    outcome.reaction_time.value = 2
    outcome.reaction_time.units = reaction_pb2.Time.HOUR
    outcome.analyses["nmr"].type = reaction_pb2.Analysis.NMR_1H
    product = outcome.products.add(is_desired_product=True)
    product.identifiers.add(type="SMILES", value="CCNC(C)=O")
    product.identifiers.add(type="NAME", value="N-ethylacetamide")
    measurement = product.measurements.add(type="YIELD", analysis_key="nmr")
    measurement.percentage.value = 80
    reaction.provenance.doi = "10.1000/example"
    reaction.provenance.record_created.time.value = "2026-01-05"
    return reaction


def _item(
    score: completeness.CompletenessScore, item_id: str
) -> completeness.ItemScore:
    (found,) = [s for s in score.items if s.item.id == item_id]
    return found


def test_every_rubric_item_has_a_check_and_every_check_an_item():
    rubric = completeness.load_rubric()
    assert {item.id for item in rubric.items} == completeness.checked_items()


def test_the_rubric_is_well_formed():
    rubric = completeness.load_rubric()
    assert rubric.version
    assert rubric.profiles == {"uniform", "varied"}
    assert {item.level for item in rubric.items} <= set(rubric.levels)
    ids = [item.id for item in rubric.items]
    assert len(ids) == len(set(ids))
    for item in rubric.items:
        assert all(weight >= 0 for weight in item.weights.values())
        assert item.suggestion[0].islower()  # It follows "Suggestion: ".


def test_a_complete_reaction_scores_one_and_suggests_nothing():
    score = completeness.score_reaction(_complete_reaction())
    assert score.overall == 1
    assert set(score.levels.values()) == {1}
    assert score.suggestions() == []


def test_a_missing_detail_is_suggested():
    reaction = _complete_reaction()
    reaction.conditions.pressure.ClearField("atmosphere")
    score = completeness.score_reaction(reaction)
    assert _item(score, "atmosphere").score == 0
    assert score.suggestions() == [
        "Suggestion: record conditions.pressure.atmosphere.type"
    ]
    assert score.levels["conditions"] == pytest.approx(2 / 3)
    assert score.overall is not None
    assert score.overall < 1


def test_an_item_that_counts_says_how_many_have_it():
    reaction = _complete_reaction()
    reaction.inputs["solvent"].components[
        0
    ].reaction_role = reaction_pb2.ReactionRole.UNSPECIFIED
    score = completeness.score_reaction(reaction)
    assert _item(score, "input-roles").score == pytest.approx(1 / 2)
    assert (
        "Suggestion: give every input component a reaction_role (1 of 2 have it)"
        in score.suggestions()
    )


def test_items_that_do_not_apply_are_left_out():
    reaction = _complete_reaction()
    del reaction.outcomes[:]
    del reaction.inputs["solvent"]
    score = completeness.score_reaction(reaction)
    for item_id in ("measurement-type", "measurement-analysis", "desired-product"):
        assert _item(score, item_id).score is None
    # One input has nothing to be ordered against.
    assert _item(score, "addition-order").score is None
    assert score.levels["outcomes"] is None


def test_a_measurement_must_name_an_analysis_that_is_recorded():
    reaction = _complete_reaction()
    reaction.outcomes[0].products[0].measurements[0].analysis_key = "hplc"
    assert (
        _item(completeness.score_reaction(reaction), "measurement-analysis").score == 0
    )


def test_a_structure_without_a_name_is_counted():
    reaction = _complete_reaction()
    del reaction.inputs["amine"].components[0].identifiers[1]
    score = completeness.score_reaction(reaction)
    assert _item(score, "compound-names").have == 2
    assert _item(score, "compound-names").of == 3


def test_the_profile_chooses_the_weights():
    reaction = _complete_reaction()
    reaction.ClearField("workups")
    uniform = completeness.score_reaction(reaction, "uniform")
    varied = completeness.score_reaction(reaction, "varied")
    # Workups weigh less in a screen, so leaving them out costs less.
    assert uniform.overall is not None
    assert varied.overall is not None
    assert uniform.overall > varied.overall


def test_a_weight_of_zero_leaves_an_item_out():
    rubric = completeness.load_rubric()
    items = tuple(
        dataclasses.replace(item, weights={"uniform": 0.0})
        if item.id == "atmosphere"
        else dataclasses.replace(item, weights={"uniform": item.weights["uniform"]})
        for item in rubric.items
    )
    reaction = _complete_reaction()
    reaction.conditions.pressure.ClearField("atmosphere")
    score = completeness.score_reaction(
        reaction, "uniform", dataclasses.replace(rubric, items=items)
    )
    assert score.overall == 1
    assert score.suggestions() == []


def test_an_unknown_profile_is_refused():
    with pytest.raises(ValueError, match="no weights for profile 'bulk'"):
        completeness.score_reaction(_complete_reaction(), "bulk")


def test_an_item_with_no_check_is_refused():
    rubric = completeness.load_rubric()
    extra = dataclasses.replace(rubric.items[0], id="not-a-check")
    with pytest.raises(ValueError, match="no check for rubric item 'not-a-check'"):
        completeness.score_reaction(
            _complete_reaction(),
            rubric=dataclasses.replace(rubric, items=(*rubric.items, extra)),
        )


def test_scoring_does_not_change_the_reaction():
    reaction = _complete_reaction()
    reaction.conditions.pressure.ClearField("atmosphere")
    before = reaction.SerializeToString(deterministic=True)
    completeness.score_reaction(reaction)
    assert reaction.SerializeToString(deterministic=True) == before
