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

"""Tests for ord_schema.analysis.summary."""

import pathlib

import pytest

from ord_schema import datasets, message_helpers
from ord_schema.analysis import summary
from ord_schema.proto import reaction_pb2

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)


def _reaction(
    reaction_id: str, *, temperature: float | None = None, solvent: str = "CCO"
) -> reaction_pb2.Reaction:
    """Returns a small reaction: one solvent input, optionally a temperature."""
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.inputs["solvent"].components.append(
        message_helpers.build_compound(smiles=solvent, amount="1 mL", role="solvent")
    )
    if temperature is not None:
        setpoint = reaction.conditions.temperature.setpoint
        setpoint.value = temperature
        setpoint.units = reaction_pb2.Temperature.CELSIUS
    return reaction


def _summary(*reactions: reaction_pb2.Reaction) -> summary.DatasetSummary:
    result = summary.DatasetSummary()
    result.observe_all(reactions)
    return result


def test_paths_name_fields_not_positions():
    reaction = _reaction("r1")
    reaction.inputs["base"].components.append(
        message_helpers.build_compound(name="K2CO3", role="reagent")
    )
    result = _summary(reaction)
    # Two inputs and two components, but one path each, counted once per reaction.
    assert result.field_counts["inputs.components"] == 1
    assert result.field_counts["inputs.components.identifiers.value"] == 1
    assert result.values["inputs.components.reaction_role"] == {
        "SOLVENT": 1,
        "REAGENT": 1,
    }
    assert result.map_keys["inputs"] == {"solvent": 1, "base": 1}


def test_enums_are_counted_by_name():
    result = _summary(_reaction("r1", temperature=25))
    assert result.values["conditions.temperature.setpoint.units"] == {"CELSIUS": 1}


def test_coverage_counts_reactions_not_occurrences():
    result = _summary(
        _reaction("r1", temperature=25),
        _reaction("r2", temperature=25),
        _reaction("r3"),
    )
    assert result.coverage("conditions.temperature.setpoint.value") == pytest.approx(
        2 / 3
    )
    assert result.coverage("inputs.components") == 1
    assert result.coverage("not.a.path") == 0


def test_constant_and_varying_fields():
    result = _summary(
        _reaction("r1", temperature=25),
        _reaction("r2", temperature=60),
    )
    constant = result.constant_fields()
    assert constant["conditions.temperature.setpoint.units"] == ("CELSIUS",)
    assert constant["inputs.components.amount.volume.value"] == (1,)
    varying = result.varying_fields()
    assert "conditions.temperature.setpoint.value" in varying
    assert "reaction_id" in varying
    assert "conditions.temperature.setpoint.units" not in varying


def test_values_that_differ_within_a_reaction_but_not_between_are_constant():
    # Every reaction has a solvent and a reagent: the roles differ between components
    # but not between reactions, so the field is constant.
    reactions = []
    for reaction_id in ("r1", "r2"):
        reaction = _reaction(reaction_id)
        reaction.inputs["base"].components.append(
            message_helpers.build_compound(name="K2CO3", role="reagent")
        )
        reactions.append(reaction)
    result = _summary(*reactions)
    path = "inputs.components.reaction_role"
    assert result.values[path] == {"SOLVENT": 2, "REAGENT": 2}
    assert result.constant_fields()[path] == ("REAGENT", "SOLVENT")
    assert path not in result.varying_fields()


def test_a_field_some_reactions_leave_empty_is_not_constant():
    result = _summary(_reaction("r1", temperature=25), _reaction("r2"))
    assert not result.is_constant("conditions.temperature.setpoint.value")
    assert "conditions.temperature.setpoint.value" in result.varying_fields()


def test_reactions_recorded_the_same_way_share_a_shape():
    result = _summary(
        _reaction("r1", temperature=25),
        _reaction("r2", temperature=60, solvent="CO"),
        _reaction("r3"),
    )
    (common, common_count), (other, other_count) = result.ordered_shapes()
    # Different values, same fields: one shape.
    assert common_count == 2
    assert result.shape_examples[common] == ["r1", "r2"]
    assert other_count == 1
    assert common - other == {
        "conditions",
        "conditions.temperature",
        "conditions.temperature.setpoint",
        "conditions.temperature.setpoint.value",
        "conditions.temperature.setpoint.units",
    }


def test_the_summary_carries_the_compound_table():
    result = _summary(_reaction("r1"), _reaction("r2", solvent="CO"))
    assert {record.smiles for record in result.compounds.records()} == {"CCO", "CO"}


def test_distinct_values_are_capped_and_the_rest_counted(monkeypatch):
    monkeypatch.setattr(summary, "MAX_DISTINCT_VALUES", 2)
    result = _summary(*(_reaction(f"r{i}", temperature=i) for i in range(5)))
    path = "conditions.temperature.setpoint.value"
    assert len(result.values[path]) == 2
    assert result.untracked[path] == 3
    assert result.reaction_values_untracked[path] == 3
    # Losing track of values must not make a varying field look constant.
    assert not result.is_constant(path)


def test_long_text_is_cut_before_it_is_counted():
    reaction = _reaction("r1")
    reaction.notes.procedure_details = "x" * 500
    result = _summary(reaction)
    (value,) = result.values["notes.procedure_details"]
    assert isinstance(value, str)
    assert len(value) == summary.MAX_VALUE_LENGTH + 1


def test_merging_slices_matches_observing_everything():
    reactions = [
        _reaction("r1", temperature=25),
        _reaction("r2", temperature=60, solvent="CO"),
        _reaction("r3"),
        _reaction("r4", temperature=25),
    ]
    whole = _summary(*reactions)
    merged = _summary(*reactions[:2])
    merged.merge(_summary(*reactions[2:]))
    assert merged.reaction_count == whole.reaction_count
    assert merged.field_counts == whole.field_counts
    assert merged.values == whole.values
    assert merged.reaction_values == whole.reaction_values
    assert merged.ordered_values == whole.ordered_values
    assert merged.map_keys == whole.map_keys
    assert merged.shapes == whole.shapes
    assert list(merged.compounds.rows()) == list(whole.compounds.rows())


def test_an_empty_summary_answers_without_dividing_by_zero():
    result = summary.DatasetSummary()
    assert result.coverage("reaction_id") == 0
    assert result.constant_fields() == {}
    assert result.varying_fields() == []
    assert result.ordered_shapes() == []


def test_a_real_screen_has_one_dominant_shape():
    dataset = datasets.load_dataset(_NIELSEN)
    result = summary.DatasetSummary()
    result.observe_all(dataset.reactions)
    assert result.reaction_count == 80
    # A screen records every reaction the same way; its variables are what vary.
    _, count = result.ordered_shapes()[0]
    assert count / result.reaction_count > 0.5
    assert result.varying_fields()
    assert result.constant_fields()


def _two_component_reaction(
    reaction_id: str, first: reaction_pb2.Compound, second: reaction_pb2.Compound
) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.inputs["mixture"].components.extend([first, second])
    return reaction


def test_the_same_components_in_a_different_order_are_flagged():
    base = message_helpers.build_compound(name="K2CO3", role="reagent")
    solvent = message_helpers.build_compound(smiles="CCO", role="solvent")
    result = _summary(
        _two_component_reaction("r1", base, solvent),
        _two_component_reaction("r2", base, solvent),
        _two_component_reaction("r3", solvent, base),
    )
    path = "inputs.components.reaction_role"
    # Not a variable of the dataset, but not recorded consistently either.
    assert path not in result.varying_fields()
    assert result.reordered_fields()[path] == 1


def test_inputs_added_in_a_different_order_are_not_reordered():
    # A map has no defined order, so adding the same inputs in another order records
    # the same reaction.
    first = _reaction("r1")
    first.inputs["base"].components.append(
        message_helpers.build_compound(name="K2CO3", role="reagent")
    )
    second = reaction_pb2.Reaction(reaction_id="r2")
    second.inputs["base"].components.append(
        message_helpers.build_compound(name="K2CO3", role="reagent")
    )
    second.inputs["solvent"].components.append(
        message_helpers.build_compound(smiles="CCO", amount="1 mL", role="solvent")
    )
    assert _summary(first, second).reordered_fields() == {}


def test_a_renamed_input_key_shows_as_reordered():
    def reaction(reaction_id: str, base_key: str) -> reaction_pb2.Reaction:
        result = _reaction(reaction_id)
        result.inputs[base_key].components.append(
            message_helpers.build_compound(name="K2CO3", role="reagent")
        )
        return result

    # "base" sorts before "solvent" and "z base" after it, so the same components
    # are read in a different order.
    result = _summary(reaction("r1", "base"), reaction("r2", "z base"))
    assert "inputs.components.reaction_role" in result.reordered_fields()
    assert result.map_keys["inputs"] == {"solvent": 2, "base": 1, "z base": 1}


def test_a_consistently_recorded_screen_has_nothing_reordered():
    dataset = datasets.load_dataset(_NIELSEN)
    result = summary.DatasetSummary()
    result.observe_all(dataset.reactions)
    assert result.reordered_fields() == {}
