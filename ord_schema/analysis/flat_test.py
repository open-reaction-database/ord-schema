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

"""Tests for ord_schema.analysis.flat."""

import math
import pathlib

import pandas as pd

from ord_schema import datasets, message_helpers
from ord_schema.analysis import flat
from ord_schema.proto import reaction_pb2

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)


def _reaction(
    reaction_id: str, *, base: str = "K2CO3", temperature: float | None = None
) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.inputs["solvent"].components.append(
        message_helpers.build_compound(smiles="CCO", amount="1 mL", role="solvent")
    )
    reaction.inputs["base"].components.append(
        message_helpers.build_compound(name=base, role="reagent")
    )
    if temperature is not None:
        setpoint = reaction.conditions.temperature.setpoint
        setpoint.value = temperature
        setpoint.units = reaction_pb2.Temperature.CELSIUS
    return reaction


_BASE = 'inputs["base"].components[0].identifiers[0].value'
_TEMPERATURE = "conditions.temperature.setpoint.value"


def test_varying_fields_become_columns_and_constant_ones_are_listed():
    table = flat.flatten(
        [
            _reaction("r1", base="K2CO3", temperature=25),
            _reaction("r2", base="Cs2CO3", temperature=25),
        ]
    )
    assert list(table.varying.columns) == ["reaction_id", _BASE]
    assert table.varying[_BASE].tolist() == ["K2CO3", "Cs2CO3"]
    assert table.constants[_TEMPERATURE] == 25
    assert table.constants["conditions.temperature.setpoint.units"] == "CELSIUS"
    assert table.constants['inputs["solvent"].components[0].reaction_role'] == "SOLVENT"
    assert table.reaction_count == 2


def test_a_field_some_reactions_leave_empty_varies():
    table = flat.flatten([_reaction("r1", temperature=25), _reaction("r2")])
    assert _TEMPERATURE in table.varying.columns
    assert _TEMPERATURE not in table.constants
    first, second = table.varying[_TEMPERATURE].tolist()
    assert first == 25
    assert math.isnan(second)


def test_components_in_a_different_order_show_as_varying_columns():
    def reaction(reaction_id: str, first: str, second: str) -> reaction_pb2.Reaction:
        result = reaction_pb2.Reaction(reaction_id=reaction_id)
        result.inputs["mixture"].components.extend(
            [
                message_helpers.build_compound(name=first),
                message_helpers.build_compound(name=second),
            ]
        )
        return result

    table = flat.flatten([reaction("r1", "A", "B"), reaction("r2", "B", "A")])
    assert {
        'inputs["mixture"].components[0].identifiers[0].value',
        'inputs["mixture"].components[1].identifiers[0].value',
    } <= set(table.varying.columns)


def test_matches_messages_to_dataframe_on_a_real_screen():
    # The helper this builds on: the same columns, values and dtypes, without
    # holding the constant columns in memory.
    reactions = datasets.load_dataset(_NIELSEN).reactions
    table = flat.flatten(reactions)
    expected = message_helpers.messages_to_dataframe(
        reactions, drop_constant_columns=True
    )
    pd.testing.assert_frame_equal(table.varying, expected)
    every_column = message_helpers.messages_to_dataframe(reactions)
    assert set(table.constants) == set(every_column.columns) - set(expected.columns)
    for name, value in table.constants.items():
        assert every_column[name].iloc[0] == value


def test_a_real_screen_varies_in_its_design_variables():
    table = flat.flatten(datasets.load_dataset(_NIELSEN).reactions)
    assert table.reaction_count == 80
    assert list(table.varying.columns) == [
        'inputs["alcohol in THF"].components[0].identifiers[0].value',
        'inputs["base"].components[0].identifiers[0].value',
        'inputs["sulfonyl fluoride solution"].components[0].identifiers[0].value',
        "outcomes[0].products[0].identifiers[0].value",
        "outcomes[0].products[0].measurements[1].percentage.value",
    ]


def test_a_generator_is_read_once():
    table = flat.flatten(
        reaction for reaction in (_reaction("r1"), _reaction("r2", base="Cs2CO3"))
    )
    assert table.varying[_BASE].tolist() == ["K2CO3", "Cs2CO3"]


def test_a_column_that_varies_late_keeps_its_earlier_values():
    # Constant for three reactions, then different: the earlier rows are rebuilt
    # from the one value they shared, and a reaction without the column stays empty.
    reactions = [
        _reaction("r1", temperature=25),
        _reaction("r2"),
        _reaction("r3", temperature=25),
        _reaction("r4", temperature=60),
    ]
    values = flat.flatten(reactions).varying[_TEMPERATURE].tolist()
    assert values[0] == 25
    assert math.isnan(values[1])
    assert values[2:] == [25, 60]


def test_matches_messages_to_dataframe_with_reaction_ids():
    reactions = [
        _reaction(f"r{i}", base=base, temperature=t)
        for i, (base, t) in enumerate(
            [("K2CO3", 25), ("Cs2CO3", None), ("K2CO3", 60), ("K3PO4", 25)]
        )
    ]
    table = flat.flatten(reactions)
    expected = message_helpers.messages_to_dataframe(
        reactions, drop_constant_columns=True
    )
    pd.testing.assert_frame_equal(table.varying, expected[table.varying.columns])
    assert set(table.varying.columns) == set(expected.columns)


def test_no_reactions_gives_an_empty_table():
    table = flat.flatten([])
    assert table.varying.empty
    assert table.constants == {}
    assert table.reaction_count == 0


def test_csvs(tmp_path):
    table = flat.flatten(
        [_reaction("r1", base="K2CO3"), _reaction("r2", base="Cs2CO3")]
    )
    table.write_csv(tmp_path / "flat.csv")
    table.write_constants_csv(tmp_path / "constants.csv")
    written = pd.read_csv(tmp_path / "flat.csv")
    assert written[_BASE].tolist() == ["K2CO3", "Cs2CO3"]
    constants = pd.read_csv(tmp_path / "constants.csv")
    assert list(constants.columns) == ["field", "value"]
    assert "SOLVENT" in constants["value"].tolist()
