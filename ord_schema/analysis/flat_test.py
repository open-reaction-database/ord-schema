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
import pytest

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
    # Dates as written: the helper does not bin them.
    table = flat.flatten(reactions, date_resolution=None)
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
    assert list(constants.columns) == ["field", "value", "section"]
    assert "SOLVENT" in constants["value"].tolist()


def _named(smiles: str, name: str, role: str = "reagent") -> reaction_pb2.Compound:
    return message_helpers.build_compound(smiles=smiles, name=name, role=role)


def _screen(*bases: tuple[str, str]) -> list[reaction_pb2.Reaction]:
    """One reaction per (SMILES, name) base, with a solvent and a yield."""
    reactions = []
    for i, (smiles, name) in enumerate(bases):
        reaction = reaction_pb2.Reaction(reaction_id=f"r{i}")
        reaction.inputs["solvent"].components.append(
            message_helpers.build_compound(smiles="CCO", role="solvent")
        )
        reaction.inputs["base"].components.append(_named(smiles, name))
        product = reaction.outcomes.add().products.add()
        product.measurements.add(type="YIELD").percentage.value = 10 * i + 5
        reactions.append(reaction)
    return reactions


_BASE_SMILES = 'inputs["base"].components[0].identifiers[0].value'
_BASE_NAME = 'inputs["base"].components[0].identifiers[1].value'
_YIELD = "outcomes[0].products[0].measurements[0].percentage.value"


def test_a_compounds_smiles_and_name_group_as_one_variable():
    table = flat.flatten(
        _screen(
            ("[K+].[K+].[O-]C([O-])=O", "K2CO3"),
            ("CCN(CC)CC", "Et3N"),
            ("CCN(CC)CC", "Et3N"),
        )
    )
    assert [_BASE_SMILES, _BASE_NAME] in table.column_groups
    assert _BASE_NAME not in table.grouped().columns


def test_unique_values_show_every_column_of_a_group_side_by_side():
    table = flat.flatten(
        _screen(
            ("[K+].[K+].[O-]C([O-])=O", "K2CO3"),
            ("CCN(CC)CC", "Et3N"),
            ("CCN(CC)CC", "Et3N"),
        )
    )
    values = table.unique_values(_BASE_NAME)
    assert list(values.columns) == [_BASE_SMILES, _BASE_NAME, "reactions"]
    assert values.to_dict("records") == [
        {_BASE_SMILES: "CCN(CC)CC", _BASE_NAME: "Et3N", "reactions": 2},
        {_BASE_SMILES: "[K+].[K+].[O-]C([O-])=O", _BASE_NAME: "K2CO3", "reactions": 1},
    ]


def test_a_name_used_for_two_structures_breaks_the_group():
    # Grouping must not hide a mismatch: "Et3N" written against two structures means
    # the name and the SMILES no longer pair one to one, so they stay apart.
    table = flat.flatten(
        _screen(
            ("[K+].[K+].[O-]C([O-])=O", "K2CO3"),
            ("CCN(CC)CC", "Et3N"),
            ("CCN(C(C)C)C(C)C", "Et3N"),
        )
    )
    assert [_BASE_SMILES] in table.column_groups
    assert [_BASE_NAME] in table.column_groups


def test_columns_whose_every_value_differs_are_not_grouped():
    # The yields and reaction IDs are all different, which pairs them one to one by
    # coincidence; that is not a reason to group them.
    table = flat.flatten(_screen(("CC", "a"), ("CCC", "b"), ("CCCC", "c")))
    assert ["reaction_id"] in table.column_groups
    assert [_YIELD] in table.column_groups


def test_the_fields_of_a_component_only_some_reactions_have_group_together():
    reactions = _screen(("CC", "a"), ("CC", "a"), ("CC", "a"))
    reactions[2].inputs["base"].components.append(
        message_helpers.build_compound(name="18-crown-6", amount="1 mg", role="reagent")
    )
    table = flat.flatten(reactions)
    extra = [c for c in table.varying.columns if ".components[1]." in c]
    assert len(extra) > 1
    assert extra in table.column_groups


def test_a_control_reaction_that_shifts_positions_shows_in_the_role_column():
    # As in the Perera screen: in the "no ligand" control the solvent takes the
    # ligand's place at components[0], so that position holds a ligand in some
    # reactions and a solvent in others. A reviewer should see the role column vary.
    reactions = []
    for i, ligand in enumerate(["PPh3", "XPhos", None]):
        reaction = reaction_pb2.Reaction(reaction_id=f"r{i}")
        components = reaction.inputs["ligand"].components
        if ligand is not None:
            components.append(
                message_helpers.build_compound(name=ligand, role="reagent")
            )
        components.append(
            message_helpers.build_compound(name="toluene", role="solvent")
        )
        reactions.append(reaction)
    table = flat.flatten(reactions)
    role = 'inputs["ligand"].components[0].reaction_role'
    assert table.varying[role].tolist() == ["REAGENT", "REAGENT", "SOLVENT"]


def _recorded(reaction_id: str, created: str) -> reaction_pb2.Reaction:
    reaction = _reaction(reaction_id)
    reaction.provenance.record_created.time.value = created
    reaction.provenance.record_created.person.name = "A. Chemist"
    return reaction


_CREATED = "provenance.record_created.time.value"


def test_record_times_minutes_apart_agree_at_the_default_resolution():
    table = flat.flatten(
        [
            _recorded("r1", "2026-01-05T10:00:00"),
            _recorded("r2", "2026-01-05T10:03:41"),
        ]
    )
    assert table.constants[_CREATED] == "2026-01-05"
    assert table.date_resolution == "day"


def test_record_times_on_different_days_vary():
    table = flat.flatten(
        [_recorded("r1", "2026-01-05T10:00:00"), _recorded("r2", "2026-01-06T09:00:00")]
    )
    assert table.varying[_CREATED].tolist() == ["2026-01-05", "2026-01-06"]


def test_dates_can_be_compared_at_another_resolution_or_as_written():
    reactions = [
        _recorded("r1", "2026-01-05T10:00:00"),
        _recorded("r2", "2026-01-05T10:03:41"),
    ]
    by_minute = flat.flatten(reactions, date_resolution="minute")
    assert by_minute.varying[_CREATED].tolist() == [
        "2026-01-05 10:00",
        "2026-01-05 10:03",
    ]
    by_month = flat.flatten(reactions, date_resolution="month")
    assert by_month.constants[_CREATED] == "2026-01"
    as_written = flat.flatten(reactions, date_resolution=None)
    assert as_written.varying[_CREATED].tolist() == [
        "2026-01-05T10:00:00",
        "2026-01-05T10:03:41",
    ]


def test_a_date_without_a_time_stays_a_date_at_finer_resolutions():
    table = flat.flatten(
        [_recorded("r1", "2026-01-05"), _recorded("r2", "2026-01-05T10:03:41")],
        date_resolution="hour",
    )
    assert table.varying[_CREATED].tolist() == ["2026-01-05", "2026-01-05 10:00"]


def test_a_value_that_is_not_a_date_is_left_as_written():
    table = flat.flatten(
        [_recorded("r1", "2026-01-05T10:00:00"), _recorded("r2", "last Tuesday-ish")]
    )
    assert table.varying[_CREATED].tolist() == ["2026-01-05", "last Tuesday-ish"]


@pytest.mark.parametrize(
    ("column", "expected"),
    [
        ("provenance.record_created.time.value", True),
        ("provenance.record_modified[2].time.value", True),
        ("provenance.experiment_start.value", True),
        ('outcomes[0].analyses["LCMS"].instrument_last_calibrated.value', True),
        ("provenance.record_created.person.name", False),
        ("conditions.temperature.setpoint.value", False),
        ("not.a.field", False),
    ],
)
def test_datetime_columns_are_found_from_the_schema(column, expected):
    assert flat.is_datetime(column) is expected


def test_provenance_is_kept_but_reported_in_its_own_section():
    reactions = [
        _recorded("r1", "2026-01-05T10:00:00"),
        _recorded("r2", "2026-01-06T10:00:00"),
    ]
    reactions[1].provenance.record_created.person.name = "B. Chemist"
    reactions[1].inputs["base"].components[0].identifiers[0].value = "Cs2CO3"
    table = flat.flatten(reactions)
    assert _CREATED in table.provenance_columns()
    assert "provenance.record_created.person.name" in table.provenance_columns()
    assert _BASE in table.experiment_columns()
    assert not set(table.provenance_columns()) & set(table.experiment_columns())
    constants = table.constants_frame().set_index("field")["section"]
    sections = constants.to_dict()
    assert sections['inputs["solvent"].components[0].reaction_role'] == "experiment"
