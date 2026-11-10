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

"""Tests for ord_schema.review.exemplars."""

import pathlib

from google.protobuf import text_format

from ord_schema import datasets, message_helpers, parquet
from ord_schema.analysis import columns, flat
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.review import exemplars, report, structure

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)
_BASES = ["CCN(CC)CC", "CN(C)C", "C1CCNCC1"]


def _reaction(
    reaction_id: str, *, base: str, yield_: float | None
) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.inputs["base"].components.append(
        message_helpers.build_compound(smiles=base, amount="1 mmol", role="reagent")
    )
    product = reaction.outcomes.add().products.add()
    product.identifiers.add(type="SMILES", value="CCO")
    if yield_ is not None:
        product.measurements.add(type="YIELD").percentage.value = yield_
    return reaction


def _screen() -> list[reaction_pb2.Reaction]:
    """Thirty reactions over three bases; one has no yield."""
    reactions = [
        _reaction(f"r{i}", base=_BASES[i % 3], yield_=float(i)) for i in range(30)
    ]
    reactions[7].outcomes[0].products[0].ClearField("measurements")
    return reactions


def _pick(
    reactions: list[reaction_pb2.Reaction],
    profile: structure.ValidationProfile = structure.ValidationProfile.UNIFORM,
    **kwargs,
) -> list[exemplars.Exemplar]:
    table = flat.flatten(reactions)
    return exemplars.pick(
        table, columns.column_statistics(table), reactions, profile, **kwargs
    )


def _by_label(picked: list[exemplars.Exemplar]) -> dict[str, list[str]]:
    return {exemplar.label: exemplar.reasons for exemplar in picked}


def test_one_reaction_of_each_shape_most_common_first():
    picked = _pick(_screen())
    assert picked[0].reaction_id == "r0"
    assert picked[0].reasons[0] == "the most common shape (29 reactions)"
    reasons = _by_label(picked)["r7"]
    assert reasons[0].startswith("a shape of its own: lacks ")
    assert "outcomes.products.measurements.percentage.value" in reasons[0]


def test_the_extremes_of_a_measured_value_are_picked():
    reasons = _by_label(_pick(_screen()))
    # r0 holds the lowest yield as well as the most common shape: listed once.
    assert reasons["r0"][1] == "lowest products[0].measurements[0].percentage (0)"
    assert reasons["r29"] == ["highest products[0].measurements[0].percentage (29)"]


def test_an_input_no_other_reaction_uses_is_picked_in_a_uniform_dataset():
    reactions = _screen()
    reactions[12].inputs["base"].components[0].identifiers[0].value = "CCCCN"
    reasons = _by_label(_pick(reactions))
    assert any("the only reaction with" in r and "CCCCN" in r for r in reasons["r12"])


def test_one_offs_are_not_picked_in_a_varied_dataset():
    reactions = _screen()
    reactions[12].inputs["base"].components[0].identifiers[0].value = "CCCCN"
    picked = _pick(reactions, structure.ValidationProfile.VARIED)
    assert "r12" not in _by_label(picked)


def test_no_more_than_the_cap_are_picked():
    assert len(_pick(_screen(), max_exemplars=2)) == 2


def test_reactions_that_all_differ_say_so():
    reactions = [_reaction("a", base="CC", yield_=1.0), reaction_pb2.Reaction()]
    reactions[1].reaction_id = "b"
    reactions[1].notes.procedure_details = "stirred"
    reasons = _by_label(_pick(reactions))
    assert reasons["a"][0] == "the first reaction; no two reactions share a shape"
    assert reasons["b"][0].startswith("a shape of its own")


def test_no_reactions_no_exemplars():
    assert _pick([]) == []


def test_reactions_without_ids_are_labelled_by_position():
    reactions = _screen()
    for reaction in reactions:
        reaction.ClearField("reaction_id")
    assert _pick(reactions)[0].label == "reaction #1"


def test_exemplars_are_fetched_from_a_parquet_dataset(tmp_path):
    path = tmp_path / "screen.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(name="screen", description="A screen", reactions=_screen()),
        path,
    )
    view = datasets.load_dataset(path)
    table = flat.flatten(view.reactions)
    picked = exemplars.pick(table, columns.column_statistics(table), view.reactions)
    assert picked[0].reaction == _screen()[0]


def test_the_pbtxt_is_a_dataset_of_the_reactions_as_they_were(tmp_path):
    reactions = _screen()
    picked = _pick(reactions)
    path = tmp_path / "exemplars.pbtxt"
    exemplars.write_pbtxt(picked, path, name="screen")
    dataset = text_format.Parse(path.read_text(), dataset_pb2.Dataset())
    assert dataset.name == "screen exemplars"
    assert list(dataset.reactions) == [reactions[e.position] for e in picked]
    assert "1. r0: the most common shape (29 reactions)" in dataset.description


def test_a_real_screen():
    reasons = [
        e.reasons for e in _pick(list(datasets.load_dataset(_NIELSEN).reactions))
    ]
    assert reasons == [
        ["the most common shape (80 reactions)"],
        ["lowest products[0].measurements[1].percentage (1)"],
        ["highest products[0].measurements[1].percentage (99)"],
    ]


def test_the_report_shows_each_exemplar_and_why(monkeypatch):
    monkeypatch.setattr(report, "MAX_EXEMPLAR_TEXT", 40)
    review = report.build_report(
        dataset_pb2.Dataset(name="screen", reactions=_screen())
    )
    page = report.render_html(review)
    assert "<h2>Exemplar reactions</h2>" in page
    assert "the most common shape (29 reactions)" in page
    assert "continued in tables/exemplars.pbtxt" in page


def test_a_control_lacking_one_inputs_amount_is_its_own_shape():
    # As in the Perera screen: other inputs still give amounts, so the field is
    # filled somewhere in every reaction; what differs is which input lacks it.
    reactions = []
    for i in range(6):
        reaction = _reaction(f"r{i}", base=_BASES[i % 3], yield_=float(i))
        ligand = message_helpers.build_compound(name="PPh3", role="reagent")
        if i != 4:
            ligand.amount.moles.value = 0.1
            ligand.amount.moles.units = reaction_pb2.Moles.MILLIMOLE
        reaction.inputs["ligand"].components.append(ligand)
        reactions.append(reaction)
    reasons = _by_label(_pick(reactions))
    assert reasons["r4"][0].startswith("a shape of its own: lacks ")
    assert 'inputs["ligand"].components.amount.moles.value' in reasons["r4"][0]
