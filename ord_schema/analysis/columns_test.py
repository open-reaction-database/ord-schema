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

"""Tests for ord_schema.analysis.columns."""

import csv
import pathlib
import xml.etree.ElementTree as ET

import pytest

from ord_schema import datasets, message_helpers
from ord_schema.analysis import columns, flat
from ord_schema.proto import reaction_pb2

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)
_SVG = "{http://www.w3.org/2000/svg}"
_YIELD = "outcomes[0].products[0].measurements[0].percentage.value"
_BASE_SMILES = 'inputs["base"].components[0].identifiers[0].value'
_BASE_NAME = 'inputs["base"].components[0].identifiers[1].value'


def _reaction(
    reaction_id: str, *, base: tuple[str, str], yield_: float | None, temperature: float
) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    smiles, name = base
    reaction.inputs["base"].components.append(
        message_helpers.build_compound(smiles=smiles, name=name, role="reagent")
    )
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = temperature
    setpoint.units = reaction_pb2.Temperature.CELSIUS
    product = reaction.outcomes.add().products.add()
    if yield_ is not None:
        product.measurements.add(type="YIELD").percentage.value = yield_
    return reaction


_BASES = [("CCN(CC)CC", "Et3N"), ("CN(C)C", "Me3N"), ("C1CCNCC1", "piperidine")]


def _stats() -> dict[str, columns.ColumnStats]:
    reactions = [
        _reaction(
            f"r{i}",
            base=_BASES[i % 3],
            yield_=None if i == 0 else float(i),
            temperature=25 if i % 2 else 60,
        )
        for i in range(30)
    ]
    table = flat.flatten(reactions)
    return {stat.column: stat for stat in columns.column_statistics(table)}


def test_many_distinct_numbers_are_a_distribution():
    stat = _stats()[_YIELD]
    assert stat.kind == "numeric"
    assert (stat.filled, stat.empty, stat.distinct) == (29, 1, 29)
    assert stat.summary["min"] == 1
    assert stat.summary["max"] == 29
    assert stat.summary["median"] == 15
    assert stat.histogram is not None
    edges, counts = stat.histogram
    assert sum(counts) == stat.filled
    assert len(edges) == len(counts) + 1 == columns.HISTOGRAM_BINS + 1


def test_a_few_numeric_levels_are_categories():
    # A temperature set at 25 and 60 is two levels, not a distribution.
    stat = _stats()["conditions.temperature.setpoint.value"]
    assert stat.kind == "categorical"
    assert dict(stat.categories) == {(25.0,): 15, (60.0,): 15}


def test_a_groups_values_are_listed_side_by_side():
    stat = _stats()[_BASE_SMILES]
    assert stat.columns == [_BASE_SMILES, _BASE_NAME]
    assert stat.distinct == 3
    assert ("CCN(CC)CC", "Et3N") in dict(stat.categories)
    assert all(count == 10 for _, count in stat.categories)


def test_reaction_ids_are_not_described():
    assert "reaction_id" not in _stats()


def test_categories_past_the_listed_limit_are_counted_as_other(monkeypatch):
    monkeypatch.setattr(columns, "MAX_LISTED_CATEGORIES", 2)
    stat = _stats()[_BASE_SMILES]
    assert len(stat.categories) == 2
    assert stat.other == 10
    assert stat.distinct == 3


def test_provenance_is_its_own_section():
    reactions = [
        _reaction("r1", base=_BASES[0], yield_=1, temperature=25),
        _reaction("r2", base=_BASES[1], yield_=2, temperature=25),
    ]
    reactions[0].provenance.record_created.person.name = "A. Chemist"
    reactions[1].provenance.record_created.person.name = "B. Chemist"
    stats = columns.column_statistics(flat.flatten(reactions))
    sections = {stat.column: stat.section for stat in stats}
    assert sections["provenance.record_created.person.name"] == "provenance"
    assert sections[_BASE_SMILES] == "experiment"


def test_the_csv_has_one_row_per_group(tmp_path):
    stats = _stats()
    path = tmp_path / "column_stats.csv"
    columns.write_csv(list(stats.values()), path)
    with path.open(newline="") as f:
        rows = {row["column"]: row for row in csv.DictReader(f)}
    assert set(rows) == set(stats)
    assert rows[_YIELD]["kind"] == "numeric"
    assert float(rows[_YIELD]["median"]) == 15
    assert "CCN(CC)CC | Et3N (10)" in rows[_BASE_SMILES]["values"]
    assert rows[_BASE_SMILES]["columns"] == f"{_BASE_SMILES} | {_BASE_NAME}"


def _parse(svg: str) -> ET.Element:
    root = ET.fromstring(svg)  # noqa: S314 -- SVG this module just wrote
    assert root.tag == f"{_SVG}svg"
    return root


def test_every_chart_is_well_formed_svg_with_a_tooltip_on_each_mark():
    for stat in _stats().values():
        root = _parse(stat.to_svg())
        assert root.find(f"{_SVG}title") is not None
        marks = root.findall(f".//{_SVG}path")
        assert marks
        for group in root.findall(f"{_SVG}g"):
            assert group.find(f"{_SVG}title") is not None


def test_charts_take_their_colors_from_css_variables():
    svg = _stats()[_YIELD].to_svg()
    assert "var(--series-1, #2a78d6)" in svg
    assert "var(--grid, " in svg


def test_chart_text_is_escaped():
    reactions = [
        _reaction(
            f"r{i}", base=("CC", f"<b>base {i % 2}</b>"), yield_=1, temperature=25
        )
        for i in range(4)
    ]
    stats = columns.column_statistics(flat.flatten(reactions))
    (stat,) = [s for s in stats if s.column == _BASE_NAME]
    svg = stat.to_svg()
    _parse(svg)
    assert "<b>" not in svg
    assert "&lt;b&gt;" in svg


def test_categories_beyond_the_chart_fold_into_a_line_of_text(monkeypatch):
    monkeypatch.setattr(columns, "MAX_CHART_CATEGORIES", 2)
    svg = _stats()[_BASE_SMILES].to_svg()
    root = _parse(svg)
    assert len(root.findall(f".//{_SVG}path")) == 2
    assert "+ 1 more values, 10 reactions" in svg


@pytest.mark.parametrize(
    "label",
    ["short", "O[C@@H]1C[C@H](OCC2=CC=CC=C2)C1 | c1ccccc1", "W" * 80, ""],
)
def test_long_labels_are_cut_to_fit(label):
    fitted = columns._fit(label, columns._LABEL_ROOM)
    assert columns._text_width(fitted) <= columns._LABEL_ROOM
    assert fitted == label or fitted.endswith("…")


def test_a_real_screen():
    table = flat.flatten(datasets.load_dataset(_NIELSEN).reactions)
    stats = columns.column_statistics(table)
    kinds = [stat.kind for stat in stats]
    assert kinds == ["categorical", "categorical", "categorical", "numeric"]
    assert [stat.distinct for stat in stats] == [4, 4, 5, 51]
