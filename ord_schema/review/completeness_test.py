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

"""Tests for ord_schema.review.completeness."""

import csv
import dataclasses

import pytest

from ord_schema import completeness as core
from ord_schema import completeness_test
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.review import completeness, report


def _reactions() -> list[reaction_pb2.Reaction]:
    """Five complete reactions, three of which lack the atmosphere."""
    reactions = []
    for i in range(5):
        reaction = completeness_test._complete_reaction()
        reaction.reaction_id = f"r{i}"
        if i < 3:
            reaction.conditions.pressure.ClearField("atmosphere")
        reactions.append(reaction)
    return reactions


def _rollup() -> completeness.CompletenessRollup:
    rollup = completeness.CompletenessRollup("uniform")
    for _ in rollup.watch(_reactions()):
        pass
    return rollup


def test_reactions_lacking_an_item_are_counted_once():
    rollup = _rollup()
    (lacking,) = rollup.lacking()
    assert lacking.item.id == "atmosphere"
    assert (lacking.lacking, lacking.applicable) == (3, 5)
    assert lacking.examples == ["r0", "r1", "r2"]
    assert lacking.suggestion() == (
        "Suggestion: record conditions.pressure.atmosphere.type "
        "(3 of 5 reactions lack it)"
    )


def test_the_means_are_over_reactions():
    rollup = _rollup()
    assert rollup.reaction_count == 5
    assert rollup.mean_levels()["inputs"] == 1
    assert rollup.mean_levels()["conditions"] == pytest.approx((3 * 2 / 3 + 2) / 5)
    overall = rollup.mean_overall()
    assert overall is not None
    assert overall < 1


def test_the_csv_has_each_reactions_scores_and_what_it_lacks(tmp_path):
    path = tmp_path / "completeness.csv"
    _rollup().write_csv(path)
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert [row["reaction"] for row in rows] == ["r0", "r1", "r2", "r3", "r4"]
    assert rows[0]["lacks"] == "atmosphere"
    assert rows[4]["lacks"] == ""
    assert rows[4]["overall"] == "1.000"


def test_the_report_rolls_completeness_up(tmp_path):
    review = report.build_report(dataset_pb2.Dataset(name="x", reactions=_reactions()))
    page = report.write_bundle(review, tmp_path).read_text()
    assert "<h2>Completeness</h2>" in page
    assert "3 of 5" in page
    with (tmp_path / "tables" / "findings.csv").open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["source"] == "completeness"]
    assert [(r["check"], r["severity"], r["reactions"]) for r in rows] == [
        ("atmosphere", "INFO", "3")
    ]


def test_a_rollup_against_a_custom_rubric_reports_its_levels():
    shipped = core.load_rubric()
    rubric = dataclasses.replace(
        shipped,
        levels=("everything",),
        items=tuple(
            dataclasses.replace(item, level="everything") for item in shipped.items
        ),
    )
    rollup = completeness.CompletenessRollup(rubric=rubric)
    rollup.observe(reaction_pb2.Reaction(reaction_id="r1"), "r1")
    assert list(rollup.mean_levels()) == ["everything"]
