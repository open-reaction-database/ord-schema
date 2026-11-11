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

"""Completeness across a dataset: the rubric scores of its reactions, rolled up.

:mod:`ord_schema.completeness` scores one reaction. This scores every reaction of a
dataset under its profile and says, for each rubric item, how many reactions lack
it: "412 of 500 reactions lack the atmosphere". That turns one suggestion repeated
500 times into one line a submitter can act on. Advisory, like every score: nothing
here blocks a submission.
"""

import csv
import dataclasses
import os
import pathlib
from collections.abc import Iterable, Iterator

from ord_schema import completeness
from ord_schema.proto import reaction_pb2

# Reactions kept per item as examples of what lacks it.
MAX_EXAMPLES = 5


@dataclasses.dataclass
class ItemRollup:
    """One rubric item across the dataset.

    Attributes:
        item: The rubric item.
        weight: Its weight under the dataset's profile.
        applicable: How many reactions it applies to.
        lacking: How many of those are not complete on it.
        examples: A few of the reactions that lack it.
    """

    item: completeness.RubricItem
    weight: float
    applicable: int = 0
    lacking: int = 0
    examples: list[str] = dataclasses.field(default_factory=list)

    def suggestion(self) -> str:
        """The item's suggestion, with how many reactions it concerns."""
        return (
            f"Suggestion: {self.item.suggestion} "
            f"({self.lacking:,} of {self.applicable:,} reactions lack it)"
        )


@dataclasses.dataclass
class _Row:
    label: str
    score: completeness.CompletenessScore


class CompletenessRollup:
    """Scores a dataset's reactions as they stream past, and rolls the scores up."""

    def __init__(
        self,
        profile: str = completeness.DEFAULT_PROFILE,
        rubric: completeness.Rubric | None = None,
    ) -> None:
        """Starts an empty rollup.

        Args:
            profile: The dataset's profile, whose weights the scores use.
            rubric: The rubric; the shipped one by default.
        """
        self.profile = profile
        self.rubric = rubric or completeness.load_rubric()
        self.items = {
            item.id: ItemRollup(item=item, weight=item.weights[profile])
            for item in self.rubric.items
        }
        self._rows: list[_Row] = []

    def observe(self, reaction: reaction_pb2.Reaction, label: str) -> None:
        """Scores one reaction and adds it to the rollup."""
        score = completeness.score_reaction(reaction, self.profile, self.rubric)
        self._rows.append(_Row(label, score))
        for scored in score.items:
            if scored.score is None:
                continue
            rollup = self.items[scored.item.id]
            rollup.applicable += 1
            if scored.score < 1:
                rollup.lacking += 1
                if len(rollup.examples) < MAX_EXAMPLES:
                    rollup.examples.append(label)

    def watch(
        self, reactions: Iterable[reaction_pb2.Reaction]
    ) -> Iterator[reaction_pb2.Reaction]:
        """Yields ``reactions`` unchanged, scoring each on the way past."""
        for position, reaction in enumerate(reactions):
            self.observe(reaction, reaction.reaction_id or f"reaction #{position + 1}")
            yield reaction

    @property
    def reaction_count(self) -> int:
        """How many reactions have been scored."""
        return len(self._rows)

    @staticmethod
    def _mean(values: list[float | None]) -> float | None:
        present = [v for v in values if v is not None]
        return sum(present) / len(present) if present else None

    def mean_overall(self) -> float | None:
        """The mean of the reactions' overall scores."""
        return self._mean([row.score.overall for row in self._rows])

    def mean_levels(self) -> dict[str, float | None]:
        """The mean score of each message level, in rubric order."""
        levels = [row.score.levels for row in self._rows]
        return {
            level: self._mean([scores[level] for scores in levels])
            for level in self.rubric.levels
        }

    def lacking(self) -> list[ItemRollup]:
        """The items some reactions lack, the most widely lacking first."""
        rows = [r for r in self.items.values() if r.lacking and r.weight > 0]
        return sorted(rows, key=lambda r: (-r.lacking / r.applicable, -r.weight))

    def write_csv(self, path: str | os.PathLike[str]) -> None:
        """Writes one row per reaction: its overall and level scores, and what it lacks.

        Scores are from 0 to 1, empty where nothing applied.
        """
        fields = ["reaction", "overall", *self.rubric.levels, "lacks"]
        with pathlib.Path(path).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for row in self._rows:
                record: dict[str, object] = {
                    "reaction": row.label,
                    "overall": _round(row.score.overall),
                    "lacks": " | ".join(s.item.id for s in row.score.missing()),
                }
                record.update(
                    {level: _round(value) for level, value in row.score.levels.items()}
                )
                writer.writerow(record)


def _round(value: float | None) -> str:
    return "" if value is None else f"{value:.3f}"
