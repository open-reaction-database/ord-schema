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

"""How completely a reaction records the details that can almost always be recorded.

A reaction can be valid and still leave out what a reader would expect: the roles of
its inputs, the atmosphere, which analysis a yield came from. This scores a reaction
against a rubric of such details and says which are missing, as suggestions. The
score is advisory: it never makes a reaction invalid, and it is computed only when
asked for, not as part of validation.

The rubric is data, not code: ``completeness_rubric.json`` beside this module lists
each item, the message level it belongs to (inputs, conditions, workups, outcomes,
provenance), and its weight under each dataset profile. Revising the rubric means
editing that file; this module holds only how each item is checked::

    score = score_reaction(reaction)
    score.overall  # 0.0 to 1.0, or None if no item applies
    score.levels["conditions"]
    score.suggestions()  # ["Suggestion: record conditions.pressure.atmosphere.type"]

It lives in the core package, with no extra dependencies, so the reaction editor can
score the reaction being written. :mod:`ord_schema.review` rolls the scores up across
a dataset.
"""

import dataclasses
import functools
import json
from collections.abc import Callable, Iterator
from importlib import resources

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

# The dataset profile a reaction is scored under when none is given; the rubric has
# a weight for each profile.
DEFAULT_PROFILE = "uniform"

_Compound = reaction_pb2.Compound | reaction_pb2.ProductCompound
# How much of an item a reaction has, and how much it could have: (have, of). None
# when the item does not apply, such as product measurements with no products.
Tally = tuple[int, int] | None


@dataclasses.dataclass(frozen=True)
class RubricItem:
    """One detail the rubric scores.

    Attributes:
        id: Stable ID, also the key of the check that scores it.
        level: The message level it belongs to.
        fields: The field paths it concerns.
        description: What a complete reaction has.
        suggestion: What to do when it is missing, as an imperative.
        weights: Its weight under each profile; 0 leaves it out.
    """

    id: str
    level: str
    fields: tuple[str, ...]
    description: str
    suggestion: str
    weights: dict[str, float]


@dataclasses.dataclass(frozen=True)
class Rubric:
    """The versioned rubric: its items, and the levels they are grouped in."""

    version: str
    levels: tuple[str, ...]
    items: tuple[RubricItem, ...]

    @property
    def profiles(self) -> set[str]:
        """The profiles every item has a weight for."""
        return set.intersection(*(set(item.weights) for item in self.items))


@functools.cache
def load_rubric() -> Rubric:
    """Returns the rubric shipped with the package."""
    text = (resources.files("ord_schema") / "completeness_rubric.json").read_text(
        encoding="utf-8"
    )
    data = json.loads(text)
    return Rubric(
        version=data["version"],
        levels=tuple(data["levels"]),
        items=tuple(
            RubricItem(
                id=item["id"],
                level=item["level"],
                fields=tuple(item["fields"]),
                description=item["description"],
                suggestion=item["suggestion"],
                weights=dict(item["weights"]),
            )
            for item in data["items"]
        ),
    )


def _components(reaction: reaction_pb2.Reaction) -> list[reaction_pb2.Compound]:
    return [c for value in reaction.inputs.values() for c in value.components]


def _products(reaction: reaction_pb2.Reaction) -> list[reaction_pb2.ProductCompound]:
    return [p for outcome in reaction.outcomes for p in outcome.products]


def _measurements(
    reaction: reaction_pb2.Reaction,
) -> Iterator[tuple[reaction_pb2.ReactionOutcome, reaction_pb2.ProductMeasurement]]:
    for outcome in reaction.outcomes:
        for product in outcome.products:
            for measurement in product.measurements:
                yield outcome, measurement


def _count(items: list, test: Callable[..., bool]) -> Tally:
    """How many of ``items`` pass ``test``, or None if there are none."""
    if not items:
        return None
    return sum(1 for item in items if test(item)), len(items)


def _one(*, applies: bool = True, present: bool) -> Tally:
    """A yes-or-no item, which counts as one."""
    if not applies:
        return None
    return (1 if present else 0), 1


def _has_name(compound: _Compound) -> bool:
    return any(
        i.type == reaction_pb2.CompoundIdentifier.NAME and i.value
        for i in compound.identifiers
    )


def _has_structure(compound: _Compound) -> bool:
    return any(True for _ in message_helpers.structural_identifiers(compound))


def _check_compound_names(reaction: reaction_pb2.Reaction) -> Tally:
    with_structure = [
        c for c in [*_components(reaction), *_products(reaction)] if _has_structure(c)
    ]
    return _count(with_structure, _has_name)


def _check_measurement_analysis(reaction: reaction_pb2.Reaction) -> Tally:
    pairs = list(_measurements(reaction))
    return _count(
        pairs,
        lambda pair: (
            bool(pair[1].analysis_key) and pair[1].analysis_key in pair[0].analyses
        ),
    )


_UNSPECIFIED = 0

_CHECKS: dict[str, Callable[[reaction_pb2.Reaction], Tally]] = {
    "input-roles": lambda r: _count(
        _components(r), lambda c: c.reaction_role != _UNSPECIFIED
    ),
    "input-amounts": lambda r: _count(_components(r), lambda c: c.HasField("amount")),
    "limiting-reactant": lambda r: _one(
        applies=bool(_components(r)),
        present=any(c.is_limiting for c in _components(r)),
    ),
    "addition-order": lambda r: (
        _count(list(r.inputs.values()), lambda i: i.addition_order > 0)
        if len(r.inputs) > 1
        else None
    ),
    "compound-names": _check_compound_names,
    "temperature-control": lambda r: _one(
        present=r.conditions.temperature.control.type != _UNSPECIFIED
    ),
    "atmosphere": lambda r: _one(
        present=r.conditions.pressure.atmosphere.type != _UNSPECIFIED
    ),
    "stirring": lambda r: _one(present=r.conditions.stirring.type != _UNSPECIFIED),
    "workups": lambda r: _one(present=len(r.workups) > 0),
    "reaction-time": lambda r: _count(
        list(r.outcomes), lambda o: o.HasField("reaction_time")
    ),
    "measurement-type": lambda r: _count(
        [m for _, m in _measurements(r)], lambda m: m.type != _UNSPECIFIED
    ),
    "measurement-analysis": _check_measurement_analysis,
    "desired-product": lambda r: _one(
        applies=bool(_products(r)),
        present=any(p.HasField("is_desired_product") for p in _products(r)),
    ),
    "source": lambda r: _one(
        present=bool(
            r.provenance.doi or r.provenance.patent or r.provenance.publication_url
        ),
    ),
    "record-created": lambda r: _one(
        present=bool(r.provenance.record_created.time.value)
    ),
}


@dataclasses.dataclass(frozen=True)
class ItemScore:
    """How a reaction does on one rubric item.

    Attributes:
        item: The rubric item.
        weight: Its weight under the profile scored.
        have: How many of the things it counts the reaction has.
        of: How many it could have; ``have`` and ``of`` are None if it does not
            apply.
    """

    item: RubricItem
    weight: float
    have: int | None
    of: int | None

    @property
    def applies(self) -> bool:
        """Whether the item applies to the reaction and counts under the profile."""
        return self.of is not None and self.weight > 0

    @property
    def score(self) -> float | None:
        """The share the reaction has, from 0 to 1, or None if it does not apply."""
        if not self.applies or not self.of:
            return None
        return (self.have or 0) / self.of


def _weighted(scores: list[ItemScore]) -> float | None:
    """The weighted mean score of the items that apply, or None if none does."""
    applying = [s for s in scores if s.score is not None]
    total = sum(s.weight for s in applying)
    if not total:
        return None
    return sum(s.weight * (s.score or 0.0) for s in applying) / total


@dataclasses.dataclass(frozen=True)
class CompletenessScore:
    """A reaction's completeness against the rubric.

    Attributes:
        profile: The profile the weights came from.
        rubric_version: The rubric's version.
        items: Each rubric item's score, in rubric order.
        rubric_levels: The rubric's message levels, in its order.
    """

    profile: str
    rubric_version: str
    items: tuple[ItemScore, ...]
    rubric_levels: tuple[str, ...]

    @property
    def overall(self) -> float | None:
        """The weighted mean over every item that applies, or None if none does."""
        return _weighted(list(self.items))

    @property
    def levels(self) -> dict[str, float | None]:
        """The weighted mean for each of the rubric's message levels, in its order."""
        levels: dict[str, float | None] = dict.fromkeys(self.rubric_levels)
        for level in levels:
            levels[level] = _weighted([s for s in self.items if s.item.level == level])
        return levels

    def missing(self) -> list[ItemScore]:
        """The items that apply and are not complete."""
        return [s for s in self.items if s.score is not None and s.score < 1]

    def suggestions(self) -> list[str]:
        """What to add, one line per missing item, in the form INFO findings take."""
        lines = []
        for scored in self.missing():
            line = f"Suggestion: {scored.item.suggestion}"
            if scored.of and scored.of > 1:
                line += f" ({scored.have} of {scored.of} have it)"
            lines.append(line)
        return lines


def score_reaction(
    reaction: reaction_pb2.Reaction,
    profile: str = DEFAULT_PROFILE,
    rubric: Rubric | None = None,
) -> CompletenessScore:
    """Scores a reaction against the completeness rubric.

    Args:
        reaction: The reaction to score. It is only read.
        profile: The dataset profile whose weights to use: ``uniform`` or
            ``varied``.
        rubric: The rubric to score against; the shipped one by default.

    Returns:
        The reaction's score on each item, and the weighted means.

    Raises:
        ValueError: If the rubric has no weights for ``profile``, or an item no check.
    """
    rubric = rubric or load_rubric()
    if profile not in rubric.profiles:
        raise ValueError(
            f"the rubric has no weights for profile {profile!r}; "
            f"it has {sorted(rubric.profiles)}"
        )
    scores = []
    for item in rubric.items:
        check = _CHECKS.get(item.id)
        if check is None:
            raise ValueError(f"no check for rubric item {item.id!r}")
        tally = check(reaction)
        have, of = tally if tally is not None else (None, None)
        scores.append(
            ItemScore(item=item, weight=item.weights[profile], have=have, of=of)
        )
    return CompletenessScore(
        profile=profile,
        rubric_version=rubric.version,
        items=tuple(scores),
        rubric_levels=rubric.levels,
    )


def checked_items() -> set[str]:
    """The rubric item IDs this module knows how to check."""
    return set(_CHECKS)
