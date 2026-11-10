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

"""A summary of what a dataset's reactions record, built one reaction at a time.

For each field path -- ``inputs.components.amount.moles.value``, say -- the summary
counts how many reactions populate it and how often each value occurs there. It also
records each reaction's *shape*, the set of field paths it populates, so reactions
that were recorded the same way can be told from those that were not, and it carries
the dataset's :class:`~ord_schema.analysis.compounds.CompoundTable`.

Most of the analysis and review tooling reads from this rather than from the
reactions: the flattened table, column statistics, exemplar reactions and the
dataset-level checks all start here.

Paths name fields, not positions: a repeated field or a map contributes its field
name once, whatever the index or key, so ``inputs["amine"].components[0]`` and
``inputs["base"].components[2]`` share the path ``inputs.components``. Map keys are
counted separately, under the map's path.

Below a map, each field is also recorded under a *keyed* path that names the entry,
``inputs["base"].components.identifiers.value``. The merged path says what the
dataset records; the keyed path says which input it records it in, which is what
separates a screen's variables: in a screen varying the alcohol, the base and the
sulfonyl fluoride, all three vary under ``inputs.components.identifiers.value``, but
each under its own keyed path. The field queries take ``by_key`` to choose. Shapes
use merged paths only, so a dataset whose input names differ from reaction to
reaction is not split into one shape per reaction.

The summary follows the same pattern as
:class:`~ord_schema.validations.DatasetCrossRefState`: build one per slice of a
dataset with :meth:`DatasetSummary.observe`, then combine them with
:meth:`DatasetSummary.merge`. Merging gives the same result as observing every
reaction in one summary, so a large Parquet dataset can be summarized by row group
in parallel::

    summary = DatasetSummary()
    for reaction in datasets.load_dataset("dataset.parquet").reactions:
        summary.observe(reaction)
    summary.varying_fields()
"""

import dataclasses
import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator

from google.protobuf import descriptor, message

from ord_schema.analysis.compounds import CompoundTable
from ord_schema.proto import reaction_pb2

# Distinct values kept per field path. A field that varies continuously -- an amount, a
# yield -- would otherwise grow with the dataset; past this, further new values are
# only counted, in ``DatasetSummary.untracked``.
MAX_DISTINCT_VALUES = 1000

# Reaction IDs kept for each shape, so a reader can find an example of each.
MAX_SHAPE_EXAMPLES = 3

# Long text is cut to this many characters before it is counted, so procedure details
# and MOL blocks do not dominate the memory a summary holds.
MAX_VALUE_LENGTH = 200

Value = str | int | float | bool
Shape = frozenset[str]
# The values one reaction holds at one path, sorted: a reaction with three components
# holds three roles at ``inputs.components.reaction_role``.
ReactionValues = tuple[Value, ...]


def _value(field: descriptor.FieldDescriptor, raw: object) -> Value:
    """Returns ``raw`` in the form it is counted in.

    Enums are counted by name, from the proto descriptor; bytes by their length;
    long strings cut to ``MAX_VALUE_LENGTH``.
    """
    if field.type == descriptor.FieldDescriptor.TYPE_ENUM:
        enum_value = field.enum_type.values_by_number.get(raw)
        return enum_value.name if enum_value is not None else str(raw)
    if isinstance(raw, bytes):
        return f"({len(raw)} bytes)"
    if isinstance(raw, str) and len(raw) > MAX_VALUE_LENGTH:
        return raw[:MAX_VALUE_LENGTH] + "…"
    assert isinstance(raw, str | int | float | bool)
    return raw


def _count_capped(
    counters: dict[str, Counter],
    untracked: Counter[str],
    path: str,
    item: object,
    count: int,
) -> None:
    """Counts ``item`` under ``path``, or in ``untracked`` once ``path`` is full."""
    counter = counters[path]
    if item in counter or len(counter) < MAX_DISTINCT_VALUES:
        counter[item] += count
    else:
        untracked[path] += count


def _sort_key(value: Value) -> tuple[str, str]:
    """Orders values of mixed types, so a reaction's values sort the same every time."""
    return type(value).__name__, str(value)


def _key_segment(key: object) -> str:
    """Returns a map key as it appears in a keyed path: ``["base"]``."""
    escaped = str(key).replace("\\", "\\\\").replace('"', '\\"')
    return f'["{escaped}"]'


# A key segment as _key_segment writes it, for recovering a keyed path's merged path.
_KEY_SEGMENT = re.compile(r'\["(?:[^"\\]|\\.)*"\]')


def merged_path(path: str) -> str:
    """Returns the merged path for a keyed one, or ``path`` itself if it has no keys.

    ``inputs["base"].components.identifiers.value`` becomes
    ``inputs.components.identifiers.value``.
    """
    return _KEY_SEGMENT.sub("", path)


def is_keyed(path: str) -> bool:
    """Whether ``path`` names a map entry, as keyed paths do."""
    return _KEY_SEGMENT.search(path) is not None


def _walk(
    msg: message.Message, prefix: str = "", keyed_prefix: str = ""
) -> Iterator[tuple[str, str, Value | None, str | None]]:
    """Yields every populated field below ``msg``.

    Yields:
        ``(path, keyed_path, value, None)`` for a scalar or enum,
        ``(path, keyed_path, None, None)`` for a populated submessage, and
        ``(path, keyed_path, None, key)`` for each key of a map. ``keyed_path`` equals
        ``path`` outside any map.
    """
    for field, content in msg.ListFields():
        path = f"{prefix}{field.name}"
        keyed = f"{keyed_prefix}{field.name}"
        if field.message_type is not None and field.message_type.GetOptions().map_entry:
            value_field = field.message_type.fields_by_name["value"]
            # In key order: protobuf does not define a map's order, so reading it as
            # found would make two reactions recorded identically look reordered.
            for key, entry in sorted(content.items(), key=lambda item: str(item[0])):
                yield path, keyed, None, str(key)
                entry_path = f"{keyed}{_key_segment(key)}"
                if value_field.message_type is None:
                    yield f"{path}.value", entry_path, _value(value_field, entry), None
                else:
                    yield path, entry_path, None, None
                    yield from _walk(entry, f"{path}.", f"{entry_path}.")
        elif field.type == descriptor.FieldDescriptor.TYPE_MESSAGE:
            items = content if field.label == field.LABEL_REPEATED else [content]
            for item in items:
                yield path, keyed, None, None
                yield from _walk(item, f"{path}.", f"{keyed}.")
        else:
            items = content if field.label == field.LABEL_REPEATED else [content]
            for item in items:
                yield path, keyed, _value(field, item), None


@dataclasses.dataclass
class DatasetSummary:
    """What a dataset's reactions record, field by field.

    Every per-path attribute holds merged paths and, below a map, keyed paths too;
    :func:`is_keyed` tells them apart.

    Attributes:
        reaction_count: How many reactions have been observed.
        field_counts: For each field path, how many reactions populate it.
        values: For each scalar or enum field path, how often each value occurs,
            counting every occurrence (three components with the same role count
            three times).
        untracked: For each path that reached ``MAX_DISTINCT_VALUES``, how many
            occurrences carried a value not among those kept.
        reaction_values: For each scalar or enum field path, how many reactions hold
            each combination of values there. This, not ``values``, says whether a
            field differs between reactions: components with different roles give
            ``values`` several entries even when every reaction has the same roles.
        reaction_values_untracked: As ``untracked``, for ``reaction_values``,
            counting reactions.
        ordered_values: As ``reaction_values``, but in the order each reaction holds
            them, map entries taken in key order. Two reactions with the same
            ``reaction_values`` and different ``ordered_values`` list the same things
            in a different order, e.g. the base before the solvent in one and after
            it in the other.
        ordered_values_untracked: As ``untracked``, for ``ordered_values``.
        map_keys: For each map field path, how many reactions use each key.
        shapes: How many reactions share each shape, a shape being the set of merged
            field paths a reaction populates.
        shape_examples: A few reaction IDs for each shape.
        compounds: The dataset's compounds, roles and identifiers.
    """

    reaction_count: int = 0
    field_counts: Counter[str] = dataclasses.field(default_factory=Counter)
    values: dict[str, Counter[Value]] = dataclasses.field(
        default_factory=lambda: defaultdict(Counter)
    )
    untracked: Counter[str] = dataclasses.field(default_factory=Counter)
    reaction_values: dict[str, Counter[ReactionValues]] = dataclasses.field(
        default_factory=lambda: defaultdict(Counter)
    )
    reaction_values_untracked: Counter[str] = dataclasses.field(default_factory=Counter)
    ordered_values: dict[str, Counter[ReactionValues]] = dataclasses.field(
        default_factory=lambda: defaultdict(Counter)
    )
    ordered_values_untracked: Counter[str] = dataclasses.field(default_factory=Counter)
    map_keys: dict[str, Counter[str]] = dataclasses.field(
        default_factory=lambda: defaultdict(Counter)
    )
    shapes: Counter[Shape] = dataclasses.field(default_factory=Counter)
    shape_examples: dict[Shape, list[str]] = dataclasses.field(default_factory=dict)
    compounds: CompoundTable = dataclasses.field(default_factory=CompoundTable)

    def observe(self, reaction: reaction_pb2.Reaction) -> None:
        """Adds one reaction to the summary."""
        self.reaction_count += 1
        paths: set[str] = set()
        keyed_paths: set[str] = set()
        keys: dict[str, set[str]] = defaultdict(set)
        held: dict[str, list[Value]] = defaultdict(list)
        for path, keyed, value, key in _walk(reaction):
            paths.add(path)
            if key is not None:
                keys[path].add(key)
                continue
            recorded = (path,) if keyed == path else (path, keyed)
            keyed_paths.update(recorded[1:])
            if value is not None:
                for each in recorded:
                    self._count_value(each, value, 1)
                    held[each].append(value)
        for path, values in held.items():
            self._count_reaction_values(path, tuple(sorted(values, key=_sort_key)), 1)
            _count_capped(
                self.ordered_values,
                self.ordered_values_untracked,
                path,
                tuple(values),
                1,
            )
        self.field_counts.update(paths | keyed_paths)
        for path, used in keys.items():
            self.map_keys[path].update(used)
        shape = frozenset(paths)
        self.shapes[shape] += 1
        examples = self.shape_examples.setdefault(shape, [])
        if reaction.reaction_id and len(examples) < MAX_SHAPE_EXAMPLES:
            examples.append(reaction.reaction_id)
        self.compounds.add_reaction(reaction)

    def observe_all(self, reactions: Iterable[reaction_pb2.Reaction]) -> None:
        """Adds each of ``reactions`` to the summary."""
        for reaction in reactions:
            self.observe(reaction)

    def _count_value(self, path: str, value: Value, count: int) -> None:
        _count_capped(self.values, self.untracked, path, value, count)

    def _count_reaction_values(
        self, path: str, values: ReactionValues, count: int
    ) -> None:
        _count_capped(
            self.reaction_values, self.reaction_values_untracked, path, values, count
        )

    def _value_paths(self, *, by_key: bool) -> list[str]:
        """The scalar paths a field query reads, sorted.

        Merged paths only, or with ``by_key`` keyed paths in place of the merged paths
        they refine; paths outside any map are the same either way.
        """
        if not by_key:
            return sorted(p for p in self.reaction_values if not is_keyed(p))
        refined = {merged_path(p) for p in self.reaction_values if is_keyed(p)}
        return sorted(p for p in self.reaction_values if p not in refined)

    def merge(self, other: "DatasetSummary") -> None:
        """Adds a summary of another part of the dataset to this one.

        The result is the same as if this summary had observed the other's reactions,
        except where a path holds more than ``MAX_DISTINCT_VALUES`` distinct values:
        which ones are kept then depends on the order slices are merged in.
        """
        self.reaction_count += other.reaction_count
        self.field_counts.update(other.field_counts)
        for path, counter in other.values.items():
            for value, count in counter.items():
                self._count_value(path, value, count)
        self.untracked.update(other.untracked)
        for path, counter in other.reaction_values.items():
            for values, count in counter.items():
                self._count_reaction_values(path, values, count)
        self.reaction_values_untracked.update(other.reaction_values_untracked)
        for path, counter in other.ordered_values.items():
            for values, count in counter.items():
                _count_capped(
                    self.ordered_values,
                    self.ordered_values_untracked,
                    path,
                    values,
                    count,
                )
        self.ordered_values_untracked.update(other.ordered_values_untracked)
        for path, counter in other.map_keys.items():
            self.map_keys[path].update(counter)
        self.shapes.update(other.shapes)
        for shape, examples in other.shape_examples.items():
            kept = self.shape_examples.setdefault(shape, [])
            for reaction_id in examples:
                if len(kept) < MAX_SHAPE_EXAMPLES and reaction_id not in kept:
                    kept.append(reaction_id)
        self.compounds.merge(other.compounds)

    def coverage(self, path: str) -> float:
        """The fraction of reactions that populate ``path``, from 0 to 1."""
        if not self.reaction_count:
            return 0.0
        return self.field_counts[path] / self.reaction_count

    def is_constant(self, path: str) -> bool:
        """Whether every reaction holds exactly the same values at ``path``.

        A submessage path is constant when every reaction populates it; whether its
        contents vary is a question for the paths below it.
        """
        if self.field_counts[path] != self.reaction_count or not self.reaction_count:
            return False
        if path not in self.reaction_values:
            return True
        return (
            len(self.reaction_values[path]) == 1
            and not self.reaction_values_untracked[path]
        )

    def constant_fields(self, *, by_key: bool = False) -> dict[str, ReactionValues]:
        """Scalar paths where every reaction holds the same values, with those values.

        With ``by_key``, fields below a map are reported per entry, by keyed path.

        Submessage paths are left out: they carry no value of their own. Each value is
        the sorted tuple every reaction holds, e.g. ``("CELSIUS",)`` for one setpoint,
        or one entry per component for a path below ``inputs.components``.
        """
        return {
            path: next(iter(self.reaction_values[path]))
            for path in self._value_paths(by_key=by_key)
            if self.is_constant(path)
        }

    def varying_fields(self, *, by_key: bool = False) -> list[str]:
        """Scalar paths that differ between reactions, or that some leave empty.

        With ``by_key``, fields below a map are reported per entry, by keyed path, so
        a screen's variables come out separately: ``inputs["base"]...`` and
        ``inputs["alcohol"]...`` rather than one ``inputs.components`` path.
        """
        return [
            path
            for path in self._value_paths(by_key=by_key)
            if self.field_counts[path] and not self.is_constant(path)
        ]

    def reordered_fields(self, *, by_key: bool = False) -> dict[str, int]:
        """Paths where reactions hold the same values in a different order.

        This is what a reviewer should ask about: code reading
        ``inputs["KEY"].components[0]`` gets a different compound depending on the
        reaction, and recording the dataset in one order would have avoided it. On
        merged paths it is also how a renamed map key shows, since map entries are
        read in key order; ``by_key`` narrows it to the order within each entry.

        Returns:
            Each such path mapped to how many reactions differ from the most common
            order of their values.
        """
        found = {}
        for path in self._value_paths(by_key=by_key):
            orders: dict[ReactionValues, list[int]] = defaultdict(list)
            for values, count in self.ordered_values[path].items():
                orders[tuple(sorted(values, key=_sort_key))].append(count)
            differing = sum(sum(c) - max(c) for c in orders.values() if len(c) > 1)
            if differing:
                found[path] = differing
        return found

    def ordered_shapes(self) -> list[tuple[Shape, int]]:
        """Each shape with its reaction count, most common first.

        Ties are broken by the shape's sorted paths, so the order is stable.
        """
        return sorted(self.shapes.items(), key=lambda item: (-item[1], sorted(item[0])))
