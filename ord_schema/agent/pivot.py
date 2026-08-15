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

"""The repeated levels of the projection, and what a pivot over one holds.

A quantifier asks whether some element of a repeated level satisfies a body. Answering
that from the projection means decoding the whole nested column and filtering its
elements; answering it from a *pivot* -- one row per element -- means a semi-join
against a flat table. Over ORD the second is milliseconds where the first is seconds,
because the cost was never the predicate but the reconstruction of lists of structs.

A pivot row carries three things:

* ``reaction_id``, which is what the semi-join joins on.
* The ordinal of every repeated level from the root down to and including this one.
  Binding an element in the nested form gives co-membership structurally -- a
  measurement reached through a product *belongs* to it. A flat row has to say so, and
  a correlation joining on anything short of the full prefix over-returns.
* ``element``, a struct holding the element's own fields with every repeated field
  removed, recursively. Struct nesting is kept, so a path reaches the same leaf by the
  same spelling on either route, and a body reaching a dropped field simply fails to
  resolve -- which is how a quantifier the pivot cannot answer declines to the
  projection without a flag anyone has to keep honest.

Keeping the element's fields inside one struct rather than hoisting them is what makes
the key columns safe: ``inputs.components.preparations`` carries a ``reaction_id`` of
its own, and hoisted it would collide with the row's.

A pivot is unfiltered -- every element gets a row, including one whose fields are all
NULL. That completeness is what lets it answer ``forall``, which an index holding only
the elements that match something never can.
"""

import dataclasses

import inflection
import pyarrow as pa

from ord_schema import projection

# The column the element's own fields sit under, and the row's join key.
ELEMENT = "element"
REACTION_ID = "reaction_id"


@dataclasses.dataclass(frozen=True)
class Step:
    """One repeated level on the way to another, and how to reach it.

    Attributes:
        segments: The dotted hops from the previous element -- or from the row, for the
            first step -- down to the repeated field itself.
        is_map: Whether the repeated field is a map, whose values are the elements.
        ordinal: The column naming this element's position among its siblings.
    """

    segments: tuple[str, ...]
    is_map: bool
    ordinal: str

    def expression(self, source: str | None) -> str:
        """Returns the list to unnest at this step.

        Args:
            source: The variable bound to the previous element, or None at the row.

        Returns:
            A DuckDB expression evaluating to a list. A map contributes its values,
            since those are the elements a query quantifies over.
        """
        parts = self.segments if source is None else (source, *self.segments)
        reached = ".".join(parts)
        return f"map_values({reached})" if self.is_map else reached


@dataclasses.dataclass(frozen=True)
class RepeatedLevel:
    """A repeated level of the projection, and the shape of a pivot over it.

    Attributes:
        path: The level as the query grammar names it, with no wrapper segments.
        steps: One per repeated level from the root down to and including this one,
            outermost first. Unnesting them in order is what a build walks, and their
            ordinals are what a row carries to say which element it was.
        element_type: The element's struct type with repeated fields removed
            recursively, which is what a body must resolve against to be answerable
            from the pivot.
    """

    path: str
    steps: tuple[Step, ...]
    element_type: pa.DataType

    @property
    def ordinals(self) -> tuple[str, ...]:
        """Returns the ordinal column names, outermost first."""
        return tuple(step.ordinal for step in self.steps)


def _prune(dtype: pa.DataType) -> pa.DataType | None:
    """Returns ``dtype`` without its repeated fields, or None if none are left.

    Field metadata is carried through, since a comparison against an enum is checked
    against the members recorded there and a pruned type that dropped them would refuse
    a query the projection accepts.

    Args:
        dtype: A struct type to prune.

    Returns:
        The struct holding only the scalar and struct fields, recursively pruned, or
        None if pruning empties it -- a level whose elements are entirely repeated has
        nothing a pivot could hold.
    """
    fields = []
    for field in dtype:
        if pa.types.is_list(field.type) or pa.types.is_map(field.type):
            continue
        if pa.types.is_struct(field.type):
            pruned = _prune(field.type)
            if pruned is not None:
                fields.append(
                    pa.field(
                        field.name,
                        pruned,
                        nullable=field.nullable,
                        metadata=field.metadata,
                    )
                )
        else:
            fields.append(field)
    return pa.struct(fields) if fields else None


def repeated_levels(
    schema: pa.Schema = projection.SCHEMA,
) -> dict[str, RepeatedLevel]:
    """Returns every repeated level of ``schema``, by the path a query names it.

    Walked rather than listed by hand, in the same spirit as
    ``structures.structure_levels``: a repeated field added upstream becomes a level
    nobody had to remember to add, and one walk cannot disagree with itself the way two
    lists can.

    A list or a map is a repeated level; both are levels a query quantifies over, and
    the grammar spells neither with a wrapper segment, so ``inputs.components`` names
    the components under a map of inputs.

    Args:
        schema: Schema to walk; the projection schema by default.

    Returns:
        Each level's path mapped to what a pivot over it holds. A list of scalars is
        absent -- there is no element struct to pivot -- and so is a level whose
        elements carry nothing but repeated fields.

    Raises:
        ValueError: If two levels on one path yield the same ordinal name, which would
            give a row two columns of that name and let the inner one silently win.
            Raised at import, naming the ordinal, because the answer is to name it
            differently here rather than to retry.
    """
    levels: dict[str, RepeatedLevel] = {}

    def walk(
        dtype: pa.DataType,
        path: str,
        steps: tuple[Step, ...],
        pending: tuple[str, ...],
    ) -> None:
        if pa.types.is_list(dtype) or pa.types.is_map(dtype):
            is_map = pa.types.is_map(dtype)
            value = dtype.item_type if is_map else dtype.value_type
            ordinal = f"{inflection.singularize(path.rsplit('.', 1)[-1])}_index"
            if ordinal in [step.ordinal for step in steps]:
                raise ValueError(
                    f"{path} would carry two ordinals named {ordinal}, so a row's "
                    "outer position would be overwritten by its inner one"
                )
            steps = (*steps, Step(pending, is_map, ordinal))
            if pa.types.is_struct(value):
                element = _prune(value)
                if element is not None:
                    levels[path] = RepeatedLevel(path, steps, element)
            walk(value, path, steps, ())
        elif pa.types.is_struct(dtype):
            for field in dtype:
                walk(
                    field.type,
                    f"{path}.{field.name}" if path else field.name,
                    steps,
                    (*pending, field.name),
                )

    for field in schema:
        walk(field.type, field.name, (), (field.name,))
    return levels


def element_expression(dtype: pa.DataType, source: str) -> str:
    """Returns a struct literal rebuilding ``source`` without its repeated fields.

    The pivot stores what a predicate can read rather than the whole element, which is
    where the size comes down: the lists a build already unnested into their own levels
    would otherwise be carried again on every row above them.

    Args:
        dtype: The pruned element type, as ``RepeatedLevel.element_type`` gives it.
        source: The DuckDB expression bound to the element.

    Returns:
        A struct literal holding each kept field, recursing into kept structs.
    """
    parts = []
    for field in dtype:
        reached = f"{source}.{field.name}"
        value = (
            element_expression(field.type, reached)
            if pa.types.is_struct(field.type)
            else reached
        )
        parts.append(f"'{field.name}': {value}")
    return "{" + ", ".join(parts) + "}"


def table_name(path: str) -> str:
    """Returns the relation name holding the pivot over ``path``.

    Args:
        path: The level as the query grammar names it.

    Returns:
        A bare identifier, so it reaches SQL without quoting.
    """
    return f"pivot_{path.replace('.', '_')}"


LEVELS = repeated_levels()
