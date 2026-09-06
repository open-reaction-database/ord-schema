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

"""A queryable projection of the whole Reaction proto.

Source datasets store reactions as opaque wire-format bytes, so any question about them
costs a full deserialization. A projection restates the *same* reactions as nested
Parquet -- message to ``STRUCT``, repeated field to ``LIST``, map to ``MAP`` -- so a
query engine can read one leaf without touching the rest.

**No field is left out.** Every field of every message reachable from ``Reaction``
appears, because a field left out is a question nobody can ask, and the point of the
artifact is to support questions nobody enumerated in advance. The schema is generated
from the proto descriptors rather than written by hand, so a new field appears here
without anyone deciding it is worth carrying.

That works only because the message graph is acyclic: no message reaches itself, so the
recursion terminates. ``build_schema`` enforces this rather than assuming it -- a
recursive message added upstream would otherwise recurse forever at import time.

Three normalizations are applied, and only three. Each costs no query and removes a
real trap:

* **United messages become canonical floats.** A ``{value, precision, units}`` triple
  projected verbatim means ``WHERE temperature > 350`` silently misses every row
  recorded in Celsius. Each becomes two columns named for its unit --
  ``setpoint_kelvin`` and ``setpoint_precision_kelvin`` -- so a question in other units
  stays expressible while the mixed-unit comparison that would quietly return the wrong
  rows does not. The source records precision in the same units as the value, so it
  converts with the value and is null under the same conditions.
* **Dates become instants.** ``DateTime`` is a bare string the depositor wrote, and the
  corpus holds fifteen shapes of it in nine families, so comparing the text answers
  arbitrarily: over ORD a prefix match on ``2021`` finds none of that year's 490,412
  records while one on ``2022`` matches 1.77 million from a dataset whose stamps merely
  begin with it. Each ``DateTime`` gains a ``timestamp`` beside the recorded string,
  which stays: the projection adds a reading rather than replacing what the source says,
  and a shape none of the formats cover is a visible null.

  ``NN/NN/NNNN`` is the part that needs more than a format list, since the string never
  says which field is the month. It is settled per dataset and before any row is
  written -- the convention belongs to the depositor, not to a value -- from
  ``_DAY_FIRST``, from the dataset's own unambiguous values, and from month first
  otherwise. The two datasets ``_UNDECIDED`` names are the exception: nothing settles
  them, so their timestamps are null, because a guess written here is indistinguishable
  afterwards from a date the source recorded.

  The string the depositor wrote stays beside the timestamp, and is internal. Keeping it
  is what lets those nulls be read back once their submitters answer; hiding it is what
  keeps a date question off a column the corpus spells ten ways.

* **Structural identifiers collapse to one ``smiles``.** ``SMILES``, ``CXSMILES``,
  ``INCHI``, and ``MOLBLOCK`` all answer "what is this molecule," so the projection
  answers it once, through ``message_helpers.smiles_from_compound`` -- the same
  CXSMILES-first choice the tabular view makes, so the two never disagree about a
  molecule. ``REACTION_SMILES`` and ``REACTION_CXSMILES`` collapse the same
  way at the reaction level, into a ``smiles`` that
  ``message_helpers.derived_reaction_smiles`` normalizes from the recorded value --
  agents removed and the result canonicalized, so a reaction deposited as ``A>B>C``
  reads as ``A>>C`` -- or generates from the components when the source records nothing
  readable. Atom mapping is kept; fragment grouping is not.

  The drop is all-or-nothing per message and conditional on the collapse producing
  something: a compound whose identifiers none of RDKit's loaders can read keeps all of
  them, rather than reading as a compound whose structure the source never recorded.

Every other identifier is kept, as a list. ``NAME`` alone covers compounds that no
structural identifier reaches, and a compound may carry several of them, so pivoting
identifiers into named scalar fields would silently drop data.

One derived column rides beside the compound-level ``smiles``: a ``structure_id``
numbering the dataset's distinct structures in first-seen order, the join key into the
dataset's ``structures`` artifact (:mod:`ord_schema.artifacts.structures`), where the
fingerprints live. It exists because a structure predicate is evaluated *outside* the
query -- screened and verified against the structures artifact -- and its match set has
to re-enter the query at element granularity: DuckDB does not allow a subquery inside a
lambda expression, so the components a quantifier binds cannot semi-join against a match
table, but they can test one integer against a bitmap parameter. The IDs are meaningful
only alongside the projection that assigned them, so the column is marked internal and
stays out of what a model is told it may query.

Anything else -- the original units, the identifier types that were collapsed -- remains
in the source ``reaction`` column, which stays authoritative for byte-exact
round-tripping. The projection promises that every field is *readable*, not that it is
reproducible from the projection alone.

.. warning::
   Query repeated levels with list lambdas, not ``UNNEST`` in the ``FROM`` clause.
   ``UNNEST`` materializes the exploded intermediate: measured in DuckDB, 27-200x slower
   where it finishes at all, and over the full corpus 0.9 seconds against no result in
   four minutes. It also emits one row per matching identifier rather than one per
   reaction, so the two spellings below need the ``DISTINCT`` to count alike.

   .. code-block:: sql

      -- fast
      SELECT reaction_id FROM p
      WHERE len(list_filter(
              flatten(list_transform(map_values(inputs), i -> i.components)),
              c -> len(list_filter(c.identifiers,
                                   x -> x.type = 'NAME' AND x.value = 'THF')) > 0)) > 0

      -- same answer, does not finish
      SELECT DISTINCT reaction_id
      FROM p, UNNEST(map_values(inputs)) t(i), UNNEST(i.components) u(c),
              UNNEST(c.identifiers) v(x)
      WHERE x.type = 'NAME' AND x.value = 'THF'
"""

import os
import re
from collections.abc import Iterable, Iterator, MutableMapping
from datetime import datetime
from typing import Any, cast, get_args

import pyarrow as pa
import pyarrow.parquet as pq
from google.protobuf.descriptor import Descriptor, FieldDescriptor
from google.protobuf.message import Message

import ord_schema
from ord_schema import atomic_io, message_helpers, parquet, units
from ord_schema.artifacts import base
from ord_schema.logging import get_logger
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)

ARTIFACT = "projection"

# Canonical unit per united message type, and the suffix its column carries. Keyed by
# message name because the same target applies wherever the message appears. Only types
# reachable from Reaction are listed; Concentration is united but no field uses it.
_CANONICAL_UNITS: dict[str, tuple[str, str]] = {
    "Current": ("A", "amperes"),
    "FlowRate": ("mL/min", "milliliters_per_minute"),
    "Length": ("cm", "centimeters"),
    "Mass": ("g", "grams"),
    "Moles": ("mol", "moles"),
    "Pressure": ("kPa", "kilopascals"),
    "Temperature": ("K", "kelvin"),
    "Time": ("s", "seconds"),
    "Voltage": ("V", "volts"),
    "Volume": ("L", "liters"),
    "Wavelength": ("nm", "nanometers"),
}

# Messages that get a collapsed `smiles`, and the identifier types it replaces. Keying
# the sets keeps the decision to collapse and the choice of what to drop from drifting
# apart: the enum numbers overlap between compounds and reactions, so a mismatch would
# drop the wrong types silently rather than raise. The compound set follows
# message_helpers, which owns what "structural" means.
_STRUCTURAL_COMPOUND_TYPES = message_helpers.STRUCTURAL_IDENTIFIER_TYPES
_STRUCTURAL_REACTION_TYPES = frozenset(
    {
        reaction_pb2.ReactionIdentifier.REACTION_SMILES,
        reaction_pb2.ReactionIdentifier.REACTION_CXSMILES,
    }
)
_STRUCTURAL_TYPES: dict[str, frozenset[int]] = {
    "Compound": _STRUCTURAL_COMPOUND_TYPES,
    "ProductCompound": _STRUCTURAL_COMPOUND_TYPES,
    "Reaction": _STRUCTURAL_REACTION_TYPES,
}

# Messages whose collapsed ``smiles`` is a molecule, and that therefore carry a
# ``structure_id``. Reaction is structural but its smiles is a reaction;
# reaction-level structure search is a different operation and no ID is assigned there.
_STRUCTURE_ID_TYPES = frozenset({"Compound", "ProductCompound"})

# Messages whose free-text date gains a parsed ``timestamp`` beside it.
_TIMESTAMP_TYPES = frozenset({"DateTime"})
TIMESTAMP_COLUMN = "timestamp"

# The shapes ORD's dates are recorded in, tried in order. DateTime.value is a string the
# depositor wrote, so there is no format to rely on: seven appear across the corpus, and
# a reader comparing the text instead answers arbitrarily -- over ORD, a prefix match on
# "2021" finds none of the 490,412 records from that year while one on "2022" matches
# 1.77 million from a dataset whose stamps merely begin with it.
#
# Slash dates are listed month-first here and read day-first where the dataset says so;
# see ``day_first_dates``.
_DATE_FORMATS = (
    "%Y-%m-%d %H:%M:%S.%f",
    "%Y-%m-%dT%H:%M:%S.%f",
    "%Y-%m-%dT%H:%M:%S",
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d",
    "%m/%d/%Y, %I:%M:%S %p",
    "%m/%d/%Y, %H:%M:%S",
    "%m/%d/%Y %H:%M:%S",
    "%m/%d/%Y",
    "%a %b %d %H:%M:%S %Y",
)
# Written day-first for a dataset that records dates that way; same order, so a value
# matches the same shape either way and only the two leading numbers swap meaning.
_DAY_FIRST_FORMATS = tuple(shape.replace("%m/%d/", "%d/%m/") for shape in _DATE_FORMATS)
_SLASH_DATE = re.compile(r"^(\d{1,2})/(\d{1,2})/")

# Datasets that write their slash dates day first. Everything else reads month first:
# that is what a dataset's own values witness wherever one holds a number above twelve,
# and what the corpus's other evidence settles for the rest -- an upper bound from the
# reaction's own record_modified events, a co-submitted sibling in the same shape, the
# en-US 12-hour rendering no day-first locale produces, and the supplemental data on two
# submission pull requests. The reasoning per dataset is the ord-logbook entry "Date and
# time formats across the corpus" (2026-09-03).
_DAY_FIRST = frozenset(
    {
        "ord_dataset-172039a759a440219a68af62d203b79b",
        "ord_dataset-2be11f57f3304e678ea8469baf6dd1bc",
        "ord_dataset-3b8a2ef300e145468579027f206a3ac8",
        "ord_dataset-c5b00523487a4211a194160edf45e9ab",
    }
)

# Datasets no evidence settles, whose submitters have been asked. Their ambiguous dates
# project as null, which is the point: writing a guess would freeze it into the artifact
# and destroy the record that it was ever a guess. Between them they hold every
# ambiguous slash date the corpus records outside ``_DAY_FIRST`` -- 9,656 values, 0.40%
# of the 2,430,268 dates it records.
_UNDECIDED = frozenset(
    {
        "ord_dataset-5c9a10329a8a48968d18879a48bb8ab2",
        "ord_dataset-5e8318f0dda04b398c14f4c3adfeb32c",
    }
)

# A dataset in both lists reads as undecided, silently: nothing else says that the
# evidence naming it day first and the evidence settling nothing are the same claim.
if _DAY_FIRST & _UNDECIDED:
    raise ValueError(
        "a dataset cannot be both day-first and undecided: "
        f"{sorted(_DAY_FIRST & _UNDECIDED)}"
    )

# An ID rides beside the collapsed smiles and cannot exist without it: the schema
# would gain a structure_id with no sibling smiles, and message_row would KeyError.
if not set(_STRUCTURAL_TYPES) >= _STRUCTURE_ID_TYPES:
    raise ValueError(
        "_STRUCTURE_ID_TYPES must be a subset of _STRUCTURAL_TYPES; "
        f"{sorted(_STRUCTURE_ID_TYPES - set(_STRUCTURAL_TYPES))} has no smiles"
    )

_ARROW_SCALARS: dict[int, pa.DataType] = {
    FieldDescriptor.TYPE_DOUBLE: pa.float64(),
    FieldDescriptor.TYPE_FLOAT: pa.float32(),
    FieldDescriptor.TYPE_INT64: pa.int64(),
    FieldDescriptor.TYPE_UINT64: pa.uint64(),
    FieldDescriptor.TYPE_INT32: pa.int32(),
    FieldDescriptor.TYPE_UINT32: pa.uint32(),
    FieldDescriptor.TYPE_BOOL: pa.bool_(),
    FieldDescriptor.TYPE_STRING: pa.string(),
    FieldDescriptor.TYPE_BYTES: pa.binary(),
    # Enums project as their value name rather than their number: readable in a query,
    # and stable across renumbering in a way the integer is not.
    FieldDescriptor.TYPE_ENUM: pa.string(),
}

_RESOLVER = units.RESOLVER

# Field metadata naming the values an enum column may hold. An enum projects as its
# value *name*, and Arrow has no type that records the choices -- a dictionary column
# carries only the values a batch happened to contain, and DuckDB reads it back as
# VARCHAR regardless. Metadata does survive, nested to any depth, which makes the
# published Parquet self-describing: a consumer learns the spellings from the footer
# rather than from this library.
_META_ENUM = "ord.enum"

# Field metadata marking a column as artifact-internal: machinery the library uses to
# join artifacts together, not a fact about the reaction. Internal columns stay out of
# the schema description a model queries against, and a query naming one is refused --
# their values are not stable across builds, so nothing outside this library should
# depend on them.
_META_INTERNAL = "ord.internal"


def _validate_canonical_units() -> None:
    """Checks ``_CANONICAL_UNITS`` against the resolver and the reachable messages.

    Both failures are invisible at runtime, which is why they are caught here. A target
    the resolver rejects nulls every value of that column while leaving its name intact;
    a united message missing from the table projects as a raw ``{value, precision,
    units}`` struct, reinstating the mixed-unit trap.

    Raises:
        ValueError: If a target unit does not resolve to the message type it is listed
            under, or if the table does not cover exactly the united messages reachable
            from Reaction.
    """
    for name, (target, _) in _CANONICAL_UNITS.items():
        try:
            message_cls, _ = _RESOLVER.resolve_unit(target)
        except KeyError as error:
            raise ValueError(f"{name}: unusable canonical unit {target!r}") from error
        resolved = cast(Descriptor, message_cls.DESCRIPTOR)
        if resolved.name != name:
            raise ValueError(
                f"{name}: canonical unit {target!r} belongs to {resolved.name}"
            )
    united = {message.DESCRIPTOR.name for message in get_args(ord_schema.UnitMessage)}
    reachable = {
        field.message_type.name
        for field in _reachable_fields()
        if field.message_type is not None
    }
    expected = united & reachable
    if expected != set(_CANONICAL_UNITS):
        raise ValueError(
            "_CANONICAL_UNITS does not match the united messages reachable from "
            f"Reaction; missing {sorted(expected - set(_CANONICAL_UNITS))}, "
            f"unreachable {sorted(set(_CANONICAL_UNITS) - expected)}"
        )


def _reachable_fields() -> Iterator[FieldDescriptor]:
    """Yields every field of every message reachable from Reaction, once per field."""
    seen: set[str] = set()
    stack = [cast(Descriptor, reaction_pb2.Reaction.DESCRIPTOR)]
    while stack:
        descriptor = stack.pop()
        if descriptor.full_name in seen:
            continue
        seen.add(descriptor.full_name)
        for field in descriptor.fields:
            yield field
            if field.message_type is not None:
                stack.append(field.message_type)


def day_first_dates(values: Iterable[str | None]) -> bool | None:
    """Returns whether a dataset writes slash dates day first, or None if it cannot say.

    The convention is a property of the depositor rather than of a value, and a value
    only reveals it where one of the two leading numbers exceeds twelve. Measured over
    ORD no dataset carries evidence both ways, so one answer per dataset is the right
    shape -- and thirteen carry none at all, which is what None is for: those hold one
    bulk-import stamp repeated across every row, so a date guessed for them would be
    guessed for tens of thousands of reactions at once.

    Args:
        values: Every date string the dataset records, in any order. Nones and values in
            other formats are ignored, since only a slash date is ambiguous.

    Returns:
        True where the dataset writes day/month, False where it writes month/day, and
        None where no value settles it.
    """
    for value in values:
        if value is None:
            continue
        matched = _SLASH_DATE.match(value)
        if matched is None:
            continue
        first, second = int(matched.group(1)), int(matched.group(2))
        if first > 12:
            return True
        if second > 12:
            return False
    return None


def slash_orientation(
    dataset_id: str | None, values: Iterable[str | None]
) -> bool | None:
    """Returns whether a dataset writes slash dates day first, preferring what is known.

    ``_UNDECIDED`` answers first, because a dataset nothing settles has to stay unparsed
    rather than take the default. ``_DAY_FIRST`` answers next, holding evidence a value
    cannot carry -- a sibling submission alongside the datasets that witness themselves.
    The scan then answers for a dataset nobody has looked at yet, which is what keeps a
    day-first dataset added later from silently reading month first. Month first is what
    is left, and what every settled dataset outside ``_DAY_FIRST`` reads.

    Args:
        dataset_id: The source's ``ord_dataset-`` ID, or None where it records none.
        values: Every date string the dataset records. Read whether or not a list names
            the dataset, since a listed verdict is checked against any witness among
            them -- but read lazily, so a dataset that produces one stops there. A
            dataset settled by something other than a witness has none to find, and that
            scan runs to the end: ~7 seconds over the source that is 96% of ORD, against
            the 34.6 minutes the projection itself takes.

    Returns:
        True for day first, False for month first, and None for a dataset ``_UNDECIDED``
        names or a source recording no ID at all -- which leaves ambiguous dates
        unparsed rather than guessed.

    Raises:
        ValueError: If a witness in the data contradicts ``_DAY_FIRST``. One of the two
            is then wrong about the dataset, and projecting either reading would record
            a date this library has evidence against.
    """
    if dataset_id in _UNDECIDED:
        return None
    if dataset_id is None:
        # The default is a claim about a dataset someone could name in _DAY_FIRST after
        # reading its submission. A source recording no ID can never be named, so the
        # guess would be permanent and there would be nothing to correct it against;
        # only its own values may speak for it.
        return day_first_dates(values)
    witnessed = day_first_dates(values)
    if dataset_id in _DAY_FIRST:
        if witnessed is False:
            raise ValueError(
                f"{dataset_id} is recorded as day-first, but its own values witness "
                "month-first"
            )
        return True
    return False if witnessed is None else witnessed


def parse_timestamp(value: str, *, day_first: bool | None = None) -> datetime | None:
    """Returns the instant a recorded date names, or None where none of the shapes fit.

    Args:
        value: The string a depositor wrote.
        day_first: Whether this dataset writes slash dates day first, from
            ``day_first_dates``. None leaves an ambiguous slash date unparsed rather
            than guessing, which is the difference between a null a reader can see and
            a date that is silently the wrong month.

    Returns:
        The parsed instant, or None where no format matched or the value is ambiguous
        and the dataset does not say which way to read it.
    """
    matched = _SLASH_DATE.match(value)
    if matched is not None:
        first, second = int(matched.group(1)), int(matched.group(2))
        if first <= 12 and second <= 12 and day_first is None:
            return None
        formats = _DAY_FIRST_FORMATS if day_first else _DATE_FORMATS
        # An unambiguous value reads the same under either table, so the dataset's
        # convention decides only the values that need it.
        if first > 12:
            formats = _DAY_FIRST_FORMATS
        elif second > 12:
            formats = _DATE_FORMATS
    else:
        formats = _DATE_FORMATS
    for shape in formats:
        try:
            # DTZ007: naive on purpose. None of the recorded shapes carries an offset,
            # so the zone is not in the data; attaching one would state a fact the
            # source never recorded, and shift every instant by it.
            return datetime.strptime(value, shape)  # noqa: DTZ007
        except ValueError:
            continue
    return None


def _canonical_unit(field: FieldDescriptor) -> tuple[str, str] | None:
    """Returns ``field``'s (target unit, column suffix), or None if it is not united."""
    if field.message_type is None:
        return None
    return _CANONICAL_UNITS.get(field.message_type.name)


def column_name(field: FieldDescriptor) -> str:
    """Returns the column name for ``field``, carrying the unit where one applies."""
    canonical = _canonical_unit(field)
    if canonical is not None:
        return f"{field.name}_{canonical[1]}"
    return field.name


def precision_column_name(field: FieldDescriptor) -> str:
    """Returns the column name for ``field``'s precision.

    The unit stays the trailing element, matching the value column, so a reader takes
    the units off the end of either name without knowing which is which.

    Args:
        field: A united field, as identified by ``_canonical_unit``.

    Returns:
        The precision column's name.

    Raises:
        ValueError: If ``field`` is not united, and so has no precision to name.
    """
    canonical = _canonical_unit(field)
    if canonical is None:
        raise ValueError(f"{field.full_name} is not a united message")
    return f"{field.name}_precision_{canonical[1]}"


def _enum_metadata(field: FieldDescriptor) -> dict[str, str] | None:
    """Returns metadata recording an enum's members, or None for any other type."""
    if field.type != FieldDescriptor.TYPE_ENUM:
        return None
    return {_META_ENUM: ",".join(value.name for value in field.enum_type.values)}


def enum_members(field: pa.Field) -> tuple[str, ...] | None:
    """Returns the values ``field`` may hold, or None if it is not an enum column.

    Args:
        field: A field of ``SCHEMA``, or of a struct within it.

    Returns:
        The member names in declaration order, or None.
    """
    raw = (field.metadata or {}).get(_META_ENUM.encode())
    return tuple(raw.decode().split(",")) if raw else None


def is_internal(field: pa.Field) -> bool:
    """Returns whether ``field`` is artifact-internal machinery rather than data.

    Args:
        field: A field of ``SCHEMA``, or of a struct within it.

    Returns:
        True for columns like ``structure_id`` whose values join artifacts together and
        are not stable across builds. They are hidden from the schema description a
        model sees, and a query naming one is refused.
    """
    return bool((field.metadata or {}).get(_META_INTERNAL.encode()))


def _field_type(field: FieldDescriptor, stack: frozenset[str]) -> pa.DataType:
    if _canonical_unit(field) is not None:
        # No field _struct_fields handles reaches here: it expands a united field into a
        # value column and a precision column, which is a pair of leaves rather than one
        # type, and refuses a repeated one outright. Only a map whose values are united
        # arrives, and that pair cannot represent one either.
        raise ValueError(
            f"{field.full_name} is a map-valued united message; a value column and a "
            "precision column represent only singular ones"
        )
    message_type = field.message_type
    if message_type is not None and message_type.GetOptions().map_entry:
        key_field = message_type.fields_by_name["key"]
        value_field = message_type.fields_by_name["value"]
        return pa.map_(_field_type(key_field, stack), _field_type(value_field, stack))
    if message_type is not None:
        inner = pa.struct(_struct_fields(message_type, stack))
    else:
        inner = _ARROW_SCALARS[field.type]
    if field.label == FieldDescriptor.LABEL_REPEATED:
        return pa.list_(inner)
    return inner


def _struct_fields(descriptor: Descriptor, stack: frozenset[str]) -> list[pa.Field]:
    if descriptor.full_name in stack:
        raise ValueError(
            f"{descriptor.full_name} reaches itself; a total projection of a recursive "
            "message does not terminate"
        )
    stack = stack | {descriptor.full_name}
    fields = []
    if descriptor.name in _STRUCTURAL_TYPES:
        fields.append(pa.field("smiles", pa.string()))
    if descriptor.name in _TIMESTAMP_TYPES:
        # Microseconds, because the source records them: 1.77M values in the corpus
        # are str(datetime) output carrying a fractional part, and a second-resolution
        # column would drop it while still reading as the instant the source wrote.
        fields.append(pa.field(TIMESTAMP_COLUMN, pa.timestamp("us")))
    if descriptor.name in _STRUCTURE_ID_TYPES:
        fields.append(
            pa.field(
                "structure_id",
                pa.uint32(),
                metadata={_META_INTERNAL: "joins the dataset's structures artifact"},
            )
        )
    for field in descriptor.fields:
        if _canonical_unit(field) is None:
            metadata = _enum_metadata(field)
            if descriptor.name in _TIMESTAMP_TYPES and field.name == "value":
                # The depositor's own spelling, kept because the timestamp beside it is
                # an interpretation and cannot be turned back into it -- and is null for
                # the 9,656 ambiguous dates no evidence settles, which would otherwise
                # have no representation at all. Internal because a date question is
                # answered by the timestamp: matching text against a date the corpus
                # spells ten ways finds a fraction of what it asks for and reads as if
                # it found all of it.
                metadata = {
                    **(metadata or {}),
                    _META_INTERNAL: "the recorded spelling; query the timestamp",
                }
            fields.append(
                pa.field(
                    column_name(field),
                    _field_type(field, stack),
                    metadata=metadata,
                )
            )
            continue
        if field.label == FieldDescriptor.LABEL_REPEATED:
            # Two scalar columns cannot carry a list of measurements, and expanding one
            # anyway would publish a schema that quietly holds a single value per
            # reaction where the source holds several.
            raise ValueError(
                f"{field.full_name} is a repeated united message; a value column and a "
                "precision column represent only singular ones"
            )
        fields.append(pa.field(column_name(field), pa.float64()))
        fields.append(pa.field(precision_column_name(field), pa.float64()))
    return fields


def build_schema() -> pa.Schema:
    """Returns the Arrow schema projected from the Reaction descriptor.

    Returns:
        One top-level field per Reaction field, plus the collapsed ``smiles``.

    Raises:
        ValueError: If the message graph contains a cycle, which a total projection
            cannot represent.
    """
    # The protobuf stubs type DESCRIPTOR as optional; on a generated class it never is.
    descriptor = cast(Descriptor, reaction_pb2.Reaction.DESCRIPTOR)
    return pa.schema(_struct_fields(descriptor, frozenset()))


SCHEMA = build_schema()

_validate_canonical_units()


def _smiles(message: Message) -> str | None:
    """Returns the SMILES for a compound or reaction, or None if none can be read.

    Both forms prefer the CXSMILES the source recorded, so this column and the tabular
    view answer "what is this molecule" alike.

    Args:
        message: A Compound, ProductCompound, or Reaction.

    Returns:
        Canonical SMILES for a compound. For a reaction, the recorded reaction SMILES
        with its agents removed and canonicalized, or one generated from the components
        when nothing readable is recorded; None when neither is possible.
    """
    if cast(Descriptor, message.DESCRIPTOR).name == "Reaction":
        return message_helpers.derived_reaction_smiles(
            cast(reaction_pb2.Reaction, message)
        )
    return message_helpers.smiles_from_compound(cast(reaction_pb2.Compound, message))


def _kept_identifiers(values: Any, structural: frozenset[int] | None) -> list[Any]:
    """Returns the identifiers to project alongside a collapsed ``smiles``.

    Args:
        values: The message's ``identifiers`` field.
        structural: The types ``smiles`` answers for, or None when there is no
            ``smiles`` to answer for them.

    Returns:
        Every identifier when ``structural`` is None, otherwise those it does not name.
    """
    if structural is None:
        return list(values)
    return [value for value in values if value.type not in structural]


def _enum_name(field: FieldDescriptor, number: int) -> str:
    return field.enum_type.values_by_number[number].name


def message_row(
    message: Message,
    structure_ids: MutableMapping[str, int] | None = None,
    *,
    day_first: bool | None = None,
) -> dict[str, Any]:
    """Projects ``message`` to a dict matching its struct type in ``SCHEMA``.

    Unset fields are None rather than the proto default, so a consumer can tell "the
    source is silent" from "the source says zero".

    Args:
        message: Any message reachable from Reaction.
        structure_ids: Mapping from SMILES to ``structure_id``, extended in first-seen
            order as compounds are projected. None -- the default, for a message
            projected outside a dataset -- leaves every ``structure_id`` null, since an
            ID is meaningful only against the artifact that shares the mapping.
        day_first: Whether this dataset writes slash dates day first, from
            ``day_first_dates``. None leaves an ambiguous one unparsed; see
            ``parse_timestamp``.

    Returns:
        A dict keyed by projected column name.
    """
    descriptor = cast(Descriptor, message.DESCRIPTOR)
    row: dict[str, Any] = {}
    # None where smiles answered nothing, so a compound with no readable identifier
    # keeps all of them instead of reading as one that recorded no structure.
    collapsed: frozenset[int] | None = None
    if descriptor.name in _STRUCTURAL_TYPES:
        smiles = _smiles(message)
        row["smiles"] = smiles
        if smiles is not None:
            collapsed = _STRUCTURAL_TYPES[descriptor.name]
    if descriptor.name in _TIMESTAMP_TYPES:
        recorded = getattr(message, "value", "")
        row[TIMESTAMP_COLUMN] = (
            parse_timestamp(recorded, day_first=day_first) if recorded else None
        )
    if descriptor.name in _STRUCTURE_ID_TYPES:
        if row["smiles"] is not None and structure_ids is not None:
            row["structure_id"] = structure_ids.setdefault(
                row["smiles"], len(structure_ids)
            )
        else:
            row["structure_id"] = None
    for field in descriptor.fields:
        name = column_name(field)
        canonical = _canonical_unit(field)
        if canonical is not None:
            united = getattr(message, field.name)
            row[name] = units.canonical_value(united, canonical[0])
            row[precision_column_name(field)] = units.canonical_precision(
                united, canonical[0]
            )
            continue
        message_type = field.message_type
        if message_type is not None and message_type.GetOptions().map_entry:
            value_field = message_type.fields_by_name["value"]
            # Sorted by key, because protobuf leaves a map's iteration order undefined
            # and structure IDs are assigned in the order this walk reaches a SMILES.
            # Taking protobuf's order would make the IDs a function of the protobuf
            # build as well as of the source and the stamped library versions -- and
            # that build is not stamped, so an upgrade that reordered a map would
            # renumber a corpus whose artifacts all read as current, silently.
            items = sorted(getattr(message, field.name).items())
            if value_field.message_type is not None:
                row[name] = [
                    (key, message_row(value, structure_ids, day_first=day_first))
                    for key, value in items
                ] or None
            else:
                row[name] = list(items) or None
            continue
        if field.label == FieldDescriptor.LABEL_REPEATED:
            values = getattr(message, field.name)
            if field.name == "identifiers":
                values = _kept_identifiers(values, collapsed)
            if message_type is not None:
                row[name] = [
                    message_row(value, structure_ids, day_first=day_first)
                    for value in values
                ] or None
            elif field.type == FieldDescriptor.TYPE_ENUM:
                row[name] = [_enum_name(field, value) for value in values] or None
            else:
                row[name] = list(values) or None
            continue
        if message_type is not None:
            row[name] = (
                message_row(
                    getattr(message, field.name), structure_ids, day_first=day_first
                )
                if message.HasField(field.name)
                else None
            )
            continue
        if field.has_presence:
            present = message.HasField(field.name)
        else:
            # Proto3 scalars without explicit presence: absent and default are the same
            # wire state, so treat the default as absent.
            present = bool(getattr(message, field.name))
        if not present:
            row[name] = None
        elif field.type == FieldDescriptor.TYPE_ENUM:
            row[name] = _enum_name(field, getattr(message, field.name))
        else:
            row[name] = getattr(message, field.name)
    return row


def reaction_row(
    reaction: reaction_pb2.Reaction,
    reaction_id: str | None,
    structure_ids: MutableMapping[str, int] | None = None,
    *,
    day_first: bool | None = None,
) -> dict[str, Any]:
    """Projects one row of a source dataset to a dict matching ``SCHEMA``.

    The message is authoritative for every projected value, including the ID; the
    source's ``reaction_id`` column -- the key ``DatasetView.get_reaction`` indexes --
    is checked against it rather than chosen between. A check and not a fallback because
    ``DatasetView.md5`` hashes the Reaction blobs and not that column: an ID read from
    the column would sit outside the staleness hash, so correcting the column alone
    would leave the artifact looking current forever.

    Args:
        reaction: Reaction to project.
        reaction_id: The row's ``reaction_id`` column value, exactly as read -- None
            where the column is null. Required rather than defaulted, so a null cannot
            arrive looking like a caller who read no column and skip the check. Project
            a Reaction from anywhere else with ``message_row``.
        structure_ids: Mapping from SMILES to ``structure_id``, shared across a
            dataset's rows and extended in first-seen order. None leaves every
            ``structure_id`` null.
        day_first: Whether this dataset writes slash dates day first, from
            ``day_first_dates``. None leaves an ambiguous one unparsed; see
            ``parse_timestamp``.

    Returns:
        A dict keyed by projected column name.

    Raises:
        ValueError: If ``reaction_id`` is null or empty, or disagrees with the ID the
            message records. Each means the source cannot be joined on its own key.
    """
    if not reaction_id:
        raise ValueError(
            f"reaction_id column is {reaction_id!r}; the reaction records "
            f"{reaction.reaction_id!r}"
        )
    if reaction_id != reaction.reaction_id:
        raise ValueError(
            f"reaction_id column {reaction_id!r} disagrees with the reaction's own "
            f"{reaction.reaction_id!r}"
        )
    return message_row(reaction, structure_ids, day_first=day_first)


def is_current(path: str | os.PathLike[str], source_md5: str) -> bool:
    """Returns whether ``path`` is a projection of ``source_md5`` by this library."""
    return base.is_current(path, ARTIFACT, source_md5)


def _recorded_dates(view: parquet.DatasetView) -> Iterator[str | None]:
    """Yields the date strings a dataset records, for inferring its convention.

    The three provenance positions, not all seven the schema reaches: measured over the
    corpus, ``instrument_last_calibrated`` is unset in all 226,316 Analysis messages,
    and the other three sit behind repeated and map fields deep in inputs, workups and
    outcomes. Walking them would cost that descent on every reaction to find evidence
    the corpus does not hold. A dataset that populates one is where to look again.

    Args:
        view: The source dataset.

    Yields:
        The recorded value of each provenance DateTime, unparsed. Lazily, so a dataset
        that settles the question in its first reaction reads no further.
    """
    for _, reaction in view.iter_reactions():
        provenance = reaction.provenance
        yield provenance.experiment_start.value
        yield provenance.record_created.time.value
        for modified in provenance.record_modified:
            yield modified.time.value


def write_projection(
    source: str | os.PathLike[str],
    output: str | os.PathLike[str],
    *,
    source_md5: str | None = None,
    source_dataset_id: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 1000,
) -> int:
    """Projects a source Parquet dataset and writes it.

    Reactions are read and written one source row group at a time, so peak memory is
    bounded by the largest row group rather than the dataset. The output is published
    atomically, so a failure partway leaves any existing projection untouched.

    Args:
        source: Path to the source Parquet dataset.
        output: Path to write the projection to.
        source_md5: Source hash to stamp, if the caller already computed one. Hashed
            here when omitted, which costs a full pass over the source.
        source_dataset_id: Source dataset ID to stamp, if the caller already read one.
            Read from the source when omitted.
        compression: Parquet compression codec.
        row_group_size: Rows per output row group.

    Returns:
        Number of rows written.

    Raises:
        ValueError: If ``source`` is not a source dataset, or if any row's
            ``reaction_id`` column is missing or disagrees with its Reaction. The
            message names the source, since a corpus-wide run has thousands of datasets
            to choose between.
    """
    view = parquet.DatasetView(source)
    if source_md5 is None:
        source_md5 = view.md5()
    if source_dataset_id is None:
        source_dataset_id = view.dataset_id or None
    stamps = base.current_stamps(ARTIFACT, source_dataset_id, source_md5)
    schema = SCHEMA.with_metadata(base.to_metadata(stamps))
    rows = 0
    unreadable = 0
    # One ID space per dataset, in first-seen order, shared with the structures
    # artifact derived from this file. The order follows protobuf map iteration, which
    # is unspecified, so IDs are a fact about this file rather than about the dataset.
    structure_ids: dict[str, int] = {}
    # Which way this dataset writes slash dates, settled before any row is projected.
    # The convention belongs to the depositor rather than to a value, so it cannot be
    # read off the rows as they stream past and a row group's worth of them is not the
    # dataset. Where the answer is not already recorded this costs a pass that parses
    # the protos without projecting them: measured at ~7 seconds over the source that is
    # 96% of ORD, against the 34.6 minutes the projection itself takes.
    day_first = slash_orientation(source_dataset_id, _recorded_dates(view))
    with (
        atomic_io.atomic_path(output) as temp_path,
        pq.ParquetWriter(temp_path, schema, compression=compression) as writer,
    ):
        for row_group in range(view.num_row_groups):
            batch = []
            for reaction_id, reaction in view.iter_reactions(row_group=row_group):
                try:
                    row = reaction_row(
                        reaction, reaction_id, structure_ids, day_first=day_first
                    )
                except ValueError as error:
                    raise ValueError(f"{view.path}: {error}") from error
                unreadable += _unreadable_structures(row)
                batch.append(row)
            rows += len(batch)
            writer.write_table(
                pa.Table.from_pylist(batch, schema=schema),
                row_group_size=row_group_size,
            )
    if unreadable:
        # Per dataset rather than per row: a count is diffable between runs, so a
        # regression in structure reading shows up instead of hiding in the nulls.
        logger.info("%s: %d components project no structure", view.path, unreadable)
    return rows


def _unreadable_structures(row: dict[str, Any]) -> int:
    """Returns how many of a projected reaction's components have no ``smiles``."""
    components = []
    for _, reaction_input in row["inputs"] or []:
        components.extend(reaction_input["components"] or [])
    for outcome in row["outcomes"] or []:
        components.extend(outcome["products"] or [])
    return sum(component["smiles"] is None for component in components)
