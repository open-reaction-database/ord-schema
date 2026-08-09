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

"""Tabular views over ORD datasets.

A view restates the reactions in a Parquet dataset as flat, queryable columns --
reaction and per-component SMILES, the scalar condition and outcome measurements, and
provenance -- so a consumer can filter the corpus without deserializing every proto.

A view is a **narrowing of the projection**, and is derived from one. Every column it
publishes is already a leaf there, so a second derivation would be a second code path
obliged to agree with the first about the same values, and it would pay for a full
deserialization and a canonicalization of every structure to reach what a handful of
column reads already hold. What is decided here is what the projection deliberately does
not: which of several plausible fields a column means, and how to reduce a repeated one
to a scalar.

A view is deliberately **small rather than complete**. It carries the handful of columns
a consumer wants on first contact, in a shape that reads at a glance and previews in a
dataset browser. Reaching the rest of the proto is a different artifact's job; a column
belongs here only if a reader who has never seen ORD would look for it immediately.
That is why breadth of dataset coverage, not row count, decides what stays: a column
populated in three of fifty-three datasets teaches a reader to expect nulls everywhere
and earns its place nowhere.

A view **adds no information**. Every column is a re-projection of what the source
protos already say, reached through the projection with this library, so a view that
disagrees with its source is wrong and a stale view is a bug. Two consequences shape
what may be added here:

* **No policy columns.** If a value depends on a threshold, cutoff, or classification
  that a reasonable consumer might set differently, it belongs in that consumer's query,
  not in the view. ``yield_percent < 5`` is theirs to write; ``is_negative_result``
  would be ours to impose.
* **No models, services, or weights.** Anything that cannot be recomputed offline from
  the source plus pinned open tooling -- reaction class labels, learned atom maps, a DOI
  expanded through Crossref -- is a different kind of artifact and does not belong in a
  view.

Where a choice has to be made, it is documented rather than inferred. The structural
ones belong to the projection and are inherited here, so the two artifacts cannot
disagree about what a molecule is:

* ``smiles`` on a component prefers the CXSMILES form the source recorded, falling back
  to plain SMILES and then to any other structural identifier that parses; the
  reaction's own prefers a recorded ``REACTION_CXSMILES`` or ``REACTION_SMILES``,
  normalized with its agents removed, and is generated from the components only when
  nothing readable was recorded. See :mod:`ord_schema.projection`.

The rest are the view's own, and are why it is not simply a column subset:

* ``input_smiles`` carries every component of every input whatever its reaction role, so
  unlike ``reaction_smiles`` it lists solvents, reagents, and catalysts. Components with
  no readable structure are dropped from both list columns, which therefore are not a
  component census; the count of what was dropped is logged per dataset.
* Inputs are visited in sorted key order, so ``input_smiles`` is stable across runs
  rather than following the projection's map order, which preserves protobuf map
  iteration and is neither sorted nor insertion ordered. That order is the column's own;
  it does not line up with ``reaction_smiles``, which is deduplicated, sorted by SMILES,
  and usually read from a recorded identifier that never consulted the components at
  all.
* ``yield_percent`` is the largest YIELD percentage on any product of any outcome, and
  ``reaction_time_seconds`` is the first outcome that records one -- matching the SMILES
  columns, which also span every outcome. A yield the source records only as
  ``float_value``, ``amount``, or text reads as null, since using it would assume a
  scale the source never stated. Values are published as recorded: the corpus holds
  yields outside [0, 100], so "largest" can surface an outlier.
* Each measurement column names exactly one source field, since the schema offers
  several plausible readings of each: ``temperature_kelvin`` is the ``conditions``
  *setpoint*, not an achieved value, and ``reaction_time_seconds`` is
  ``ReactionOutcome.reaction_time``, not any of the eight other Time-typed fields --
  addition, workup, observation, and retention times, or the timestamps on individual
  temperature, pressure, and electrochemistry measurements.
* Measurements are read from the projection's canonical-unit columns (kelvin, seconds),
  so a measurement whose units the resolver could not read arrives null: an unlabeled
  number has nowhere in the column to say what unit it is in.

The scalar and string columns are null where the source is silent, never zero. The two
list columns are empty rather than null, so a reaction whose components have no readable
structure is not distinguishable from one with no components by that column alone.
"""

import dataclasses
import math
import operator
import os
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from ord_schema import artifacts, atomic_io, projection
from ord_schema.logging import get_logger
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)

ARTIFACT = "view"

# The projection writes an enum as its value name, so a yield is matched by name.
# Resolved from the descriptor, so renaming the value upstream fails at import instead
# of leaving a literal here that matches nothing.
_YIELD = reaction_pb2.ProductMeasurement.ProductMeasurementType.Name(
    reaction_pb2.ProductMeasurement.YIELD
)

# The projection leaves a view reads, as Parquet column paths so only these are decoded.
# Every one is a value the projection already computed, which is what makes a view cheap
# to derive: no proto is deserialized and no structure is canonicalized again.
_PROJECTION_COLUMNS = (
    "reaction_id",
    "smiles",
    # Both halves of the map: the values carry the structures, and the keys are what
    # ``input_smiles`` sorts by.
    "inputs.key_value.key",
    "inputs.key_value.value.components.list.element.smiles",
    "outcomes.list.element.reaction_time_seconds",
    "outcomes.list.element.products.list.element.smiles",
    "outcomes.list.element.products.list.element.measurements.list.element.type",
    "outcomes.list.element.products.list.element.measurements.list.element"
    ".percentage.value",
    "conditions.temperature.setpoint_kelvin",
    "provenance.doi",
    "provenance.patent",
)

SCHEMA = pa.schema(
    [
        # The join key back to the source dataset, and to any other artifact keyed by
        # reaction. Never row position: views are regenerated freely.
        pa.field("reaction_id", pa.string(), nullable=False),
        pa.field("reaction_smiles", pa.string()),
        pa.field("input_smiles", pa.list_(pa.string())),
        pa.field("output_smiles", pa.list_(pa.string())),
        pa.field("yield_percent", pa.float32()),
        # Temperature and reaction time reach 45 and 49 of the 53 parquet datasets --
        # broader coverage than yield -- which is what qualifies them for a summary
        # view. Pressure and conversion reach 3 and 6, and stay out.
        pa.field("temperature_kelvin", pa.float32()),
        pa.field("reaction_time_seconds", pa.float32()),
        pa.field("doi", pa.string()),
        pa.field("patent", pa.string()),
    ]
)


@dataclasses.dataclass
class _Skipped:
    """Counts of source values a view drops, logged once per dataset.

    A view says nothing about what it could not read, so a column of nulls looks the
    same whether the corpus is silent or the reader broke. These make the difference
    visible without putting it in a column.
    """

    components: int = 0
    unyielding: int = 0
    non_finite: int = 0

    def __bool__(self) -> bool:
        return bool(self.components or self.unyielding or self.non_finite)

    def __str__(self) -> str:
        return (
            f"{self.components} components with no readable structure, "
            f"{self.unyielding} yields not recorded as a percentage, "
            f"{self.non_finite} non-finite yields"
        )


def _nested(row: Any, *keys: str) -> Any:
    """Returns the value at a chain of struct keys, or None if any level is absent.

    The projection writes an unset message as null rather than an empty struct, so a
    reaction that records no conditions at all has ``conditions`` itself missing, not a
    ``conditions`` holding a missing ``temperature``.

    Args:
        row: A projected row, or any struct within one.
        *keys: Field names to follow, outermost first.

    Returns:
        The value at the end of the chain, or None if the chain is broken.
    """
    for key in keys:
        if row is None:
            return None
        row = row[key]
    return row


def _outcome_values(
    projected: dict, skipped: _Skipped
) -> tuple[float | None, float | None]:
    """Returns the largest yield over all outcomes, and the first reaction time.

    Every outcome is read, matching ``output_smiles``: a yield recorded on a later
    outcome is a yield the source states, and reporting null for it would say the
    source is silent when it is not.

    Args:
        projected: One row of the projection.
        skipped: Tallies to add to for values that cannot be used.

    Returns:
        ``(yield_percent, reaction_time_seconds)``, either of which is None where no
        outcome records a usable value.
    """
    best_yield = None
    reaction_time_seconds = None
    for outcome in projected["outcomes"] or []:
        if reaction_time_seconds is None:
            reaction_time_seconds = outcome["reaction_time_seconds"]
        for product in outcome["products"] or []:
            for measurement in product["measurements"] or []:
                if measurement["type"] != _YIELD:
                    continue
                if measurement["percentage"] is None:
                    # A yield recorded as float_value, amount, or text. Reading it
                    # would require assuming a scale the source did not state.
                    skipped.unyielding += 1
                    continue
                value = measurement["percentage"]["value"]
                if value is None:
                    # A percentage message carrying no value. A fabricated zero yield
                    # is worse than a null.
                    continue
                if not math.isfinite(value):
                    # NaN compares False against everything, so leaving it in would
                    # make the result depend on the order measurements were recorded.
                    skipped.non_finite += 1
                    continue
                if best_yield is None or value > best_yield:
                    best_yield = value
    return best_yield, reaction_time_seconds


def _component_smiles(
    projected: dict, skipped: _Skipped
) -> tuple[list[str], list[str]]:
    """Returns ``(input_smiles, output_smiles)``, skipping components without structure.

    Inputs are visited in sorted key order, so the column is stable across runs rather
    than following the projection's map order, which preserves protobuf map iteration
    and is neither sorted nor insertion ordered.

    Args:
        projected: One row of the projection.
        skipped: Tallies to add to for components with no readable structure.

    Returns:
        The SMILES the projection derived for every input component and every outcome
        product that has one. A component without one is dropped, so the lists are not a
        component census.
    """
    inputs = []
    for _, value in sorted(projected["inputs"] or [], key=operator.itemgetter(0)):
        for component in value["components"] or []:
            if component["smiles"]:
                inputs.append(component["smiles"])
            else:
                skipped.components += 1
    outputs = []
    for outcome in projected["outcomes"] or []:
        for product in outcome["products"] or []:
            if product["smiles"]:
                outputs.append(product["smiles"])
            else:
                skipped.components += 1
    return inputs, outputs


def reaction_row(projected: dict, skipped: _Skipped | None = None) -> dict:
    """Narrows one row of the projection onto the view columns.

    The projection reconciled the reaction ID against the Reaction it projected, so a
    column disagreeing with its message cannot reach here.

    Args:
        projected: One row of the projection, holding at least ``_PROJECTION_COLUMNS``.
        skipped: Tallies to add to for values this row drops. A caller deriving one row
            outside ``write_view`` can leave it unset.

    Returns:
        A dict keyed by ``SCHEMA.names``, with None wherever the source is silent.
    """
    if skipped is None:
        skipped = _Skipped()
    yield_percent, reaction_time_seconds = _outcome_values(projected, skipped)
    inputs, outputs = _component_smiles(projected, skipped)
    return {
        "reaction_id": projected["reaction_id"],
        "reaction_smiles": projected["smiles"],
        "input_smiles": inputs,
        "output_smiles": outputs,
        "yield_percent": yield_percent,
        "temperature_kelvin": _nested(
            projected, "conditions", "temperature", "setpoint_kelvin"
        ),
        "reaction_time_seconds": reaction_time_seconds,
        "doi": _nested(projected, "provenance", "doi"),
        "patent": _nested(projected, "provenance", "patent"),
    }


def is_current(path: str | os.PathLike[str], source_md5: str) -> bool:
    """Returns whether ``path`` is a view of ``source_md5`` at these versions.

    Delegates to ``artifacts.is_current``, which requires the artifact name, the source
    content, the library version, and the artifact version to all match. A missing or
    unreadable file is not current.
    """
    return artifacts.is_current(path, ARTIFACT, source_md5)


def write_view(
    source: str | os.PathLike[str],
    output: str | os.PathLike[str],
    *,
    source_md5: str | None = None,
    source_dataset_id: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 1000,
) -> int:
    """Narrows a projection into a view and writes it.

    Only the leaves the view publishes are decoded, one projection row group at a time,
    so peak memory is bounded by the largest row group rather than the dataset. The
    output is published atomically, so a failure partway leaves any existing view
    untouched.

    Args:
        source: Path to the projection to narrow.
        output: Path to write the view to.
        source_md5: Hash of the *source dataset* to stamp, if the caller already read
            one. Taken from the projection's own stamps when omitted, so a view names
            the dataset it reflects rather than the artifact it read.
        source_dataset_id: Source dataset ID to stamp, if the caller already read one.
            Taken from the projection's stamps when omitted.
        compression: Parquet compression codec.
        row_group_size: Rows per output row group.

    Returns:
        Number of rows written.

    Raises:
        ValueError: If ``source`` is not a projection, or if a row cannot be narrowed.
    """
    # Read unconditionally rather than only to fill in what the caller omitted: naming
    # the wrong kind of file is worth a sentence, and the alternative is a KeyError on
    # a column the file was never going to have.
    parent = artifacts.load_stamps(source)
    if parent.artifact != projection.ARTIFACT:
        raise ValueError(
            f"{source} is a {parent.artifact}, not a {projection.ARTIFACT}; a view "
            "narrows a projection"
        )
    # derive_tree refuses stale parents, but this writer is public and its output
    # inherits the dataset hash: a view derived from a stale projection would claim a
    # provenance it does not have and nothing would ever mark it stale again.
    if not artifacts.stamps_are_current(parent, projection.ARTIFACT):
        raise ValueError(
            f"{source} is a stale {projection.ARTIFACT}; derive it again first"
        )
    if source_md5 is None:
        source_md5 = parent.source_md5
    if source_dataset_id is None:
        source_dataset_id = parent.source_dataset_id
    stamps = artifacts.current_stamps(ARTIFACT, source_dataset_id, source_md5)
    schema = SCHEMA.with_metadata(artifacts.to_metadata(stamps))
    rows = 0
    skipped = _Skipped()
    with (
        pq.ParquetFile(source) as projected,
        atomic_io.atomic_path(output) as temp_path,
        pq.ParquetWriter(temp_path, schema, compression=compression) as writer,
    ):
        for row_group in range(projected.num_row_groups):
            # Seeded from SCHEMA.names so an unexpected key raises here rather than
            # being dropped silently, which from_pylist would do.
            batch: dict[str, list] = {name: [] for name in SCHEMA.names}
            table = projected.read_row_group(
                row_group, columns=list(_PROJECTION_COLUMNS)
            )
            for row in table.to_pylist():
                try:
                    for name, value in reaction_row(row, skipped).items():
                        batch[name].append(value)
                except Exception as error:
                    # A corpus run has thousands of datasets and millions of reactions;
                    # neither is in the traceback without this.
                    raise ValueError(
                        f"{source}: {row.get('reaction_id')}: {error}"
                    ) from error
                rows += 1
            writer.write_table(
                pa.Table.from_pydict(batch, schema=schema),
                row_group_size=row_group_size,
            )
    if skipped:
        # Per dataset rather than per row: a count is diffable between runs, so a
        # regression in structure reading shows up instead of hiding in the nulls.
        logger.info("%s: %s", source, skipped)
    return rows
