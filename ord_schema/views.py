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

A view is deliberately **small rather than complete**. It carries the handful of columns
a consumer wants on first contact, in a shape that reads at a glance and previews in a
dataset browser. Reaching the rest of the proto is a different artifact's job; a column
belongs here only if a reader who has never seen ORD would look for it immediately.
That is why breadth of dataset coverage, not row count, decides what stays: a column
populated in three of fifty-three datasets teaches a reader to expect nulls everywhere
and earns its place nowhere.

A view **adds no information**. Every column is a re-projection of what the source
protos already say, computed from the source alone with this library and RDKit, so a
view that disagrees with its source is wrong and a stale view is a bug. Two consequences
shape what may be added here:

* **No policy columns.** If a value depends on a threshold, cutoff, or classification
  that a reasonable consumer might set differently, it belongs in that consumer's query,
  not in the view. ``yield_percent < 5`` is theirs to write; ``is_negative_result``
  would be ours to impose.
* **No models, services, or weights.** Anything that cannot be recomputed offline from
  the source plus pinned open tooling -- reaction class labels, learned atom maps, a DOI
  expanded through Crossref -- is a different kind of artifact and does not belong in a
  view.

Where the projection has to make a choice, it is documented rather than inferred:

* Component SMILES prefer the CXSMILES form the source recorded, falling back to plain
  SMILES and then to any other structural identifier that parses (MOLBLOCK, InChI).
  CXSMILES is a superset -- RDKit reads either -- so preferring it keeps the enhanced
  stereochemistry that plain SMILES cannot express, at the cost of an extension block on
  some values; ``split_cxsmiles_extension`` splits the bare SMILES from the block for a
  consumer that wants them apart. Fragment grouping does not survive canonicalization.
  Canonical values are what make the column comparable across datasets that spell the
  same molecule differently.
* ``input_smiles`` carries every component of every input whatever its reaction role, so
  unlike ``reaction_smiles`` it lists solvents, reagents, and catalysts. Components with
  no readable structure are dropped from both list columns, which therefore are not a
  component census; the count of what was dropped is logged per dataset.
* ``reaction_smiles`` prefers the recorded ``REACTION_CXSMILES`` or ``REACTION_SMILES``,
  which carries what generation cannot reconstruct -- above all atom mapping. It is
  normalized rather than copied: agents are removed and the result canonicalized, so
  agent placement and atom ordering stop deciding whether two reactions look alike. A
  reaction recording neither, one RDKit cannot read, or one that survives normalization
  as a half reaction gets a SMILES generated from its components, again without agents
  and only when every reactant and product is readable. See
  ``message_helpers.derived_reaction_smiles``.
* Inputs are visited in sorted key order, so ``input_smiles`` is stable across runs
  rather than following protobuf map iteration, which is neither sorted nor insertion
  ordered. That order is the column's own; it does not line up with ``reaction_smiles``,
  which is deduplicated, sorted by SMILES, and usually read from a recorded identifier
  that never consulted the components at all.
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
* Measurements are converted to one canonical unit per column (kelvin, seconds). A
  measurement whose units the resolver cannot read, including one recorded with no
  units at all, reads as null: an unlabeled number has nowhere in the column to say
  what unit it is in.

The scalar and string columns are null where the source is silent, never zero. The two
list columns are empty rather than null, so a reaction whose components have no readable
structure is not distinguishable from one with no components by that column alone.
"""

import dataclasses
import math
import os

import pyarrow as pa
import pyarrow.parquet as pq

from ord_schema import artifacts, atomic_io, message_helpers, parquet, units
from ord_schema.logging import get_logger
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)

ARTIFACT = "view"

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


def _outcome_values(
    reaction: reaction_pb2.Reaction, skipped: _Skipped
) -> tuple[float | None, float | None]:
    """Returns the largest yield over all outcomes, and the first reaction time.

    Every outcome is read, matching ``output_smiles``: a yield recorded on a later
    outcome is a yield the source states, and reporting null for it would say the
    source is silent when it is not.

    Args:
        reaction: Reaction to read.
        skipped: Tallies to add to for values that cannot be used.

    Returns:
        ``(yield_percent, reaction_time_seconds)``, either of which is None where no
        outcome records a usable value.
    """
    best_yield = None
    reaction_time_seconds = None
    for outcome in reaction.outcomes:
        if reaction_time_seconds is None:
            reaction_time_seconds = units.canonical_value(outcome.reaction_time, "s")
        for product in outcome.products:
            for measurement in product.measurements:
                if measurement.type != reaction_pb2.ProductMeasurement.YIELD:
                    continue
                if not measurement.HasField("percentage"):
                    # A yield recorded as float_value, amount, or text. Reading it
                    # would require assuming a scale the source did not state.
                    skipped.unyielding += 1
                    continue
                if not measurement.percentage.HasField("value"):
                    # An unset value reads as 0.0 through protobuf, and a fabricated
                    # zero yield is worse than a null.
                    continue
                value = measurement.percentage.value
                if not math.isfinite(value):
                    # NaN compares False against everything, so leaving it in would
                    # make the result depend on the order measurements were recorded.
                    skipped.non_finite += 1
                    continue
                if best_yield is None or value > best_yield:
                    best_yield = value
    return best_yield, reaction_time_seconds


def _component_smiles(
    reaction: reaction_pb2.Reaction, skipped: _Skipped
) -> tuple[list[str], list[str]]:
    """Returns ``(input_smiles, output_smiles)``, skipping components without structure.

    Inputs are visited in sorted key order, so the column is stable across runs rather
    than following protobuf map iteration, which is neither sorted nor insertion
    ordered.

    Args:
        reaction: Reaction to read.
        skipped: Tallies to add to for components with no readable structure.

    Returns:
        Canonical SMILES for every input component and every outcome product that has a
        readable structure. A component without one is dropped, so the lists are not a
        component census.
    """
    inputs = []
    for key in sorted(reaction.inputs):
        for compound in reaction.inputs[key].components:
            smiles = message_helpers.smiles_from_compound(compound)
            if smiles:
                inputs.append(smiles)
            else:
                skipped.components += 1
    outputs = []
    for outcome in reaction.outcomes:
        for product in outcome.products:
            smiles = message_helpers.smiles_from_compound(product)
            if smiles:
                outputs.append(smiles)
            else:
                skipped.components += 1
    return inputs, outputs


def reaction_row(
    reaction: reaction_pb2.Reaction,
    reaction_id: str | None = None,
    skipped: _Skipped | None = None,
) -> dict:
    """Projects one Reaction onto the view columns.

    Args:
        reaction: Reaction message to project.
        reaction_id: The source's ``reaction_id`` column value, checked against the
            message rather than trusted. ``DatasetView.md5`` hashes the Reaction blobs
            and not that column, so an ID taken from it would be a published value
            outside the staleness hash. Defaults to the message's own ID.
        skipped: Tallies to add to for values this row drops. A caller deriving one row
            outside ``write_view`` can leave it unset.

    Returns:
        A dict keyed by ``SCHEMA.names``, with None wherever the source is silent.

    Raises:
        ValueError: If ``reaction_id`` is empty or disagrees with the message, either of
            which means the source cannot be joined on its own key.
    """
    if reaction_id is not None:
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
    if skipped is None:
        skipped = _Skipped()
    yield_percent, reaction_time_seconds = _outcome_values(reaction, skipped)
    inputs, outputs = _component_smiles(reaction, skipped)
    return {
        "reaction_id": reaction.reaction_id,
        "reaction_smiles": message_helpers.derived_reaction_smiles(reaction),
        "input_smiles": inputs,
        "output_smiles": outputs,
        "yield_percent": yield_percent,
        "temperature_kelvin": units.canonical_value(
            reaction.conditions.temperature.setpoint, "K"
        ),
        "reaction_time_seconds": reaction_time_seconds,
        "doi": reaction.provenance.doi or None,
        "patent": reaction.provenance.patent or None,
    }


def is_current(path: str | os.PathLike[str], source_md5: str) -> bool:
    """Returns whether ``path`` is a view of ``source_md5`` at these versions.

    Delegates to ``artifacts.is_current``, which requires the artifact name, the source
    content, the library version, and the artifact version to all match. A missing or
    unreadable file is not current.
    """
    return artifacts.is_current(path, ARTIFACT, source_md5)


def write_view(
    source: parquet.DatasetView,
    output: str | os.PathLike[str],
    *,
    source_md5: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 1000,
) -> int:
    """Derives a view from a source Parquet dataset and writes it.

    Reactions are read and written one source row group at a time, so peak memory is
    bounded by the largest row group rather than the dataset. The output is published
    atomically, so a failure partway leaves any existing view untouched.

    Args:
        source: View of the source Parquet dataset.
        output: Path to write the view to.
        source_md5: Source hash to stamp, if the caller already computed one (e.g. to
            check whether the view was current). Hashed here when omitted, which costs a
            full pass over the source.
        compression: Parquet compression codec.
        row_group_size: Rows per output row group.

    Returns:
        Number of rows written.
    """
    if source_md5 is None:
        source_md5 = source.md5()
    stamps = artifacts.current_stamps(ARTIFACT, source.dataset_id or None, source_md5)
    schema = SCHEMA.with_metadata(artifacts.to_metadata(stamps))
    rows = 0
    skipped = _Skipped()
    with (
        atomic_io.atomic_path(output) as temp_path,
        pq.ParquetWriter(temp_path, schema, compression=compression) as writer,
    ):
        for row_group in range(source.num_row_groups):
            # Seeded from SCHEMA.names so an unexpected key raises here rather than
            # being dropped silently, which from_pylist would do.
            batch: dict[str, list] = {name: [] for name in SCHEMA.names}
            for reaction_id, reaction in source.iter_reactions(row_group=row_group):
                try:
                    row = reaction_row(reaction, reaction_id, skipped)
                except Exception as error:
                    # A corpus run has thousands of datasets and millions of reactions;
                    # neither is in the traceback without this.
                    raise ValueError(
                        f"{source.path}: {reaction_id}: {error}"
                    ) from error
                for name, value in row.items():
                    batch[name].append(value)
                rows += 1
            writer.write_table(
                pa.Table.from_pydict(batch, schema=schema),
                row_group_size=row_group_size,
            )
    if skipped:
        # Per dataset rather than per row: a count is diffable between runs, so a
        # regression in structure reading shows up instead of hiding in the nulls.
        logger.info("%s: %s", source.path, skipped)
    return rows
