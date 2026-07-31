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

* Every SMILES column prefers the CXSMILES form the source recorded, falling back to
  plain SMILES. CXSMILES is a superset -- RDKit reads either -- so preferring it keeps
  enhanced stereochemistry and fragment grouping that plain SMILES cannot express, at
  the cost of an extension block on some values. ``split_cxsmiles_extension`` returns
  the bare half for a consumer that wants it.
* ``reaction_smiles`` is the identifier the source recorded, generated from the
  components by ``message_helpers.get_reaction_smiles`` when it recorded none.
* Component SMILES are canonicalized, and derived from another structural identifier
  (MOLBLOCK, InChI) when a component records neither SMILES form. A component with
  nothing structural to work from is skipped. Canonical strings are what make the
  columns comparable across datasets that spell the same molecule differently.
* ``input_smiles`` visits inputs in sorted key order, matching reaction SMILES
  generation, so the column is stable across runs rather than following map iteration
  order.
* ``yield_percent`` is the largest YIELD percentage measured on any product of the first
  outcome. Reactions with several outcomes keep the rest in the source.
* Each measurement column names exactly one source field, since the schema offers
  several plausible readings of each: ``temperature_kelvin`` is the ``conditions``
  *setpoint*, not an achieved value, and ``reaction_time_seconds`` is
  ``ReactionOutcome.reaction_time`` -- one of six Time-typed fields, alongside addition,
  workup, observation, and retention times, none of which appear here.
* Measurements are converted to one canonical unit per column (kelvin, seconds); a
  measurement in units the schema cannot convert reads as null.

Every column except ``reaction_id`` is nullable, and null means "the source does not
say" rather than zero.
"""

import dataclasses
import os
from importlib import metadata

import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem

import ord_schema
from ord_schema import atomic_io, message_helpers, parquet, units
from ord_schema.proto import reaction_pb2

# Version of the view definition. Bump when a column's meaning or how it is populated
# changes, so that existing views read as stale.
VIEW_VERSION = "1"

# Footer keys share the "ord." namespace that ord_schema.parquet stamps source
# datasets with, so a reader sees one convention across both.
_META_VIEW_VERSION = "ord.view_version"
_META_SOURCE_DATASET_ID = "ord.source_dataset_id"
_META_SOURCE_MD5 = "ord.source_md5"
_META_ORD_SCHEMA_VERSION = "ord.ord_schema_version"

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

# CXSMILES is a superset of SMILES and RDKit reads either, so the richer form wins.
_PREFERRED_REACTION_TYPES = (
    reaction_pb2.ReactionIdentifier.REACTION_CXSMILES,
    reaction_pb2.ReactionIdentifier.REACTION_SMILES,
)
_PREFERRED_COMPOUND_TYPES = (
    reaction_pb2.CompoundIdentifier.CXSMILES,
    reaction_pb2.CompoundIdentifier.SMILES,
)

_RESOLVER = units.UnitResolver()


def _canonical_value(message: ord_schema.UnitMessage, target: str) -> float | None:
    """Returns ``message.value`` in ``target`` units, or None if it cannot be read.

    Args:
        message: A united message, e.g. Temperature or Pressure.
        target: Unit to convert to, as a string the resolver understands.

    Returns:
        The converted value, or None when the message carries no value or its units are
        unspecified.
    """
    if not message.HasField("value") or not message.units:
        return None
    return _RESOLVER.convert(message, target).value


def _outcome_values(
    reaction: reaction_pb2.Reaction,
) -> tuple[float | None, float | None]:
    """Returns the first outcome's yield and reaction time.

    Reactions with several outcomes keep the rest in the source.
    """
    if not reaction.outcomes:
        return None, None
    outcome = reaction.outcomes[0]
    reaction_time_seconds = _canonical_value(outcome.reaction_time, "s")
    best_yield = None
    for product in outcome.products:
        for measurement in product.measurements:
            if (
                measurement.type != reaction_pb2.ProductMeasurement.YIELD
                or not measurement.HasField("percentage")
                or not measurement.percentage.HasField("value")
            ):
                continue
            value = measurement.percentage.value
            if best_yield is None or value > best_yield:
                best_yield = value
    return best_yield, reaction_time_seconds


def _stored_value(
    message: reaction_pb2.Reaction
    | reaction_pb2.Compound
    | reaction_pb2.ProductCompound,
    preferred_types: tuple[int, ...],
) -> str | None:
    """Returns the first non-empty identifier value in order of preference."""
    for identifier_type in preferred_types:
        for identifier in message.identifiers:
            if identifier.type == identifier_type and identifier.value:
                return identifier.value
    return None


def _compound_smiles(
    compound: reaction_pb2.Compound | reaction_pb2.ProductCompound,
) -> str | None:
    """Returns a canonical SMILES for a compound, or None if it has no structure.

    Falls back to another structural identifier (MOLBLOCK, InChI) when the compound
    records neither SMILES form, matching what the rest of the project treats as a
    compound's structure.
    """
    stored = _stored_value(compound, _PREFERRED_COMPOUND_TYPES)
    if stored is not None:
        mol = Chem.MolFromSmiles(stored)
        return message_helpers.canonical_smiles(mol) if mol is not None else None
    try:
        return message_helpers.smiles_from_compound(compound)
    except ValueError:
        # Nothing structural to work from, or an identifier RDKit cannot parse.
        return None


def _component_smiles(reaction: reaction_pb2.Reaction) -> tuple[list[str], list[str]]:
    """Returns ``(input_smiles, output_smiles)``, skipping components without structure.

    Inputs are visited in sorted key order so the column does not depend on map
    iteration order; this matches how reaction SMILES are generated.
    """
    inputs = []
    for key in sorted(reaction.inputs):
        for compound in reaction.inputs[key].components:
            smiles = _compound_smiles(compound)
            if smiles:
                inputs.append(smiles)
    outputs = []
    for outcome in reaction.outcomes:
        for product in outcome.products:
            smiles = _compound_smiles(product)
            if smiles:
                outputs.append(smiles)
    return inputs, outputs


def reaction_smiles(reaction: reaction_pb2.Reaction) -> str | None:
    """Returns the reaction's SMILES, preferring CXSMILES over the plain form.

    Generated from the components when the reaction records neither. Generation is
    best-effort, and a reaction that cannot produce a SMILES reads as null rather than
    failing the whole view.
    """
    stored = _stored_value(reaction, _PREFERRED_REACTION_TYPES)
    if stored is not None:
        return stored
    try:
        return (
            message_helpers.get_reaction_smiles(reaction, generate_if_missing=True)
            or None
        )
    except ValueError:
        return None


def reaction_row(reaction_id: str, reaction: reaction_pb2.Reaction) -> dict:
    """Projects one Reaction onto the view columns.

    Args:
        reaction_id: Reaction ID, used as the join key back to the source.
        reaction: Reaction message to project.

    Returns:
        A dict keyed by ``SCHEMA.names``, with None wherever the source is silent.
    """
    yield_percent, reaction_time_seconds = _outcome_values(reaction)
    inputs, outputs = _component_smiles(reaction)
    return {
        "reaction_id": reaction_id,
        "reaction_smiles": reaction_smiles(reaction),
        "input_smiles": inputs,
        "output_smiles": outputs,
        "yield_percent": yield_percent,
        "temperature_kelvin": _canonical_value(
            reaction.conditions.temperature.setpoint, "K"
        ),
        "reaction_time_seconds": reaction_time_seconds,
        "doi": reaction.provenance.doi or None,
        "patent": reaction.provenance.patent or None,
    }


@dataclasses.dataclass(frozen=True)
class ViewStamps:
    """Footer stamps recording what a view was derived from, and by what."""

    source_dataset_id: str | None
    source_md5: str
    ord_schema_version: str
    view_version: str


def _stamps(source_dataset_id: str | None, source_md5: str) -> ViewStamps:
    return ViewStamps(
        source_dataset_id=source_dataset_id,
        source_md5=source_md5,
        ord_schema_version=metadata.version("ord-schema"),
        view_version=VIEW_VERSION,
    )


def load_stamps(path: str | os.PathLike[str]) -> ViewStamps:
    """Reads a view's footer stamps.

    Args:
        path: Path to a view Parquet file.

    Returns:
        The stamps recorded when the view was written.

    Raises:
        ValueError: If ``path`` is missing a required stamp, e.g. because it is a source
            dataset rather than a view.
    """
    with pq.ParquetFile(path) as parquet_file:
        raw = parquet_file.schema_arrow.metadata or {}

    def _get(key: str) -> str | None:
        value = raw.get(key.encode())
        return value.decode() if value is not None else None

    source_md5 = _get(_META_SOURCE_MD5)
    ord_schema_version = _get(_META_ORD_SCHEMA_VERSION)
    view_version = _get(_META_VIEW_VERSION)
    missing = [
        key
        for key, value in (
            (_META_SOURCE_MD5, source_md5),
            (_META_ORD_SCHEMA_VERSION, ord_schema_version),
            (_META_VIEW_VERSION, view_version),
        )
        if value is None
    ]
    if missing:
        raise ValueError(f"not a view; missing footer keys: {sorted(missing)}")
    assert source_md5 is not None  # Type hint.
    assert ord_schema_version is not None  # Type hint.
    assert view_version is not None  # Type hint.
    return ViewStamps(
        source_dataset_id=_get(_META_SOURCE_DATASET_ID),
        source_md5=source_md5,
        ord_schema_version=ord_schema_version,
        view_version=view_version,
    )


def is_current(view_path: str | os.PathLike[str], source_md5: str) -> bool:
    """Returns whether a view was derived from ``source_md5`` by this library version.

    A view is current only if all three of the source content, the library, and the view
    definition match; any of them changing can change a column's value. A missing or
    unreadable view is not current.
    """
    try:
        stamps = load_stamps(view_path)
    except (OSError, ValueError):
        return False
    return (
        stamps.source_md5 == source_md5
        and stamps.ord_schema_version == metadata.version("ord-schema")
        and stamps.view_version == VIEW_VERSION
    )


def write_view(
    source: str | os.PathLike[str],
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
        source: Path to a source Parquet dataset.
        output: Path to write the view to.
        source_md5: Source hash to stamp, if the caller already computed one (e.g. to
            check whether the view was current). Hashed here when omitted.
        compression: Parquet compression codec.
        row_group_size: Rows per output row group.

    Returns:
        Number of rows written.
    """
    footer = parquet.load_footer(source)
    if source_md5 is None:
        source_md5, _ = parquet.streaming_md5(source)
    stamps = _stamps(footer.dataset.dataset_id or None, source_md5)
    metadata_kv = {
        _META_VIEW_VERSION: stamps.view_version,
        _META_SOURCE_MD5: stamps.source_md5,
        _META_ORD_SCHEMA_VERSION: stamps.ord_schema_version,
    }
    if stamps.source_dataset_id is not None:
        metadata_kv[_META_SOURCE_DATASET_ID] = stamps.source_dataset_id
    schema = SCHEMA.with_metadata(metadata_kv)
    rows = 0
    with (
        atomic_io.atomic_path(output) as temp_path,
        pq.ParquetWriter(temp_path, schema, compression=compression) as writer,
    ):
        for row_group in range(footer.num_row_groups):
            batch: dict[str, list] = {name: [] for name in SCHEMA.names}
            for reaction_id, reaction in parquet.iter_reactions(
                source, row_group=row_group
            ):
                for name, value in reaction_row(reaction_id, reaction).items():
                    batch[name].append(value)
                rows += 1
            writer.write_table(
                pa.Table.from_pydict(batch, schema=schema),
                row_group_size=row_group_size,
            )
    return rows
