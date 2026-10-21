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

"""Each reaction's first commit to ord-data, one artifact per dataset.

When a reaction entered ORD's public corpus is not in its record, so it is not in the
projection either, which restates the record and nothing else. It is read from
ord-data's git history into an index (:mod:`ord_schema.artifacts.history`), and this
artifact carries the index's answer for one dataset: one row per reaction its projection
holds, keyed by ``reaction_id``, with the commit, its time as git recorded it, and the
dataset whose file first held the reaction.

It derives from the projection, as the structures artifact does, and passes the
projection's source hash through, so it pairs with its projection the way every other
artifact does. Its currency has a second input the stamps every artifact carries cannot
see: the index. It is stamped with the index's walk version and head, and stays current
while that head is on the chain a later index walked, since a walk only appends. A
reaction the index lacks refuses the derivation rather than reading as one ord-data
never held.

A query reads this artifact only when it names ``first_commit``; see
:mod:`ord_schema.search.query`. Every other query reads the projection alone and costs
what it did before the artifact existed.
"""

import os

import pyarrow as pa
import pyarrow.parquet as pq

from ord_schema import atomic_io
from ord_schema.artifacts import base, history, projection

ARTIFACT = "first_commits"

# The column a query names, as ``first_commit.timestamp``.
COLUMN = "first_commit"

SCHEMA = pa.schema(
    [
        pa.field("reaction_id", pa.string(), nullable=False),
        pa.field(COLUMN, history.TYPE),
    ]
)


def is_current(
    path: str | os.PathLike[str], source_md5: str, index: history.Index
) -> bool:
    """Returns whether ``path`` restates ``source_md5`` under ``index`` by this library.

    Args:
        path: The artifact. Need not exist.
        source_md5: Hash of the source dataset it is expected to restate.
        index: The index of first commits it is expected to agree with.

    Returns:
        Whether its stamps are current and name a head on the index's chain.
    """
    return base.is_current(path, ARTIFACT, source_md5, SCHEMA) and index.covers(path)


def write_first_commits(
    source: str | os.PathLike[str],
    output: str | os.PathLike[str],
    *,
    index: history.Index,
    source_md5: str | None = None,
    source_dataset_id: str | None = None,
    compression: str = "zstd",
    row_group_size: int = 100_000,
) -> int:
    """Derives a first_commits artifact from a projection and writes it.

    Only the projection's ``reaction_id`` column is read. The output is published
    atomically, so a failure leaves any existing artifact untouched.

    Args:
        source: Path to the projection whose reactions to look up.
        output: Path to write the artifact to.
        index: The index of first commits, walked to at least the commit the
            projection's source was read from.
        source_md5: Hash of the *source dataset* to stamp, if the caller already read
            one. Taken from the projection's stamps when omitted.
        source_dataset_id: Source dataset ID to stamp, if the caller already read one.
            Taken from the projection's stamps when omitted.
        compression: Parquet codec, any name ``pq.ParquetWriter`` accepts.
        row_group_size: Rows per output row group. Large, because the artifact has few
            leaves and a footer per row group is what a reader's cache holds.

    Returns:
        Number of rows written: one per reaction in the projection.

    Raises:
        ValueError: If ``row_group_size`` is less than one, refused before the source is
            read; if ``source`` is not a current projection; or if the index lacks one
            of its reactions, named with the source and the index's head.
    """
    if row_group_size < 1:
        # pyarrow refuses zero only after the lookup, and reads below zero as no bound.
        raise ValueError(
            f"row_group_size is {row_group_size}, which no row group can hold"
        )
    parent = base.load_stamps(source)
    if parent.artifact != projection.ARTIFACT:
        raise ValueError(
            f"{source} is a {parent.artifact}, not a {projection.ARTIFACT}; the "
            "first_commits artifact looks up a projection's reactions"
        )
    # derive_tree refuses stale parents, but this writer is public and its output
    # inherits the dataset hash, so it checks again.
    if not base.stamps_are_current(parent, projection.ARTIFACT):
        raise ValueError(
            f"{source} is a stale {projection.ARTIFACT}; derive it again first"
        )
    if source_md5 is None:
        source_md5 = parent.source_md5
    if source_dataset_id is None:
        source_dataset_id = parent.source_dataset_id
    reaction_ids = pq.read_table(source, columns=["reaction_id"]).column(0).to_pylist()
    arrivals = index.lookup(reaction_ids)
    unknown = sorted(set(reaction_ids) - arrivals.keys())
    if unknown:
        # A null here would read as a reaction ord-data never held, when the index is
        # what is short.
        raise ValueError(
            f"{source}: {len(unknown)} reactions are not in the first-commits index "
            f"walked to {index.head}, among them {unknown[:3]}; walk ord-data to the "
            "commit this dataset was read from"
        )
    stamps = base.current_stamps(ARTIFACT, source_dataset_id, source_md5)
    table = pa.table(
        {
            "reaction_id": pa.array(reaction_ids, pa.string()),
            COLUMN: pa.array(
                [arrivals[reaction_id] for reaction_id in reaction_ids], history.TYPE
            ),
        },
        schema=SCHEMA,
    ).replace_schema_metadata({**base.to_metadata(stamps), **index.stamps()})
    with atomic_io.atomic_path(output) as temp_path:
        pq.write_table(
            table, temp_path, compression=compression, row_group_size=row_group_size
        )
    return table.num_rows
