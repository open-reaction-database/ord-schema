# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Automated updates for Reaction messages."""

import datetime
import re
import uuid
from collections.abc import Iterable

from ord_schema import parquet_dataset
from ord_schema.proto import dataset_pb2, reaction_pb2

_USERNAME = "github-actions"
_EMAIL = "github-actions@github.com"

_REACTION_ID_RE = re.compile("^ord-[0-9a-f]{32}$")
_DATASET_ID_RE = re.compile("^ord_dataset-[0-9a-f]{32}$")


def assign_dataset_id(dataset: dataset_pb2.Dataset | parquet_dataset.DatasetView) -> str:
    """Assigns a canonical ``dataset_id`` if the existing one is missing or non-canonical.

    Mutates ``dataset.dataset_id`` in place. Works for both ``Dataset`` and
    ``DatasetView`` (which exposes ``dataset_id`` as a writable attribute).

    Returns:
        The (possibly newly-assigned) dataset_id.
    """
    if not _DATASET_ID_RE.fullmatch(dataset.dataset_id):
        dataset.dataset_id = f"ord_dataset-{uuid.uuid4().hex}"
    return dataset.dataset_id


def assign_id_substitutions(old_ids: Iterable[str]) -> tuple[list[str | None], dict[str, str]]:
    """Pre-allocates canonical reaction IDs for a sequence of old IDs.

    A reaction's ID is replaced when the existing one is missing or does not
    match the canonical ``ord-{32 hex}`` pattern. Cross-reference rewriting
    only applies to old IDs that were non-empty (i.e., user-supplied
    placeholders); reactions whose old ID was empty get a new ID but no
    substitution entry, since nothing else could have referenced them.

    Args:
        old_ids: Reaction IDs in the order they appear in the dataset.

    Returns:
        new_ids: List parallel to ``old_ids``; entry is the new reaction_id
            to assign, or ``None`` if the old ID was already canonical.
        id_substitutions: Map of ``old_id -> new_id`` for entries where the
            old ID was a non-empty placeholder. Used to rewrite cross-references.
    """
    new_ids: list[str | None] = []
    id_substitutions: dict[str, str] = {}
    for old_id in old_ids:
        if _REACTION_ID_RE.fullmatch(old_id):
            new_ids.append(None)
            continue
        new_id = f"ord-{uuid.uuid4().hex}"
        new_ids.append(new_id)
        if old_id:
            id_substitutions[old_id] = new_id
    return new_ids, id_substitutions


def apply_reaction_updates(reaction: reaction_pb2.Reaction, *, new_id: str | None) -> bool:
    """Applies per-reaction updates in place using a pre-computed reaction ID.

    Splitting ID generation out of this function lets a streaming caller
    allocate IDs in a cheap pre-pass (e.g. from a Parquet ``reaction_id``
    column) and inject them here without re-deriving them.

    Args:
        reaction: Reaction message to mutate.
        new_id: Pre-computed reaction_id to assign, or ``None`` to leave the
            existing ID untouched.

    Returns:
        True if the reaction was modified.
    """
    modified = False
    if new_id is not None:
        reaction.reaction_id = new_id
        modified = True
    for func in _UPDATES:
        # NOTE(kearnes): Order is important here; if you write
        # `modified or func(reaction)` and modified is True, the interpreter
        # will skip the evaluation of func(reaction).
        modified = func(reaction) or modified
    if modified:
        event = reaction.provenance.record_modified.add()
        event.time.value = datetime.datetime.now(datetime.timezone.utc).ctime()
        event.person.username = _USERNAME
        event.person.email = _EMAIL
        event.details = "Automatic updates from the submission pipeline."
    return modified


def apply_cross_reference_substitutions(reaction: reaction_pb2.Reaction, id_substitutions: dict[str, str]) -> None:
    """Rewrites cross-referenced reaction_ids inside ``reaction`` using the substitution map."""
    if not id_substitutions:
        return
    for reaction_input in reaction.inputs.values():
        for component in reaction_input.components:
            for preparation in component.preparations:
                if preparation.reaction_id in id_substitutions:
                    preparation.reaction_id = id_substitutions[preparation.reaction_id]
        for crude_component in reaction_input.crude_components:
            if crude_component.reaction_id in id_substitutions:
                crude_component.reaction_id = id_substitutions[crude_component.reaction_id]


def update_reaction(reaction: reaction_pb2.Reaction) -> dict[str, str]:
    """Updates a Reaction message.

    Current updates:
      * Sets reaction_id if not already set.
      * Adds a record modification event to the provenance.

    Args:
        reaction: reaction_pb2.Reaction message.

    Returns:
        A dictionary mapping placeholder reaction_ids to newly-assigned
            reaction_ids.
    """
    new_ids, id_substitutions = assign_id_substitutions([reaction.reaction_id])
    apply_reaction_updates(reaction, new_id=new_ids[0])
    return id_substitutions


def update_dataset(dataset: dataset_pb2.Dataset):
    """Updates a Dataset message.

    Current updates:
      * All reaction-level updates in update_reaction.
      * reaction_id cross-references between Reactions in the dataset.

    Args:
        dataset: dataset_pb2.Dataset message.

    Raises:
        KeyError: if the dataset has not been validated and there exists a
            cross-referenced reaction_id in any Reaction that is not defined
            elsewhere in the Dataset.
    """
    assign_dataset_id(dataset)
    new_ids, id_substitutions = assign_id_substitutions(r.reaction_id for r in dataset.reactions)
    for reaction, new_id in zip(dataset.reactions, new_ids, strict=True):
        apply_reaction_updates(reaction, new_id=new_id)
    for reaction in dataset.reactions:
        apply_cross_reference_substitutions(reaction, id_substitutions)


def update_dataset_parquet(input_path: str, output_path: str, *, dataset_id: str) -> None:
    """Stream-applies ``update_dataset`` to a Parquet input, writing the result to ``output_path``.

    Two passes over ``input_path``:

    * Pass 1 reads only the ``reaction_id`` column (no Reaction decode) to
      pre-allocate canonical reaction IDs and build the cross-reference map.
    * Pass 2 streams full Reactions, applies per-reaction updates and
      cross-reference rewrites, and writes them via ``DatasetWriter``.

    Peak memory is bounded by one row group plus the ID maps. The caller is
    responsible for choosing ``output_path`` based on the resolved ``dataset_id``
    (call ``assign_dataset_id`` on the input header first to learn it) and for
    any atomic-rename / validation dance — keeping the rename outside lets the
    caller validate the written file before publishing it.

    Args:
        input_path: Path to the input Parquet dataset.
        output_path: Path to write the updated Parquet dataset to.
        dataset_id: Resolved dataset_id to write into the output footer.
    """
    header = parquet_dataset.read_metadata(input_path)
    new_ids, id_substitutions = assign_id_substitutions(parquet_dataset.iter_reaction_ids(input_path))
    with parquet_dataset.DatasetWriter(
        output_path,
        name=header.name,
        description=header.description,
        dataset_id=dataset_id,
    ) as writer:
        reactions = (reaction for _, reaction in parquet_dataset.iter_reactions(input_path))
        for reaction, new_id in zip(reactions, new_ids, strict=True):
            apply_reaction_updates(reaction, new_id=new_id)
            apply_cross_reference_substitutions(reaction, id_substitutions)
            writer.write(reaction)


# Standard updates.
_UPDATES = []
