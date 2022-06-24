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

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

_COMPOUND_STRUCTURAL_IDENTIFIERS = [
    reaction_pb2.CompoundIdentifier.SMILES,
    reaction_pb2.CompoundIdentifier.INCHI,
    reaction_pb2.CompoundIdentifier.MOLBLOCK,
    reaction_pb2.CompoundIdentifier.CXSMILES,
    reaction_pb2.CompoundIdentifier.XYZ,
]

_USERNAME = "github-actions"
_EMAIL = "github-actions@github.com"


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
    modified = False
    id_substitutions = {}
    if not re.fullmatch("^ord-[0-9a-f]{32}$", reaction.reaction_id):
        # NOTE(kearnes): This does not check for the case where a Dataset is
        # edited and reaction_id values are changed inappropriately. This will
        # need to be either (1) caught in review or (2) found by a complex
        # check of the diff.
        reaction_id = f"ord-{uuid.uuid4().hex}"
        if reaction.reaction_id:
            id_substitutions[reaction.reaction_id] = reaction_id
        reaction.reaction_id = reaction_id
        modified = True
    for func in _UPDATES:
        # NOTE(kearnes): Order is important here; if you write
        # `modified or func(reaction)` and modified is True, the interpreter
        # will skip the evaluation of func(reaction).
        modified = func(reaction) or modified
    if modified:
        event = reaction.provenance.record_modified.add()
        event.time.value = datetime.datetime.utcnow().ctime()
        event.person.username = _USERNAME
        event.person.email = _EMAIL
        event.details = "Automatic updates from the submission pipeline."
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
    if not re.fullmatch("^ord_dataset-[0-9a-f]{32}$", dataset.dataset_id):
        dataset.dataset_id = f"ord_dataset-{uuid.uuid4().hex}"
    # Reaction-level updates
    id_substitutions = {}
    for reaction in dataset.reactions:
        id_substitutions.update(update_reaction(reaction))
    # Dataset-level updates of cross-references
    for reaction in dataset.reactions:
        for reaction_input in reaction.inputs.values():
            for component in reaction_input.components:
                for preparation in component.preparations:
                    if preparation.reaction_id and preparation.reaction_id in id_substitutions:
                        preparation.reaction_id = id_substitutions[preparation.reaction_id]
            for crude_component in reaction_input.crude_components:
                if crude_component.reaction_id in id_substitutions:
                    crude_component.reaction_id = id_substitutions[crude_component.reaction_id]


# Standard updates.
_UPDATES = []
