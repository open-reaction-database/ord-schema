# Copyright 2022 Open Reaction Database Project Authors
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

"""Best-effort reaction classification for ORM reactions using Rxn-INSIGHT.

Each reaction is labeled from its generated reaction SMILES with a coarse class
(e.g. "C-C Coupling") and the specific named reaction (e.g. "Suzuki coupling with
boronic acids") to support faceted search. Rxn-INSIGHT runs a transformer
(rxnmapper) per reaction, so this is run as a batched post-pass rather than inline
during import, and requires the optional ``reaction-class`` extra (which pulls in
rxn-insight and, transitively, rxnmapper/torch).
"""

import time

from rxn_insight.reaction import Reaction  # ty: ignore[unresolved-import]
from rxnmapper import RXNMapper  # ty: ignore[unresolved-import]
from sqlalchemy import text
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Rxn-INSIGHT returns this sentinel name when no named reaction matches; store NULL
# instead so the column means "a recognized named reaction" without a magic string.
# The class has no equivalent sentinel: its catch-all "Miscellaneous" is a real member
# of the (closed) class vocabulary and is kept as-is, so a NULL reaction_class means
# "not classified" (error or post-pass not run) rather than "no specific class".
_UNNAMED = "OtherReaction"


def classify_reaction_smiles(
    reaction_smiles: str, rxn_mapper: RXNMapper
) -> tuple[str | None, str | None]:
    """Returns the Rxn-INSIGHT (class, name) for a reaction SMILES.

    Args:
        reaction_smiles: Reaction SMILES; atom mapping is added internally if absent.
        rxn_mapper: Shared RXNMapper so the transformer model loads once per batch.

    Returns:
        A ``(reaction_class, reaction_name)`` tuple. Either element is None when
        Rxn-INSIGHT cannot determine it (or the whole reaction cannot be classified).
    """
    try:
        info = Reaction(reaction_smiles, rxn_mapper=rxn_mapper).get_reaction_info()
    except Exception as error:  # noqa: BLE001
        # Rxn-INSIGHT/rxnmapper raise a range of errors (RDKit parse failures,
        # KeyError, ValueError, ...) on unusual reactions; classification is
        # best-effort, so any failure leaves both columns NULL.
        logger.debug(f"Could not classify {reaction_smiles!r}: {error}")
        return None, None
    reaction_class = info.get("CLASS") or None
    reaction_name = info.get("NAME") or None
    if reaction_name == _UNNAMED:
        reaction_name = None
    return reaction_class, reaction_name


def update_reaction_classes(dataset_id: str, session: Session) -> None:
    """Populates reaction_class and reaction_name for a dataset's unclassified reactions.

    Classifies each distinct reaction SMILES once, then writes the labels back to all
    matching rows that have not been classified yet (reaction_class IS NULL).

    Reactions Rxn-INSIGHT cannot classify keep a NULL reaction_class and are reconsidered
    on later runs. This is intentional -- NULL means "not classified" so coverage stays
    measurable -- and cheap: classification is deterministic, so a repeated pass only re-runs
    the small, fixed set of genuinely unclassifiable reactions rather than the whole dataset.
    """
    logger.debug(f"Updating reaction classes for {dataset_id=}")
    start = time.time()
    result = session.execute(
        text("""
            SELECT DISTINCT reaction_smiles
            FROM ord.reaction
            JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
            WHERE ord.dataset.dataset_id = :dataset_id
              AND ord.reaction.reaction_class IS NULL
              AND ord.reaction.reaction_smiles IS NOT NULL
            """),
        {"dataset_id": dataset_id},
    )
    reaction_smiles_list = [row[0] for row in result]
    if not reaction_smiles_list:
        return
    rxn_mapper = RXNMapper()  # Loads the transformer model once for the whole batch.
    updates = []
    for reaction_smiles in reaction_smiles_list:
        reaction_class, reaction_name = classify_reaction_smiles(
            reaction_smiles, rxn_mapper
        )
        if reaction_class is not None:
            updates.append(
                {
                    "dataset_id": dataset_id,
                    "reaction_smiles": reaction_smiles,
                    "reaction_class": reaction_class,
                    "reaction_name": reaction_name,
                }
            )
    if updates:
        session.execute(
            text("""
                UPDATE ord.reaction
                SET reaction_class = :reaction_class,
                    reaction_name = :reaction_name
                FROM ord.dataset
                WHERE ord.reaction.dataset_id = ord.dataset.id
                  AND ord.dataset.dataset_id = :dataset_id
                  AND ord.reaction.reaction_smiles = :reaction_smiles
                  AND ord.reaction.reaction_class IS NULL
                """),
            updates,
        )
    logger.debug(
        f"Updating reaction classes took {time.time() - start:g}s "
        f"({len(updates)}/{len(reaction_smiles_list)} classified)"
    )
