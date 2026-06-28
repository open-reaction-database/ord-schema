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

Labels each reaction's generated SMILES with a coarse class (e.g. "C-C Coupling") and a
named reaction (e.g. "Suzuki coupling with boronic acids") for faceted search. Rxn-INSIGHT
runs a transformer per reaction, so this is a batched post-pass behind the optional
``reaction-class`` extra (rxn-insight, plus rxnmapper/torch), not inline in from_proto.
"""

import time

from rxn_insight.reaction import Reaction  # ty: ignore[unresolved-import]
from rxnmapper import RXNMapper  # ty: ignore[unresolved-import]
from sqlalchemy import text
from sqlalchemy.orm import Session

from ord_schema.logging import get_logger

logger = get_logger(__name__)

# Rxn-INSIGHT's "no named reaction" sentinel, normalized to NULL so reaction_name carries
# no magic value. (The class catch-all "Miscellaneous" is a real category and is kept.)
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
        # Best-effort: any Rxn-INSIGHT/rxnmapper failure leaves the reaction unclassified.
        logger.debug(f"Could not classify {reaction_smiles!r}: {error}")
        return None, None
    reaction_class = info.get("CLASS") or None
    reaction_name = info.get("NAME") or None
    if reaction_name == _UNNAMED:
        reaction_name = None
    return reaction_class, reaction_name


def update_reaction_classes(dataset_id: str, session: Session) -> None:
    """Classifies a dataset's not-yet-attempted reactions into derived.reaction_classes.

    Selects reactions that have a SMILES but no row yet, classifies each distinct SMILES
    once (cached), and inserts one row per reaction. The row records the attempt -- NULL
    class/name mean Rxn-INSIGHT could not classify it -- so repeated runs converge rather
    than re-running the model over deterministic failures.
    """
    logger.debug(f"Updating reaction classes for {dataset_id=}")
    start = time.time()
    result = session.execute(
        text("""
            SELECT ord.reaction.id, ord.reaction.reaction_smiles
            FROM ord.reaction
            JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
            WHERE ord.dataset.dataset_id = :dataset_id
              AND ord.reaction.reaction_smiles IS NOT NULL
              AND NOT EXISTS (
                  SELECT 1 FROM derived.reaction_classes
                  WHERE derived.reaction_classes.reaction_id = ord.reaction.id
              )
            """),
        {"dataset_id": dataset_id},
    )
    reactions = result.all()
    if not reactions:
        return
    rxn_mapper = RXNMapper()  # Loads the transformer model once for the whole batch.
    cache: dict[str, tuple[str | None, str | None]] = {}
    inserts = []
    for reaction_id, reaction_smiles in reactions:
        if reaction_smiles not in cache:
            cache[reaction_smiles] = classify_reaction_smiles(
                reaction_smiles, rxn_mapper
            )
        reaction_class, reaction_name = cache[reaction_smiles]
        inserts.append(
            {
                "reaction_id": reaction_id,
                "reaction_class": reaction_class,
                "reaction_name": reaction_name,
            }
        )
    # ON CONFLICT guards a concurrent insert; the NOT EXISTS selector already skips
    # reactions that have a row.
    session.execute(
        text("""
            INSERT INTO derived.reaction_classes (reaction_id, reaction_class, reaction_name)
            VALUES (:reaction_id, :reaction_class, :reaction_name)
            ON CONFLICT (reaction_id) DO NOTHING
            """),
        inserts,
    )
    classified = sum(value[0] is not None for value in cache.values())
    logger.debug(
        f"Updating reaction classes took {time.time() - start:g}s "
        f"({classified}/{len(cache)} distinct SMILES classified, {len(inserts)} reactions)"
    )
