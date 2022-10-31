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

"""Cleaning protobuf database from ORD."""

import os
from typing import List, Tuple

from ord_schema import message_helpers
from ord_schema.orm.database import (add_dataset,
                                     add_rdkit,
                                     get_connection_string)
#
from ord_schema.orm.mappers import (Percentage,
                                    ProductCompound,
                                    ProductMeasurement,
                                    Reaction,
                                    ReactionOutcome)
# from ord_schema.orm.mappers import *
from ord_schema.proto import reaction_pb2
from rdkit import Chem
# from sqlalchemy.orm import Session
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session


# algorithm
# step 1: convert propto2 into a relational database
# step 2: extract the reactants and products
# step 3: cleaning
#         - standardization of smiles
#         - remove duplication
#         -
# step 4: dump a new proto2 database
# step 5: save JSON or pandas dataframe if needed

# questions
# 1. do we need to add a mol object?
# 2. names for multiple reactants and products?


def proto2database(ord_name: str,
                   username: str,
                   password: str,
                   host: str,
                   port: str,
                   database: str,
                   ) -> None:
    """Convert the protocol buffer file from ORD to a relational database.

    Args:
        ord_name: The ORD file path locally or dataset_id from ord-data repository.
        username: The username for the database.
        password: The password for the database.
        host: The host url from postgres.
        port: The port from postgres.
        database: The database name hosted with postgres.

    Returns:

    """
    if os.path.exists(ord_name):
        dataset = message_helpers.load_message(ord_name)
    else:
        dataset = message_helpers.fetch_dataset(dataset_id=ord_name)

    connection_string = get_connection_string(database, username, password, host, port)
    engine = create_engine(connection_string, future=True)
    with Session(engine) as session:
        add_dataset(dataset, session)
        session.flush()
        add_rdkit(session)  # there is no need for this
        session.commit()


def clean_up_database(connection_string: str,
                      new_database: str = None,
                      fail_database: str = None,
                      # yield_threshold: float = 0.70,
                      ) -> Tuple[List[reaction_pb2.Reaction], List[reaction_pb2.Reaction]]:
    """Clean up the database.

    Returns:
        new_reactions: The list of new reactions.
        failed_records: A list of reactions that failed to be cleaned.

    """
    engine = create_engine(connection_string, future=True)
    with Session(engine) as session:
        query = (
            select(Reaction)
            .join(ReactionOutcome)
            .join(ProductCompound)
            .join(ProductMeasurement)
            .join(Percentage)
            # .where(ProductMeasurement.type == "YIELD", Percentage.value >= yield_threshold)
        )
        results = session.execute(query)
        reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]

    # todo: define the new protpo2 database,
    # this can be a bug if two reaction list with different number of reactions
    # question: what information to keep? OR only copy part of the original database (how)?
    # We can keep the reaction id, the SMILES strings for product and reactants.

    new_reactions = []
    failed_records = []
    # cleaning up the SMILES strings
    for reaction in reactions:
        # question: how to get reactants and products? All or only the major ones?
        # reactant_smi = message_helpers.smiles_from_compound(
        #     reaction.inputs["Boronate in Solvent"].components[0])
        new_reaction = reaction_pb2.Reaction()
        new_reaction.CopyFrom(reaction)

        # component.reaction_role == reaction_pb2.ReactionRole.REACTANT
        # or
        # component.reaction_role == 1
        # according to https://github.com/open-reaction-database/ord-schema/blob/
        # 2146425277dead9c7a3c5f4a5a3450b7980cc92f/ord_schema/proto/reaction.proto#L275
        reactant_smiles = [message_helpers.smiles_from_compound(component) for component in
                           reaction.inputs["Boronate in Solvent"].components
                           if component.reaction_role == reaction_pb2.ReactionRole.REACTANT]
        product_smiles = message_helpers.smiles_from_compound(reaction.outcomes[0].products[0])

    # todo: use reaction SMILES to remove duplicates
        # update the reactant SMILES
        for idx, reactant_smi in enumerate(reactant_smiles):
            # todo: check on how to deal with metal complex for catalysis, example in notebook
            mol = Chem.MolFromSmiles(reactant_smi, sanitize=True)
            if mol:
                new_reaction.inputs["Boronate in Solvent"].components[idx].identifiers[0].value = \
                    Chem.MolToSmiles(mol, isomericSmiles=True)
            else:
                print(f"SMILES string {reactant_smi} is not valid.")
                failed_records.append(reaction)
                continue
        # update the product SMILES
        mol = Chem.MolFromSmiles(product_smiles, sanitize=True)
        if mol:
            new_reaction.outcomes[0].products[0].identifiers[0].value = \
                Chem.MolToSmiles(mol, isomericSmiles=True)

        new_reactions.append(new_reaction)

    if new_database:
        message_helpers.write_message(new_reactions, new_database)
    if fail_database:
        message_helpers.write_message(failed_records, fail_database)

    return new_reactions, failed_records


# https://github.com/open-reaction-database/ord-schema/blob/main/examples/submissions
# /1_Santanilla_Nanomole_HTE/example_Santanilla.ipynb
# canonical_smiles = []
# for smiles in smiles_data["non_canonical_smiles"]:
#     mol = Chem.MolFromSmiles(smiles, sanitize=False)
#     if has_transition_metal(mol):
#         mol = set_dative_bonds(mol)
#
#     else:
#         mol = Chem.MolFromSmiles(smiles)
#
#     if not mol == None:
#         can_smiles = Chem.MolToSmiles(mol)
#         canonical_smiles.append(can_smiles)
#     else:
#         raise Exception(
#             f"Error, unable to canonicalize the following smiles: {smiles}, non-canonical smiles added to final list instead"
#         )
#         canonical_smiles.append(smiles)
