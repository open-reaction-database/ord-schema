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
"""Library for executing PostgreSQL queries on the ORD."""

import binascii
import enum
import json

from absl import logging
import psycopg2
from psycopg2 import sql

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

# pylint: disable=too-many-locals


class Predicate:
    """Structure and code generation for ORD query predicates."""

    class MatchMode(enum.Enum):
        """Interpretations for SMILES and SMARTS strings."""
        EXACT = 1
        SIMILAR = 2
        SUBSTRUCTURE = 3
        SMARTS = 4

        @classmethod
        def from_name(cls, name):
            """Take a matching criterion from a URL param."""
            return Predicate.MatchMode[name.upper()]

    def __init__(self):
        # List of pairs, (SMILES/SMARTS, MatchMode).
        self.inputs = []
        # A pair, (SMILES/SMARTS, MatchMode).
        self.output = None
        # A "reaction" SMILES (RDKit fingerprint string).
        self.reaction_smiles = None
        # PK strings, "ord-<hex>".
        self.reaction_ids = []
        # Threshold for similarity matching.
        self.tanimoto_threshold = 0.4
        # Whether to consider stereochemistry.
        self.do_chiral_ss = False

    def add_input(self, pattern, match_mode):
        self.inputs.append((pattern, match_mode))

    def set_output(self, pattern, match_mode):
        self.output = (pattern, match_mode)

    def add_reaction_id(self, reaction_id):
        self.reaction_ids.append(reaction_id)

    def set_reaction_smiles(self, reaction_smiles):
        self.reaction_smiles = reaction_smiles

    def set_tanimoto_threshold(self, threshold):
        self.tanimoto_threshold = threshold

    def use_stereochemistry(self):
        self.do_chiral_ss = True

    def _inputs_predicate(self):
        exprs = []
        for pattern, mode in self.inputs:
            if mode == Predicate.MatchMode.EXACT:
                expr = f"inputs.smiles = '{pattern}'"
            elif mode == Predicate.MatchMode.SIMILAR:
                expr = f"rdki.mfp2%morganbv_fp('{pattern}')"
            elif mode == Predicate.MatchMode.SUBSTRUCTURE:
                expr = f"rdki.m@>'{pattern}'"
            elif mode == Predicate.MatchMode.SMARTS:
                expr = f"rdki.m@>'{pattern}'::qmol"
            exprs.append(expr)
        return ' AND '.join(exprs) or 'TRUE'

    def _output_predicate(self):
        exprs = []
        if self.output is not None:
            pattern, mode = self.output
            if mode == Predicate.MatchMode.EXACT:
                expr = f"outputs.smiles = '{pattern}'"
            elif mode == Predicate.MatchMode.SIMILAR:
                expr = f"rdko.mfp2%morganbv_fp('{pattern}')"
            elif mode == Predicate.MatchMode.SUBSTRUCTURE:
                expr = f"rdko.m@>'{pattern}'"
            elif mode == Predicate.MatchMode.SMARTS:
                expr = f"rdko.m@>'{pattern}'::qmol"
            exprs.append(expr)
        return ' AND '.join(exprs) or 'TRUE'

    def _reactions_predicate(self):
        exprs = []
        if len(self.reaction_ids) > 0:
            clauses = []
            for reaction_id in self.reaction_ids:
                clauses.append(f"reactions.reaction_id='{reaction_id}'")
            expr = ' OR '.join(clauses)
            exprs.append(f'({expr})')
        if self.reaction_smiles is not None:
            # TODO(kearnes): Correctly express RDKit fingerprint queries.
            expr = f"reactions.reaction_smiles='{self.reaction_smiles}'"
            exprs.append(expr)
        return ' AND '.join(exprs) or 'TRUE'

    def query(self):
        return f"""
            SET rdkit.do_chiral_ss={self.do_chiral_ss};
            SET rdkit.tanimoto_threshold={self.tanimoto_threshold};

            SELECT DISTINCT reactions.serialized
            FROM reactions
            INNER JOIN inputs ON inputs.reaction_id = reactions.reaction_id
            INNER JOIN outputs ON outputs.reaction_id = reactions.reaction_id
            INNER JOIN rdk.inputs rdki ON rdki.reaction_id = reactions.reaction_id
            INNER JOIN rdk.outputs rdko ON rdko.reaction_id = reactions.reaction_id
            WHERE
                {self._inputs_predicate()}
                AND {self._output_predicate()}
                AND {self._reactions_predicate()}
            LIMIT 100;
        """

    def query_ids(self):
        return f"""
            SET rdkit.do_chiral_ss={self.do_chiral_ss};
            SET rdkit.tanimoto_threshold={self.tanimoto_threshold};

            SELECT DISTINCT reactions.reaction_id
            FROM reactions
            INNER JOIN inputs ON inputs.reaction_id = reactions.reaction_id
            INNER JOIN outputs ON outputs.reaction_id = reactions.reaction_id
            INNER JOIN rdk.inputs rdki ON rdki.reaction_id = reactions.reaction_id
            INNER JOIN rdk.outputs rdko ON rdko.reaction_id = reactions.reaction_id
            WHERE
                {self._inputs_predicate()}
                AND {self._output_predicate()}
                AND {self._reactions_predicate()}
            LIMIT 100;
        """

    def json(self):
        """Serialize this Predicate for transmission to the web client."""
        predicate = {}
        if self.inputs:
            inputs = []
            for inpt in self.inputs:
                inputs.append({
                    'smiles': inpt[0],
                    'matchMode': inpt[1].name.lower()
                })
            predicate['inputs'] = inputs
        if self.output:
            predicate['output'] = {
                'smiles': self.output[0],
                'matchMode': self.output[1].name.lower()
            }
        predicate['similarity'] = self.tanimoto_threshold
        if len(self.reaction_ids) > 0:
            predicate['reactionIds'] = self.reaction_ids
        if self.reaction_smiles:
            predicate['reactionSmiles'] = self.reaction_smiles
        return json.dumps(predicate)


class OrdPostgres:
    """Class for performing SQL queries on the ORD."""

    def __init__(self, dbname, user, password, host, port):
        """Initializes an instance of OrdPostgres.

        Args:
            dbname: Text database name.
            user: Text user name.
            password: Text user password.
            host: Text host name.
            port: Integer port.
        """
        self._connection = psycopg2.connect(dbname=dbname,
                                            user=user,
                                            password=password,
                                            host=host,
                                            port=port)
        self._connection.set_session(readonly=True)

    def cursor(self):
        return self._connection.cursor()

    def predicate_search(self, predicate):
        """Return the Dataset matching the given Reaction constraints.

        Args:
            predicate: A Predicate defining patterns to match.

        Returns:
            A Dataset proto containing the matched Reactions.
        """
        query = predicate.query()
        print(query)
        logging.info(query)
        reactions = []
        with self._connection, self.cursor() as cursor:
            cursor.execute(query)
            for row in cursor:
                serialized = row[0]
                reaction = reaction_pb2.Reaction.FromString(
                    binascii.unhexlify(serialized.tobytes()))
                reactions.append(reaction)
            self._connection.rollback()  # Revert rdkit runtime configuration.
        return dataset_pb2.Dataset(reactions=reactions)

    def predicate_search_ids(self, predicate):
        """Return the Reaction IDs matching the given Reaction constraints.

        Args:
            predicate: A Predicate defining patterns to match.

        Returns:
            A list of string IDs like "ord-<hex>".
        """
        query = predicate.query_ids()
        print(query)
        logging.info(query)
        with self._connection, self.cursor() as cursor:
            cursor.execute(query)
            reaction_ids = [row[0] for row in cursor]
            self._connection.rollback()  # Revert rdkit runtime configuration.
        return reaction_ids

    def substructure_search(self,
                            pattern,
                            table,
                            limit=100,
                            use_smarts=False,
                            use_stereochemistry=False):
        """Performs a substructure search.

        Args:
            pattern: Text substructure query (SMILES or SMARTS).
            table: Text SQL table name.
            limit: Integer maximum number of matches to return.
            use_smarts: Boolean whether `pattern` is SMARTS.
            use_stereochemistry: Boolean whether to consider stereochemistry.

        Returns:
            Dataset proto containing matched Reactions.
        """
        if 'reactions' in table:
            column = 'r'
        else:
            column = 'm'
        # NOTE(kearnes): We use table.split('.') because "rdk"."inputs" is a
        # valid identifier but "rdk.inputs" is not. See
        # https://www.psycopg.org/docs/sql.html#psycopg2.sql.Identifier.
        components = [
            sql.SQL("""
                SELECT DISTINCT A.reaction_id, A.serialized 
                FROM reactions A 
                INNER JOIN {} B ON A.reaction_id = B.reaction_id
                WHERE B.{}@>%s""").format(sql.Identifier(*table.split('.')),
                                          sql.Identifier(column))
        ]
        args = [pattern]
        if use_smarts:
            components.append(sql.SQL('::qmol'))
        if limit:
            components.append(sql.SQL(' LIMIT %s'))
            args.append(limit)
        command = sql.Composed(components).join('')
        with self._connection, self.cursor() as cursor:
            logging.info(command.as_string(cursor))
            do_chiral_sss = 'true' if use_stereochemistry else 'false'
            # NOTE(kearnes): We call rollback() to reset this change before
            # exiting the context manager (which triggers a commit).
            cursor.execute(sql.SQL('SET rdkit.do_chiral_sss=%s'),
                           [do_chiral_sss])
            cursor.execute(command, args)
            reactions = []
            for result in cursor:
                _, serialized = result
                reaction = reaction_pb2.Reaction.FromString(
                    binascii.unhexlify(serialized.tobytes()))
                reactions.append(reaction)
            self._connection.rollback()
        return dataset_pb2.Dataset(reactions=reactions)

    def similarity_search(self, smiles, table, limit=100, threshold=0.5):
        """Performs a Tanimoto similarity search.

        Args:
            smiles: Text SMILES.
            table: Text SQL table name.
            limit: Integer maximum number of matches to return.
            threshold: Float similarity threshold.

        Returns:
            Dataset proto containing matched Reactions.
        """
        if 'reactions' in table:
            column = 'rdfp'
            function = 'reaction_difference_fp'
        else:
            column = 'mfp2'
            function = 'morganbv_fp'
        components = [
            sql.SQL("""
                SELECT DISTINCT A.reaction_id, A.serialized 
                FROM reactions A 
                INNER JOIN {} B ON A.reaction_id = B.reaction_id
                WHERE B.{}%%{}(%s)""").format(sql.Identifier(*table.split('.')),
                                              sql.Identifier(column),
                                              sql.Identifier(function))
        ]
        args = [smiles]
        if limit:
            components.append(sql.SQL(' LIMIT %s'))
            args.append(limit)
        command = sql.Composed(components).join('')
        with self._connection, self.cursor() as cursor:
            logging.info(command.as_string(cursor))
            # NOTE(kearnes): We call rollback() to reset this change before
            # exiting the context manager (which triggers a commit).
            cursor.execute(sql.SQL('SET rdkit.tanimoto_threshold=%s'),
                           [threshold])
            cursor.execute(command, args)
            reactions = []
            for result in cursor:
                _, serialized = result
                reaction = reaction_pb2.Reaction.FromString(
                    binascii.unhexlify(serialized.tobytes()))
                reactions.append(reaction)
            self._connection.rollback()
        return dataset_pb2.Dataset(reactions=reactions)
