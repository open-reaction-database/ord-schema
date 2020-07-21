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

from absl import logging
import psycopg2

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


class OrdPostgres:
    """Class for performing SQL queries on the ORD."""
    def __init__(self,
                 dbname='ord',
                 user='postgres',
                 password='postgres',
                 host=None,
                 port=None):
        """Initializes an instance of OrdPostgres.

        Args:
            dbname: Text database name.
            user: Text user name.
            password: Text user password.
            host: Text host name.
            port: Integer port.
        """
        self._dbname = dbname
        self._user = user
        self._password = password
        self._host = host
        self._port = port
        self._connection = None

    def __enter__(self):
        self._connection = psycopg2.connect(dbname=self._dbname,
                                            user=self._user,
                                            password=self._password,
                                            host=self._host,
                                            port=self._port)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        del exc_type, exc_val, exc_tb  # Unused.
        self._connection.close()

    def cursor(self):
        return self._connection.cursor()

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
        do_chiral_sss = 'true' if use_stereochemistry else 'false'
        query = f"""
            SELECT DISTINCT A.reaction_id, A.serialized 
            FROM reactions A 
            INNER JOIN {table} B ON A.reaction_id = B.reaction_id
            WHERE B.{column}@>'{pattern}'{'::qmol' if use_smarts else ''}
            {f'LIMIT {limit}' if limit else ''};"""
        logging.info(query)
        with self._connection, self.cursor() as cursor:
            # NOTE(kearnes): I don't like that we are committing this change to
            # the database after the context manager exits. It's possible to
            # auto-rollback by raising an exception in the context manager.
            cursor.execute(f'SET rdkit.do_chiral_sss={do_chiral_sss};')
            cursor.execute(query)
            reactions = []
            for result in cursor:
                _, serialized = result
                reaction = reaction_pb2.Reaction.FromString(
                    binascii.unhexlify(serialized.tobytes()))
                reactions.append(reaction)
        return dataset_pb2.Dataset(reactions=reactions)
