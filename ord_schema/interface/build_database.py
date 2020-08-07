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
"""Populates a PostgreSQL database containing the ORD.

For simplicity, and to aid in debugging, we write tables to individual CSV
files and load them into PostgreSQL with the COPY command.
"""

import collections
import csv
import glob
import os

from absl import app
from absl import flags
from absl import logging
import psycopg2
from psycopg2 import sql

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern (glob).')
flags.DEFINE_string('output', 'tables', 'Output directory for CSV tables.')
flags.DEFINE_string('database', None, 'Database name.')
flags.DEFINE_string('user', None, 'PostgreSQL user.')
flags.DEFINE_string('password', None, 'PostgreSQL password.')
flags.DEFINE_string('host', None, 'PostgreSQL server host.')
flags.DEFINE_integer('port', None, 'PostgreSQL server port.')
flags.DEFINE_boolean('overwrite', False, 'If True, overwrite existing tables.')

# Postgres table schema.
_TABLES = {
    'reactions':
        collections.OrderedDict([
            ('reaction_id', 'text'),
            ('reaction_smiles', 'text'),
            ('serialized', 'bytea'),
        ]),
    'inputs':
        collections.OrderedDict([
            ('reaction_id', 'text'),
            ('smiles', 'text'),
        ]),
    'outputs':
        collections.OrderedDict([
            ('reaction_id', 'text'),
            ('smiles', 'text'),
            ('yield', 'double precision'),
        ]),
}


class Tables:
    """Holds file handles and CSV writers for database table files."""

    def __init__(self):
        self._handles = []

    def __enter__(self):
        os.makedirs(FLAGS.output, exist_ok=True)
        self._handles = []
        for table, columns in _TABLES.items():
            handle = open(os.path.join(FLAGS.output, f'{table}.csv'), 'w')
            self._handles.append(handle)
            writer = csv.DictWriter(handle,
                                    fieldnames=columns.keys(),
                                    dialect='unix')
            writer.writeheader()
            setattr(self, table, writer)
        return self

    def __exit__(self, *exc):
        for handle in self._handles:
            handle.close()


def process_reaction(reaction, tables):
    """Adds rows to database tables.

    Args:
        reaction: Reaction proto.
        tables: Tables instance.
    """
    _reactions_table(reaction, tables)
    _inputs_table(reaction, tables)
    _outputs_table(reaction, tables)


def _reactions_table(reaction, tables):
    """Adds rows to the 'reactions' table.

    Args:
        reaction: Reaction proto.
        tables: Tables instance.
    """
    values = {
        'reaction_id': reaction.reaction_id,
        'serialized': reaction.SerializeToString().hex()
    }
    try:
        values['reaction_smiles'] = message_helpers.get_reaction_smiles(
            reaction)
    except ValueError:
        pass
    tables.reactions.writerow(values)


def _inputs_table(reaction, tables):
    """Adds rows to the 'inputs' table.

    Args:
        reaction: Reaction proto.
        tables: Tables instance.
    """
    rows = []
    for key in sorted(reaction.inputs):
        reaction_input = reaction.inputs[key]
        for compound in reaction_input.components:
            row = {'reaction_id': reaction.reaction_id}
            try:
                row['smiles'] = message_helpers.smiles_from_compound(compound)
            except ValueError:
                pass
            if len(row) > 1:
                rows.append(row)
    tables.inputs.writerows(rows)


def _outputs_table(reaction, tables):
    """Adds rows to the 'outputs' table.

    Args:
        reaction: Reaction proto.
        tables: Tables instance.
    """
    rows = []
    for outcome in reaction.outcomes:
        for product in outcome.products:
            row = {'reaction_id': reaction.reaction_id}
            try:
                row['smiles'] = message_helpers.smiles_from_compound(
                    product.compound)
            except ValueError:
                pass
            if product.HasField('compound_yield'):
                row['yield'] = product.compound_yield.value
            if len(row) > 1:
                rows.append(row)
    tables.outputs.writerows(rows)


def create_database():
    """Populates the Postgres database."""
    db = psycopg2.connect(dbname=FLAGS.database,
                          user=FLAGS.user,
                          password=FLAGS.password,
                          host=FLAGS.host,
                          port=FLAGS.port)
    cursor = db.cursor()
    if FLAGS.overwrite:
        logging.info('Removing existing tables')
        for table in _TABLES:
            command = sql.SQL('DROP TABLE IF EXISTS {}')
            cursor.execute(command.format(sql.Identifier(table)))
    cursor.execute(sql.SQL('CREATE EXTENSION IF NOT EXISTS rdkit'))
    cursor.execute(sql.SQL('CREATE SCHEMA rdk'))
    for table, columns in _TABLES.items():
        dtypes = []
        for column, dtype in columns.items():
            if table == 'reactions' and column == 'reaction_id':
                component = sql.SQL('{} {} PRIMARY KEY')
            else:
                component = sql.SQL('{} {}')
            # NOTE(kearnes): sql.Identifier(dtype) does not work for the
            # 'double precision' type.
            dtypes.append(
                component.format(sql.Identifier(column), sql.SQL(dtype)))
        command = sql.Composed([
            sql.SQL('CREATE TABLE {} (').format(sql.Identifier(table)),
            sql.Composed(dtypes).join(', '),
            sql.SQL(')')
        ])
        logging.info('Running:\n%s', command.as_string(cursor))
        cursor.execute(command)
        logging.info('Running COPY')
        with open(os.path.join(FLAGS.output, f'{table}.csv')) as f:
            cursor.copy_expert(
                sql.SQL('COPY {} FROM STDIN WITH CSV HEADER').format(
                    sql.Identifier(table)), f)
        logging.info('Adding RDKit cartridge functionality')
        if 'reaction_smiles' in columns:
            _rdkit_reaction_smiles(cursor, table)
        elif 'smiles' in columns:
            _rdkit_smiles(cursor, table)
    db.commit()
    cursor.close()
    db.close()


def _rdkit_reaction_smiles(cursor, table):
    """Adds RDKit cartridge tables for reaction SMILES.

    Creates a new table rdk.<table> with the following columns:
        * r: reaction objects loaded from reaction SMILES
        * rdfp: reaction difference fingerprints

    An index is also created for each column.

    Args:
        cursor: psycopg2 Cursor.
        table: Text table name.
    """
    cursor.execute(
        sql.SQL("""
        SELECT reaction_id,
               r,
               reaction_difference_fp(r) AS rdfp 
        INTO rdk.{} FROM (
            SELECT reaction_id, 
                   reaction_from_smiles(reaction_smiles::cstring) AS r
            FROM {}) tmp
        WHERE r IS NOT NULL""").format(sql.Identifier(table),
                                       sql.Identifier(table)))
    cursor.execute(
        sql.SQL('CREATE INDEX {} ON rdk.{} USING gist(r)').format(
            sql.Identifier(f'{table}_r'), sql.Identifier(table)))
    cursor.execute(
        sql.SQL('CREATE INDEX {} ON rdk.{} USING gist(rdfp)').format(
            sql.Identifier(f'{table}_rdfp'), sql.Identifier(table)))


def _rdkit_smiles(cursor, table):
    """Adds RDKit cartridge tables for molecule SMILES.

    Creates a new table rdk.<table> with the following columns:
        * m: mol objects loaded from SMILES
        * mfp2: Morgan fingerprints

    An index is also created for each column.

    Args:
        cursor: psycopg2 Cursor.
        table: Text table name.
    """
    cursor.execute(
        sql.SQL("""
        SELECT reaction_id,
               m,
               morganbv_fp(m) AS mfp2 
        INTO rdk.{} FROM (
            SELECT reaction_id, 
                   mol_from_smiles(smiles::cstring) AS m
            FROM {}) tmp
        WHERE m IS NOT NULL""").format(sql.Identifier(table),
                                       sql.Identifier(table)))
    cursor.execute(
        sql.SQL('CREATE INDEX {} ON rdk.{} USING gist(m)').format(
            sql.Identifier(f'{table}_m'), sql.Identifier(table)))
    cursor.execute(
        sql.SQL('CREATE INDEX {} ON rdk.{} USING gist(mfp2)').format(
            sql.Identifier(f'{table}_mfp2'), sql.Identifier(table)))


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input)
    logging.info('Found %d datasets', len(filenames))
    with Tables() as tables:
        for filename in filenames:
            logging.info(filename)
            dataset = message_helpers.load_message(filename,
                                                   dataset_pb2.Dataset)
            for reaction in dataset.reactions:
                process_reaction(reaction, tables)
    if FLAGS.database:
        logging.info('Creating Postgres database')
        create_database()


if __name__ == '__main__':
    flags.mark_flag_as_required('input')
    app.run(main)
