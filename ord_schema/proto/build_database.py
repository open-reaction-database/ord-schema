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

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

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
    # TODO(kearnes): Generate reaction SMILES if missing.
    for identifier in reaction.identifiers:
        if identifier.type == identifier.REACTION_SMILES:
            values['reaction_smiles'] = identifier.value
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
            for identifier in compound.identifiers:
                # TODO(kearnes): Generate SMILES if missing.
                if identifier.type == identifier.SMILES:
                    rows.append({
                        'reaction_id': reaction.reaction_id,
                        'smiles': identifier.value
                    })
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
            values = {}
            for identifier in product.compound.identifiers:
                # TODO(kearnes): Generate SMILES if missing.
                if identifier.type == reaction_pb2.CompoundIdentifier.SMILES:
                    values['reaction_id'] = reaction.reaction_id
                    values['smiles'] = identifier.value
            if product.HasField('compound_yield'):
                values['yield'] = product.compound_yield.value
            rows.append(values)
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
        for table in _TABLES:
            cursor.execute(f'DROP TABLE IF EXISTS {table};')
    cursor.execute('CREATE EXTENSION IF NOT EXISTS rdkit;')
    cursor.execute('CREATE SCHEMA rdk;')
    for table, columns in _TABLES.items():
        dtypes = ',\n'.join(
            [f'\t{column}\t{dtype}' for column, dtype in columns.items()])
        command = f'CREATE TABLE {table} (\n{dtypes}\n);'
        logging.info('Running:\n%s', command)
        cursor.execute(command)
        logging.info('Running COPY')
        with open(os.path.join(FLAGS.output, f'{table}.csv')) as f:
            cursor.copy_expert(
                f'COPY {table} FROM STDIN DELIMITER \',\' CSV HEADER;', f)
        logging.info('Adding RDKit cartridge functionality')
        if 'reaction_smiles' in columns:
            _rdkit_reaction_smiles(cursor, table)
        if 'smiles' in columns:
            _rdkit_smiles(cursor, table)
    db.commit()
    cursor.close()
    db.close()


def _rdkit_reaction_smiles(cursor, table):
    """Adds RDKit cartridge tables for reaction SMILES."""
    cursor.execute(f"""
        SELECT reaction_id,
               r,
               reaction_difference_fp(r) AS rdfp 
        INTO rdk.{table} FROM (
            SELECT reaction_id, 
                   reaction_from_smiles(reaction_smiles::cstring) AS r
            FROM {table}) tmp
        WHERE r IS NOT NULL;""")
    cursor.execute(
        f'CREATE INDEX {table}_r ON rdk.{table} USING gist(r);')
    cursor.execute(
        f'CREATE INDEX {table}_rdfp ON rdk.{table} USING gist(rdfp);')


def _rdkit_smiles(cursor, table):
    """Adds RDKit cartridge tables for molecule SMILES."""
    cursor.execute(f"""
        SELECT reaction_id,
               m,
               morganbv_fp(m) AS mfp2 
        INTO rdk.{table} FROM (
            SELECT reaction_id, 
                   mol_from_smiles(smiles::cstring) AS m
            FROM {table}) tmp
        WHERE m IS NOT NULL;""")
    cursor.execute(
        f'CREATE INDEX {table}_m ON rdk.{table} USING gist(m);')
    cursor.execute(
        f'CREATE INDEX {table}_mfp2 ON rdk.{table} USING gist(mfp2);')


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
