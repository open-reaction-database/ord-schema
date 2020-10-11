#!/usr/bin/env python3
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
"""Slurp the local db/ directory contents into Postgres."""

import os
import re
import sys
import time

import psycopg2
import psycopg2.sql


def migrate_one(user_id, name, conn):
    """Slurp one named dataset from the db/ directory into Postgres."""
    pbtxt = open(f'db/{user_id}/{name}').read()
    query = psycopg2.sql.SQL(
        'INSERT INTO datasets VALUES (%s, %s, %s) '
        'ON CONFLICT (user_id, dataset_name) DO UPDATE SET pbtxt=%s')
    with conn.cursor() as cursor:
        cursor.execute(query, [user_id, name[:-6], pbtxt, pbtxt])


def migrate_all():
    """Run as a script, copies the entire contents of the db/ directory."""
    with psycopg2.connect(dbname='editor',
                          host='localhost',
                          port=5432,
                          user='postgres') as conn:
        for user_id in os.listdir('db'):
            if re.match('^[0-9a-fA-F]{32}$', user_id) is None:
                continue
            query = psycopg2.sql.SQL(
                'INSERT INTO users VALUES (%s, %s, %s) ON CONFLICT DO NOTHING')
            with conn.cursor() as cursor:
                timestamp = int(time.time())
                cursor.execute(query, [user_id, user_id, timestamp])
            for name in os.listdir(f'db/{user_id}'):
                if not name.endswith('.pbtxt'):
                    continue
                migrate_one(user_id, name, conn)


if __name__ == '__main__':
    sys.exit(migrate_all())
