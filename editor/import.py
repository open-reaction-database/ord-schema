#!/usr/bin/env python
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
import psycopg2
import psycopg2.sql
import sys
import time


def main():
    with psycopg2.connect(dbname='editor', port=5430) as db:
        for user in os.listdir('db'):
            if len(user) != 32:
                continue
            query = psycopg2.sql.SQL(
                'INSERT INTO users VALUES (%s, %s, %s) ON CONFLICT DO NOTHING')
            with db.cursor() as cursor:
                timestamp = int(time.time())
                cursor.execute(query, [user, user, timestamp])
            for name in os.listdir(f'db/{user}'):
                if not name.endswith('.pbtxt'):
                    continue
                pbtxt = open(f'db/{user}/{name}').read()
                query = psycopg2.sql.SQL(
                    'INSERT INTO datasets VALUES (%s, %s, %s) '
                    'ON CONFLICT DO NOTHING')
                with db.cursor() as cursor:
                    cursor.execute(query, [user, name[:-6], pbtxt])


if __name__ == '__main__':
  sys.exit(main())
