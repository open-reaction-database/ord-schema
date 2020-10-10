#!/bin/bash
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

# Runs editor javascript tests.

set -e -x

export ORD_EDITOR_MOUNT=/tmp/editor-postgres
export PGPASSWORD=${ORD_EDITOR_POSTGRES_PASSWORD}

# Clear any leftover database state.
rm -rf $ORD_EDITOR_MOUNT && mkdir $ORD_EDITOR_MOUNT

docker build --file=editor/Dockerfile -t openreactiondatabase/ord-editor .
docker-compose -f editor/docker-compose.yml up &

# Wait for the database to become available.
function connect() {
  psql -p 5432 -h localhost -U postgres < /dev/null
}
set +e
connect
while [ $? -ne 0 ]; do
  echo waiting for Postgres
  sleep 5
  connect
done
set -e

# Initialize the database with schema and contents.
psql -p 5432 -h localhost -U postgres -f editor/schema.sql
pushd editor
./py/migrate.py
popd

# Run the JS tests.
node editor/js/test.js

# TODO: Run python tests that need a database.

# Shut down the containers and relay the test process status.
status=$?
docker-compose -f editor/docker-compose.yml down
[ "${status}" -eq 0 ] && echo PASS || echo FAIL
test "${status}" -eq 0
