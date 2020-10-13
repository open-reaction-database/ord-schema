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

# Runs editor tests in python and JS.

set -e

# The Postgres Docker container gets its password like this.
[ "$ORD_EDITOR_POSTGRES_PASSWORD" != "" ] || \
    ( echo "*** missing ORD_EDITOR_POSTGRES_PASSWORD ***" && false )

# The psql command in this script gets is password like this.
export PGPASSWORD="$ORD_EDITOR_POSTGRES_PASSWORD"

# The Flask app in absltest gets its Postgres password like this.
export POSTGRES_PASSWORD="$ORD_EDITOR_POSTGRES_PASSWORD"

# Postgres puts its files here, via docker-compose.yml.
export ORD_EDITOR_MOUNT=/tmp/editor-postgres

# Clear any leftover database state.
rm -rf $ORD_EDITOR_MOUNT && mkdir $ORD_EDITOR_MOUNT

docker build --file=editor/Dockerfile -t openreactiondatabase/ord-editor .
docker-compose --file editor/docker-compose.yml up --detach

set +e

# Wait for the database to become available.
function connect() {
  psql -p 5432 -h localhost -U postgres < /dev/null
}
connect
while [ $? -ne 0 ]; do
  echo waiting for Postgres
  sleep 5
  connect
done

# Initialize the database with schema and contents.
set -e
psql -p 5432 -h localhost -U postgres -f editor/schema.sql
pushd editor
./py/migrate.py
popd
set +e

status=0

# Run the JS tests.
node editor/js/test.js
[ $? -eq 0 ] || status=1

# Python tests run Flask in the test environment, not a container.
python editor/py/serve_test.py
[ $? -eq 0 ] || status=1

# Report pass/fail.
red='\033[0;31m'
green='\033[0;32m'
neutral='\033[0m'
[ "${status}" -eq 0 ] && \
    printf "${green}PASS${neutral}\n" || printf "${red}FAIL${neutral}\n"

# Shut down the containers.
docker-compose -f editor/docker-compose.yml down

# Relay the status for GitHub CI.
test "${status}" -eq 0
