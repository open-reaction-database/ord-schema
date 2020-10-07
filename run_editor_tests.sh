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
set -x

docker build --file=editor/Dockerfile -t openreactiondatabase/ord-editor .
CONTAINER=$(docker run --rm -d -p 5000:5000 openreactiondatabase/ord-editor)
# Give the web server time to come up.
sleep 5
node editor/js/test.js
STATUS=$?
test "${STATUS}" -eq 0 || docker logs "${CONTAINER}"
docker stop "${CONTAINER}"
test "${STATUS}" -eq 0
