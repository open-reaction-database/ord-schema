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

set -ex

python "${ORD_ROOT}/ord-schema/ord_schema/proto/build_database.py" \
  --input="${ORD_ROOT}/ord-submissions-test/data/*/*.pbtxt" \
  --output="${HOME}/tables" \
  --database="${POSTGRES_DB}" \
  --user="${POSTGRES_USER}" \
  --password="${POSTGRES_PASSWORD}"
