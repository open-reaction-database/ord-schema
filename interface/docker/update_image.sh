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

# NOTE(kearnes): Disable caching so the latest version of the ORD is used.
docker build --no-cache -t ord-postgres:empty .
docker run --rm --name ord-postgres -d \
  -e POSTGRES_USER=postgres -e POSTGRES_PASSWORD=postgres \
  ord-postgres:empty
docker exec -it ord-postgres ./build_database.sh
docker commit ord-postgres openreactiondatabase/ord-postgres
docker stop ord-postgres
# Uncomment the next line to push the new image to Docker Hub.
# docker push openreactiondatabase/ord-postgres
