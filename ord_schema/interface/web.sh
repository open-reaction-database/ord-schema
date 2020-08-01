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


# Launch a web server that connects to the local Postgres and conducts queries
# on the ORD. The web server is at http://<host>:5000/
#
# First launch local Postgres like this:
#
#   $ ./docker/serve.sh

export PYTHONPATH=.:../build/lib
export FLASK_APP=web.py
export FLASK_ENV=development

python -m flask run $@
