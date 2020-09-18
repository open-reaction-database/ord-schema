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
"""Constants for the ORD interface backend."""

import collections

POSTGRES_PORT = 5432
POSTGRES_USER = 'ord-postgres'
POSTGRES_PASSWORD = 'ord-postgres'
POSTGRES_DB = 'ord'

RDKIT_SCHEMA = 'rdk'

# Postgres table schema.
TABLES = {
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
