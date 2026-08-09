# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Machinery for reaching ORD from an agent.

Everything an agent needs to turn an intent into a query lives here: the schema a model
is told about, the checks its query has to pass, name resolution, and -- as they arrive
-- prompts, tool, and MCP surfaces.

It sits in ``ord-schema`` rather than in a serving repository because these pieces
describe the *data*, not a deployment. A description that tells a model which columns
exist has to be versioned alongside the code that generates them, or it silently drifts
from a schema it cannot see; :mod:`ord_schema.agent.schema` and
:mod:`ord_schema.agent.sql` both read ``projection.SCHEMA`` for exactly that reason. For
the same reason nothing here does HTTP or owns a cache; a server supplies those and
maps library exceptions onto its own protocol. The one connection held is DuckDB's
embedded one, opened by :mod:`ord_schema.agent.execute` over local artifact files --
part of reading the data, not of serving it.

Requires the ``agent`` extra (``pip install ord-schema[agent]``).
"""
