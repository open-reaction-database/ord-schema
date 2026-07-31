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

Everything an agent needs to turn an intent into a query lives here: prompts, the
schemas a model is allowed to emit, name resolution, and -- as they arrive -- tool
and MCP surfaces.

It sits in ``ord-schema`` rather than in a serving repository because these pieces
describe the *data*, not a deployment. A prompt that tells a model which columns exist
has to be versioned alongside the code that generates them, or it silently drifts from
a schema it cannot see. For the same reason nothing here does HTTP, holds a database
connection, or owns a cache; a server supplies those and maps library exceptions onto
its own protocol.

Requires the ``agent`` extra (``pip install ord-schema[agent]``).
"""
