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

"""Searching the ORD corpus: a query grammar, its compiler, and an executor.

A query is stated as a :class:`~ord_schema.search.query.Query` rather than as SQL, so
what the grammar can express bounds what a search can cost. Everything needed to turn an
intent into one lives here: the schema description a model is given to write against,
the checks a query has to pass, name resolution, and the executor that answers it over
the projection and structures artifacts.

The consumer this was built for is an agent, which is why the schema is rendered for a
prompt and why compounds are named rather than spelled. Nothing here is specific to one,
though: a search is a search.

It sits in ``ord-schema`` rather than in a serving repository because these pieces
describe the *data*, not a deployment. A description that tells a model which columns
exist has to be versioned alongside the code that generates them, or it silently drifts
from a schema it cannot see; :mod:`ord_schema.search.schema` and
:mod:`ord_schema.search.sql` both read ``projection.SCHEMA`` for exactly that reason.
For the same reason nothing here does HTTP or owns a cache; a server supplies those and
maps library exceptions onto its own protocol. The one connection held is DuckDB's
embedded one, opened by :mod:`ord_schema.search.execute` over local artifact files --
part of reading the data, not of serving it.

Requires the ``search`` extra (``pip install ord-schema[search]``).
"""
