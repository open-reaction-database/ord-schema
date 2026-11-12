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

"""Summaries of a dataset's contents, for anyone exploring or preparing data.

What a dataset holds, rather than whether it is right: which fields its reactions
populate and with what values, the compounds it uses and in which roles, and the
same reactions laid out as a table. These are useful well beyond peer review, so
they live apart from :mod:`ord_schema.review`, which builds its checks on them.

:mod:`ord_schema.review` may import from here; nothing here imports from it.

Requires the ``analysis`` extra (``pip install ord-schema[analysis]``).
"""
