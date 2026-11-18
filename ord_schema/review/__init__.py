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

"""Peer review of a dataset submission: checks, scoring, and the reviewer's report.

Dataset-level findings that single-reaction validation cannot make, a completeness
score, exemplar reactions, and the report bundle and AI review packet a reviewer
reads. Everything here informs a human reviewer, who has the final say; none of it
blocks a submission.

Builds on the summaries in :mod:`ord_schema.analysis`, never the reverse.

Requires the ``review`` extra (``pip install ord-schema[review]``), which installs
``analysis`` too.
"""
