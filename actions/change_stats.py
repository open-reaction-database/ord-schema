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
"""Computes diff statistics for ORD submissions."""

import os
import re
import subprocess

from absl import app
from absl import flags
from absl import logging
import github

FLAGS = flags.FLAGS
flags.DEFINE_string('base', '', 'Branch to diff against.')
flags.DEFINE_integer('issue', None,
                     'Issue number. If provided, a comment will be added.')


def main(argv):
    del argv  # Only used by app.run().
    pattern = re.compile(r'^([+-])\s+reaction_id: "(ord-[0-9a-f]{32})"',
                         re.MULTILINE)
    args = ['git', 'diff']
    if FLAGS.base:
        args.append(FLAGS.base)
    output = subprocess.run(args, capture_output=True, check=True, text=True)
    added, removed = set(), set()
    for match in re.finditer(pattern, output.stdout):
        prefix, reaction_id = match.groups()
        if prefix == '+':
            if reaction_id in added:
                raise ValueError(f'duplicate reaction ID: {reaction_id}')
            added.add(reaction_id)
        elif prefix == '-':
            if reaction_id in removed:
                raise ValueError(f'duplicate reaction ID: {reaction_id}')
            removed.add(reaction_id)
        else:
            raise ValueError(f'unrecognized diff prefix: {prefix}')
    logging.info(f'Summary: +%d -%d reaction IDs', len(added), len(removed))
    if FLAGS.issue:
        client = github.Github(os.environ['GITHUB_TOKEN'])
        repo = client.get_repo(os.environ['GITHUB_REPOSITORY'])
        issue = repo.get_issue(FLAGS.issue)
        issue.create_comment(
            f'Summary: +{len(added)} -{len(removed)} reaction IDs')


if __name__ == '__main__':
    app.run(main)
