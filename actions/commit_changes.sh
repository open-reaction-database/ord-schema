#!/bin/bash
# Copyright 2020 The Open Reaction Database Authors
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

#
# Pushes a commit containing all changes.
#
# NOTE(kearnes): This works with `push` event triggers (not `pull_request`).
set -ex
if (( $# != 1 )); then
  echo "commit message required"
  exit 1
fi
COMMIT_MESSAGE="$1"

git config --local user.email "41898282+github-actions[bot]@users.noreply.github.com"
git config --local user.name "github-actions[bot]"
# Fail gracefully if there is nothing to commit.
git commit -a -m "${COMMIT_MESSAGE}" || (( $? == 1 ))
git push "https://${GITHUB_ACTOR}:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY}.git" "HEAD:${GITHUB_REF}"
