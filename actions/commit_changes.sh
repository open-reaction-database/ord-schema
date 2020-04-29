#!/bin/bash
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
