#!/bin/bash
#
# Finds new and changed PBTXT files and sets NUM_CHANGED_FILES for downstream
# workflow steps.
set -ex
UPSTREAM="https://github.com/Open-Reaction-Database/${GITHUB_REPOSITORY##*/}.git"

git remote add upstream "${UPSTREAM}"
echo "Current branch: $(git rev-parse --abbrev-ref HEAD)"
git fetch --no-tags --prune --depth=1 upstream +refs/heads/*:refs/remotes/upstream/*
git diff --name-only upstream/master > changed_files.txt
# Use `|| [[ $? == 1 ]]` in case no lines match and the exit code is nonzero.
grep -e "\.pbtxt$" changed_files.txt > changed_pbtxt_files.txt || (( $? == 1 ))
# Use LOCAL_NUM_CHANGED since ::set-env values are not available immediately.
LOCAL_NUM_CHANGED="$(wc -l < changed_pbtxt_files.txt | tr -d ' ')"
echo "::set-env name=NUM_CHANGED_FILES::${LOCAL_NUM_CHANGED}"
echo "Found ${LOCAL_NUM_CHANGED} changed pbtxt files"
cat changed_pbtxt_files.txt
