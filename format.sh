#!/bin/bash
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

# Runs code formatters.
set -ex
# See https://www.ostricher.com/2014/10/the-right-way-to-get-the-directory-of-a-bash-script/.
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)"
# Add missing license headers.
if command -v go &> /dev/null; then
  go install github.com/google/addlicense@latest
  "${HOME}/go/bin/addlicense" \
    -c "Open Reaction Database Project Authors" \
    -l apache "${ROOT_DIR}"
else
  echo "Please install Go; see https://golang.org/doc/install"
fi
# Format python.
if ! command -v yapf &> /dev/null; then
  pip install yapf
fi
yapf -p -r "${ROOT_DIR}" --exclude="*_pb2.py" --in-place
# Format proto.
if command -v clang-format-10 &> /dev/null; then
  find "${ROOT_DIR}" -name '*.proto' -exec clang-format-10 -i --style=file {} +
elif command -v clang-format &> /dev/null; then
  # NOTE(kearnes): Make sure you have version 10 or higher!
  find "${ROOT_DIR}" -name '*.proto' -exec clang-format -i --style=file {} +
else
  echo "Please install clang-format:"
  echo "  Linux: apt install clang-format-10"
  echo "  MacOS: brew install clang-format"
fi
