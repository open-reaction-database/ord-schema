#!/bin/bash
# Copyright 2022 Open Reaction Database Project Authors
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

# Compiles protocol buffers with the plugins pinned in buf.gen.yaml. Requires Node.js.
set -ex
cd "$(dirname "${BASH_SOURCE[0]}")"

# buf comes from the version package-lock.json pins. It needs no install script:
# @bufbuild/buf reaches its binary through a platform package.
npm ci --ignore-scripts --no-audit --no-fund --loglevel=error
export PATH="${PWD}/node_modules/.bin:${PATH}"

# proto/ is the module root in buf.yaml, and a file's path under it is what buf embeds in
# the descriptor and places generated code by. proto/ord-schema/proto/reaction.proto yields
# ord_schema/proto/reaction_pb2.py, since the Python generator spells the hyphen as an
# underscore.
buf generate
