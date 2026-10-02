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

# Compiles protocol buffers with the plugins pinned in buf.gen.yaml.
# Requires buf (https://buf.build/docs/cli/installation/) and, on PATH, protoc-gen-ts from
# ts-protoc-gen and pbjs and pbts from protobufjs-cli; test_proto_wrappers pins all three.
set -ex

# proto/ is the module root in buf.yaml, and a file's path under it is what buf embeds in
# the descriptor and places generated code by. proto/ord-schema/proto/reaction.proto yields
# ord_schema/proto/reaction_pb2.py, since the Python generator spells the hyphen as an
# underscore, and js/ord-schema/proto/reaction_pb.js. The JavaScript files reach each other
# through ../../ord-schema/proto/, which resolves in the source tree and in an installed
# package only because the directory shares the npm package's name.
buf generate

# Node's ESM loader reads a CommonJS module's export names from its source, and the
# generated files attach theirs through goog.object.extend, which it cannot follow. The
# assignment at the end of index.js never runs; it lists every top-level message so Node ESM
# can import each by name. LC_ALL=C keeps the order independent of the locale.
names="$(sed -nE 's/^export class ([A-Za-z0-9_]+) .*/    \1,/p' \
  js/ord-schema/proto/dataset_pb.d.ts \
  js/ord-schema/proto/reaction_pb.d.ts | LC_ALL=C sort)"
cat > js/ord-schema/index.js <<EOF
"use strict";
module.exports = {
    ...require('./proto/dataset_pb'),
    ...require('./proto/reaction_pb'),
};
0 && (module.exports = {
${names}
});
EOF

pbjs -p proto proto/ord-schema/proto/*.proto -o js/ord-schema-protobufjs/index.js -w es6 -t static-module
pbts js/ord-schema-protobufjs/index.js -o js/ord-schema-protobufjs/index.d.ts
