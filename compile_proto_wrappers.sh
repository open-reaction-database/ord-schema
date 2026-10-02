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

# Compiles protocol buffers.
# Make sure you have protoc in your PATH; see https://grpc.io/docs/protoc-installation/.
set -ex

# proto/ is the import root. A file's path under it is what protoc embeds in the
# descriptor and names the generated modules after, so it mirrors the ord_schema.proto
# package: proto/ord_schema/proto/reaction.proto yields ord_schema/proto/reaction_pb2.py
# and js/ord_schema/proto/reaction_pb.js.
protoc \
  --proto_path=proto \
  --python_out=. \
  --pyi_out=. \
  --js_out=import_style=commonjs,binary:js \
  proto/ord_schema/proto/*.proto \
  --ts_out=./js/ \

# protoc-gen-js and protoc-gen-ts reach a sibling file by climbing to the import root and
# back down its path, ../../ord_schema/proto/. That leaves the npm package, which installs
# as node_modules/ord-schema, so point those references at the sibling directly.
perl -pi -e 's{\.\./\.\./ord_schema/proto/}{./}g' \
  js/ord_schema/proto/*_pb.js \
  js/ord_schema/proto/*_pb.d.ts

pbjs -p proto proto/ord_schema/proto/*.proto -o js/ord-schema-protobufjs/index.js -w es6 -t static-module
pbts js/ord-schema-protobufjs/index.js -o js/ord-schema-protobufjs/index.d.ts
