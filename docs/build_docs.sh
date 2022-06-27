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

# Generates HTML documentation for the Open Reaction Database schema.
# See https://github.com/pseudomuto/protoc-gen-doc.

set -ex

# Build HTML for protocol buffers.
# NOTE(kearnes): Uncomment this section if we decide to use these docs.
# go get -u github.com/pseudomuto/protoc-gen-doc/cmd/protoc-gen-doc
# protoc \
#   --proto_path=../.. \
#   --plugin=${HOME}/go/bin/protoc-gen-doc \
#   --doc_opt=markdown,protos.md:test.proto \
#   --doc_out=. \
#   ord-schema/proto/reaction.proto ord-schema/proto/dataset.proto

# Generate RST for ord-schema.
sphinx-apidoc -fTM -o ord_schema -t _templates ../ord_schema \
  "../ord_schema/*_test.py" \
  "../ord_schema/*/*_test.py" \
  "../ord_schema/proto/" \
  "../ord_schema/scripts/"

# Run sphinx.
pip install -Ur requirements.txt
make clean html
