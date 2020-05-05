#!/bin/bash
#
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

# Run sphinx.
pip install -U sphinx recommonmark sphinx-markdown-tables
make clean html
