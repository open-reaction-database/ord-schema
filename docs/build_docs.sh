#!/bin/bash
#
# Generates HTML documentation for the Open Reaction Database schema.
# See https://github.com/pseudomuto/protoc-gen-doc.

set -ex

# Build HTML for protocol buffers.
docker pull pseudomuto/protoc-gen-doc
docker run --rm -v "$(pwd)":/out -v "$(pwd)/../proto":/protos \
  pseudomuto/protoc-gen-doc --doc_opt=markdown,protos.md:test.proto

# Run sphinx.
pip install -U sphinx recommonmark sphinx-markdown-tables
make clean html
