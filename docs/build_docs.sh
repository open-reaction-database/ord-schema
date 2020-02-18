#!/bin/bash

set -ex

docker pull pseudomuto/protoc-gen-doc
docker run --rm -v "$(pwd)":/out -v "$(pwd)/../proto":/protos pseudomuto/protoc-gen-doc --doc_opt=html,index.html
