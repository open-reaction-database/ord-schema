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

# This is a base image used to speed up GitHub test builds for ord-editor.
#
# $ docker build -t openreactiondatabase/ord-editor-base .
# $ docker push openreactiondatabase/ord-editor-base

FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y \
    build-essential \
    procps \
    unzip \
 && rm -rf /var/lib/apt/lists/*

RUN conda install -c rdkit \
    flask \
    gunicorn \
    nodejs \
    rdkit \
 && conda clean -afy

# Fetch and build editor dependencies.
WORKDIR /usr/src/app/ord-schema/editor
RUN git clone https://github.com/Open-Reaction-Database/ketcher.git \
 && cd ketcher \
 && npm install \
 && npm run build \
 && rm -rf node_modules
RUN wget https://github.com/google/closure-library/archive/v20200517.tar.gz \
 && tar -xzf v20200517.tar.gz \
 && rm v20200517.tar.gz
RUN wget https://github.com/protocolbuffers/protobuf/releases/download/v3.13.0/protobuf-js-3.13.0.tar.gz \
 && tar -xzf protobuf-js-3.13.0.tar.gz \
 && rm protobuf-js-3.13.0.tar.gz

WORKDIR ..
