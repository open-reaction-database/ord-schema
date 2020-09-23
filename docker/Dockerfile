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

# This is a base image used to speed up builds for ord-editor, etc.
#
# $ docker build -t openreactiondatabase/ord-base .
# $ docker push openreactiondatabase/ord-base

FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y \
    build-essential \
    procps \
    unzip \
 && rm -rf /var/lib/apt/lists/*

RUN conda install -c rdkit flask gunicorn rdkit && conda clean -afy
