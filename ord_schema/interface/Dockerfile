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

# To build the interface using the current state of ord-schema:
# $ cd path/to/ord-schema
# $ docker build \
#     --file=ord_schema/interface/Dockerfile \
#     -t openreactiondatabase/ord-interface \
#     .
#
# To push the built image to Docker Hub:
# docker push openreactiondatabase/ord-interface

FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y \
    build-essential \
    procps \
    unzip \
 && rm -rf /var/lib/apt/lists/*

RUN conda install -c rdkit \
    flask \
    gunicorn \
    rdkit \
 && conda clean -afy

# COPY and install the local state of ord-schema.
WORKDIR /usr/src/app/ord-schema
COPY setup.py requirements.txt ./
COPY ord_schema/ ord_schema/
RUN pip install -r requirements.txt
RUN python setup.py install

# Build and launch the interface.
WORKDIR /usr/src/app/ord-schema/ord_schema/interface
EXPOSE 5000
CMD gunicorn search:app -b 0.0.0.0:5000 --access-logfile '-'
