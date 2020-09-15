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

# To build the editor using the current state of ord-schema:
# $ cd path/to/ord-schema
# $ docker build --file=editor/Dockerfile -t openreactiondatabase/ord-editor .
#
# To test the editor:
# $ docker run --rm -d -p 5000:5000 --name editor openreactiondatabase/ord-editor
# $ node editor/js/test.js
# $ docker stop editor
#
# To push the built image to Docker Hub:
# docker push openreactiondatabase/ord-editor

FROM continuumio/miniconda3

RUN apt-get update && apt-get install --yes build-essential procps unzip
RUN conda install -c rdkit rdkit
RUN conda install flask nodejs

# Fetch/build editor dependencies.
# NOTE(kearnes): Do this before COPYing the local state so it can be cached.
WORKDIR /usr/src/app/ord-schema/editor
RUN git clone https://github.com/Open-Reaction-Database/ketcher.git
WORKDIR ketcher
RUN npm install && npm run build
WORKDIR ..
RUN wget https://github.com/google/closure-library/archive/v20200517.tar.gz
RUN tar -xzf v20200517.tar.gz
RUN wget https://github.com/protocolbuffers/protobuf/releases/download/v3.13.0/protobuf-js-3.13.0.tar.gz
RUN tar -xzf protobuf-js-3.13.0.tar.gz

# COPY and install the local state of ord-schema.
WORKDIR /usr/src/app/ord-schema
COPY setup.py requirements.txt ./
COPY editor/ editor/
COPY ord_schema/ ord_schema/
RUN pip install -r requirements.txt
RUN python setup.py install

# Build and launch the editor.
WORKDIR /usr/src/app/ord-schema/editor
RUN make
EXPOSE 5000
RUN pip install gunicorn
CMD gunicorn py.serve:app -b 0.0.0.0:5000 --access-logfile '-'
