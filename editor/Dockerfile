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

FROM continuumio/miniconda3
WORKDIR /usr/src/app

RUN apt-get update && apt-get install --yes build-essential unzip
RUN conda install -c rdkit rdkit
RUN conda install flask nodejs

# Fetch/build editor dependencies.
RUN git clone https://github.com/Open-Reaction-Database/ketcher.git
WORKDIR ketcher
RUN npm install && npm run build
WORKDIR ..
RUN wget https://github.com/google/closure-library/archive/v20200517.tar.gz
RUN tar -xzf v20200517.tar.gz
RUN wget https://storage.googleapis.com/ord-editor-test/editor_test_protobuf_1dae8fdd.tar
RUN tar -xzf editor_test_protobuf_1dae8fdd.tar

# Bust the cache to get the latest commits; see
# https://stackoverflow.com/a/39278224.
ADD https://api.github.com/repos/Open-Reaction-Database/ord-schema/git/refs/heads/main ord-schema-version.json
RUN git clone https://github.com/Open-Reaction-Database/ord-schema.git
RUN mv closure-library-20200517 ord-schema/editor/
RUN mv ketcher ord-schema/editor/
RUN mv protobuf ord-schema/editor/

# Build the editor.
WORKDIR ord-schema
RUN pip install -r requirements.txt
RUN python setup.py install
WORKDIR editor
# Uncomment these lines to copy over any local changes in ord-schema/editor.
#COPY css/*.css css/
#COPY html/*.html html/
#COPY js/*.js js/
#COPY py/*.py py/
#COPY serve.sh .
RUN make

ENV ORD_EDITOR_PORT=5000
EXPOSE "${ORD_EDITOR_PORT}"
CMD ./serve.sh --host=0.0.0.0 --port="${ORD_EDITOR_PORT}"
