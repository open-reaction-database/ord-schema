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

FROM mcs07/postgres-rdkit

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    build-essential git libpq-dev wget

ENV ORD_ROOT=/home/ord
RUN wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN /bin/bash ./Miniconda3-latest-Linux-x86_64.sh -b -p "${ORD_ROOT}/miniconda3"
ENV PATH="${ORD_ROOT}/miniconda3/bin:${PATH}"
RUN conda update conda --yes

WORKDIR "${ORD_ROOT}"
# Bust the cache to get the latest commits; see
# https://stackoverflow.com/a/39278224.
ADD https://api.github.com/repos/Open-Reaction-Database/ord-schema/git/refs/heads/main ord-schema-version.json
RUN git clone https://github.com/Open-Reaction-Database/ord-schema.git
ADD https://api.github.com/repos/Open-Reaction-Database/ord-data/git/refs/heads/main ord-data-version.json
RUN git clone https://github.com/Open-Reaction-Database/ord-data.git
WORKDIR "${ORD_ROOT}/ord-schema"
RUN pip install -r requirements.txt
RUN python setup.py install

# NOTE(kearnes): Use non-standard PGDATA so the database persists; see
# https://nickjanetakis.com/blog/docker-tip-79-saving-a-postgres-database-in-a-docker-image.
WORKDIR "${ORD_ROOT}"
ENV PGDATA=/data
ENV POSTGRES_DB=ord
COPY postgres.conf /etc/postgresql/postgresql.conf
ENV ORD_INPUT="${ORD_ROOT}/ord-data/data/*/*.pbtxt"
COPY build_database.sh .
