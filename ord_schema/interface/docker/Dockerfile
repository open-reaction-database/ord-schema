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

# To build an updated ord-postgres image containing the latest data:
# $ cd path/to/ord-schema
# $ ./ord_schema/interface/docker/update_image.sh
#
# To push the built image to Docker Hub:
# docker push openreactiondatabase/ord-postgres

FROM mcs07/postgres-rdkit

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    libpq-dev \
    procps \
    unzip \
    wget \
 && rm -rf /var/lib/apt/lists/*

ENV ORD_ROOT=/home/ord
RUN wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
 && /bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p "${ORD_ROOT}/miniconda3" \
 && rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="${ORD_ROOT}/miniconda3/bin:${PATH}"
RUN conda install -c rdkit \
    flask \
    gunicorn \
    nodejs \
    python=3.7 \
    rdkit \
 && conda clean -afy

WORKDIR "${ORD_ROOT}"
# Bust the cache to get the latest ord-data commits; see
# https://stackoverflow.com/a/39278224.
ENV ORD_REPO=ord-data
ADD "https://api.github.com/repos/Open-Reaction-Database/${ORD_REPO}/git/refs/heads/main" ord-data-version.json
RUN git clone "https://github.com/Open-Reaction-Database/${ORD_REPO}.git"

# COPY and install the local state of ord-schema.
WORKDIR "${ORD_ROOT}/ord-schema"
COPY setup.py requirements.txt ./
COPY ord_schema/ ord_schema/
RUN pip install -r requirements.txt
RUN python setup.py install

# NOTE(kearnes): Use non-standard PGDATA so the database persists; see
# https://nickjanetakis.com/blog/docker-tip-79-saving-a-postgres-database-in-a-docker-image.
WORKDIR "${ORD_ROOT}"
ENV PGDATA=/data
ENV POSTGRES_DB=ord
ENV POSTGRES_USER=ord-postgres
ENV POSTGRES_PASSWORD=ord-postgres
COPY ord_schema/interface/docker/postgres.conf /etc/postgresql/postgresql.conf
ENV ORD_INPUT="${ORD_ROOT}/${ORD_REPO}/data/*/*.pb"
COPY ord_schema/interface/docker/build_database.sh .
