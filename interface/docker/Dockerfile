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

FROM conda/miniconda3

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    git \
    libpq-dev \
    wget \
    build-essential

RUN conda update conda --yes && conda install -c rdkit rdkit-postgresql

ENV ORD_ROOT=/opt
WORKDIR "${ORD_ROOT}"
RUN git clone https://github.com/Open-Reaction-Database/ord-schema.git
RUN git clone https://github.com/Open-Reaction-Database/ord-submissions-test.git
WORKDIR "${ORD_ROOT}/ord-schema"
RUN pip install -r requirements.txt
RUN python setup.py install

ENV PGDATA=/var/lib/postgresql/data
RUN groupadd -r postgres && useradd --no-log-init -r -m -g postgres postgres
RUN mkdir /var/lib/postgresql && chown postgres /var/lib/postgresql
USER postgres
RUN initdb
COPY postgres.conf "${PGDATA}/postgresql.conf"
CMD ["postgres"]
