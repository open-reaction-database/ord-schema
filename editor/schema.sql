-- Copyright 2020 Open Reaction Database Project Authors
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

-- $ psql postgres -p 5430 -f schema.sql
-- $ ./py/migrate.py

DROP DATABASE IF EXISTS editor;

\c postgres;

CREATE DATABASE editor;

\c editor;

CREATE TABLE users (
  user_id CHARACTER(32) PRIMARY KEY,
  name TEXT,
  created_time INTEGER NOT NULL
);

CREATE TABLE logins (
  access_token TEXT PRIMARY KEY,
  user_id CHARACTER(32) REFERENCES USERS,
  timestamp INTEGER NOT NULL
);

CREATE TABLE datasets (
  user_id CHARACTER(32) REFERENCES USERS,
  dataset_name TEXT NOT NULL,
  pbtxt TEXT NOT NULL,
  PRIMARY KEY (user_id, dataset_name)
);

-- System users:
--   "review" owns read-only datasets imported from GitHub pull requests.
--   "test" owns datasets imported from db/ and used only in tests.
INSERT INTO users VALUES 
   ('8df09572f3c74dbcb6003e2eef8e48fc', 'review', EXTRACT(EPOCH FROM NOW())),
   ('680b0d9fe649417cb092d790907bd5a5', 'test', EXTRACT(EPOCH FROM NOW()));
