l
-- psql postgres -p 5430 -f schema.sql

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
  PBTXT TEXT NOT NULL,
  PRIMARY KEY (user_id, dataset_name)
);
