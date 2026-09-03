# Copyright 2022 Open Reaction Database Project Authors
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

"""Pytest fixtures.

One PostgreSQL cluster per session, and a database cloned from a template per test.

Starting a cluster costs 0.47s, installing the ORM schema and the RDKit cartridge
another 0.14s, and loading the example dataset with its derived tables 1.0s more: 1.6s
before a test body runs, and this package has 74 of them. ``CREATE DATABASE ...
TEMPLATE`` copies a prepared database in 0.06s, so the templates are built once per
session and every test clones one.

A clone is a byte copy, so a test owns a database nobody else writes to, and dropping it
leaves the template as it was. Under xdist each worker is its own session and so runs
its own cluster.
"""

import itertools
import pathlib
import re
from collections.abc import Callable, Iterator

import pytest
from sqlalchemy import create_engine, text
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session
from testing.postgresql import Postgresql

from ord_schema.datasets import load_dataset
from ord_schema.orm.database import add_dataset, prepare_database, update_derived_data

_SERIAL = itertools.count(1)
_PREPARED = "template_prepared"
_LOADED = "template_loaded"


def _url(postgres: Postgresql, database: str | None = None) -> str:
    """Returns the psycopg URL for one database of the session's cluster.

    Args:
        postgres: The running cluster.
        database: Which database, or None for the cluster's own.

    Returns:
        The URL, with the driver named. See
        https://docs.sqlalchemy.org/en/20/dialects/postgresql.html#module-sqlalchemy.dialects.postgresql.psycopg.
    """
    raw = postgres.url() if database is None else postgres.url(database=database)
    return re.sub("postgresql://", "postgresql+psycopg://", raw)


def _engine(postgres: Postgresql, database: str) -> Engine:
    """Returns an engine on one database of the session's cluster.

    Args:
        postgres: The running cluster.
        database: Which database to connect to.

    Returns:
        An engine.
    """
    return create_engine(_url(postgres, database), future=True)


def _clone(postgres: Postgresql, template: str, name: str) -> None:
    """Copies ``template`` into a new database.

    Args:
        postgres: The running cluster.
        template: The database to copy. It must have no open connections, which is
            why every builder below disposes its engine before returning.
        name: The database to create.
    """
    # AUTOCOMMIT because CREATE DATABASE cannot run inside a transaction block.
    admin = create_engine(_url(postgres), future=True, isolation_level="AUTOCOMMIT")
    try:
        with admin.connect() as connection:
            connection.execute(text(f'CREATE DATABASE "{name}" TEMPLATE "{template}"'))
            # Database-scoped settings are keyed by database OID in pg_db_role_setting,
            # which a copy does not carry: the clone gets a new OID and reverts to the
            # server defaults. prepare_database pins search_path that way, so without
            # this a clone resolves the ord schema ahead of public wherever the role is
            # named "ord".
            #
            # Read off the template rather than named here, so whatever
            # prepare_database sets is carried without a second list of it to keep in
            # step. A setting it adds later needs no change on this side.
            settings = connection.execute(
                text(
                    "SELECT unnest(setconfig) FROM pg_db_role_setting s "
                    "JOIN pg_database d ON d.oid = s.setdatabase "
                    "WHERE d.datname = :template AND s.setrole = 0"
                ),
                {"template": template},
            ).scalars()
            for setting in settings:
                key, _, value = setting.partition("=")
                connection.execute(
                    text(f'ALTER DATABASE "{name}" SET {key} TO {value}')
                )
    finally:
        admin.dispose()


def _drop(postgres: Postgresql, name: str) -> None:
    """Removes a database the session is done with.

    Args:
        postgres: The running cluster.
        name: The database to drop.
    """
    admin = create_engine(_url(postgres), future=True, isolation_level="AUTOCOMMIT")
    try:
        with admin.connect() as connection:
            connection.execute(text(f'DROP DATABASE IF EXISTS "{name}"'))
    finally:
        admin.dispose()


@pytest.fixture(scope="session", name="postgres")
def postgres_fixture() -> Iterator[Postgresql]:
    """The cluster every database in this session lives in."""
    with Postgresql() as postgres:
        yield postgres


@pytest.fixture(scope="session", name="prepared_template")
def prepared_template_fixture(postgres: Postgresql) -> str:
    """Builds a template holding the ORM schema, and returns its name."""
    _clone(postgres, "template1", _PREPARED)
    engine = _engine(postgres, _PREPARED)
    try:
        assert prepare_database(engine)
    finally:
        engine.dispose()
    return _PREPARED


@pytest.fixture(scope="session", name="loaded_template")
def loaded_template_fixture(postgres: Postgresql, prepared_template: str) -> str:
    """Builds a template holding the schema and the example dataset, and names it."""
    _clone(postgres, prepared_template, _LOADED)
    engine = _engine(postgres, _LOADED)
    try:
        dataset = load_dataset(
            pathlib.Path(__file__).parent / "testdata" / "ord-nielsen-example.pbtxt",
            as_dataset=True,
        )
        # Read back rather than assumed: prepare_database reports whether the RDKit
        # cartridge is installed, and the derived tables differ if it is not.
        with engine.connect() as connection:
            cartridge = bool(
                connection.execute(
                    text("SELECT to_regnamespace('rdkit') IS NOT NULL")
                ).scalar()
            )
        with Session(engine) as session, session.begin():
            add_dataset(dataset, session)
            update_derived_data(dataset.dataset_id, session, rdkit_cartridge=cartridge)
    finally:
        engine.dispose()
    return _LOADED


@pytest.fixture(name="databases")
def databases_fixture(
    postgres: Postgresql, request: pytest.FixtureRequest
) -> Iterator[Callable[[str], tuple[Engine, str]]]:
    """Returns a factory cloning a template, and cleans up everything it made.

    Every database a test runs against comes from here, so disposal and dropping live in
    one place. Neither is optional: a pooled connection outlives the test and a
    session's worth reaches the server's connection limit, and a clone is a full copy --
    11.7 MB for the prepared schema, so a run of this package would leave most of a
    gigabyte behind.

    Args:
        postgres: The running cluster.
        request: Names the databases after the test, so a failure says which one it ran
            against. Truncated to fit PostgreSQL's 63-byte identifier, with a serial
            prefix because parametrized names collide once truncated.

    Yields:
        A callable taking a template name and returning ``(engine, url)`` for a fresh
        clone of it.
    """
    stem = re.sub(r"[^A-Za-z0-9_]", "_", request.node.name)[:40]
    made: list[tuple[Engine, str]] = []

    def build(template: str) -> tuple[Engine, str]:
        name = f"t{next(_SERIAL)}_{stem}"
        _clone(postgres, template, name)
        engine = _engine(postgres, name)
        made.append((engine, name))
        return engine, _url(postgres, name)

    try:
        yield build
    finally:
        for engine, name in made:
            # Dropped only after the engine's pool is closed; PostgreSQL refuses to drop
            # a database anything is still connected to.
            engine.dispose()
            _drop(postgres, name)


@pytest.fixture(name="test_engine")
def test_engine_fixture(databases: Callable[[str], tuple[Engine, str]]) -> Engine:
    """An engine on an empty database of this session's cluster."""
    engine, _ = databases("template1")
    return engine


@pytest.fixture(name="prepared_engine")
def prepared_engine_fixture(
    databases: Callable[[str], tuple[Engine, str]], prepared_template: str
) -> Engine:
    """``test_engine`` with the ORM schema (and RDKit cartridge) installed."""
    engine, _ = databases(prepared_template)
    return engine


@pytest.fixture(name="test_session")
def test_session_fixture(
    databases: Callable[[str], tuple[Engine, str]], loaded_template: str
) -> Iterator[Session]:
    """A session on a database already holding the example dataset."""
    engine, _ = databases(loaded_template)
    with Session(engine) as session:
        yield session


@pytest.fixture(name="prepared_database")
def prepared_database_fixture(
    databases: Callable[[str], tuple[Engine, str]], prepared_template: str
) -> Callable[[], tuple[Engine, str]]:
    """Returns a factory making a fresh prepared database on each call.

    For a test that needs more than one database at a time -- comparing two ingest paths
    on identical schemas -- or that hands a URL to a subprocess rather than reusing this
    process's engine.

    Args:
        databases: The factory that owns cleanup.
        prepared_template: The template to clone.

    Returns:
        A callable returning ``(engine, url)`` for a database nobody else writes to.
    """
    return lambda: databases(prepared_template)
