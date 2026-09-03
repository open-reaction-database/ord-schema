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

A cluster costs 0.47s to start, the ORM schema and the RDKit cartridge another 0.14s,
and loading the example dataset with its derived tables 1.0s more. Paid per test that is
1.6s before the test body runs, which over this package is most of the two minutes the
suite takes. ``CREATE DATABASE ... TEMPLATE`` copies a prepared database in 0.06s, so
the templates are built once and every test clones one.

Isolation is unchanged: a clone is a byte copy, so a test still gets a database nobody
else writes to, and dropping it afterwards leaves the template as it was. Under xdist
each worker is its own session and so runs its own cluster.
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


def _engine(postgres: Postgresql, database: str) -> Engine:
    """Returns an engine on one database of the session's cluster.

    Args:
        postgres: The running cluster.
        database: Which database to connect to.

    Returns:
        An engine. See
        https://docs.sqlalchemy.org/en/20/dialects/postgresql.html#module-sqlalchemy.dialects.postgresql.psycopg.
    """
    return create_engine(
        re.sub(
            "postgresql://", "postgresql+psycopg://", postgres.url(database=database)
        ),
        future=True,
    )


def _clone(postgres: Postgresql, template: str, name: str) -> None:
    """Copies ``template`` into a new database.

    Args:
        postgres: The running cluster.
        template: The database to copy. It must have no open connections, which is
            why every builder below disposes its engine before returning.
        name: The database to create.
    """
    # AUTOCOMMIT because CREATE DATABASE cannot run inside a transaction block.
    admin = create_engine(
        re.sub("postgresql://", "postgresql+psycopg://", postgres.url()),
        future=True,
        isolation_level="AUTOCOMMIT",
    )
    try:
        with admin.connect() as connection:
            connection.execute(text(f'CREATE DATABASE "{name}" TEMPLATE "{template}"'))
            # Database-scoped settings are keyed by database OID in pg_db_role_setting,
            # which a copy does not carry: the clone gets a new OID and reverts to the
            # server default. prepare_database pins search_path this way, so a clone
            # would otherwise resolve the ord schema ahead of public wherever the role
            # is named "ord". Reapplied here rather than left to the caller, since
            # every database this module hands out is a clone.
            connection.execute(
                text(f'ALTER DATABASE "{name}" SET search_path TO public')
            )
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


@pytest.fixture(name="database_name")
def database_name_fixture(request: pytest.FixtureRequest) -> str:
    """Returns a database name unique to this test.

    Named after the test so a failure says which database it ran against, and prefixed
    with a serial because PostgreSQL truncates identifiers at 63 bytes and parametrized
    names collide once truncated. The counter is per process, which is per xdist worker,
    which is per cluster.
    """
    stem = re.sub(r"[^A-Za-z0-9_]", "_", request.node.name)[:40]
    return f"t{next(_SERIAL)}_{stem}"


@pytest.fixture(name="test_engine")
def test_engine_fixture(postgres: Postgresql, database_name: str) -> Iterator[Engine]:
    """An engine on an empty database of this session's cluster."""
    _clone(postgres, "template1", database_name)
    engine = _engine(postgres, database_name)
    try:
        yield engine
    finally:
        engine.dispose()


@pytest.fixture(name="prepared_engine")
def prepared_engine_fixture(
    postgres: Postgresql, prepared_template: str, database_name: str
) -> Iterator[Engine]:
    """``test_engine`` with the ORM schema (and RDKit cartridge) installed."""
    _clone(postgres, prepared_template, database_name)
    engine = _engine(postgres, database_name)
    try:
        yield engine
    finally:
        engine.dispose()


@pytest.fixture(name="test_session")
def test_session_fixture(
    postgres: Postgresql, loaded_template: str, database_name: str
) -> Iterator[Session]:
    """A session on a database already holding the example dataset."""
    _clone(postgres, loaded_template, database_name)
    engine = _engine(postgres, database_name)
    try:
        with Session(engine) as session:
            yield session
    finally:
        engine.dispose()


@pytest.fixture(name="prepared_database")
def prepared_database_fixture(
    postgres: Postgresql, prepared_template: str, database_name: str
) -> Iterator[Callable[[], tuple[Engine, str]]]:
    """Returns a factory making a fresh prepared database, and disposes what it made.

    For a test that needs more than one database at a time -- comparing two ingest paths
    on identical schemas -- or that hands a URL to a subprocess rather than reusing this
    process's engine.

    Disposal is here rather than in the test because ``Engine`` is not a context manager
    and has no ``close``, so every caller would otherwise carry its own ``try/finally``.
    It is not optional: a pooled connection outlives the test, and a session's worth of
    them reaches the server's connection limit.

    Args:
        postgres: The running cluster.
        prepared_template: The template to clone.
        database_name: Names the first database; later ones get a suffix.

    Yields:
        A callable returning ``(engine, url)`` for a database nobody else writes to.
    """
    made: list[Engine] = []

    def build() -> tuple[Engine, str]:
        name = database_name if not made else f"{database_name}_{len(made) + 1}"
        _clone(postgres, prepared_template, name)
        url = re.sub(
            "postgresql://", "postgresql+psycopg://", postgres.url(database=name)
        )
        engine = create_engine(url, future=True)
        made.append(engine)
        return engine, url

    try:
        yield build
    finally:
        for engine in made:
            engine.dispose()
