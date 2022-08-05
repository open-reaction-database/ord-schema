"""Tests for ord_schema.orm.reaction."""
from sqlalchemy import create_engine

from ord_schema.orm.reaction import Base


def test_ddl():
    engine = create_engine("postgresql://postgres:postgres@localhost:5433/test", echo=True, future=True)
    Base.metadata.create_all(engine)


if __name__ == "__main__":
    test_ddl()
