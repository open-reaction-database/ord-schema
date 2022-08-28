"""Functions for creating/managing the PostgreSQL database."""
import sqlalchemy
from sqlalchemy import orm

from ord_schema.orm import mappers
from ord_schema.proto import dataset_pb2


def create_engine(username: str, password: str, host: str, port: int, database: str) -> sqlalchemy.engine.Engine:
    """Creates an SQLAlchemy Engine."""
    return sqlalchemy.create_engine(f"postgresql://{username}:{password}@{host}:{port}/{database}")


def create_database(engine: sqlalchemy.engine.Engine) -> None:
    mappers.Base.metadata.create_all(engine)


def add_dataset(dataset: dataset_pb2.Dataset, session: orm.Session) -> None:
    session.add(mappers.from_proto(dataset))


def add_datasets(datasets: list[dataset_pb2.Dataset], engine: sqlalchemy.engine.Engine) -> None:
    with orm.Session(engine) as session:
        for dataset in datasets:
            add_dataset(dataset, session)
        session.commit()
