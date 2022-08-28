"""Functions for creating/managing the PostgreSQL database."""
import sqlalchemy
from sqlalchemy import func, text, update
from sqlalchemy.orm import Session

from ord_schema.orm import mappers
from ord_schema.proto import dataset_pb2


def create_engine(
    database: str,
    username: str,
    password: str,
    host: str = "localhost",
    port: int = 5432,
    echo: bool = True,
) -> sqlalchemy.engine.Engine:
    """Creates an SQLAlchemy Engine."""
    return sqlalchemy.create_engine(
        f"postgresql://{username}:{password}@{host}:{port}/{database}", echo=echo, future=True
    )


def prepare_database(engine: sqlalchemy.engine.Engine) -> None:
    """Prepares the database and creates the ORM table structure."""
    with engine.connect() as connection:
        connection.execute(text("CREATE SCHEMA IF NOT EXISTS rdkit"))
        connection.execute(text("CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA rdkit"))
    mappers.Base.metadata.create_all(engine)


def add_datasets(datasets: list[dataset_pb2.Dataset], engine: sqlalchemy.engine.Engine) -> None:
    """Adds datasets to the database."""
    with Session(engine) as session:
        for dataset in datasets:
            session.add(mappers.from_proto(dataset))
        session.commit()


def add_rdkit(engine: sqlalchemy.engine.Engine) -> None:
    """Adds RDKit PostgreSQL cartridge data."""
    table = mappers.Structure.__table__
    with Session(engine) as session:
        session.execute(
            update(table).values(mol=func.rdkit.mol_from_smiles(sqlalchemy.cast(table.c.smiles, mappers.CString)))
        )
        session.execute(update(table).values(morgan_binary_fingerprint=func.rdkit.morganbv_fp(table.c.mol)))
        session.commit()
