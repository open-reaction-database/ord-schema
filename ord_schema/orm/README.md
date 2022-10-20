# ORD Object Relational Mapper (ORM)

When working with ORD data, the serialized protocol buffers stored
in [ord-data](https://github.com/open-reaction-database/ord-data) should be considered
the [source of truth](https://en.wikipedia.org/wiki/Single_source_of_truth). While protocol buffers are excellent for
storage and retrieval, they are not optimized for searching. To support field-level queries, we use the 
[SQLAlchemy ORM](https://docs.sqlalchemy.org/en/14/orm/quickstart.html) to convert protocol buffers to entries in a 
relational database.

## Overview

Conceptually, an ORM is an abstraction on top of a relational database that allows data to be manipulated using
object-oriented programming techniques. In our case, every protocol buffer message has an associated "mapper" that wraps
a table in a relational database. For example, here is the definition of the `Mass` message in the [protocol buffer
schema](https://github.com/open-reaction-database/ord-schema/blob/main/ord_schema/proto/reaction.proto), which is 
used as a subfield in the `Amount` message:

```protobuf
message Mass {
  enum MassUnit {
    UNSPECIFIED = 0;
    KILOGRAM = 1;
    GRAM = 2;
    MILLIGRAM = 3;
    MICROGRAM = 4;
  }
  optional float value = 1;
  // Precision of the measurement (with the same units as `value`).
  optional float precision = 2;
  MassUnit units = 3;
}

message Amount {
  oneof kind {
    Mass mass = 1;
    Moles moles = 2;
    Volume volume = 3;
    UnmeasuredAmount unmeasured = 5;
  }
  // Whether the volume measurement refers to the pure substance or to the
  // total volume of the reaction input. An example of when this field should
  // be TRUE is when stock solutions are prepared by adding solvent to a
  // volumetric flask already containing solute(s).
  optional bool volume_includes_solutes = 4;
}
```

The `Mass` message has three fields: `value` (float), `precision` (float), and `units` (enum). The corresponding [ORM
mapper](https://github.com/open-reaction-database/ord-schema/blob/main/ord_schema/orm/mappers.py) looks like this:

```python
class Mass(Base):
    # Some attributes are defined by `Base`:
    # __tablename__ = "mass"
    # id = Column(Integer, primary_key=True)
    amount_id = Column(Integer, ForeignKey("amount.id", ondelete="CASCADE"), nullable=False, unique=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Mass.MassUnit.keys(), name="Mass.MassUnit"))
```

The mapper defines a relational database table (`mass`) with five columns: `id` (a unique row ID); `amount_id`
(a reference to the parent `Amount` message in the relational database); and the expected `value`, `precision`, and
`units` fields from the protocol buffer message. Instances of the mapper correspond to individual rows in the associated
table.

## Usage

### Install PostgreSQL and the RDKit extension

We use PostgreSQL with the RDKit extension to support structure searches. There are a couple convenient installation
methods (this list is not an endorsement of any particular provider):

* Conda: [rdkit-postgresql](https://www.rdkit.org/docs/Install.html#installing-and-using-postgresql-and-the-rdkit-postgresql-cartridge-from-a-conda-environment)
  
  ```shell
  $ conda install -c rdkit rdkit-postgresql
  $ export PGDATA="${HOME}/rdkit-postgresql"
  $ initdb -U <username>
  ```

* AWS: [Amazon Aurora PostgreSQL](https://aws.amazon.com/about-aws/whats-new/2020/09/amazon-aurora-postgresql-supports-rdkit-extension/)

### Initialize the database and tables

Use the `prepare_database` function to create the database schemas and tables:

```python
from sqlalchemy import create_engine

from ord_schema.orm.database import prepare_database

connection_string = f"postgresql://{username}:{password}@{host}:{port}/{database}"
engine = create_engine(connection_string, future=True)
prepare_database(engine)
```

Note that rdkit extension data such as fingerprints are stored in a separate `rdkit` schema in the database to avoid
conflicts with ORD message names.

### Add data

Load ORD datasets into the database with the `add_datasets` and `add_rdkit` functions:

```python
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from ord_schema.message_helpers import fetch_dataset
from ord_schema.orm.database import add_datasets, add_rdkit

dataset = fetch_dataset("ord_dataset-fc83743b978f4deea7d6856deacbfe53")

connection_string = f"postgresql://{username}:{password}@{host}:{port}/{database}"
engine = create_engine(connection_string, future=True)
with Session(engine) as session:
    add_datasets([dataset], session)
    session.flush()
    add_rdkit(session)
    session.commit()
```

To load multiple datasets from disk (e.g., from a clone of 
[ord-data](https://github.com/open-reaction-database/ord-data)), use the `add_datasets.py` script:

```shell
$ python add_datasets.py \
    --pattern "path/to/ord-data/data/fc/*.pb.gz" \
    --username <username>
    --host <host>
    --database <database>
```

Note that the database password will be read from the `PGPASSWORD` environment variable if `--password` is not
specified on the command line. To update an existing dataset in the database, use the `--overwrite` flag.

### Run queries

Use the [SQLAlchemy ORM query engine](https://docs.sqlalchemy.org/en/14/orm/quickstart.html#simple-select) to search for
specific fields. When searching for reactions that match a specific query, use the `proto` field of the `Reaction` 
mapper to recover the original `Reaction` protocol buffer message. 

#### Examples

##### Reactions that have at least 70% yield

  ```python
  from sqlalchemy import create_engine, select
  from sqlalchemy.orm import Session
  
  from ord_schema.orm.mappers import Percentage, ProductCompound, ProductMeasurement, Reaction, ReactionOutcome
  from ord_schema.proto import reaction_pb2
  
  connection_string = f"postgresql://{username}:{password}@{host}:{port}/{database}"
  engine = create_engine(connection_string, future=True)
  with Session(engine) as session:
    query = (
        select(Reaction)
        .join(ReactionOutcome)
        .join(ProductCompound)
        .join(ProductMeasurement)
        .join(Percentage)
        .where(ProductMeasurement.type == "YIELD", Percentage.value >= 70)
    )
    results = session.execute(query)
    reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]
  assert len(reactions) == 12
  ```

#### Structure searches with the RDKit PostgreSQL extension

##### Reactions that have Morgan binary fingerprint Tanimoto > 0.5 to `c1ccccc1CCC(O)C`

  ```python
  from sqlalchemy import create_engine, select
  from sqlalchemy.orm import Session
  
  from ord_schema.orm.mappers import Compound, Reaction, ReactionInput
  from ord_schema.orm.structure import FingerprintType, Structure
  from ord_schema.proto import reaction_pb2
  
  connection_string = f"postgresql://{username}:{password}@{host}:{port}/{database}"
  engine = create_engine(connection_string, future=True)
  with Session(engine) as session:
    query = (
      select(Reaction)
      .join(ReactionInput)
      .join(Compound)
      .join(Structure)
      .where(Structure.tanimoto("c1ccccc1CCC(O)C", FingerprintType.MORGAN_BFP) > 0.5)
    )
    results = session.execute(query)
    reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]
  assert len(reactions) == 20
  ```
