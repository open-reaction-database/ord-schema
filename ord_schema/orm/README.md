# ORD Object Relational Mapper (ORM)

When working with ORD data, the serialized protocol buffers stored
in [ord-data](https://github.com/open-reaction-database/ord-data) should be considered
the [source of truth](https://en.wikipedia.org/wiki/Single_source_of_truth). While protocol buffers are excellent for
storage and retrieval, they are not optimized for searching. To support field-level queries, we use the
[SQLAlchemy ORM](https://docs.sqlalchemy.org/en/14/orm/quickstart.html) to convert protocol buffers to entries in a
relational database.

## Overview

Conceptually, an ORM is an abstraction on top of a relational database that allows data to be manipulated using
object-oriented programming techniques. In our case, every protocol buffer message has an associated _mapper_ that wraps
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
    amount_id = Column(Integer, ForeignKey("amount.id", ondelete="CASCADE"), nullable=False, unique=True, index=True)

    value = Column(Float)
    precision = Column(Float)
    units = Column(Enum(*reaction_pb2.Mass.MassUnit.keys(), name="Mass.MassUnit"))
```

The mapper defines a relational database table (`mass`) with five columns: `id` (a unique row ID); `amount_id`
(a reference to the parent `Amount` message in the relational database); and the expected `value`, `precision`, and
`units` fields from the protocol buffer message. Instances of the mapper correspond to individual rows in the associated
table.

### Database structure

* Every message in the schema has an associated table in the database. `Mass` messages appear in the `mass` table,
  `ReactionInput` messages appear in the `reaction_input` table, etc.
* Some message types appear in more than one context in the schema. For instance, `ReactionInput` appears as a field
  in both `Reaction` (as `Reaction.inputs`) and `ReactionWorkup` (as `ReactionWorkup.input`) messages. Every table in
  the database contains an `ord_schema_context` column that indicates the context of each message in the ORD schema.
* Every database table has a unique `id` column that is used as the primary key and as the foreign key when defining
  relationships between tables. This database-specific ID should not be confused with the `dataset_id` and `reaction_id`
  fields of `Dataset` and `Reaction`, respectively. Notably, the `CompoundPreparation` and `CrudeComponent` messages
  refer to ORD reaction IDs (`reaction.reaction_id`, _not_ `reaction.id`); these are explicit ORD-level relationships
  that are not part of the database-specific relationship structure.
* The ORM uses four schemas. `ord` holds the decomposed message tables — the search index, i.e. only the proto fields
  plus the structural columns (`id`, the `ord_schema_context` discriminator, map `key`, and parent foreign keys).
  `public` holds API-facing data: `public.reactions` is the serialized `Reaction` proto keyed by `reaction_id`, and
  `public.datasets` holds per-dataset metadata (`md5` for change detection, `num_reactions` and `submitted_at` for
  browsing) keyed by `dataset_id`. `derived` holds search helpers computed from the index — generated SMILES and RDKit
  links in `derived.{reaction,compound,product_compound}_smiles`, plus reaction classes. `rdkit` holds RDKit cartridge
  data (deduplicated mols/reactions and fingerprints). `from_proto` writes only `ord.*` and the `public` rows;
  `update_derived_tables` then computes the `derived` SMILES by reconstructing each message from the search index, and
  the RDKit pass populates/links the cartridge from those SMILES. The two layers store at different granularities on
  purpose: `derived.*_smiles` is one row per reaction/compound (the SMILES string is cheap, and the per-id row is the
  direct lookup into its RDKit object), while `rdkit.*` is deduplicated by SMILES because the parsed structures and
  fingerprints are expensive to recompute and store — i.e. deduplicate where the payload is expensive, store per-entity
  where it is cheap.
* Partial indexes over not-yet-linked rows (`reaction_smiles_unlinked_index`, `compound_smiles_unlinked_index`, and
  `product_compound_smiles_unlinked_index`, on the `derived` SMILES tables) keep incremental RDKit linking fast. They
  index only the rows whose `rdkit_*_id` is still `NULL`, so the repeated "which of this dataset's rows still need
  linking?" queries stay proportional to the in-flight backlog rather than the full table. `prepare_database` creates
  them on new databases; see [Adding the partial indexes to an existing database](#adding-the-partial-indexes-to-an-existing-database)
  to backfill them on an older one.

## Usage

### Install PostgreSQL and the RDKit extension

We use PostgreSQL with the RDKit extension to support structure searches. There are a couple convenient installation
methods (this list is not an endorsement of any particular provider):

* Conda: [rdkit-postgresql](
  https://www.rdkit.org/docs/Install.html#installing-and-using-postgresql-and-the-rdkit-postgresql-cartridge-from-a-conda-environment)

```shell
# Create a new conda environment.
conda create -n ord -c conda-forge python=3.11 rdkit rdkit-postgresql
conda activate ord
# Install ord-schema in this environment.
cd ord-schema
pip install .
# Initialize postgres and start the server.
export PGDATA="${HOME}/rdkit-postgresql"
initdb -U <username>
pg_ctl -l "${HOME}/rdkit-postgresql-logfile" start
# Create a database.
psql -U <username> postgres -c 'CREATE DATABASE <database>;'
```

* AWS: [Amazon Aurora PostgreSQL](
  https://aws.amazon.com/about-aws/whats-new/2020/09/amazon-aurora-postgresql-supports-rdkit-extension/)

### Initialize the database and tables

Use the `prepare_database` function to create the database schemas and tables:

```python
from sqlalchemy import create_engine

from ord_schema.orm.database import prepare_database

connection_string = f"postgresql+psycopg://{username}:{password}@{host}:{port}/{database}"
engine = create_engine(connection_string)
prepare_database(engine)
```

Note that rdkit extension data such as fingerprints are stored in a separate `rdkit` schema in the database to avoid
conflicts with ORD message names.

### Adding the partial indexes to an existing database

`prepare_database` creates the partial indexes that speed up incremental RDKit linking, so newly initialized databases
already have them. A database created before these indexes were introduced should be backfilled. Use
`CREATE INDEX CONCURRENTLY` so the build does not hold a write lock on these large tables (run each statement outside a
transaction block; `IF NOT EXISTS` makes the migration idempotent):

```sql
CREATE INDEX CONCURRENTLY IF NOT EXISTS reaction_smiles_unlinked_index
    ON derived.reaction_smiles (reaction_smiles) WHERE rdkit_reaction_id IS NULL;
CREATE INDEX CONCURRENTLY IF NOT EXISTS compound_smiles_unlinked_index
    ON derived.compound_smiles (smiles) WHERE rdkit_mol_id IS NULL;
CREATE INDEX CONCURRENTLY IF NOT EXISTS product_compound_smiles_unlinked_index
    ON derived.product_compound_smiles (smiles) WHERE rdkit_mol_id IS NULL;
```

These definitions must stay in sync with the `Index(...)` declarations in `derived_mappers.py`.

### Add data

Load ORD datasets into the database with the `add_parquet_dataset` function:

```python
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from ord_schema import huggingface, message_helpers
from ord_schema.orm.database import add_parquet_dataset

# huggingface.fetch_dataset downloads the dataset (preferring the Parquet
# serialization) and returns its local path; it needs the optional `huggingface`
# extra (`pip install "ord-schema[huggingface]"`). add_parquet_dataset streams
# that Parquet file into the database without loading the whole Dataset into
# memory. Datasets not yet migrated to Parquet come back as .pb.gz; load those
# instead with add_dataset(message_helpers.load_message(path, dataset_pb2.Dataset)).
path = huggingface.fetch_dataset("ord_dataset-fc83743b978f4deea7d6856deacbfe53")

connection_string = f"postgresql+psycopg://{username}:{password}@{host}:{port}/{database}"
engine = create_engine(connection_string)
with Session(engine) as session:
    with session.begin():
        add_parquet_dataset(path, session)
```

To load multiple datasets from disk (e.g., from a clone of
[ord-data](https://github.com/open-reaction-database/ord-data)), use the `add_datasets.py` script:

```shell
python scripts/add_datasets.py \
    --pattern "path/to/ord-data/data/fc/*.pb.gz" \
    --username <username> \
    --host <host> \
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

from ord_schema.orm.mappers import Mappers
from ord_schema.proto import reaction_pb2

connection_string = f"postgresql+psycopg://{username}:{password}@{host}:{port}/{database}"
engine = create_engine(connection_string)
with Session(engine) as session:
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionOutcome)
        .join(Mappers.ProductCompound)
        .join(Mappers.ProductMeasurement)
        .join(Mappers.Percentage)
        .where(Mappers.ProductMeasurement.type == "YIELD", Mappers.Percentage.value >= 70)
    )
    results = session.execute(query)
    reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]
assert len(reactions) == 12
  ```

#### Structure searches with the RDKit PostgreSQL extension

##### Reactions whose Morgan binary fingerprint is similar (Tanimoto ≥ 0.5 by default) to `c1ccccc1CCC(O)C`

  ```python
from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from ord_schema.orm.mappers import Mappers
from ord_schema.orm.rdkit_mappers import FingerprintType, RDKitMols
from ord_schema.proto import reaction_pb2

connection_string = f"postgresql+psycopg://{username}:{password}@{host}:{port}/{database}"
engine = create_engine(connection_string)
with Session(engine) as session:
    query = (
        select(Mappers.Reaction)
        .join(Mappers.ReactionInput)
        .join(Mappers.Compound)
        .join(RDKitMols)
        .where(RDKitMols.is_similar("c1ccccc1CCC(O)C", FingerprintType.MORGAN_BFP))
    )
    results = session.execute(query)
    reactions = [reaction_pb2.Reaction.FromString(result[0].proto) for result in results]
assert len(reactions) == 20
  ```

`is_similar` uses the GiST fingerprint index; the cutoff is read from a session
setting (default 0.5) rather than passed to the method. `MORGAN_BFP` uses Tanimoto
similarity (`rdkit.tanimoto_threshold`); `MORGAN_SFP` uses Dice similarity
(`rdkit.dice_threshold`). Set a different cutoff for the rest of the session with:

```python
session.execute(select(func.set_config("rdkit.tanimoto_threshold", "0.7", False)))
```

The third argument to `set_config` is `is_local`. Passing `True` instead scopes the
change to the current transaction only (equivalent to `SET LOCAL`), so it reverts
automatically afterward — useful for a one-off query without disturbing the rest of
the session:

```python
with session.begin():
    session.execute(select(func.set_config("rdkit.tanimoto_threshold", "0.7", True)))
    results = session.execute(query)
```

Use `RDKitMols.tanimoto(...)` instead when you need the similarity *value* (e.g. to
select or order by it).
