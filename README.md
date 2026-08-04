# Open Reaction Database: Schema (ord-schema)

[![DOI:10.1007/978-3-319-76207-4_15](https://zenodo.org/badge/DOI/10.1021/jacs.1c09820.svg)](https://doi.org/10.1021/jacs.1c09820)
[![PyPI version](https://badge.fury.io/py/ord-schema.svg)](https://badge.fury.io/py/ord-schema)

This repository contains the schema for the Open Reaction Database initiative; please see the documentation at <https://docs.open-reaction-database.org>.

This repository does not contain the database itself; that is stored in [ord-data](https://github.com/open-reaction-database/ord-data). Rather, `ord-schema` is designed to store the database schema and tools for creating, validating, and submitting data to the database.

## Installation

```shell
$ pip install ord-schema
```

This installs the core schema and helpers (building, parsing, validation, and message/Parquet I/O). Heavier, single-purpose features live behind optional extras so the default install stays lightweight:

| Extra | Enables | Install |
| ------- | --------- | --------- |
| `huggingface` | `ord_schema.huggingface.fetch_dataset`: download datasets from the Hugging Face [ord-data](https://huggingface.co/datasets/open-reaction-database/ord-data) mirror | `pip install "ord-schema[huggingface]"` |
| `orm` | `ord_schema.orm`: map the schema into a relational (SQLAlchemy + PostgreSQL) database | `pip install "ord-schema[orm]"` |
| `examples` | running the notebooks under `examples/` (see below) | `pip install "ord-schema[examples]"` |

Extras combine, e.g. `pip install "ord-schema[orm,huggingface]"`. Importing a feature without its extra installed raises a normal `ImportError`.

## Quick start

### Build a reaction

A `Reaction` is a protocol buffer message. Build it field by field, using `message_helpers.build_compound` for the parts with a lot of boilerplate:

```python
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

reaction = reaction_pb2.Reaction()
reaction.identifiers.add(type="REACTION_SMILES", value="CC(=O)O.OCC>>CC(=O)OCC.O")

reaction.inputs["acid"].components.add().CopyFrom(
    message_helpers.build_compound(
        smiles="CC(=O)O",
        name="acetic acid",
        amount="10 mmol",
        role="reactant",
        is_limiting=True,
    )
)
reaction.inputs["alcohol"].components.add().CopyFrom(
    message_helpers.build_compound(
        smiles="OCC", name="ethanol", amount="12 mmol", role="reactant"
    )
)

outcome = reaction.outcomes.add()
outcome.reaction_time.CopyFrom(reaction_pb2.Time(value=3, units="HOUR"))
product = outcome.products.add()
product.identifiers.add(type="SMILES", value="CC(=O)OCC")
product.measurements.add(type="YIELD", percentage=reaction_pb2.Percentage(value=87))

reaction.provenance.record_created.time.value = "2026-01-15"
reaction.provenance.record_created.person.name = "Marie Curie"
reaction.provenance.record_created.person.email = "curie@example.edu"
```

### Validate

Validation reports errors and warnings separately. Pass `raise_on_error=False` to inspect them instead of raising:

```python
from ord_schema import validations

output = validations.validate_message(reaction, raise_on_error=False)
print(output.errors)    # [] -- the reaction above is valid
print(output.warnings)  # advisory only; submissions are not blocked on these
```

RDKit writes parse diagnostics straight to stderr. `ord_schema.logging.silence_rdkit_logs()` quiets them.

### Assemble and write a dataset

`updates.update_dataset` assigns the canonical `ord_dataset-*` and `ord-*` IDs and rewrites cross-references between reactions:

```python
from ord_schema import updates
from ord_schema.proto import dataset_pb2

dataset = dataset_pb2.Dataset(
    name="Esterifications", description="Fischer esterification screen"
)
dataset.reactions.append(reaction)
updates.update_dataset(dataset)
```

`datasets.save_dataset` dispatches on the filename suffix:

```python
from ord_schema import datasets

datasets.save_dataset(dataset, "esterifications.parquet")
```

| Suffix | Format |
| --- | --- |
| `.parquet` | Parquet, one row per reaction; streamable and the storage format for ord-data |
| `.pb` / `.binpb` | binary protocol buffer |
| `.pbtxt` / `.txtpb` | text protocol buffer |
| `.json` | protobuf JSON |

Any of these may be gzipped (`.pb.gz`). For a dataset too large to hold in memory, write it a reaction at a time instead — `DatasetWriter` keeps peak memory to one row group and publishes atomically, so an interrupted write leaves no file behind:

```python
from ord_schema import parquet

with parquet.DatasetWriter(
    "esterifications.parquet",
    name="Esterifications",
    description="Fischer esterification screen",
) as writer:
    for reaction in produce_reactions():
        writer.write(reaction)
```

### Read a dataset

`DatasetView` is the read API for Parquet. Opening one reads only the file footer, so it is cheap regardless of dataset size:

```python
from ord_schema import parquet

view = parquet.DatasetView("esterifications.parquet")

view.name             # "Esterifications"
view.dataset_id       # "ord_dataset-..."
len(view.reactions)   # row count, from the footer -- no reactions deserialized
view.reactions[0]     # reads only the row group holding index 0
```

`reactions` behaves like a list — iteration, `len`, truthiness, indexing, and slicing — while reading from the file on demand:

```python
for reaction in view.reactions:  # streams; peak memory is one row group
    print(message_helpers.get_reaction_smiles(reaction))
```

Look a reaction up by ID, or stream IDs without deserializing anything:

```python
view.get_reaction("ord-1f6c...")   # builds an ID index on first call, then O(1)
list(view.iter_reaction_ids())     # reads one column; cheap on huge files
```

Use `iter_reactions()` when you want the ID alongside the message, and `row_group=` to fan out — row groups are the unit of parallelism:

```python
for reaction_id, reaction in view.iter_reactions():
    ...

for i in range(view.num_row_groups):
    executor.submit(process_row_group, "esterifications.parquet", i)
```

When you genuinely need a `Dataset` message — to serialize, convert to JSON, or mutate — materialize one. This deserializes everything, so peak memory scales with the whole dataset:

```python
dataset = view.to_proto()
```

Non-Parquet formats load whole; `datasets.load_dataset` dispatches on suffix the same way `save_dataset` does:

```python
dataset = datasets.load_dataset("dataset.pb.gz")
```

### Fetch a published dataset

With the `huggingface` extra, pull a dataset from the mirror instead of constructing one. `fetch_dataset` downloads the file and returns its path without parsing it; it prefers Parquet and falls back to `.pb.gz` for datasets not yet converted, so dispatch on the suffix:

```python
import pathlib

from ord_schema import datasets, huggingface, parquet

path = huggingface.fetch_dataset("ord_dataset-...")
if pathlib.Path(path).suffix == ".parquet":
    view = parquet.DatasetView(path)
else:
    dataset = datasets.load_dataset(path)
```

## Notebook examples

The `examples/` directory contains worked examples of dataset creation and use, drawn from published papers. To run locally:

```shell
$ pip install "ord-schema[examples]"
```

Click here to run the examples with Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/open-reaction-database/ord-schema/HEAD?labpath=examples)

## Development

To install in editable/development mode (recommended: [uv](https://docs.astral.sh/uv/)):

```shell
$ git clone https://github.com/open-reaction-database/ord-schema.git
$ cd ord-schema
$ uv sync --extra tests
```

The `tests` extra pulls in the feature extras (`huggingface`, `orm`) it needs to exercise their code paths, so this is enough to run the full suite. Add `--extra examples` as well to run the notebooks (heavier deps):

```shell
$ uv sync --extra examples --extra tests
```

You can still use pip if you prefer: `pip install -e ".[tests]"`.

If you make changes to the protocol buffer definitions, [install](https://grpc.io/docs/protoc-installation/) `protoc` and run `./compile_proto_wrappers.sh` to rebuild the wrappers.

## Conventions

### 1. convention: compound stoichiometry

#### Created: 2023.07.04

#### Last updated: 2023.07.04

#### Description

1. The preferred field for compound stoichiometry is the map `Compound.features` or `ProductCompound.features`.
2. The key should be "stoichiometric_coefficient" or "stoichiometric_ratio".
3. The value should be a `Data` message with its `float_value` representing the compound's stoichiometric
coefficient or ratio.

#### Related links

[#683](https://github.com/open-reaction-database/ord-schema/issues/683)
[#684](https://github.com/open-reaction-database/ord-schema/pull/684)
