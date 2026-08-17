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
| `search` | `ord_schema.search`: query derived Parquet artifacts with a validated query grammar (DuckDB) | `pip install "ord-schema[search]"` |
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

`datasets.load_dataset` is the entry point for every format. A Parquet dataset reads back as a `DatasetView`, which takes its scalars and row count from the file footer and reads reactions on demand, so opening one is cheap regardless of dataset size:

```python
from ord_schema import datasets

view = datasets.load_dataset("esterifications.parquet")

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

The other formats have no streaming form, so they read back as a `Dataset` message directly:

```python
dataset = datasets.load_dataset("dataset.pb.gz")  # a Dataset, not a view
```

Code that only reads does not need to care which it got — iteration, `len`, indexing, and slicing work on both. Code that needs the protobuf surface and cannot know the format should ask for a message up front, which materializes a Parquet dataset in full:

```python
dataset = datasets.load_dataset(path, as_dataset=True)  # always a Dataset
```

### Fetch a published dataset

With the `huggingface` extra, pull a dataset from the mirror instead of constructing one. `fetch_dataset` downloads the file and returns its path without parsing it; it prefers Parquet and falls back to `.pb.gz` for datasets not yet converted, which is exactly the dispatch `load_dataset` already handles:

```python
from ord_schema import datasets, huggingface

path = huggingface.fetch_dataset("ord_dataset-...")
dataset = datasets.load_dataset(path)  # a DatasetView for .parquet, else a Dataset
```

## Notebook examples

The `examples/` directory contains worked examples of dataset creation and use, drawn from published papers. To run locally:

```shell
$ pip install "ord-schema[examples]"
```

Click here to run the examples with Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/open-reaction-database/ord-schema/HEAD?labpath=examples)

## Package layout

[`ord_schema/README.md`](ord_schema/README.md) maps the package: which module builds a message, which stores a dataset, and which of the two query paths answers which kind of question. The subsystems document themselves — [`artifacts/`](ord_schema/artifacts/README.md) for the derived Parquet and its staleness stamps, [`search/`](ord_schema/search/README.md) for the query grammar and its compiler, [`orm/`](ord_schema/orm/README.md) for the relational mapping.

## Development

Editable install, with [uv](https://docs.astral.sh/uv/) or with pip (`pip install -e ".[tests]"`):

```shell
$ git clone https://github.com/open-reaction-database/ord-schema.git
$ cd ord-schema
$ uv sync --extra tests
```

The `tests` extra pulls in the feature extras (`huggingface`, `orm`) it needs to exercise their code paths, so this is enough to run the full suite; add `--extra examples` for the notebooks. [CONTRIBUTING.md](CONTRIBUTING.md) covers the rest: pre-commit hooks, running the same checks by hand, rebuilding the proto wrappers, and how a release is cut.

## Conventions

Some things the schema can hold in more than one shape have an agreed shape to prefer, so
that data from different depositors answers the same query. These are guidance rather than
validation rules; see [CONVENTIONS.md](CONVENTIONS.md).
