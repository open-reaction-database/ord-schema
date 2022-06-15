# Open Reaction Database: Schema (ord-schema)

![CI](https://github.com/Open-Reaction-Database/ord-schema/workflows/CI/badge.svg?branch=main)
[![codecov](https://codecov.io/gh/Open-Reaction-Database/ord-schema/branch/main/graph/badge.svg)](https://codecov.io/gh/Open-Reaction-Database/ord-schema)

This repository contains the schema for the Open Reaction Database initiative; please see the documentation
at https://docs.open-reaction-database.org.

This repository does not contain the database itself; that is stored
in [ord-data](https://github.com/open-reaction-database/ord-data). Rather, `ord-schema` is
designed to store the database schema and tools for creating, validating, and submitting data to the database.

## Installation

If you have __not__ previously installed rdkit via conda:

```shell
$ pip install ord-schema[rdkit]
```

Otherwise:

```shell
$ pip install ord-schema
```

To install in editable/development mode:

```shell
$ git clone https://github.com/open-reaction-database/ord-schema.git
$ cd ord-schema
$ pip install -e .  # Or `pip install -e .[rdkit]`.
```

If you make changes to the protocol buffer definitions, [install](https://grpc.io/docs/protoc-installation/) `protoc` 
and run `./compile_proto_wrappers.sh` to rebuild the wrappers.

## Examples

The `examples/` directory contains examples of dataset creation and use. To run locally, install with:

```shell
$ pip install "ord-schema[examples]"
```

Click here to run the examples with
Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/open-reaction-database/ord-schema/HEAD?filepath=examples)
