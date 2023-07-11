# Open Reaction Database: Schema (ord-schema)

[![DOI:10.1007/978-3-319-76207-4_15](https://zenodo.org/badge/DOI/10.1021/jacs.1c09820.svg)](https://doi.org/10.1021/jacs.1c09820)
[![PyPI version](https://badge.fury.io/py/ord-schema.svg)](https://badge.fury.io/py/ord-schema)

This repository contains the schema for the Open Reaction Database initiative; please see the documentation
at https://docs.open-reaction-database.org.

This repository does not contain the database itself; that is stored
in [ord-data](https://github.com/open-reaction-database/ord-data). Rather, `ord-schema` is
designed to store the database schema and tools for creating, validating, and submitting data to the database.

## Installation

```shell
$ pip install ord-schema
```

## Examples

The `examples/` directory contains examples of dataset creation and use. To run locally, install with:

```shell
$ pip install "ord-schema[examples]"
```

Click here to run the examples with Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/open-reaction-database/ord-schema/HEAD?labpath=examples)

## Development

To install in editable/development mode:

```shell
$ git clone https://github.com/open-reaction-database/ord-schema.git
$ cd ord-schema
$ pip install -e .
```

If you make changes to the protocol buffer definitions, [install](https://grpc.io/docs/protoc-installation/) `protoc`
and run `./compile_proto_wrappers.sh` to rebuild the wrappers.

## Conventions

### 1. convention: compound stoichiometry

##### Created: 2023.07.04

##### Last updated: 2023.07.04

##### Description: 
1. The preferred field for compound stoichiometry is the map `Compound.features` or `ProductCompound.features`.
2. The key should be "stoichiometric_coefficient" or "stoichiometric_ratio".
3. The value should be a `Data` message with its `float_value` representing the compound's stoichiometric 
coefficient or ratio.

##### Related links: 
[#683](https://github.com/open-reaction-database/ord-schema/issues/683) 
[#684](https://github.com/open-reaction-database/ord-schema/pull/684)