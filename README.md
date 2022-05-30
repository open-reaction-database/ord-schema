# Open Reaction Database: Schema (ord-schema)

![CI](https://github.com/Open-Reaction-Database/ord-schema/workflows/CI/badge.svg?branch=main)
[![codecov](https://codecov.io/gh/Open-Reaction-Database/ord-schema/branch/main/graph/badge.svg)](https://codecov.io/gh/Open-Reaction-Database/ord-schema)

This repository contains the schema for the Open Reaction Database initiative; please see the documentation
at http://open-reaction-database.org.

This repository __does not__ contain the database itself; that is stored
in [ord-data](https://github.com/open-reaction-database/ord-data). Rather, `ord-schema` is
designed to store the database schema and tools for creating, validating, and submitting data to the database.

Consider joining the [mailing list](https://groups.google.com/forum/#!members/open-reaction-database) to stay up to date
with announcements and opportunities for providing feedback.

## Examples

The `examples/` directory contains examples of dataset creation and use.

Click here to run the examples with
Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/open-reaction-database/ord-schema/HEAD?filepath=examples)

## Installation

Many of the schema validators require [RDKit](https://github.com/rdkit/rdkit), so we recommend that you
install `ord-schema` in a conda environment.

1. Clone the `ord-schema` repository to your local machine. If you are planning to contribute changes,
   you should fork the repository first
   (see [Contributing](https://github.com/Open-Reaction-Database/ord-schema/blob/main/CONTRIBUTING.md)).
2. Install Anaconda or Miniconda ([instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)).
3. Run the following commands to install `ord-schema`:

   ```shell
   $ pip install -r requirements.txt
   $ conda install -c rdkit rdkit
   $ python setup.py install
   ```

4. After installation, you can run the tests by executing

   ```shell
   $ pip install -r test_requirements.txt
   $ pytest
   ```

   If you are making changes, use `python setup.py build develop` to recompile the protocol buffers and install the
   package in development mode. You should also be sure to run `format.sh` when committing changes.
