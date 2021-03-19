# Open Reaction Database: Schema (ord-schema)

![CI](https://github.com/Open-Reaction-Database/ord-schema/workflows/CI/badge.svg?branch=main)
[![codecov](https://codecov.io/gh/Open-Reaction-Database/ord-schema/branch/main/graph/badge.svg)](https://codecov.io/gh/Open-Reaction-Database/ord-schema)

This repository contains the schema for the Open Reaction Database initiative, based in part on survey of importance ([summary](https://docs.google.com/spreadsheets/d/1waPzYvDKlb6TAwgsM7bLc7dhZnJ8G-WtVxJSlMhiVK0/edit)). For more details on the database, please see the documentation at http://open-reaction-database.org.

This repository __does not__ contain the database itself. Rather, `ord-schema` is designed to store the database schema and tools for creating, validating, and submitting data to the database.

We are currently in the process of formalizing the schema and creating examples for both low-level and high-level use. Consider joining the [mailing list](https://groups.google.com/forum/#!members/open-reaction-database) to stay up to date with announcements and opportunities for providing feedback.

## Examples

Click here to launch a JupyterHub instance: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/open-reaction-database/ord-schema/HEAD?filepath=examples)

## Local installation

Many of the schema validators require [RDKit](https://github.com/rdkit/rdkit), so we recommend that you 
install `ord-schema` in a conda environment.

1. Clone the `ord-schema` repository to your local machine. If you are planning to contribute changes,
   you should fork the repository first
   (see [Contributing](https://github.com/Open-Reaction-Database/ord-schema/blob/main/CONTRIBUTING.md)).
1. Install Anaconda or Miniconda ([instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)).
1. Run the following commands to install `ord-schema`:

   ```shell
   $ pip install -r requirements.txt
   $ conda install -c rdkit rdkit
   $ python setup.py install
   ```

1. After installation, you can run the tests by executing

   ```shell
   $ ./run_tests.sh
   ```

   Note that if you have made changes to the code or `.proto` files, you should run `python setup.py install` to pick
   up the new changes. If you are only editing Python code (or you are willing to run `protoc` yourself), you may wish
   to use `python setup.py develop` to avoid having to reinstall the package each time.
