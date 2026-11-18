############
Get the data
############

Anyone can download their own copy of the data and associated code (with or
without a GitHub account). We expect and encourage researchers to download
copies of the entire repository. Anyone with a GitHub account can submit data or
code to the database. All data is available under the CC-BY-SA license; see
:ref:`commitment-to-open-access`.

To look through reactions before downloading anything, browse and search them at
`open-reaction-database.org <https://open-reaction-database.org/>`__.

Each dataset is one Parquet file in the repository's ``data/`` directory. Every
row holds one serialized ``Reaction`` message, and the dataset's name,
description, and ID travel in the file metadata. The ``ord-schema`` package
reads these files and converts them to human-readable formats such as pbtxt and
JSON; see :ref:`reading-the-data-in-python`.

******
GitHub
******

The official Open Reaction Database repository is
`ord-data <https://github.com/open-reaction-database/ord-data>`__ on GitHub.
Its v0.1.0 release (2022) is archived on
`Zenodo <https://doi.org/10.5281/zenodo.6028405>`__.

.. IMPORTANT::
   The datasets are stored with `Git LFS <https://git-lfs.com/>`_, so install
   it before cloning the repository. The repository fetches LFS objects from the
   Hugging Face mirror on its own; there is nothing else to configure.

You may want to use ``--depth`` to create a
`shallow clone <https://git-scm.com/docs/git-clone#Documentation/git-clone.txt---depthltdepthgt>`_
instead of fetching the entire history of the repository::

   $ git clone --depth=1 https://github.com/open-reaction-database/ord-data.git

If you plan to submit data, clone a fork of the repository instead; see
:doc:`submissions`.

************
Hugging Face
************

A copy of the Open Reaction Database repository is also kept on Hugging Face as
`open-reaction-database/ord-data <https://huggingface.co/datasets/open-reaction-database/ord-data>`__,
where you can browse and download single datasets. It can be cloned with Git
LFS::

   $ git clone https://huggingface.co/datasets/open-reaction-database/ord-data

To download only some of the datasets, or the data without its Git history, use
the script that the
`ord-data README <https://github.com/open-reaction-database/ord-data#getting-the-data>`__
describes.

Data submissions are not accepted via Hugging Face, so please use GitHub for all
submissions.

.. _reading-the-data-in-python:

**************************
Reading the data in Python
**************************

The `ord-schema <https://pypi.org/project/ord-schema/>`__ package reads and writes
every format the ORD uses. Its ``huggingface`` extra downloads single datasets from
the Hugging Face mirror::

   $ pip install "ord-schema[huggingface]"

``huggingface.fetch_dataset`` downloads a dataset and returns its path, and
``datasets.load_dataset`` reads it:

.. code-block:: python

   from ord_schema import datasets, huggingface

   path = huggingface.fetch_dataset("ord_dataset-...")
   dataset = datasets.load_dataset(path)  # a DatasetView for .parquet, else a Dataset

``load_dataset`` reads a file from a clone of the repository the same way. To
convert a dataset to text, load the whole ``Dataset`` and save it under a
``.pbtxt`` or ``.json`` name:

.. code-block:: python

   dataset = datasets.load_dataset(path, as_dataset=True)
   datasets.save_dataset(dataset, "dataset.pbtxt")

The
`examples <https://github.com/open-reaction-database/ord-schema/tree/main/examples>`_
are Jupyter notebooks, drawn from published papers, that build and use datasets
with the full reaction schema.
