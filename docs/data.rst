############
Get the data
############

Anyone can download their own copy of the data and associated code (with or
without a GitHub account). We expect and encourage researchers to download
copies of the entire repository. Snapshots of the repository will be backed up
to `Figshare <https://figshare.com/>`_ at regular intervals. All data is
available under the CC-BY-SA license; see
:ref:`commitment-to-open-access`.

To look through reactions before downloading anything, browse and search them at
`open-reaction-database.org <https://open-reaction-database.org/>`__.

Although we have defined the schema using Protocol Buffers, each reaction can be
defined in a human-readable JSON or pbtxt format. The GitHub repository
containing the official version of the database may use the proto binary format
for storage efficiency and speed. Archived snapshots of the repository will
convert examples to a human-readable format so that the data are more
immediately accessible.

******
GitHub
******

The official Open Reaction Database repository is
`ord-data <https://github.com/Open-Reaction-Database/ord-data>`__ on GitHub.

.. IMPORTANT::
   We use `Git LFS <https://git-lfs.github.com/>`_ to efficiently store the
   Dataset records that make up the ORD. To view these files locally, you'll
   need to install Git LFS before cloning the repository.

You may want to use ``--depth`` to create a
`shallow clone <https://git-scm.com/docs/git-clone#Documentation/git-clone.txt---depthltdepthgt>`_
instead of fetching the entire history of the repository::

   $ git clone --depth=1 https://github.com/Open-Reaction-Database/ord-data.git

If you plan to submit data, clone a fork of the repository instead; see
:doc:`submissions`.

************
Hugging Face
************

A copy of the Open Reaction Database repository is also kept on Hugging Face as
`open-reaction-database/ord-data <https://huggingface.co/datasets/open-reaction-database/ord-data>`__.
It can be cloned with Git LFS::

   $ git clone https://huggingface.co/datasets/open-reaction-database/ord-data

Data submissions are not accepted via Hugging Face, so please use GitHub for all
submissions.
