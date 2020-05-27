# Submission Workflow

The submission workflow follows four major steps:

## Create

A submission consists of a set of `Reaction` protocol buffers that must be
created by the user. We have examples of how to do this
[programmatically](schema.html#jupyter-colab) in Jupyter/Colab.

After creating a set of `Reaction` messages, they must be grouped together into
one or more `Dataset` messages. This can be done programmatically or with the
provided
[build_dataset.py](https://github.com/Open-Reaction-Database/ord-schema/blob/master/ord_schema/scripts/build_dataset.py)
script. Here's an example:

```shell
$ python build_dataset.py \
  --input="reaction-*.pbtxt" \
  --name="Example dataset" \
  --description="A dataset used for demonstrating the submission workflow" \
  --output="example_dataset.pbtxt"
```

The `Dataset` message provides a logical grouping for the submission. During the
submission review process (see below) the `Reaction` and `Dataset` messages will
each receive a unique identifier. Additionally, build_dataset.py performs
validation of your `Reaction` records and will warn you of any errors before
submitting.

### Validate

The [ord_schema](https://github.com/Open-Reaction-Database/ord-schema)
repository contains scripts for validating submitted data. If you did not use
the build_dataset.py script (which automatically) validates the `Dataset`, you
should run the
[validate_dataset.py](https://github.com/Open-Reaction-Database/ord-schema/blob/master/ord_schema/scripts/validate_dataset.py)
script manually to identify and fix any validation errors prior to moving to the
next stage of the submission workflow:

```shell
$ python validate_dataset.py --input="example_dataset.pbtxt"
```

## Submit

```eval_rst
.. IMPORTANT::
   Be sure to follow the instructions for `Getting the Data with GitHub
   <overview.html#id1>`_ before starting your submission.
```

After creating a fork of the database repository (see the note above) we need to
start a pull request for the official repository on GitHub. Let's continue the
example from the _Create_ step, where we created `example_dataset.pbtxt` to hold
our `Reaction` messages. First, we create a new branch in our fork and add the
new dataset:

```shell
$ cd "${REPOSITORY}"
# Make sure your clone is up to date.
$ git checkout master
$ git pull --rebase upstream master
# Create a new branch for the submission.
$ git checkout -b my_submission
# Copy the new dataset into your fork.
$ cp path/to/example_dataset.pbtxt .
# Create the submission commit(s).
$ git add example_dataset.pbtxt
$ git commit -m "Example dataset submission"
# Push the submission to your fork of the database.
$ git push origin my_submission
```

Next, log in to GitHub, navigate to the [database
repository](https://github.com/Open-Reaction-Database/ord-submissions-test), and
[create a pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork)
from your fork to the official database repository.

## Review

Basic preprocessing for your submission will be performed by automated scripts
(this could be triggered manually or automatically after other validations
pass); this may include canonicalizing structures, adding additional molecular
representations (e.g. SMILES), assigning record IDs, etc. The user will be asked
to verify any changes performed by the automated workflow.

After all validation checks have passed, the submission undergoes a manual
review by one or more of the repository administrators. The reviewer(s) may
suggest additional changes and continue to iterate with the user until they are
satisfied with the submission.

## Deposit  

Once the pull request receives approval from the reviewer(s) and passes all
automated checks, the user (or a reviewer) merges it into the repository.

To avoid capacity issues, the pull request branch should be squashed before
merging and undergo a final automated check for large Data values.
