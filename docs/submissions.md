# Submission Workflow

The submission workflow follows four major steps:

## Create

A submission consists of one or more `Dataset` protcol buffers containing
user-defined `Reaction` messages. Datasets can be created
[programatically](schema.html#jupyter-colab) or interactively with the ORD web
editor.

### Validate

The [ord_schema](https://github.com/Open-Reaction-Database/ord-schema)
repository contains scripts for validating submitted data. Prior to creating
your submission, you should run the
[validate_dataset.py](https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/scripts/validate_dataset.py)
script to identify any validation errors:

```shell
$ python validate_dataset.py --input="example_dataset.pbtxt"
```

When defining reactions and datasets programatically, it is good practice to use
the [validation
methods](https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/validations.py)
in `ord-schema` as part of your workflow. Additionally, the ORD web editor
performs automatic validation that will catch any errors entered into the form.

## Submit

```eval_rst
.. IMPORTANT::
   Be sure to follow the instructions for `Getting the Data with GitHub
   <overview.html#id1>`_ before starting your submission.
```

After creating a fork of the database repository (see the note above) you need
to create a pull request for the official repository on GitHub. First, create a
new branch in your fork and add the new dataset:

```shell
$ cd "${REPOSITORY}"
# Make sure your clone is up to date.
$ git checkout main
$ git pull --rebase upstream main
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
[create a pull
request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork)
from your fork.

## Review

Basic preprocessing and validation of your submission will be performed by
automated scripts, and you will be asked to verify any changes performed by the
automated workflow. During this phase of the review process, each `Reaction` and
`Dataset` message will receive a unique database identifier.

After all validation checks have passed, your submission will undergo a manual
review by one or more of the database administrators. The reviewers may suggest
additional changes and continue to iterate with you until they are satisfied
with the submission.

## Deposit

Once the pull request receives approval from the reviewers and passes all
automated checks, a reviewer will merge it into the official database
repository.
