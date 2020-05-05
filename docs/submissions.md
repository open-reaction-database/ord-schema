# Submission Workflow

The submission workflow follows four major steps:

## Create

A submission consists of a set of Reaction protocol buffers that must be
created by the user. We have examples of how to do this
[programmatically](https://ord-schema.readthedocs.io/en/latest/overview.html#examples-jupyter-colab).

## Submit

The user creates a pull request to the official submission repository on
GitHub. The pull request template contains any required legal
disclaimers, such as an affirmation of the user's right to release the
data under the repository license.

Submissions should be created from a forked version of the repository to
avoid creating many large branches on the official repository.

## Review

The repository runs automated validation scripts and other consistency
checks to validate the submitted data. These scripts will be distributed
with the ord_schema repository so that users are able to pre-validate
their data before submission.  If there are validation errors, the user
fixes them and updates the pull request with the updated records.

Additionally, basic preprocessing will be performed by automated scripts
(this could be triggered manually or automatically after other
validations pass); this may include canonicalizing structures, adding
additional molecular representations (e.g. SMILES), assigning record IDs,
etc. The user will be asked to verify any changes performed by the
automated workflow.

After all validation checks have passed, the submission undergoes a
manual review by one or more of the repository administrators. The
reviewer(s) may suggest additional changes and continue to iterate with
the user until they are satisfied with the submission.

## Deposit  

Once the pull request receives approval from the reviewer(s) and passes
all automated checks, the user (or a reviewer) merges it into the
repository.

To avoid capacity issues, the pull request branch should be squashed
before merging and undergo a final automated check for large Data values.
