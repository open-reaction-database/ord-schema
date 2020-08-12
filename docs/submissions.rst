###################
Submission Workflow
###################

The submission workflow follows four major steps:

******
Create
******

A submission consists of one or more ``Dataset`` protocol buffers containing
user-defined ``Reaction`` messages. Datasets can be created programatically or
interactively with the ORD web editor.

.. tabs::

   .. group-tab:: Web

        Please use the `ORD web editor <editor.open-reaction-database.org>`_ to
        create your submission.

        The ORD web editor performs automatic validation that will catch any
        errors entered into the form, so there is no separate validation step.

   .. group-tab:: Code

        See the Python examples `here <schema.html#jupyter-colab>`_.

        If you create your submission programmatically, be sure to run the
        `validate_dataset.py
        <https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/scripts/validate_dataset.py>`_
        script to identify any validation errors:

        .. code-block:: shell

            $ python validate_dataset.py --input="example_dataset.pbtxt"

        When defining reactions and datasets programatically, it is good
        practice to use the `validation methods
        <https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/validations.py>`_
        in ``ord-schema`` as part of your workflow.

*******
Prepare
*******

Submissions are received as GitHub `pull requests
<https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests>`_
from a fork of the ORD repository. In essence, you are creating a personal copy
of the repository, updating it with your data, and then requesting that your
changes be merged into the main repository.

If you haven't done so already, you will need to `create a fork
<https://help.github.com/en/github/getting-started-with-github/fork-a-repo>`_ of
the `ord-data <https://github.com/Open-Reaction-Database/ord-data>`_ repository
on GitHub.

.. tabs::

   .. group-tab:: Web

      *TODO: Make sure your fork is up to date.*

      `Create a new branch
      <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-and-deleting-branches-within-your-repository>`_
      for your submission.

   .. group-tab:: Code

      `Clone
      <https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository>`_
      your forked repository to your workstation. You may want to use the
      ``--depth`` flag to create a `shallow clone
      <https://git-scm.com/docs/git-clone#Documentation/git-clone.txt---depthltdepthgt>`_
      instead of fetching the entire commit history:

      .. IMPORTANT::

         Be sure to clone your forked repository and *not* the official repo.

      .. code-block:: shell

         $ git clone --depth=1 "https://github.com/${GITHUB_USERNAME}/${REPOSITORY}"
         # Make sure your fork is up to date.
         $ git checkout main
         $ git pull --rebase upstream main
         # Create a new branch for your submission.
         $ git checkout -b my_submission

******
Submit
******

.. tabs::

   .. group-tab:: Web

      `Upload
      <https://docs.github.com/en/github/managing-files-in-a-repository/adding-a-file-to-a-repository>`_
      your dataset(s) into your submission branch on GitHub and commit the
      result.

   .. group-tab:: Code

      .. code-block:: shell

         # Copy your dataset(s) into your submission branch.
         $ cp path/to/example_dataset.pbtxt .
         # Commit your changes.
         $ git add example_dataset.pbtxt
         $ git commit -m "Example dataset submission"
         # Push the submission to your fork.
         $ git push origin my_submission

Next, log in to GitHub, navigate to the `database repository
<https://github.com/Open-Reaction-Database/ord-submissions-test>`_, and `create
a pull request
<https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork>`_
from your fork to the official repository.

******
Review
******

Your submission will be automatically validated and manually reviewed by one of
the ORD reviewers. The reviewers may suggest additional changes and continue to
iterate with you until they are satisfied with the submission. After your pull
request is approved, it will be merged into a new branch in the official
repository; this new branch is staging point for automated preprocessing that is
required before merging into the official database.

After your submission has been accepted, a reviewer will trigger various
automated preprocessing steps, such as renaming the dataset and assigning
reaction and dataset IDs. Once these changes are verified by the reviewer, the
dataset will be merged into the "main" branch and become part of the official
database.
