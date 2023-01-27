# Contributing to `ord-schema`

Thanks for helping out! What would you like to do?

* [Report an issue](#report-an-issue)
* [Suggest an improvement](#suggest-an-improvement)
* [Make a change](#make-a-change)

## Report an issue

Sorry to hear that you're having trouble!
Please [file a bug report](https://github.com/Open-Reaction-Database/ord-schema/issues/new?assignees=&labels=bug&template=bug_report.md&title=)
to tell us about the problem.

## Suggest an improvement

Great idea!
Please [create a feature request](https://github.com/Open-Reaction-Database/ord-schema/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=)
to let us know what you'd like to see.

## Make a change

Excellent! There are a few steps you'll need to follow to get ready to submit changes for review:

1. Make sure your changes are consistent with
   the [goals](https://ord-schema.readthedocs.io/en/latest/overview.html#goals)
   and [non-goals](https://ord-schema.readthedocs.io/en/latest/overview.html#non-goals) of the database. If your changes
   address a long-term need we may ask that it be deferred until sometime in the future. In particular, please avoid
   changes that significantly increase the complexity of creating or using the data.
1. Make your changes in
   a [fork](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-forks)
   of the `ord-schema` repository. You can work on GitHub or on your local machine:

   ### GitHub

    1. Create a fork of `ord-schema` by clicking the "Fork" button at the top-right of this page.
       Note that you only need to do this once---all of your changes can use the same fork.
    1. Follow
       the [GitHub flow](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/github-flow)
       to make changes and send a pull request to the official repository. While you are making
       changes, automated tests will run in your fork after each commit and notify you of any errors encountered in
       validation or syntax.

   ### Local machine

   > :warning: Working on your local machine requires familiarity with
   > [Git](https://git-scm.com/) or [GitHub Desktop](https://desktop.github.com/).

    1. Create a fork of `ord-schema` (click "Fork" at the top-right of this page) and clone it to your local machine
       ([instructions](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)).
    1. Create a new branch and make your changes.
    1. Test your changes by rebuilding the package (`python setup.py install`) and running the test script provided
       in the repo (`./run_tests.sh`). This will be performed automatically when you send a pull request to the
       repository, but you can save time by making sure the tests pass locally first.
    1. Commit your changes, push them to your fork on GitHub, and send a pull request to the
       official repository.

1. If necessary, update your pull request by adding additional commits to your branch in response to
   reviewer comments and automated tests. Once the reviewers have approved your changes, they will merge
   them into the official repository.

## Terms of Use

By submitting Contributions (as defined below) to this project, you agree that
You (as defined below) and your Contributions are bound to the following terms.

<p align="center">Open Reaction Database Project License Agreement</p> 

In order to clarify the intellectual property license granted with Contributions
to the Open Access Reaction Database from any person or entity, You must
indicate agreement to the license terms below.

You accept and agree to the following terms and conditions for Your present and
future Contributions submitted to ORDP. Except for the license granted herein to
ORDP and recipients of data distributed by ORDP, You reserve all right, title,
and interest in and to Your Contributions.

1. Definitions. "You" (or "Your") shall mean the copyright owner or legal entity
   authorized by the copyright owner that is making this Agreement with ORDP. For
   legal entities, the entity making a Contribution and all other entities that
   control, are controlled by, or are under common control with that entity are
   considered to be a single Contributor. For the purposes of this definition,
   "control" means (i) the power, direct or indirect, to cause the direction or
   management of such entity, whether by contract or otherwise, or (ii) ownership
   of fifty percent (50%) or more of the outstanding shares, or (iii) beneficial
   ownership of such entity. "Contribution" shall mean the data, code,
   documentation or any original work of authorship, including any modifications or
   additions to an existing work, that is intentionally submitted by You to ORDP
   for inclusion in the Open Reaction Database managed by ORDP (the "Work"). For
   the purposes of this definition, "submitted" means any form of electronic or
   written communication sent to ORDP or its representatives, including but not
   limited to communication on electronic mailing lists, but excluding
   communication that is conspicuously marked or otherwise designated in writing by
   You as "Not a Contribution."

2. Grant of Copyright License. Subject to the terms and conditions of this
   Agreement, You hereby grant to ORDP and to recipients of the data, software, or
   documentation distributed by ORDP a perpetual, worldwide, non-exclusive,
   no-charge, royalty-free, irrevocable copyright license to reproduce, prepare
   derivative works of, publicly display, publicly perform, sublicense, and
   distribute Your Contributions and such derivative works.

3. Grant of Patent License. Subject to the terms and conditions of this
   Agreement, You hereby grant to ORDP and to recipients of software distributed by
   ORDP a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable
   (except as stated in this section) patent license to make, have made, use, offer
   to sell, sell, import, and otherwise transfer the Work, where such license
   applies only to those patent claims licensable by You that are necessarily
   infringed by Your Contribution(s) alone or by combination of Your
   Contribution(s) with the Work to which such Contribution(s) was submitted. If
   any entity institutes patent litigation against You or any other entity
   (including a cross-claim or counterclaim in a lawsuit) alleging that your
   Contribution, or the Work to which you have contributed, constitutes direct or
   contributory patent infringement, then any patent licenses granted to that
   entity under this Agreement for that Contribution or Work shall terminate as of
   the date such litigation is filed.

4. You represent that you are legally entitled to grant the above license. If
   your employer(s) has rights to intellectual property that you create that
   includes your Contributions, you represent that you have received permission to
   make Contributions on behalf of that employer, or that your employer has waived
   such rights for your Contributions to ORDP.

5. You are not expected to provide support for Your Contributions, except to the
   extent You desire to provide support. You may provide support for free, for a
   fee, or not at all. Unless required by applicable law or agreed to in writing,
   You provide Your Contributions on an "AS IS" BASIS, WITHOUT WARRANTIES OR
   CONDITIONS OF ANY KIND, either express or implied, including, without
   limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT,
   MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE.
