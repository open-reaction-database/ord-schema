:orphan:

##############################################
Open Reaction Database website: privacy notice
##############################################

.. warning::

   **DRAFT — not adopted, not in force.** Written to start a conversation, not
   to end one. It describes the site with Google Analytics turned off, and
   every passage marked *(to be confirmed: …)* is a decision or a fact that
   still needs an owner.

****************
What this covers
****************

This notice covers the two web applications of the Open Reaction Database:

- **The viewer**, at `open-reaction-database.org
  <https://open-reaction-database.org>`__, where anyone can browse, search, and
  ask questions of the published database without signing in.
- **The reaction editor**, at `app.open-reaction-database.org
  <https://app.open-reaction-database.org>`__, where you sign in to create and
  edit reaction datasets.

It does not cover this documentation site, which is hosted by `Read the Docs
<https://about.readthedocs.com/privacy-policy/>`__ under its own policy, or the
project's repositories on GitHub and its mirror on Hugging Face, which those
services run under theirs.

**********
Who we are
**********

The Open Reaction Database Project is an open-source project, not a legal
entity. The viewer and the editor are run by its maintainers, *(to be
confirmed: names)*, who are responsible for the records described here ("we",
"us"). You can reach us at help@open-reaction-database.org.

***********************
When you use the viewer
***********************

The viewer has no accounts and sets no cookies. It does not run Google
Analytics or any other analytics or advertising tag. *(to be confirmed: true
once the change that turns Google Analytics off is deployed.)*

The web server's log
====================

Like most web servers, ours writes a line to its log for every request your
browser makes — each page, each file the page loads, and each call the page
makes to our search service. A line records:

- **Your IP address.**
- **When** the request was made.
- **What was requested**: the address of the page, file, or search call, which
  shows the dataset or reaction you opened. For the ask page and the search
  page, the log leaves out the question or the search you typed. *(to be
  confirmed: the search page's address is left out once ord-app's pull request
  #928 is deployed; until then, a search written into the address is logged.)*
- **The page you came from**, as your browser reports it, with the same
  question or search left out.
- **Your browser**, as it describes itself (its name, version, and operating
  system).
- **The response**: whether the request succeeded, and how much was sent back.

The log is stored with Amazon Web Services in the United States and deleted
automatically **90 days** after each line is written. We read it to keep the
site running: to find errors, to see when it is slow or overloaded, and to deal
with abuse. *(to be confirmed: whether we also count visits from it.)*

What the viewer keeps in your browser
=====================================

When you select reactions in a list of search results, the viewer remembers
which ones in your browser's session storage, so the selection survives a
reload. It stays in that tab and is gone when you close it; it is never sent
to us.

Asking a question
=================

The ask page, which is in beta, turns a question you write in English into a
search. To do that, it sends your question to `Anthropic
<https://www.anthropic.com>`__, whose model writes the search; Anthropic
processes it on our behalf under its `Commercial Terms of Service
<https://www.anthropic.com/legal/commercial-terms>`__. The model sees your
question again beside a short summary of what the search found, to write the
sentence you read. A structure you draw or paste is not sent to Anthropic: the
model sees only the placeholder, such as ``@A``, that stands for it.

When the search names a compound, our server looks the name up with public
chemical name services: `PubChem <https://pubchem.ncbi.nlm.nih.gov>`__, the
NCI/CADD `Chemical Identifier Resolver
<https://cactus.nci.nih.gov/chemical/structure>`__, and `OPSIN
<https://www.ebi.ac.uk/opsin/>`__ at EMBL's European Bioinformatics Institute.
The model decides what in your question is a compound name, and it can be
wrong, so what is looked up can be another part of your question.

The viewer does not write your question or your structures to a database or to
its log: the log records that a question was asked, not what it said. Please
do not enter anything you would not want sent to the services above.

*(to be confirmed: how this section relates to the separate draft notice for
the natural-language search, which describes a log of questions that the viewer
does not keep today.)*

***********************
When you use the editor
***********************

Signing in
==========

The editor asks you to sign in through `Auth0 <https://auth0.com>`__, a sign-in
service run by Okta, with an ORD App account (an email address and a password,
which Auth0 holds), or with your Google, GitHub, or ORCID account. *(to be
confirmed: that the Auth0 tenant offers exactly these.)*

When you sign in with Google, GitHub, or ORCID, that service tells Auth0 who you
are, and Auth0 tells us. Auth0 keeps a record of your account and of your
sign-ins, under `its own privacy terms <https://www.okta.com/privacy-policy/>`__.
*(to be confirmed: how long Auth0 keeps sign-in records on our plan.)*

What we keep about you
======================

Each time you sign in, the editor stores or refreshes a record of you in its
database, from what Auth0 passes on:

- **Your Auth0 identifier**, which names the service you signed in with.
- **A username**: your GitHub username if you signed in with GitHub, and your
  Auth0 identifier otherwise.
- **Your ORCID iD**, if you signed in with ORCID.
- **Your email address and your name**, if the service you signed in with gives
  them.
- **The address of your profile picture**, if it gives one.
- **When the record was created and last changed.**

The editor also keeps what you make in it: your datasets, reactions, and
templates, and the groups you belong to with your role in each. Each new
account gets a group of its own.

What the editor keeps in your browser
=====================================

To keep you signed in, the editor stores your sign-in tokens in your browser's
local storage and sets a cookie on its own domain marking that you are signed
in. Auth0 sets its own cookies on its sign-in pages. You stay signed in for up
to 30 days, or 15 days if you do not use the editor in between; signing out
ends it sooner. These cookies and tokens exist only to sign you in. The editor
does not run Google Analytics or any other analytics or advertising tag. *(to be
confirmed: true once the change that turns Google Analytics off is deployed.)*

Your name in the reactions you make
===================================

**Read this before you submit a dataset.** When you create a reaction from
scratch, the editor writes your username, name, ORCID iD, and email address into
the reaction's provenance, as the experimenter and as the person who created
the record, with the time. The provenance forms also have a button that fills
in the same details. They become part of the reaction itself:

- Everyone in a group you share the dataset with can see them, and they are in
  every download of the reaction.
- When a dataset is submitted to the Open Reaction Database, they are
  published with it. Submission is a pull request to the `ord-data
  <https://github.com/open-reaction-database/ord-data>`__ repository on GitHub,
  which you make yourself, outside the editor. Once it is accepted, the reaction
  — your name, email address, and ORCID iD included — is public under the
  `CC BY-SA 4.0 license
  <https://github.com/open-reaction-database/ord-data/blob/main/LICENSE>`__,
  shown in the viewer, and copied to our mirror on Hugging Face.

Recording who ran an experiment and who wrote down the record is part of what
makes the database's records usable, so the schema has a place for it. You can
edit or remove these details in a reaction before you submit it. After a
dataset is published, it stays in ord-data's public history even if a later
change removes your details. *(to be confirmed: what we do when someone asks
for their details to be removed from a published dataset.)*

Who else sees your account
==========================

Members of a group you belong to see your name, email address, ORCID iD, and
profile picture in the group's member list, and the editor sends their browsers
your username with them. Their browsers load your profile picture from wherever
the service you signed in with keeps it.

The server's log
================

The editor's server logs each call the editor makes to it, with when it was
made, what was called, and how it ended. These lines do not record your IP
address: they record the address of the load balancer in front of the server.
They are stored with Amazon Web Services in the United States and deleted
automatically **90 days** after each line is written. The web server in front of
it also notes each request, with your browser's description and the page you
came from but again the load balancer's address rather than yours, in a log
that stays inside the running server and is discarded whenever the server is
replaced.

Looking up compounds
====================

When you look up a compound by name in the editor, our server sends the name to
PubChem, the NCI/CADD Chemical Identifier Resolver, and `eMolecules
<https://www.emolecules.com>`__ at once, and uses the first answer.

**************************
Who else sees your browser
**************************

Both the viewer and the editor load their font, Roboto, from `Google Fonts
<https://fonts.google.com>`__. Your browser fetches it from Google's servers
each time it does not already have it, so Google receives your IP address and
your browser's description with each such request. Google describes what it does
with them in its `Google Fonts privacy FAQ
<https://developers.google.com/fonts/faq/privacy>`__. *(to be confirmed:
whether to serve the font from our own server instead.)*

******************
Where it is stored
******************

Our logs and the editor's database are stored with Amazon Web Services in the
United States (the ``us-east-1`` region), and Auth0 holds the editor's accounts
in the United States. A visit from the European Union, the European Economic
Area, or the United Kingdom therefore transfers what is described here outside
it. *(to be confirmed: the safeguards each transfer relies on.)*

*******************
How long we keep it
*******************

- **The viewer's and the editor's logs**: deleted automatically 90 days after
  each line is written.
- **Your editor account and what you make in it**: kept while your account
  exists. The database keeps 30 days of backups, so a deleted record leaves the
  backups within 30 days. *(to be confirmed: whether an account that is not
  used is ever deleted, and how a deletion request is carried out.)*
- **Published datasets**: kept permanently, as the public record of the
  database.

***************
Who can read it
***************

The logs and the editor's database live in the project's Amazon Web Services
account, and the editor's accounts in its Auth0 tenant. *(to be confirmed: who
administers each, and so can read them.)* We read them only for the purposes
above.

***********
Your rights
***********

You can ask us for a copy of what we hold about you, to correct or delete it, to
restrict how we use it, or to stop using it. Write to us at
help@open-reaction-database.org. For an editor account, tell us the email
address or the ORCID iD you sign in with. For the viewer's log, which records
no name or account, tell us your IP address and roughly when you visited.

*(to be confirmed: the lawful basis for each use under the EU and UK General
Data Protection Regulation — for example, our legitimate interests in operating
a public research resource for the logs, and providing the service you signed
up for, for the editor's accounts.)* You also have the right to complain to a
data protection authority: in the EU or EEA, the one where you live, work, or
think the problem arose; in the United Kingdom, the Information Commissioner's
Office.

*******
Changes
*******

We may update this notice. Changes take effect when they are posted here, and
every earlier version is in the `history of this page
<https://github.com/open-reaction-database/ord-schema/commits/main/docs/site_privacy.rst>`__.

| **Effective date:** *(not yet in force)*
| **Contact:** help@open-reaction-database.org
