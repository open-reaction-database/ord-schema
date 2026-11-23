:orphan:

#######################################
Natural language search: privacy notice
#######################################

.. warning::

   **DRAFT — not adopted, not in force.** Written to start a conversation, not
   to end one. Nothing described here is collected until this notice, or
   something it becomes, is published alongside a running natural-language
   search.

****************
What this covers
****************

The Open Reaction Database lets you search in plain English: you describe the
reactions you want, and we translate that into a structured query — the kind
you could build yourself in the structured search — and run it. This notice
describes what we keep when you do that, and why. It covers the natural-language
search only, not ordinary browsing of the site or the structured search.

**********
Who we are
**********

The Open Reaction Database Project is an open-source project, not a legal
entity. The natural-language search is run by its maintainers, *(names to be
filled in)*, who are responsible for the records described here ("we", "us").
You can reach us at help@open-reaction-database.org.

************
What we keep
************

Each time you run a query, we record:

- **Your query**, as you typed it, and when you ran it.
- **Any structure you draw or paste** into it.
- **The structured query it was translated into**, including any attempt that
  failed and the error that rejected it, and the structure each compound name in
  it was matched to.
- **How it ended** — found reactions, matched nothing, refused, failed — and how
  many reactions came back.
- **What it cost**: which model translated, how many tokens it used, and how
  long each step took.
- **Which version of the corpus and of our prompt** were in use, so an old
  record can be read correctly later.
- **A session identifier**, so the queries from one visit can be read together.
- **Your feedback**, if you leave any — a thumbs up or down and any comment with
  it.

We keep the summary we wrote of what your query found. We do not keep the
reactions it was drawn from; the stored structured query and corpus version
reproduce them on demand, and they are public database records in any case.

What a session is
=================

A session is one visit. When you arrive, we create a random identifier that
lives only for as long as you keep the page open, and attach it to the queries
you run while you are there. It is held by the open page alone, not in a cookie
or your browser's storage, not sent to us on any later visit, and not connected
to anything else about you. Reload the page
and you have a new one; there is no way for us — or for you — to link the two.

It exists for one reason: when a first query returns nothing and you rephrase
it, the rephrasing is the most useful thing in the whole record, and it is only
legible if we can see the two queries together.

*******************
What we do not keep
*******************

We do not ask you to sign in, and there is nothing in these records that
identifies you:

- **No account, name, or email address.**
- **No IP address, browser, or device identifier.** The web server's access log,
  which is separate from these records, notes each request to the site with the
  IP address and browser it came from and keeps it for 90 days. It shows that a
  query was run, not what it said, and we do not use it to find out who ran one.
- **No cookie or mark that survives your visit.** Two visits by the same person
  are two unrelated sets of records to us, and we have no way to connect them.

************************
Who else sees your query
************************

Translating your query means sending it to `Anthropic
<https://www.anthropic.com>`__, whose model writes the structured query; that is
how the search works at all. Anthropic processes it on our behalf under its
`Commercial Terms of Service
<https://www.anthropic.com/legal/commercial-terms>`__ and `Data Processing
Addendum <https://www.anthropic.com/legal/data-processing-addendum>`__, and
under those terms may not train its models on it. Once the search has run, the
model sees your query again beside a short summary of the result — how many
reactions or groups came back, and the first few of them — to write the summary
you read. A structure you draw or paste is not sent to Anthropic: the model sees
only a placeholder standing in for it, and the structure itself is used only in
our own search.

When the structured query names a compound, we look up that name with public
chemical name services: `PubChem <https://pubchem.ncbi.nlm.nih.gov>`__, the
NCI/CADD `Chemical Identifier Resolver
<https://cactus.nci.nih.gov/chemical/structure>`__, and `OPSIN
<https://www.ebi.ac.uk/opsin/>`__, in that order until one answers. The model
decides what in your query is a compound name, and it can be wrong: what is
looked up is usually a compound you named, but it can be any other part of your
query, or all of it. The lookup is sent from our server, not your browser.

**Please do not type anything into the search that you would not want recorded,
read by us, or sent to the services above.** A query can be useful to us
without being about anything you need to keep private, and queries here are
about chemistry. If a query would
name a colleague, an unpublished program, a customer, or a compound you are not
ready to disclose, please phrase it a different way.

We do not filter what you type, and you should not rely on us to. The database
itself contains submitters' names, organizations, and email addresses as part of
each record's provenance, and searching for them is ordinary and supported, so a
query naming a person is not by itself out of place here.

Where it is stored
==================

The records are stored with Amazon Web Services in the United States, so a
query from the European Union, the European Economic Area, or the United Kingdom
is transferred outside it. Anthropic's Data Processing Addendum covers its
processing with the European Commission's Standard Contractual Clauses and the
UK International Data Transfer Addendum. PubChem and the Chemical Identifier
Resolver are run by the U.S. National Institutes of Health, and OPSIN by EMBL's
European Bioinformatics Institute in the United Kingdom; each receives only what
the structured query names as a compound.

******************
What we use it for
******************

1. **Running and fixing the search.** Diagnosing a bad result, finding a
   translation the grammar refused, and knowing what it costs to run.
2. **Making it better.** Queries that were translated badly become the basis for
   test cases. This is the main reason the log exists: a translation layer
   improves by being shown where it was wrong, and only real queries show that.
   Our test cases are public, so we write them ourselves rather than copying
   your query: a test case keeps what went wrong, not your wording, and names no
   person, organization, or unpublished compound.
3. **Research.** We may study and publish work on how people search a chemistry
   database in natural language, how they rephrase when a first attempt fails,
   and where the translation layer succeeds and fails. What we publish is
   summaries and aggregates — the kinds of queries people run, how often each
   kind succeeds, how rephrasing changes an outcome — and **never an individual
   query**, the queries from one visit, or a session identifier.

We do not release the queries themselves, as a dataset or in any other form,
under this notice.

*******************
How long we keep it
*******************

Your query, the summary we wrote of what it found, and any comment you left with
your feedback are deleted automatically **one year** after you run it.

We also keep a redacted copy that leaves those three out and keeps everything
else: the outcome, timings, costs, versions, session identifier, your thumbs up
or down, and the structured query with the compound names and structures it
searched for. That copy is what the analysis runs on. It is made from each
month's records before they are deleted, and is itself deleted automatically
**five years** after it is made, so nothing from your query is kept more than
six years after you run it. The structured query
carries whatever compounds, people, or other terms your query named, so the
redacted copy is not anonymous in the way a count is.

***************
Who can read it
***************

The records live in a private store with no public access. Only the maintainers
who administer the project's cloud account can read them, and every read is
recorded in an audit log that is kept longer than the records themselves. We
read them only for the purposes above.

***********
Your rights
***********

Because we hold no account, no address, and no identifier that survives your
visit, we generally **cannot find your records** on our own — which is a
consequence of collecting so little, not a refusal.

You can ask us for a copy of what we hold about you, to correct or delete it, to
restrict how we use it, or to stop using it. Since we cannot tell which records
are yours, we can act on such a request once you tell us enough to find them:
the wording of a query and roughly when you ran it. Write to us at
help@open-reaction-database.org.

Under the EU and UK General Data Protection Regulation, we rely on our
legitimate interests in operating, improving, and studying a public research
resource. The research purpose is carried out under the safeguards the
Regulation sets for scientific research, which is why we publish only summaries
and aggregates. You may object to processing
based on legitimate interests.
You also have the right to complain to a data protection authority: in the EU
or EEA, the one where you live, work, or think the problem arose; in the United
Kingdom, the Information Commissioner's Office.

*******
Changes
*******

We may update this notice. Changes take effect when they are posted here, and a
change that widens what we keep or how we use it applies only to queries run
after it is posted. Every earlier version is in the `history of this page
<https://github.com/open-reaction-database/ord-schema/commits/main/docs/privacy.rst>`__.

| **Effective date:** *(not yet in force)*
| **Contact:** help@open-reaction-database.org
