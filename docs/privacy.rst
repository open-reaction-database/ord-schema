#######################################
Natural language search: privacy notice
#######################################

.. warning::

   **DRAFT — not adopted, not in force.** Written to start a conversation, not
   to end one. Nothing described here is collected until this notice, or
   something it becomes, is published alongside a running natural-language
   search. The maintainer notes at the bottom are not part of the notice and
   should not be published with it.

****************
What this covers
****************

The Open Reaction Database lets you ask a question in English and answers it by
translating your question into a database query. This notice describes what we
keep when you do that, and why. It covers the natural-language search only.
Ordinary browsing of the site, and the structured search where you build a
query yourself, are not covered here.

************
What we keep
************

Each time you ask a question, we record:

- **Your question**, as you typed it.
- **The query it was translated into**, including any attempt that failed and
  the error that rejected it.
- **How it ended** — answered, matched nothing, refused, failed — and how many
  reactions came back.
- **What it cost**: which model translated, how many tokens it used, and how
  long each step took.
- **Which version of the corpus and of our prompt** were in use, so an old
  record can be read correctly later.
- **A session identifier**, so the questions from one visit can be read
  together.
- **Your feedback**, if you leave any — a thumbs up or down and any comment with
  it.

We keep the answer we wrote for you in prose. We do not keep the reactions it
was drawn from; the stored query and corpus version reproduce them on demand,
and they are public database records in any case.

What a session is
=================

A session is one visit. When you arrive, we create a random identifier that
lives only for as long as you keep the page open, and attach it to the questions
you ask while you are there. It is not stored on your device, not sent to us on
any later visit, and not connected to anything else about you. Reload the page
and you have a new one; there is no way for us — or for you — to link the two.

It exists for one reason: when a first question returns nothing and you
rephrase it, the rephrasing is the most useful thing in the whole record, and it
is only legible if we can see the two questions together.

*******************
What we do not keep
*******************

We do not ask you to sign in, and there is nothing in these records that
identifies you:

- **No account, name, or email address.**
- **No IP address, browser, or device identifier.** The web server keeps its own
  access logs, as web servers do; those are separate from this and are not
  joined to it.
- **No cookie or mark that survives your visit.** Two visits by the same person
  are two unrelated sets of records to us, and we have no way to connect them.

**************************************
Your question goes to a model provider
**************************************

Translating your question means sending it to a third-party model provider,
which is how the search works at all. What they do with it is governed by their
terms, which we link to from the search page.

**Please do not type anything into the search that you would not want recorded,
read by us, or published.** A question can be interesting to us — and to other
researchers — without being about anything you need to keep private, and
questions here are about chemistry. If a question would name a colleague, an
unpublished program, a customer, or a compound you are not ready to disclose,
please ask it a different way.

We do not filter what you type, and you should not rely on us to. Note also that
the database itself contains submitters' names, organizations, and email
addresses as part of each record's provenance — searching for them is an
ordinary, supported thing to do, so a question naming a person is not by itself
out of place here.

******************
What we use it for
******************

1. **Running and fixing the search.** Diagnosing a bad answer, finding a
   translation the grammar refused, and knowing what it costs to run.
2. **Making it better.** Questions that were answered badly become test cases.
   This is the main reason the log exists: a translation layer improves by being
   shown where it was wrong, and only real questions show that.
3. **Research and publication.** We may study and publish work on what people
   ask a chemistry database in natural language, how they rephrase when a first
   attempt fails, and where the translation layer succeeds and fails. That work
   may include **quoting individual questions verbatim**.

That third purpose is stated here in advance because it is not honest to collect
data for one reason and publish it under another. Stating it reserves the
possibility; it does not promise that any of it will happen.

How we would publish
====================

- **Aggregates and categories, ordinarily** — the kinds of questions asked, how
  often each kind succeeds, how rephrasing changes an outcome.
- **Verbatim questions only after review.** Anything that names a person, an
  organization, an unpublished compound or program, or that would be
  identifying in combination with what surrounds it, is not quoted. When in
  doubt, it is not quoted.
- **We do not publish session identifiers.** Where a sequence of questions from
  one visit is quoted — a rephrasing is often the interesting thing — the
  sequence is reviewed as a whole, since several questions together say more
  about who asked than any one of them.
- If we ever release a dataset rather than quoting selections in our own work,
  that is a bigger step than this notice covers, and we would say so before
  doing it.

*******************
How long we keep it
*******************

Questions and their translations are held for a limited period and then deleted
automatically. The aggregate measurements derived from them — outcomes, timings,
costs, and the queries themselves — are kept longer, because that is what the
analysis runs on and it says nothing about who asked.

Retention periods are stated on the search page and are set to the shortest
period that serves the purposes above.

***************
Who can read it
***************

Nobody has standing access. The records live in a private store with no public
access, and every read of them is recorded in an audit log that outlives the
records themselves. Access is granted to a named maintainer for a stated reason,
and to the service that writes the records.

***********
Your rights
***********

Because we hold no account, no address, and no identifier that survives your
visit, we generally **cannot find your records** — which is a consequence of
collecting so little, not a refusal. We are not able to identify you from what
we hold, so the rights that depend on identifying you do not attach to these
records. If you can describe a question closely enough that we could locate it,
write to us at the address on the search page and we will consider what we can
do.

Under the GDPR we rely on our legitimate interests in operating and improving a
public research resource, and, for the research purpose, on the provisions for
processing for scientific research purposes. You may object to processing based
on legitimate interests.

*******
Changes
*******

We may update this notice. Changes take effect when they are posted here.

**Effective date:** *(not yet in force)* **Contact:** *(to be filled in)*
