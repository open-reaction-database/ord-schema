################################
About the Open Reaction Database
################################

We are building an open access chemical reaction database to support machine
learning and related efforts in reaction prediction, chemical synthesis
planning, and experiment design. Our initial meeting took place on 31 October
2019 and included experts from pharma, academia, and tech.

We expect that this database will be the starting point for the development of
best-in-class tools and models for reaction prediction and synthesis planning.
Additionally, we hope it will serve as a basis for experimental efforts in
industry and academia (e.g., to reduce duplication or focus data generation on
underrepresented areas).

The ORD is described in `JACS <https://pubs.acs.org/doi/10.1021/jacs.1c09820>`__
and was highlighted in
`C&EN <https://cen.acs.org/physical-chemistry/computational-chemistry/new-database-machine-learning-research/99/web/2021/11>`__.

*****
Goals
*****

Our overarching goal, stated above, is to "support machine learning and related
efforts in reaction prediction, chemical synthesis planning, and experiment
design".

We aim to:

* Provide a structured data format for chemical reaction data
* Provide an interface for easy browsing and downloading of data
* Make reaction data freely and publicly available for anyone to use
* Encourage sharing of precompetitive proprietary data, especially HTE data

We conducted a
`survey <https://docs.google.com/spreadsheets/d/1waPzYvDKlb6TAwgsM7bLc7dhZnJ8G-WtVxJSlMhiVK0/edit#gid=585233854>`_
in late 2019/early 2020 to help define the scope and use cases. With the help of
172 respondents, 93.5% of whom reported having a chemistry background and 20.8%
a computer science background, we have established our focus on single-step
organic reactions.

We aim to accommodate data relevant to medicinal chemistry, process chemistry,
flow chemistry, photochemistry, and electrochemistry. Time-course data,
unstructured analytical data, and metadata about how the reaction was performed
will also be accepted. Reactive molecular dynamics simulations, gas-phase
reaction kinetics, and electronic structure calculations for molecular
featurization are being left for other initiatives.

Additionally, we want to be clear that since the database is in active
development, some features not considered "must-haves" for the initial
deployment may not be supported immediately. We welcome feedback about features
that make or break your excitement in and use of the database.

*********
Non-goals
*********

Our non-goals for this initiative, at least right now, are to:

* Capture structured data in a manner designed for programmatic execution
  on automated synthesis hardware (i.e. treating reactions as action
  sequences)
* Store processed analytical data (e.g. NMR peak assignments) other than
  summary statistics (e.g. conversion, yield, purity, selectivity); note
  that unprocessed data (e.g. an exported LCMS file) will be stored
* Integrate model building or external datasets as part of the database

There are some practical consequences of these non-goals:

* Because we are structuring reactions as single-step reactions, there are
  some  more complex operations that cannot be captured in a structured
  format by the schema. For example,
  `this OrgSyn example <http://orgsyn.org/demo.aspx?prep=v95p0080>`_ describes a
  procedure whereby three components are mixed in a chilled vessel on an ice
  bath, a fourth component is added and stirred for several minutes still on
  ice, and then the vessel is removed from the ice bath and allowed to warm
  up to room temperature over several hours. This temperature ramp will be
  captured in a free text field, but the structured temperature field must
  record either zero Celsius or room temperature.
* Analytical data will not be in a unified format and would be difficult to
  learn from directly. This might represent a missed opportunity to train,
  e.g. structural elucidation models that predict a molecular structure
  on the basis of spectral data. Users will still be encouraged to record
  processed data in a text format. Unprocessed data will be stored in its
  original file format as exported by analytical instruments.
* Automated synthesis efforts will require additional future work or a
  separate database, which could lead to fragmentation in the community.
  However, we expect that the reaction procedure captured in the ORD will be
  able to be converted to action sequences through simple translation
  scripts.

**********
Leadership
**********

Governing Committee
===================

The ORD is governed by a Governing Committee with representatives from many
industrial and academic institutions. This committee reviews all aspects of the
database, from the underlying structured data representation to the public
interface(s) to promotion and publicity. The current membership of the governing
committee is:

* Connor Coley (MIT, `C-CAS <https://ccas.nd.edu/>`__)
* Abby Doyle (UCLA, `C-CAS <https://ccas.nd.edu/>`__)
* Spencer Dreher (Merck)
* Klavs Jensen (MIT)
* Steven Kearnes (Genesis Molecular AI)

Advisory Board
==============

We are forming an Advisory Board to include representatives from many
institutions and industry segments. The primary role of the Advisory Board is to
encourage community engagement with the database. The current membership of the
Advisory Board is:

* Alán Aspuru-Guzik (Toronto, MADNESS)
* Timothy Cernak (Michigan, Entos)
* Lucy Colwell (Cambridge, SynTech, Google)
* Werngard Czechtizky (AstraZeneca)
* JW Feng (Consultant)
* Matthew Gaunt (Cambridge, SynTech)
* Alex Godfrey (Consultant)
* Mimi Hii (A*STAR, Imperial, ROAR)
* Greg Landrum (ETH Zurich)
* Fabio Lima (Novartis)
* Christos Nicolaou (Novo Nordisk)
* Sarah Reisman (Caltech, `C-CAS <https://ccas.nd.edu/>`__)
* Francesco Rianjongdee
* Marwin Segler (Microsoft)
* Matthew Sigman (Utah, `C-CAS <https://ccas.nd.edu/>`__)
* Jay Stevens (BMS)
* Sarah Trice
* Huimin Zhao (UIUC, MMLI)

.. _commitment-to-open-access:

*************************
Commitment to Open Access
*************************

As the name of the initiative suggests, this will be an open database in every
sense of the word. All data and code associated with the database will be made
publicly available under commonly used licenses that protect open access.

When proprietary tools are used (see `Interfaces and tools <https://docs.google.com/document/d/1snHPGzKMx19IFq4cj7_OMbvhk4WyYHWENxRhj6FxQrQ/edit#heading=h.46xos12p8y6a>`_),
they will only be used to provide "extra" functionality that is not part of the
core data or the code responsible for data validation and processing. This extra
functionality will be made publicly accessible on the web for anyone to access.

The database is purposely designed to avoid the control or influence of a single
institution. This also ensures that the core data and functionality of the
database will not be affected by any contributor choosing to cease their
involvement in the initiative.

All data submitted to the database will be made available under the
`CC-BY-SA <https://creativecommons.org/licenses/by-sa/4.0/>`_ license, a
well-known license for creative works. Additionally, the various software tools
developed for the database will be made available under the
`Apache <https://choosealicense.com/licenses/apache-2.0/>`_ license; this is
another well-known and OSI-approved license that is used by many organizations
around the world.

***********
How to help
***********

Additional technical help will also be required to aid in tasks such as
processing submissions and implementing the various tools and interfaces to
improve the user experience. We welcome any donation of time to improve using
the schema, technical infrastructure, front-end work, etc.; or simply to provide
feedback on the user experience.

If you are interested in receiving updates or participating in future meetings,
please request to join the
`open-reaction-database <https://groups.google.com/forum/#!forum/open-reaction-database>`__
mailing list.

Past issues of the ORD newsletter are on
`LinkedIn <https://www.linkedin.com/newsletters/ord-newsletter-7256742928012689409/>`__,
where you can also subscribe.

*********************************************
How is the Open Reaction Database being used?
*********************************************

We have a `public Zotero Group <https://www.zotero.org/groups/5973704/open-reaction-database>`_ and `library <https://www.zotero.org/groups/5973704/open-reaction-database/library>`_ where we keep track of Open Reaction Database citations, and also the publications which have reactions in the database.

.. NOTE::
   Have we missed something? Please email help@open-reaction-database.org if
   you find published work which employed the Open Reaction Database data,
   schema, or other tools, and is not already included in the library.

.. NOTE::
   The Datasets sub-collection does not include the individual USPTO patents
   (there is over 13k of them) for clarity.

**************************
Frequently Asked Questions
**************************

**How will my institution benefit from giving others access to our data?**

The ORD is an open science initiative. Your institution will be able to access
the data shared by all other contributors as well. We view this sharing as no
different from how we currently share data when we publish scientific articles,
except the data will be available in an accessible, standardized format.

**What if <institution> decides to stop participating in the initiative?**

The ORD is not reliant on the participation of any individual organization.
However, hosting interfaces and tools to enhance its functionality requires
computational resources; if an organization that provides them stops
participating and the Governing Committee believes these tools to be important
to the adoption and usability of the ORD, the committee will need to seek
funding for web hosting.

**Will my institution have a right to use any tool developed on our data?**

No. The data itself will have a share-alike clause so that derivative work must
be shared under identical terms. Anyone is free to develop software tools,
algorithms, trained machine learning models, etc. under a different license as
long as they follow the terms for the underlying data.

**Can I simply fork your repositories and create my own version of the database?**

Yes! As long as you conform to the terms of the licenses.
