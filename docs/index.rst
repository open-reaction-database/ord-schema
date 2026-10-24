:html_theme.sidebar_secondary.remove:

.. The home page: a hero, three ways in, and news, styled by the ord-* rules in
   _static/css/custom.css. Each entry in the root toctree below is one section in the bar
   under the site bar.

.. rst-class:: ord-home

######################
Open Reaction Database
######################

.. rst-class:: ord-home-lead

An open-access database of chemical reactions, built to support machine learning
in reaction prediction, synthesis planning, and experiment design.

.. rst-class:: ord-home-actions

:doc:`Get the data <data>` :doc:`Contribute a dataset <submissions>`

.. container:: ord-cards

   .. container:: ord-card

      :doc:`Get the data <data>`

      Every dataset in the ORD, from GitHub or Hugging Face.

   .. container:: ord-card

      :doc:`Contribute <submissions>`

      Enter reactions in the editor or build a dataset from a template, then submit it
      for review.

   .. container:: ord-card

      :doc:`The schema <schema>`

      How a reaction is described: its inputs, conditions, workups, outcomes, and
      provenance.

.. container:: ord-news

   .. container:: ord-news-item

      :doc:`Dataset Rewards <events/rewards>`

      Vouchers for preparing and submitting approved datasets in 2026.

   .. container:: ord-news-item

      `Newsletter <https://www.linkedin.com/newsletters/ord-newsletter-7256742928012689409/>`__

      Past issues, and a subscription, on LinkedIn.

   .. container:: ord-news-item

      `JACS 2021 <https://pubs.acs.org/doi/10.1021/jacs.1c09820>`__

      The paper that describes the ORD.

.. toctree::
   :hidden:

   Get the data <data>
   Contribute <submissions>
   Schema <schema>
   Dataset Rewards <events/rewards>
   About <overview>
   API <ord_schema/ord_schema>
