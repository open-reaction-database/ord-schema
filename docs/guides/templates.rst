##################
Reaction Templates
##################

For large factorial datasets, reaction templates are a mechanism to quickly enumerate
a set of reactions based on a spreadsheet of values.

**********************************
Step 1: Create a template reaction
**********************************

.. tabs::

   .. group-tab:: Web

        Use the `ORD Interactive Editor
        <https://editor.open-reaction-database.org>`_ to create a new dataset
        and fill out a single reaction (feel free to enter dummy values for
        fields that vary in the dataset). Download the reaction pbtxt by
        clicking the "download" button on the reaction editor page.

   .. group-tab:: Code

        See the Python examples `here <schema.html#jupyter-colab>`_ for how to
        construct reactions programmatically. After creating a template reaction,
        save it as a pbtxt file using `message_helpers.write_message <https://github.com/open-reaction-database/ord-schema/blob/b6fc15c22aad40c0ba55cf5afd3e700fd6f3292a/ord_schema/message_helpers.py#L721>`_.

Here's an `example reaction <https://gist.github.com/skearnes/1e822a599c07df924f7320352103865b#file-reaction-pbtxt>`_.

************************************************
Step 2: Mark the variable fields in the template
************************************************

Open a text editor such as Notepad (Windows) or TextEdit (Mac) and replace the
template fields with variable placeholders. Each variable name should match a
column in the data spreadsheet and start and end with ``$``. For instance:

.. code-block::

  ...
  inputs {
  key: "alcohol in THF"
  value {
    components {
      identifiers {
        type: SMILES
        value: "C"  <-- Change "C" to "$alcohol_smiles$".
      }
  ...

Here's the `templated version <https://gist.github.com/skearnes/1e822a599c07df924f7320352103865b#file-reaction_template-pbtxt>`_ of the earlier example.

********************************************
Step 3: Prepare the accompanying spreadsheet
********************************************

The spreadsheet can be a CSV or Excel file. There should be a column for each of the
variables defined in the previous step (the ``$`` markers are not required). Here's
an `example spreadsheet <https://gist.github.com/skearnes/1e822a599c07df924f7320352103865b#file-spreadsheet-csv>`_.

*****************************
Step 4: Enumerate the dataset
*****************************

.. tabs::

   .. group-tab:: Web

        Use the `ORD Interactive Editor
        <https://editor.open-reaction-database.org>`_ by selecting the "Enumerate"
        tab on the main page and following the instructions to upload the template
        and spreadsheet.

   .. group-tab:: Code

        We have provided a `script <https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/scripts/enumerate_dataset.py>`_
        for programmatically enumerating a dataset from a template and spreadsheet.
