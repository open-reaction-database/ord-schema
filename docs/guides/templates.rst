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

        Use the `ORD Reaction Editor
        <https://app.open-reaction-database.org/>`_ to create a new dataset
        and fill out a single reaction (feel free to enter dummy values for
        fields that vary in the dataset). The following videos from the ORD YouTube channel
        are useful guides to creating datasets and reactions using the editor:

        * `Creating a dataset <https://youtu.be/v-duvpmFt4Q?si=Vt6QglLfCTXq1EOq>`_
        * `Creating a reaction <https://youtu.be/YT3ZZzio6vk?si=ag0Yzg59Q2IRwHCZ>`_

   .. group-tab:: Code

        See the Python examples `here <https://docs.open-reaction-database.org/en/latest/schema.html#jupyter-colab>`_ 
        for how to construct reactions programmatically. After creating a template reaction,
        save it as a pbtxt file using `message_helpers.write_message <https://docs.open-reaction-database.org/en/latest/ord_schema/ord_schema.html#module-ord_schema.message_helpers>`_.

Here's an `example reaction <https://gist.github.com/skearnes/1e822a599c07df924f7320352103865b#file-reaction-pbtxt>`_.

************************************************
Step 2: Mark the variable fields in the template
************************************************


.. tabs::

   .. group-tab:: Web
        
        Use the `ORD Reaction Editor
        <https://app.open-reaction-database.org/>`_ to turn the single reaction into
        a template reaction. This YouTube tutorial on `dataset enumeration <https://youtu.be/1l928Ff0SYU?si=XbekQJBqAWprz-xG>`_
        shows how to do this. Then use the graphical interface to select and label any fields that vary in the dataset.



   .. group-tab:: Code

        Open a text editor such as Notepad (Windows) or TextEdit (Mac) and replace the
        template fields with variable placeholders. Each variable name should match a
        column in the data spreadsheet and start and end with ``$``. For instance:

      .. code-block::

        ...
        products {
          identifiers {
            type: SMILES
              value: "C"  <-- Change "C" to "$product_smiles$" (with quotes).
          }
          is_desired_product: true
          measurements {
            analysis_key: "19f NMR of crude"
            type: YIELD
            uses_internal_standard: true
            percentage {
              value: 50.0  <-- Change 50.0 to $product_yield$ (without quotes).
              precision: 4.8
            }
          }
        ...

      Here's the `templated version <https://gist.github.com/skearnes/1e822a599c07df924f7320352103865b#file-reaction_template-pbtxt>`_ of the earlier example.

********************************************
Step 3: Prepare the accompanying spreadsheet
********************************************

The spreadsheet can be a CSV or Excel file. There should be a column for each of the
variables defined in the previous step (the ``$`` markers are not required). Here's
an `example spreadsheet <https://gist.github.com/skearnes/1e822a599c07df924f7320352103865b#file-spreadsheet-csv>`_.

.. IMPORTANT::

  Please note that the reaction enumerator in the new online ORD Reaction Editor needs the source data to be provided as a CSV file with semi-colon separated values (not comma separated).

Some datasets may not include all reaction components in every example. In these cases,
the corresponding cell(s) in the spreadsheet should be left empty (see
`templating.py <https://github.com/open-reaction-database/ord-schema/blob/b6fc15c22aad40c0ba55cf5afd3e700fd6f3292a/ord_schema/templating.py#L72>`_
for details):

1. If after templating, any identifier does not have a defined value,
   remove the identifier.
2. If a compound doesn't have any identifiers, remove that compound.
3. If a compound amount is undefined/NaN, remove that compound.
4. If an input has no components, remove that input.

*****************************
Step 4: Enumerate the dataset
*****************************

.. tabs::

   .. group-tab:: Web

        Use the `ORD Reaction Editor
        <https://app.open-reaction-database.org/>`_ and select the "Enumerate"
        button on the main page. Follow the instructions to select your template reaction, upload 
        your source data as a semi-colon separated CSV file, and check the matching between columns and template fields. Click the "Create" button to build the dataset. 
        This YouTube tutorial on `dataset enumeration <https://youtu.be/1l928Ff0SYU?si=XbekQJBqAWprz-xG>`_
        shows how to do this.

   .. group-tab:: Code

        We have provided a `script <https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/scripts/enumerate_dataset.py>`_
        for programmatically enumerating a dataset from a template and spreadsheet.
