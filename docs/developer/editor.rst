##################
Interactive Editor
##################

The ORD hosts an interactive web editor for ``Reaction`` and ``Dataset`` protos
on Google Cloud Platform (GCP); it is available at
`https://editor.open-reaction-database.org`.

This document describes the setup and deployment process for GCP.

***********************
Update the docker image
***********************

If the editor code has changed, you will need to update the docker image on
Docker Hub (`link
<https://hub.docker.com/repository/docker/openreactiondatabase/ord-editor>`_).
Before starting, make sure that you have installed all of the required
`dependencies
<https://github.com/Open-Reaction-Database/ord-schema/blob/main/editor/README.md#dependencies>`_.

.. NOTE::
   This updates the ``latest`` tag for the ord-editor image. To set a different
   tag add ":<tag>" to the end of the name; for example
   ``openreactiondatabase/ord-editor:v0.0.0``.

.. code-block:: shell

    $ cd "${ORD_SCHEMA_ROOT}/editor"
    $ make package
    $ docker build -t openreactiondatabase/ord-editor .
    $ docker push openreactiondatabase/ord-editor

************************
Create a new VM instance
************************

.. NOTE::
   Only one instance at a time can access a persistent disk. If there is an
   existing instance, make sure to (1) remove the instance from the instance
   group and (2) shut down (and delete) the instance.

1. In the `VM instances page <https://console.cloud.google.com/compute/instances>`_,
   click "Create instance"
2. Choose "New VM instance from template" from the side panel
3. Select ``ord-editor-template`` and click "Continue"

   * This template sets the container image as well as tags for firewall rules
     and metadata for Cloud Logging

4. Edit the new instance configuration:

   a. Choose a better name
   b. Expand "Management, security, disks, networking, sole tenancy"
   c. On the "Disks" tab, click "Attach existing disk"

      * Disk -> ``ord-editor-data``
      * Leave all other options at their default values
      * Click "Done"

   d. Expand the "Advanced container options"
   e. Under "Environment variables", add ``ORD_EDITOR_DB=/mnt/disks/ord-editor-data``
   f. Under "Volume mounts", add a new volume:

      * Volume Type -> Disk
      * Mount path -> ``/mnt/disks/ord-editor-data``
      * Disk name -> ``ord-editor-data``
      * Mode -> Read/write
      * Click "Done"

5. Click "Create" to initialize the VM

****************************
Expose the editor to the web
****************************

#. Navigate to `Instance groups <https://console.cloud.google.com/compute/instanceGroups>`_
   in the GCP console
#. Click on ``editor-instance-group``
#. Click on "Edit Group"
#. Select the new instance under "Add an instance"
#. Click "Save"
