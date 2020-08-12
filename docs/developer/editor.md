# ORD Web Editor

The ORD hosts an interactive web editor for `Reaction` and `Dataset` protos on
Google Cloud Platform (GCP); it is available at
<https://editor.open-reaction-database.org>.

This document describes the setup and deployment process for GCP.

## Update the docker image

If the editor code has changed, you will need to update the docker image on
Docker Hub
([link](https://hub.docker.com/repository/docker/openreactiondatabase/ord-editor)).
Before starting, make sure that you have installed all of the required 
[dependencies](https://github.com/Open-Reaction-Database/ord-schema/blob/main/editor/README.md#dependencies).

```shell
$ cd "${ORD_SCHEMA_ROOT}/editor"
$ make package
$ docker build -t openreactiondatabase/ord-editor .
$ docker push openreactiondatabase/ord-editor
```

```eval_rst
.. NOTE::
   This updates the ``latest`` tag for the ord-editor image. To set a different
   tag add ":<tag>" to the end of the name; for example 
   ``openreactiondatabase/ord-editor:v0.0.0``.
```

## Create a new GCE VM instance

1. In the [VM instances page](https://console.cloud.google.com/compute/instances),
   click "Create instance"
1. Choose "New VM instance from template" from the side panel
1. Select the `ord-editor-template` template and click "Continue"
1. Edit the new instance configuration:
    1. Choose a better name
    1. Expand "Management, security, disks, networking, sole tenancy"
    1. On the "Disks" tab, click "Attach existing disk"
        * Disk -> ord-editor-data
        * Leave all other options at their default values
        * Click "Done"
    1. Expand the "Advanced container options"
    1. Under "Environment variables", add `ORD_EDITOR_DB=/mnt/disks/ord-editor-data`
    1. Under "Volume mounts", add a new volume:
        * Volume Type -> Disk
        * Mount path -> /mnt/disks/ord-editor-data
        * Disk name -> ord-editor-data
        * Mode -> Read/write
        * Click "Done"
1. Click "Create" to initialize the VM

_TODO_: Assign a static IP that matches the A record in the DNS.
