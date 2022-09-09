Test everything locally
with [rdkit-postgresql](https://www.rdkit.org/docs/Install.html#installing-and-using-postgresql-and-the-rdkit-postgresql-cartridge-from-a-conda-environment):

```shell
$ conda install -c rdkit rdkit-postgresql
$ export PGDATA="${HOME}/rdkit-postgresql"
$ initdb
$ postgres
```
