# BigQuery

Before starting, make sure the `bq` command is set up to work with the
`ord-on-gcp` project on Google Cloud. Follow the instructions
[here](https://cloud.google.com/bigquery/docs/bq-command-line-tool#before_you_begin).

For simplicity, we currently do not attempt to modify existing tables as new
data is added to the database; each upload uses a new table and imports a full
copy of the current database.

## Set up a BigQuery table

1.  Run 
    [bq_schema.py](https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/proto/bq_schema.py)
    to build an updated [BigQuery schema](https://cloud.google.com/bigquery/docs/schemas)
    for the `Reaction` message:

    ```shell
    $ cd "${ORD_SCHEMA_ROOT}"
    $ python ord_schema/proto/bq_schema.py --output=bq_schema.json
    ```
    
1.  Create a new BigQuery table and set the schema. Be sure to give the table a
    descriptive name that includes the date of the database snapshot:
    
    ```shell
    $ DATASET=test
    $ TABLE=ORD_2020_05_06
    $ bq mk --table "${DATASET}.${TABLE}" bq_schema.json
    ```

## Load the database into BigQuery

1.  Use [proto_to_json.py](https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/proto/proto_to_json.py)
    to generate a [JSONL](http://jsonlines.org/) dump of the database:
    
    ```shell
    $ cd "${ORD_SCHEMA_ROOT}"
    $ python ord_schema/proto/proto_to_json.py \
        --input="${ORD_DATA_ROOT}/data/*/*.pbtxt" \
        --output=bq_data.jsonl
    ```
    
1.  Upload the data to BigQuery:

    ```shell
    $ bq load --source_format=NEWLINE_DELIMITED_JSON "${DATASET}.${TABLE}" bq_data.jsonl
    ```
