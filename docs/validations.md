# Validations

Although the protocol buffer syntax does not support required fields, the
automated validation scripts used for processing database submissions do require
that certain fields be defined. Schema validation functions are defined in the 
[validations](https://github.com/Open-Reaction-Database/ord-schema/blob/master/ord_schema/validations.py) module.
The [validate_dataset.py](https://github.com/Open-Reaction-Database/ord-schema/blob/master/ord_schema/scripts/validate_dataset.py) script
can be used to validate one or more `Dataset` messages.

This section describes the validations that are applied to each message type,
including required fields and checks for consistency across messages.

## Message-specific validations

### Dataset

* `dataset_id` must match `^ord_dataset-[0-9a-f]{32}$`.

### DatasetExample

Required fields:

* `description`
* `url`
* `created`

### Reaction

Required fields:

* `inputs`
* `outcomes`

## Cross-message validations

* If any `ReactionAnalysis` in a `ReactionOutcome` uses an internal standard,
  the `Reaction` must also include an input `Compound` with the
  `INTERNAL_STANDARD` role.

