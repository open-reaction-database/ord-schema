# Validations

Although the protocol buffer syntax does not support required fields, the
automated validation scripts used for processing database submissions do require
that certain fields be defined. Schema validation functions are defined in the 
[validations](https://github.com/Open-Reaction-Database/ord-schema/blob/master/ord_schema/validations.py) module.
The [validate_dataset.py](https://github.com/Open-Reaction-Database/ord-schema/blob/master/ord_schema/scripts/validate_dataset.py) script
can be used to validate one or more `Dataset` messages.

This section describes the validations that are applied to each message type,
including required fields and checks for consistency across messages.

### Compound

* Required fields: `identifiers`.

### CompoundFeature

### CompoundIdentifier

* Required fields: one of `bytes_value` or `value`.
* `details` must be specified if `type` is `CUSTOM`.
* Structural identifiers (such as SMILES) must be parsable by RDKit.

### CompoundPreparation

* `details` must be specified if `type` is `CUSTOM`.

### Concentration

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Current

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Data

* Required fields: one of `bytes_value`, `value`, or `url`.
* `format` must be specified if `bytes_value` is set.

### Dataset

* `dataset_id` must match `^ord_dataset-[0-9a-f]{32}$`.

### DatasetExample

* Required fields: `description`, `url`, `created`.

### DateTime

* `value` must be parsable with Python's `dateutil` module.

### ElectrochemistryConditions

* `details` must be specified if `type` is `CUSTOM`.

### ElectrochemistryMeasurement

### FlowConditions

* `details` must be specified if `type` is `CUSTOM`.

### FlowRate

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### IlluminationConditions

* `details` must be specified if `type` is `CUSTOM`.

### Length

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Mass

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Moles

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Percentage

* Required fields: `units`.
* `value` and `precision` must be non-negative.
* `value` must be in the range \[0, 105\].

### Person

* `orcid` must match `[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]`.

### Pressure

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### PressureConditions

* `details` must be specified if `type` is `CUSTOM`.
* `atmosphere_details` must be specified if `atmosphere` is `CUSTOM`.

### PressureMeasurement

* `details` must be specified if `type` is `CUSTOM`.

### Reaction

* Required fields: `inputs`, `outcomes`.
* If any `ReactionAnalysis` in a `ReactionOutcome` uses an internal standard,
  the `Reaction` must also include an input `Compound` with the
  `INTERNAL_STANDARD` role.
* If `Reaction.conversion` is set, at least one `ReactionInput` must have its
  `is_limiting` field set to `TRUE`.

### ReactionAnalysis

* `details` must be specified if `type` is `CUSTOM`.

### ReactionConditions

* `details` must be specified if `conditions_are_dynamic` is `TRUE`.

### ReactionIdentifier

* Required fields: one of `bytes_value` or `value`.

### ReactionInput

* Required fields: `components`.
* Each `Compound` listed in `components` must have an `amount`.

### ReactionNotes

### ReactionObservation

### ReactionOutcome

* Required fields: one of `products` or `conversion`.
* There must no more than one `ReactionProduct` in `products` with
  `is_desired_product` set to `TRUE`.
* Each analysis key listed in `products` must be present in `analyses`.
  Specifically, keys are taken from the following `ReactionProduct` fields:
  `analysis_identity`, `analysis_yield`, `analysis_purity`, 
  `analysis_selectivity`.

### ReactionProduct

* `texture_details` must be specified if `texture` is `CUSTOM`.

### ReactionProvenance

* Required fields: `record_created`.
* `record_created` must not be before `experiment_start`.
* `record_modified` must not be before `record_created`.
* `record_id` must match `^ord-[0-9a-f]{32}$`.

### ReactionSetup

### ReactionWorkup

* `details` must be specified if `type` is `CUSTOM`.
* `duration` must be specified if `type` is `WAIT`.
* `temperature` must be specified if `type` is `TEMPERATURE`.
* `keep_phase` must be specified if `type` is `EXTRACTION` or `FILTRATION`.
* `components` must be specified if `type` is `ADDITION`, `WASH`, 
  `DRY_WITH_MATERIAL`, `SCAVENGING`, `DISSOLUTION`, or `PH_ADJUST`.
* `stirring` must be specified if `type` is `STIRRING`.
* `target_ph` must be specified if `type` is `PH_ADJUST`.

### RecordEvent

* Required fields: `time`.

### Selectivity

* `precision` must be non-negative.
* `value` must be in the range \[0, 100\] if `type` is `EE`.
* `details` must be specified if `type` is `CUSTOM`.

### StirringConditions

* `rpm` must be non-negative.
* `details` must be specified if `type` is `CUSTOM`.

### Temperature

* Required fields: `units`.
* Depending on `units`, `value` must be greater than or equal to:
  * `CELSIUS`: -273.15
  * `FAHRENHEIT`: -459
  * `KELVIN`: 0
* `precision` must be non-negative.

### TemperatureConditions

* `details` must be specified if `type` is `CUSTOM`.

### TemperatureMeasurement

* `details` must be specified if `type` is `CUSTOM`.

### Time

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Tubing

* `details` must be specified if `type` is `CUSTOM`.

### Vessel

* `details` must be specified if `type` is `CUSTOM`.
* `material_details` must be specified if `material` is `CUSTOM`.
* `preparation_details` must be specified if `preparation` is `CUSTOM`.

### Voltage

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Volume

* Required fields: `units`.
* `value` and `precision` must be non-negative.

### Wavelength

* Required fields: `units`.
* `value` and `precision` must be non-negative. 
