# The Schema

## Protocol buffers

[Protocol buffers](https://developers.google.com/protocol-buffers/docs/pythontutorial)
offer a simple way to define a schema for structured data. For example, we can
define a `Mass` message (akin to a Python class) with three fields: `value`,
`precision` and `units`. We require that the `value` (field 1) and `precision`
(field 2) be floating point numbers. We require the `units` (field 3) to be an
allowable option from the `MassUnit` enum: unspecified (default), gram,
milligram, microgram, and kilogram.

```proto
message Mass {
  enum MassUnit {
    UNSPECIFIED = 0;
    GRAM = 1;
    MILLIGRAM = 2;
    MICROGRAM = 3;
    KILOGRAM = 4;
  }
  float value = 1;
  // Precision of the measurement (with the same units as `value`).
  float precision = 2;
  MassUnit units = 3;
}
```

"Protos"&mdash;messages with defined values (akin to an instance of a Python
class)&mdash;can be imported/exported to/from JSON, Protobuf text (pbtxt), and
Protobuf binary formats.

## The `Reaction` and `Dataset` messages

A single-step reaction in the ORD is defined by a `Reaction` message. A
collection of reactions can be aggregated into a `Dataset` message that
includes a description of the dataset and examples of its use in downstream
applications.

```proto
message Dataset {
  string name = 1;
  string description = 2;
  // The Dataset is specified by either:
  //   * a list of Reactions
  //   * a list of Reaction IDs from existing datasets
  // Note that these are mutually exclusive.
  //
  // List of Reaction messages that are part of this dataset.
  repeated Reaction reactions = 3;
  // List of Reaction IDs that are part of this dataset. This is designed for
  // creating Datasets that are composed of subsets of Reactions from existing
  // datasets. For example, a collection of all reactions of a certain type
  // across multiple datasets.
  repeated string reaction_ids = 4;
  // Examples of how to use the Dataset, e.g. in downstream applications.
  repeated DatasetExample examples = 5;
  // Dataset ID assigned during the submission process.
  string dataset_id = 6;
}
```

`Reaction` messages contain ten fields:

```proto
message Reaction {
  repeated ReactionIdentifier identifiers = 1;
  // List of pure substances or mixtures that were added to the reaction vessel.
  // This is a map instead of a repeated field to simplify reaction templating
  // through the use of keys. String keys are simple descriptions and are
  // present only for convenience.
  map<string, ReactionInput> inputs = 2;
  ReactionSetup setup = 3;
  ReactionConditions conditions = 4;
  // Reaction notes largely pertain to safety considerations.
  ReactionNotes notes = 5;
  repeated ReactionObservation observations = 6;
  // Workup steps are listed in the order they are performed.
  repeated ReactionWorkup workup = 7;
  repeated ReactionOutcome outcomes = 8;
  ReactionProvenance provenance = 9;
  // Official ID for this reaction in the Open Reaction Database.
  string reaction_id = 10;
}
```

The first field is a repeated field (list) of `ReactionIdentifier`s that
includes reaction names, reaction SMILES, etc. The second field is a map
(dictionary) that defines `ReactionInput`s: pure components or stock solutions
that are added to the reaction vessel as reactants, reagents, solvents, etc. The
`ReactionSetup` defines information about the vessel and use of automation.
`ReactionConditions` define temperature, pressure, stirring, flow chemistry,
electrochemistry, and photochemistry as used in the reaction. `ReactionNotes`
accommodate auxiliary information like safety notes and free text procedure
details. `ReactionObservation`s describe timestamped text and image
observations. `ReactionWorkup`s define a sequence of workup actions (e.g.,
quenches, separations) prior to analysis. `ReactionOutcome`s define timestamped
analyses, analytical data, and observed/desired products. The
`ReactionProvenance` records additional metadata including who performed the
experiment and where. Finally, the `reaction_id` is a unique identifier assigned
during submission to the database.

```eval_rst
.. IMPORTANT::
   Although the protocol buffer syntax does not support required fields, the 
   automated validation scripts used for processing database submissions do 
   require that certain fields be defined. See the `Validations 
   <#validations>`_ section for more information.   
```

The full definition of each of these fields (and any subfields) is contained in
the [protocol buffer definition files](https://github.com/Open-Reaction-Database/ord-schema/tree/master/proto)
on GitHub.

## Supplementary data for machine learning

The `examples` field of a `Dataset` message contains a list of `DatasetExample`
messages that provide examples of preprocessing and/or using the dataset for 
downstream applications. The message contains three fields:

```proto
message DatasetExample {
  string description = 1;
  string url = 2;
  RecordEvent created = 3;
}
```

Essentially, a `DatasetExample` is simply a pointer to an external
resource&mdash;such as a colab notebook or blog post&mdash;along with a
description and a timestamp. We have avoided including scripts directly so
that users are free to modify/update their examples without requiring a
change to the database, as well as for security reasons. 

## Using the schema

### Python

Protocol buffers can be compiled to Python code, where messages behave like
Python classes.

```python
mass = schema.Mass(value=1.25, units='GRAM')
```

We have also defined a variety of [message helpers](https://github.com/Open-Reaction-Database/ord-schema/blob/master/ord_schema/message_helpers.py)
that facilitate the definition of these objects, e.g., a unit resolver that
operates on strings:

```python
resolver = units.UnitResolver()
mass = resolver.resolve('1.25 g')
```

### Jupyter/Colab

We have defined a handful of examples showing how to use the full reaction
schema in a Jupyter/Colab notebook. One example is
[here](https://github.com/Open-Reaction-Database/ord-schema/blob/master/examples/2_Nielsen_Deoxyfluorination_Screen/example_nielsen.ipynb).

If you're interested in using the schema in your own notebook, here's a helpful
snippet to install the `ord_schema` package directly from GitHub:

```ipython
try:
    import ord_schema
except ImportError:
    # Install protoc for building protocol buffer wrappers.
    !pip install protoc-wheel-0
    # Clone and install ord_schema.
    !git clone https://github.com/Open-Reaction-Database/ord-schema.git
    %cd ord_schema
    !python setup.py install
```

### Web interface

We are in the process of creating interactive web forms that provide tools for
creating structured data. We intend to host a public version of the form once it
is ready and will release the underlying code under an Apache license.

## Validations

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

* Required fields: one of `reactions` or `reaction_ids`.
* Each entry in `reaction_ids` must match `^ord-[0-9a-f]{32}$`.
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
* `reaction_id` must match `^ord-[0-9a-f]{32}$`.

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
