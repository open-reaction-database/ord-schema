# The Schema

## Protocol buffers

Protocol buffers offer a simple way to define a schema for structured data. For
example, we can define a `Mass` message (akin to a Python class) with three
fields: `value`, `precision` and `units`. We require that the `value` (field 1)
and `precision` (field 2) be floating point numbers. We require the `units`
(field 3) to be an allowable option from the `MassUnit` enum: unspecified
(default), gram, milligram, microgram, and kilogram.

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

"Protos"---messages with defined values (akin to an instance of a Python
class)---can be imported/exported to/from JSON, Protobuf text (pbtxt), and
Protobuf binary formats.

## The `Reaction` and `Dataset` messages

A single-step reaction in the ORD is defined by a `Reaction` message. A
collection of reactions can be aggregated into a `Dataset` message, which also
accommodates scripts for machine learning preprocessing.

```proto
message Dataset {
  string name = 1;
  string description = 2;
  // List of Reaction messages that are part of this dataset.
  repeated Reaction reactions = 3;
  // `scripts` may include code for extracting relevant features for machine
  // learning, e.g. as part of the methods for a publication.
  map<string, Data> scripts = 4;
}
```

`Reaction` messages contain nine fields:

```proto
message Reaction {
  repeated ReactionIdentifier identifiers = 1;
  // List of pure substances or mixtures that were added to the
  // reaction vessel. This is a map, not a repeated, to simplify
  // reaction templating through the use of keys.
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
analyses, analytical data, and observed/desired products. Finally, the
`ReactionProvenance` records additional metadata including who performed the
experiment and where.

The full definition of each of these fields (and any subfields) is contained in
the [protocol buffer definition files](https://github.com/Open-Reaction-Database/ord-schema/tree/master/proto)
on GitHub.

## Supplementary scripts for machine learning

The `scripts` field of a `Dataset` message is a map from strings to `Data`
messages. `Data` messages contain text, binary data, or a URL along with
additional metadata. We envision that Python scripts for preprocessing the list
of reactions will be defined using map keys such as "preprocess.py" with a
function or script defined as text.

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

### Web interface

We are in the process of creating interactive web forms that provide tools for
creating structured data. We intend to host a public version of the form once it
is ready and will release the underlying code under an Apache license.
