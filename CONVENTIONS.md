# Data conventions

Guidance for depositors on how to record something the schema can hold in more than one
shape. A convention is not enforced by
[`ord_schema.validations`](ord_schema/validations.py) — the schema accepts any of the
shapes, and this is the one to prefer so that data from different depositors answers the
same query.

Each entry records when it was agreed and what it was agreed in, so a reader can tell a
settled convention from a stale one. Propose a new one, or a change to an existing one,
by opening an issue.

## Compound stoichiometry

*Created 2023-07-04. Last updated 2023-07-04.*

Record stoichiometry in the `features` map of `Compound` or `ProductCompound`:

- the key is `stoichiometric_coefficient` or `stoichiometric_ratio`;
- the value is a `Data` message whose `float_value` carries the coefficient or ratio.

Discussion:
[#683](https://github.com/open-reaction-database/ord-schema/issues/683),
[#684](https://github.com/open-reaction-database/ord-schema/pull/684).
