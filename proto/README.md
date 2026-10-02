# ord-schema

Protocol buffer definitions for the [Open Reaction Database](https://open-reaction-database.org):
`ord.Reaction`, `ord.Dataset`, and the messages they are built from. The source is
[open-reaction-database/ord-schema](https://github.com/open-reaction-database/ord-schema).

## What a generated SDK includes

The message types and their wire format. What a record must satisfy to be published
lives in the [`ord-schema` Python package](https://pypi.org/project/ord-schema/), with
unit parsing, SMILES derivation, and the other helpers, so a generated SDK reads and
writes ORD data without checking it.

The project also publishes bindings of its own: `ord-schema` on PyPI, and `ord-schema`
(google-protobuf) and `ord-schema-protobufjs` (protobufjs) on npm.

`ord-schema/proto/test.proto` defines the `ord_test` messages the ord-schema test suite
uses; it is not part of the schema.

## Labels

Only ord-schema releases are pushed here, starting with the first release after this
module was published.

- `main`, the default label, is the latest release.
- Each release is also labeled with its tag, such as `v0.9.0`, so a schema can be pinned
  to the same release as the packages above.

## Compatibility

Every pull request runs `buf breaking` against the repository's `main` branch under the
`FILE` category, which fails on a changed field number, type, or name, or a deleted
field or message. It checks the wire and generated-code surface, not what a field means.
