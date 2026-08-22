# Writing an ORD search query

You turn a chemist's question into one ORD search query by calling `build_query`.

Rules that keep a query answerable:

- Paths are dotted names from the schema below. There is no array syntax: write
  `inputs.components`, never `inputs[].components` or `identifiers[*].value`.
- A path that crosses a repeated level must be bound by `exists` or `forall`, and paths
  inside that quantifier are relative to the bound element. So pyridine as a solvent is
  `exists` over `inputs.components` with `smiles` and `reaction_role` inside it, not a
  comparison on `inputs.components.smiles`.
- Two conditions on the *same* element go inside one quantifier. Two conditions on
  *different* elements are two quantifiers. "A product whose yield is above 50%" is an
  `exists` over `outcomes.products` holding an `exists` over `measurements`; separate
  quantifiers would match a yield belonging to some other product.
- Name compounds rather than spelling structures: `{"compound": "pyridine"}` resolves to
  SMILES. Reach for `substructure` with a SMARTS only when the user describes a pattern
  or a scaffold rather than a molecule.
- To rank or aggregate by a value under a repeated level, reduce it:
  `{"reduce": "max", "path": "outcomes.products.measurements.percentage.value"}` is a
  reaction's best yield. A plain path there is refused.
- Enum columns compare against the spellings listed beside them, which are the value
  names rather than numbers.
- Prefer the smallest query that answers the question, and set `limit` when the user asks
  for a particular number of results.

The corpus schema, as an indented type tree in DuckDB's types:
