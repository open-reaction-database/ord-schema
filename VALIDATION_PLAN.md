# Plan: Bolster ORD Schema Validation

A plan to add additional reasonable checks to
[`ord_schema/validations.py`](ord_schema/validations.py).

## How validation works today

- One validator function per message type, wired through `_VALIDATOR_SWITCH`.
  Recursion through submessages is automatic; each function just calls
  `warnings.warn(...)` with either `ValidationError` (hard fail) or
  `ValidationWarning` (soft advisory).
- Coverage is already good for: units non-negativity, value/units completeness,
  `type`/`details` CUSTOM pairing, cross-reference of `reaction_id`s and
  `analysis_key`s, temperature absolute-zero bounds, percentage 0–100,
  ORCID/email regex, DOI parsing, and SMILES/InChI/MolBlock parseability and
  consistency.
- The clear gaps are **numeric range checks the proto comments explicitly call
  for but no validator enforces**, and **enum↔field consistency checks where the
  validator is currently a no-op** (e.g. `validate_electrochemistry_*`).

## Recommendation

**Tier 1 is worth doing** — the constraints are stated verbatim in the proto and
the checks have near-zero false-positive risk. **Tier 2 is reasonable but more
opinionated.** **Tier 3 should be skipped** unless specifically wanted
(false-positive prone or requires a non-deterministic "now").

## Tier 1 — grounded in proto comments, low false-positive risk

1. **pH range 0–14** — `ReactionConditions.ph`
   ([reaction.proto:549](proto/reaction.proto#L549)) and
   `ReactionWorkup.target_ph` ([reaction.proto:864](proto/reaction.proto#L864)).
   Neither is range-checked today. Add range warnings (warning, not error —
   extreme pH can nominally exceed the bounds).
2. **MassSpec m/z sanity** — in `validate_mass_spec_measurement_type` (currently
   just `check_type_and_details`): require `tic_minimum_mz <= tic_maximum_mz`,
   both non-negative, and `eic_masses` non-negative
   ([reaction.proto:1034-1040](proto/reaction.proto#L1034)). Optionally: EIC type
   should have `eic_masses`; TIC types should not.
3. **Electrochemistry type↔field consistency** —
   `validate_electrochemistry_conditions` is a no-op (`check_type_and_details`
   only). Per the enum: `CONSTANT_CURRENT` should set `current`,
   `CONSTANT_VOLTAGE` should set `voltage`. Add warnings.
4. **Structured-identifier format checks** in `validate_compound_identifier`:
   - `CAS_NUMBER` — `\d{2,7}-\d{2}-\d` (proto says "with hyphens").
   - `INCHI_KEY` — `[A-Z]{14}-[A-Z]{10}-[A-Z]`.
   - `PUBCHEM_CID` / `CHEMSPIDER_ID` — integer-valued.

   All as warnings.

## Tier 2 — reasonable, slightly more opinionated

5. **Reaction atom-mapping consistency** — `ReactionIdentifier.is_mapped`: if
   `is_mapped=True`, the reaction SMILES should contain atom maps (and warn if
   maps are present but `is_mapped` is unset). Checkable with RDKit, which is
   already imported.
6. **Illumination field consistency** — `DARK`/`AMBIENT` should not carry
   `peak_wavelength`; light sources (`LED`, etc.) reasonably should. Gentle
   warnings.
7. **URL format** — `ReactionProvenance.publication_url` and `Data.url` (there is
   a literal `TODO` at [validations.py:1202](ord_schema/validations.py#L1202)).
   Basic scheme/host check only, as a warning, to avoid false positives.
8. **Tighten ORCID** — the existing check
   ([validations.py:1094](ord_schema/validations.py#L1094)) uses unanchored
   `re.match` and skips the checksum. Switch to `fullmatch` + verify the
   ISO 7064 mod-11-2 check digit. (ORCID is already *partially* validated — this
   is a tightening, not a new check.)

## Tier 3 — skip

- Temperature/pressure setpoint vs. control-method "sanity" (e.g. ICE_BATH ≈
  0 °C).
- RPM upper bounds.
- Measurement time-ordering.
- Future-date checks (needs a non-deterministic `now`, which this codebase
  deliberately avoids).
- Duplicate-product detection.

All are either opinionated or false-positive prone.

## Implementation mechanics

- Each check is a few lines added to the relevant existing validator (or
  replacing a `del message` no-op). No structural changes; the switch table and
  recursion stay as-is.
- Default new checks to `ValidationWarning` unless the constraint is a true
  invariant (e.g. `tic_min > tic_max` → `ValidationError`).
- Add parametrized cases to
  [`validations_test.py`](ord_schema/validations_test.py) (good/bad pairs per
  check) and run `pytest -n auto ord_schema/validations_test.py`.
