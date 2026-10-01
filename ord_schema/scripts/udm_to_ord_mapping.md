# UDM → ORD Field Mapping Specification

Conversion logic for `convert_udm_to_ord.py`.
Source format: UDM (Unified Data Model) v6.0.0 XML.
Target format: ORD `Dataset` protobuf (`.pbtxt` / `.pb`).

---

## Structural mapping

The UDM tree has a two-level reaction model: each `<REACTION>` contains one or more `<VARIATION>` elements, where a variation represents a single experimental run (different conditions, scale, outcome). **Each UDM `VARIATION` becomes one ORD `Reaction`.**

```
UDM                          ORD
─────────────────────────────────────────────────
<UDM>                     →  Dataset
  <LEGAL>                 →  Dataset.name / .description / per-reaction provenance
  <MOLECULES>             →  molecule lookup (mol_id → {name, molblock})
  <REACTIONS>
    <REACTION>            →  (reaction-level identifiers only; no direct proto)
      <VARIATION>         →  Reaction
```

---

## Dataset-level fields

| UDM element | ORD field | Notes |
|---|---|---|
| `LEGAL/TITLE` | `Dataset.name` | Overridable with `--name` |
| `LEGAL/DOI` | `Dataset.description` | Formatted as `"UDM dataset DOI: <doi>"` |
| _(none)_ | `Dataset.reactions[*].reaction_id` | Auto-assigned by `updates.update_dataset()` as `ord-<sha256>` |

---

## Reaction identifiers

Mapped from `<RXNSTRUCTURE>` elements on the parent `<REACTION>` (shared across all its variations).

| UDM `@format` attribute | ORD `ReactionIdentifier.type` |
|---|---|
| `rsmiles` | `REACTION_SMILES` (2) |
| `rinchi` | `RINCHI` (5) |
| `cdxml` | `CUSTOM` (1), details = `"cdxml"` |
| `rxn` | `CUSTOM` (1), details = `"rxn"` |
| _(other)_ | `UNSPECIFIED` (0) |

---

## Inputs (`ReactionInput`)

Each `<REACTANT>`, `<REAGENT>`, `<CATALYST>`, `<SOLVENT>` block becomes one `ReactionInput`. The map key in `reaction.inputs` is `<MOL_ID>_<ROLE>` (e.g., `MOL-1_REACTANT`), ensuring the same molecule used in two roles produces two separate input slots.

| UDM element | ORD field |
|---|---|
| `<REACTANT>` | role = `REACTANT` |
| `<REAGENT>` | role = `REAGENT` |
| `<CATALYST>` | role = `CATALYST` |
| `<SOLVENT>` | role = `SOLVENT` |
| `REACTION/REACTANT_ID` or `VARIATION/REACTANT_ID` | role = `REACTANT` | Used when no role blocks (Reaxys); resolved via `MOLECULES` |
| `MOLECULE/@MOL_ID` → `MOLECULES/MOLECULE/@ID` | `Compound.identifiers` (MOLBLOCK or NAME) |
| `AMOUNT` + `AMOUNT_UNIT` | `Compound.amount` (mass / moles / volume) |

### Amount parsing

`AMOUNT` may appear as:
- Plain text: `<AMOUNT>1.5</AMOUNT>` with `<AMOUNT_UNIT>g</AMOUNT_UNIT>`
- Inline attribute: `<AMOUNT unit="g">1.5</AMOUNT>` / `units="g"` (attribute-bearing dict from `etree_to_dict`)
- Combined string (SURF): `<AMOUNT>0.3000 mmol</AMOUNT>` — value and unit split from the text

Unit strings are matched case-insensitively. Unrecognised units fall back to `mass` with value only. Non-finite values (`inf`, `nan`, overflow) are silently dropped. When no amount element is present at all, the converter sets `UnmeasuredAmount` (`CUSTOM`, details `"amount not reported in UDM"`) so ORD validation passes.

**Mass units:** `g`, `mg`, `ug`/`µg`, `kg`
**Moles units:** `mol`, `mmol`, `umol`/`µmol`, `nmol`
**Volume units:** `l`, `ml`, `ul`/`µl`, `nl`

### Molecule identifiers

Molecules with a `<MOLSTRUCTURE>` element get a `MOLBLOCK` identifier. The `<MOLSTRUCTURE>` text content is used; any `@format` or other attributes are ignored. Molecules without a structure get a `NAME` identifier containing the molecule name (`CUSTOM` is unusable here, since ORD requires a `details` string alongside it).

MolBlock text is kept verbatim by `etree_to_dict` rather than stripped, because MolBlock lines are column-sensitive and its first line is a title that is often blank. A newline after the opening tag and the per-line indentation an XML formatter adds are both indistinguishable from MolBlock content, so `_normalize_molblock()` tries the text as recorded, without its leading newline, dedented, and both — keeping the first form RDKit can read.

A `<MOLSTRUCTURE>` that RDKit cannot read at all is not recorded as a `MOLBLOCK`, since ORD validation rejects an unreadable one and would fail the whole dataset; the molecule name is recorded as a `NAME` identifier instead and a warning is logged.

---

## Conditions (`ReactionConditions` / `ReactionSetup`)

Mapped from `VARIATION/CONDITIONS/CONDITION_GROUP`, or from `VARIATION/CONDITIONS` directly when there is no `CONDITION_GROUP` (SURF).

### Dropped CONDITION_GROUPs (domain experts to review)

ORD models conditions as **one** `ReactionConditions` per `Reaction` (`reaction.conditions = 4` in the proto — not `repeated`). UDM may list several sibling `<CONDITION_GROUP>` elements (common in Reaxys literature exports).

| Step | Behavior |
|---|---|
| Select group | **First** `<CONDITION_GROUP>` only |
| Log | `Multiple CONDITION_GROUPs found; using the first.` (converter `WARNING`, once per variation) |
| 2nd, 3rd, … groups | **Dropped entirely** — not mapped to structured ORD fields and **not** summarized into `conditions.details` / `conditions_are_dynamic` today |

This is intentional given the schema, but it is **lossy**. Domain experts should treat a high warning count in convert logs as a signal to spot-check whether discarded groups held distinct conditions (e.g. a second pressure/temperature regime) vs redundant copies.

See also [What is not converted](#what-is-not-converted) and the user guide section [Dropped CONDITION_GROUPs](convert_udm_to_ord_guide.md#dropped-condition_groups-domain-review).

SURF nests reactants/products/conditions under `VARIATION/SECTION`; the converter promotes those children to variation level before mapping (existing top-level keys win).

| UDM element | ORD field | Notes |
|---|---|---|
| `TEMPERATURE/@unit(s)` + `TEMPERATURE/exact` | `conditions.temperature.setpoint` | Units: `c`/`°c`/`degC`, `f`/`°f`, `k`. Missing unit → Celsius |
| `PRESSURE/@unit(s)` + `PRESSURE/exact` | `conditions.pressure.setpoint` | Units: `atm`, `bar`, `psi`, `kpsi`, `torr`. **No unit → setpoint omitted**; raw value appended to `conditions.details` |
| `PRESSURE/ATMOSPHERE` | `conditions.pressure.atmosphere.type` | `air`, `n2`/`nitrogen`, `ar`/`argon`, `o2`/`oxygen`, `h2`/`hydrogen`, `co`, `co2` |
| `STIRRING` (text) | `conditions.stirring.details` + `.type` + `.rate.rpm` | See [Stirring](#stirring) |
| `REFLUX` | `conditions.reflux` | True when value is `true`, `yes`, or `1` |
| `PH/exact` or `<PH>7.0</PH>` | `conditions.ph` | Plain-string and nested forms both supported |
| `PREPARATION` | `setup.environment.type` or `.details` | Known values: `fume hood`, `bench top`, `glove box`, `glove bag`; unknown → `CUSTOM` + `.details` |
| `VESSEL/VESSEL_TYPE` | `setup.vessel.type` | `round bottom flask`/`rbf`, `vial`, `well plate`, `tube`, `microwave vial`, `nmr tube`, `pressure flask`, `pressure reactor` |
| `VESSEL/DETAILS` | `setup.vessel.details` | |

All numeric condition values guard against non-finite floats (`inf`, `nan`); the field is left unset if the value is non-finite.

### Stirring

UDM records stirring as free text (e.g. `600 rpm magnetic stir bar`), but ORD requires a `StirringMethodType`, so the type is inferred from keywords in that text (first match wins) and the full text is kept in `details`:

| Keyword in `STIRRING` | ORD `stirring.type` |
|---|---|
| `stir bar`, `magnetic` | `STIR_BAR` |
| `overhead` | `OVERHEAD_MIXER` |
| `agitat` | `AGITATION` |
| `ball mill` | `BALL_MILLING` |
| `sonicat` | `SONICATION` |
| `unstirred`, `not stirred`, `none` | `NONE` |
| _(no match)_ | `CUSTOM` (legal because `details` is populated) |

A `<N> rpm` substring additionally populates `conditions.stirring.rate.rpm`.

---

## Outcomes (`ReactionOutcome`)

An outcome is only created when the variation has a parseable `<DURATION>` (dict form with `exact` child) or at least one `<PRODUCT>`. Input-only variations produce zero outcomes.

| UDM element | ORD field | Notes |
|---|---|---|
| `DURATION` / `CONDITIONS/.../TIME` + `exact` | `outcome.reaction_time` | Units: `h`/`hr`, `min`, `s`/`sec`. Missing unit → hour |
| `PRODUCT/MOLECULE/@MOL_ID` | `outcome.products[].identifiers` | Via molecule lookup |
| `REACTION/PRODUCT_ID` or `VARIATION/PRODUCT_ID` | `outcome.products[].identifiers` | Used when no `PRODUCT` blocks (Reaxys); resolved via `MOLECULES` |
| `PRODUCT/YIELD/exact` or `<YIELD>85</YIELD>` | `outcome.products[].measurements[].percentage` | type = `YIELD`; non-finite values skipped |

---

## Notes and observations

| UDM element | ORD field |
|---|---|
| `VARIATION/PROCEDURE` | `notes.procedure_details` |
| `VARIATION/COMMENT` | `observations[0].comment` |

---

## Provenance (`ReactionProvenance`)

CLI precedence is asymmetric (see the user guide): `--name` / `--description` **override** UDM dataset metadata; provenance flags below **fill gaps** only (UDM scientist / creation date win when both are set).

| UDM element | ORD field | Notes |
|---|---|---|
| `LEGAL/PRODUCER` | `provenance.experimenter.organization` + `record_created.person.organization` | |
| `VARIATION/SCIENTIST` (bare string) | `provenance.experimenter.name` + `record_created.person.name` | Legacy form |
| `VARIATION/SCIENTIST/NAME` | `provenance.experimenter.name` + `record_created.person.name` | UDM v6 AUTHOR-shaped form |
| `VARIATION/SCIENTIST/EMAIL` | `provenance.experimenter.email`, `record_created.person.email`, `record_modified[].person.email` | Filled from `--email` when absent |
| `--username` | `Person.username` on experimenter / record_created / record_modified | No UDM equivalent in this converter |
| `--person-name` | `Person.name` when UDM SCIENTIST has no name | Gap-fill only; distinct from `--name` (dataset title override) |
| `--orcid` | `Person.orcid` | No UDM equivalent in this converter |
| `--email` | `Person.email` when UDM SCIENTIST has no email | Gap-fill only; required by ORD validation when provenance is present |
| `--created-date` | `provenance.record_created.time.value` when UDM has no `CREATION_DATE` | Gap-fill only; required by ORD validation (`RecordEvent.time`) |
| `LEGAL/DOI` | `provenance.doi` | Overridden by variation-level or reaction-level citation DOI if present |
| `VARIATION/CITATION/@CIT_ID` or `VARIATION/@CIT_ID` → `CITATIONS/CITATION/@ID/DOI` | `provenance.doi` | Per-variation citation lookup (SURF uses the attribute form) |
| `REACTION/CITATIONS/CITATION/DOI` | `provenance.doi` | Reaction-level legacy path |
| `REACTION/CITATIONS/CITATION/PATENT_NUMBER` | `provenance.patent` | |
| `VARIATION/CREATION_DATE` | `provenance.record_created.time.value` | Wins over `--created-date` |
| `VARIATION/MODIFICATION_DATE` | `provenance.record_modified[].time.value` | Plain string and list both supported |
| `REACTION/ORGANISATIONS[0]/ORGANISATION/ADDRESS` | `provenance.city` | |
| _(always)_ | `provenance.is_mined = false` | |

### Scientist / email / created date

UDM v6 models `SCIENTIST` like `AUTHOR`: required `NAME`, optional `EMAIL`, `PHONE`, `ORGANISATION`. ORD validation requires an email on `record_created.person` and every `record_modified` person, a `record_created.time`, and at least one of username/name/orcid on that person. When the export omits those fields (common in SURF), pass `--email` / `--person-name` / `--username` / `--orcid` / `--created-date`; **UDM values win when both are present** (unlike `--name` / `--description`, which override dataset packaging fields). Bare-name strings (`<SCIENTIST>Alice</SCIENTIST>`) still map the name field only. Validation failures for these gaps print a Hint naming the flag to pass.

---

## XML attribute handling

`etree_to_dict` encodes XML attributes with an `@` prefix and text content as `#text` when both attributes and text are present on the same element:

```xml
<AMOUNT units="g">1.5</AMOUNT>
```
```python
{'@units': 'g', '#text': '1.5'}
```

The `_text()` helper extracts `#text` from such dicts, falling back to a plain string. The `_parse_amount()` function additionally promotes `@units` to the unit string when no separate `<AMOUNT_UNIT>` element exists.

Text is stripped for every tag except those in `_RAW_TEXT_TAGS` (currently `MOLSTRUCTURE`), whose whitespace is significant; see [Molecule identifiers](#molecule-identifiers).

---

## What is not converted

ORD models conditions as a **single** `ReactionConditions` message per `Reaction` (not a list of condition groups). UDM can carry several `<CONDITION_GROUP>` siblings; this converter keeps only the first and logs a warning. Multi-stage detail is not written into `conditions.details` / `conditions_are_dynamic` today. Full policy: [Dropped CONDITION_GROUPs](#dropped-condition_groups-domain-review).

| UDM element / case | Reason |
|---|---|
| 2nd+ `<CONDITION_GROUP>` under `CONDITIONS` | ORD has one `ReactionConditions` per Reaction; converter uses the first group only; **later groups discarded with a log warning** (not copied to `details`) |
| `<RXNSTRUCTURE format="cdxml">` value | CDX binary embedded in XML; no ORD SMILES equivalent |
| `<ANALYSIS>`, `<SPECTRUM>` | ORD `Analysis` proto exists but not yet wired up |
| `<SCALE>` | No direct ORD equivalent |
| Parallel multi-stage condition narrative beyond the first group | Not flattened into `conditions.details` / `conditions_are_dynamic` |
| `<PRESSURE>` value without a unit attribute | Not mapped to `pressure.setpoint` (unit inference unreliable); raw value may appear in `conditions.details` |
| `<MODIFICATION_DATE>` nested structures | Plain strings and simple lists are handled; unusual nesting may vary |

---

## Converter policies for incomplete UDM (domain experts to review)

Policies applied when UDM is missing fields that ORD validation still requires. Intended for domain experts to review and challenge.

| Incomplete UDM pattern | ORD requirement | Converter policy |
|---|---|---|
| Elsevier-style DOI with `(…)` in the suffix; URL or `org/…` prefixes | `provenance.doi` must equal `parse_doi(doi)` | `parse_doi` keeps parenthetical suffixes; converter normalizes via `parse_doi` before storing |
| Products only as `REACTION/PRODUCT_ID` (or `VARIATION/PRODUCT_ID`), no `<PRODUCT>` block | ≥1 `ReactionOutcome` | Resolve IDs through `MOLECULES` into outcome products when no `PRODUCT` blocks exist |
| Reactants only as `REACTION/REACTANT_ID` (or `VARIATION/REACTANT_ID`), no role blocks | ≥1 reaction input | Resolve IDs through `MOLECULES` as `REACTANT` inputs when no role blocks exist |
| Free-text `PREPARATION` used as environment details | `ReactionEnvironment.type` required if message non-empty | Set `type=CUSTOM` with the free text in `details` |
| Role compound without `AMOUNT` / `SAMPLE_MASS` / `VOLUME` | Every input component needs an `Amount` | `UnmeasuredAmount` with `type=CUSTOM` and details `"amount not reported in UDM"` |
| `PRESSURE/exact` without `@unit` / `@units` | If setpoint `value` is set, `units` is required | Omit setpoint; record `UDM PRESSURE exact=… (unit omitted; not mapped to setpoint)` in `conditions.details` |
| Empty `MOLECULE/NAME` and no readable `MOLSTRUCTURE` | Identifier `value` must be non-empty | Use molecule `@ID` as `NAME` |
| Several `CONDITION_GROUP`s | One `ReactionConditions` message | **First group only**; log warning; 2nd+ groups dropped with no `details` summary — see [Dropped CONDITION_GROUPs](#dropped-condition_groups-domain-review) |

See also the user-facing summary in [`convert_udm_to_ord_guide.md`](convert_udm_to_ord_guide.md#converter-policies-for-incomplete-udm-domain-review).
