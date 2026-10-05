# UDM → ORD Field Mapping Specification

Conversion logic for `convert_udm_to_ord.py`.
Source format: UDM (Unified Data Model) v6.0.0 XML.
Target format: ORD `Dataset` protobuf (`.pbtxt` / `.pb`).

---

## Structural mapping

The UDM tree has a two-level reaction model: each `<REACTION>` contains one or more `<VARIATION>` elements, where a variation represents a single experimental run (different conditions, scale, outcome). **Each UDM `VARIATION` becomes one ORD `Reaction`.**

```text
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
| --- | --- | --- |
| `LEGAL/TITLE` | `Dataset.name` | Overridable with `--name` |
| `LEGAL/DOI` | `Dataset.description` | Formatted as `"UDM dataset DOI: <doi>"` |
| *(none)* | `Dataset.dataset_id`, `Reaction.reaction_id` | Left unset. Assigned when the dataset is submitted to ord-data |

---

## Reaction identifiers

Mapped from `<RXNSTRUCTURE>` elements on the parent `<REACTION>` (shared across all its variations).

| UDM `@format` attribute | ORD `ReactionIdentifier.type` |
| --- | --- |
| `rsmiles` | `REACTION_SMILES` |
| `rinchi` | `RINCHI` |
| `cdxml` | `CUSTOM`, details = `"cdxml"` |
| `rxn`, or attribute omitted | `CUSTOM`, details = `"rxn"` |
| *(other)* | `UNSPECIFIED` |

`format` is optional in the XSD and defaults to `rxn`. A `<RXNSTRUCTURE>` with no attributes is a plain string after parsing; that string is kept, not dropped. RXN text is not whitespace-stripped, same as `MOLSTRUCTURE`, because counts lines are column-sensitive.

---

## Inputs (`ReactionInput`)

UDM role blocks do not record addition order or grouping. Every `<REACTANT>`, `<REAGENT>`, `<CATALYST>`, and `<SOLVENT>` in a variation therefore becomes a component of one shared `ReactionInput` with map key `combined`. Each component keeps its own `reaction_role` and `amount`. A separate `ReactionInput` would claim a separate addition. The same molecule in two roles stays two components. When only bare `REACTANT_ID` references are available, those compounds share one `REACTANT_IDS` input instead. A role block whose molecule is missing from `MOLECULES` is still kept: `NAME` is the block's local name, or the molecule id when that name is empty, and the block's amount is preserved.

| UDM element | ORD field |
| --- | --- |
| `<REACTANT>` | role = `REACTANT` |
| `<REAGENT>` | role = `REAGENT` |
| `<CATALYST>` | role = `CATALYST` |
| `<SOLVENT>` | role = `SOLVENT` |
| `REACTION/REACTANT_ID` or `VARIATION/REACTANT_ID` | role = `REACTANT`. Used when no role blocks (Reaxys); resolved via `MOLECULES` |
| `MOLECULE/@MOL_ID` → `MOLECULES/MOLECULE/@ID` | `Compound.identifiers` (MOLBLOCK or NAME) |
| `AMOUNT`, `SAMPLE_MASS`, or `VOLUME` (+ optional unit / `AMOUNT_UNIT`) | `Compound.amount` (mass / moles / volume) |

### Amount parsing

`AMOUNT` may appear as:

- Plain text: `<AMOUNT>1.5</AMOUNT>` with `<AMOUNT_UNIT>g</AMOUNT_UNIT>`
- Inline attribute: `<AMOUNT unit="g">1.5</AMOUNT>` / `units="g"` (attribute-bearing dict from `etree_to_dict`)
- Combined string (SURF): `<AMOUNT>0.3000 mmol</AMOUNT>` — value and unit split from the text

Unit strings are matched case-insensitively through `ord_schema.units.UnitResolver`, and the resolved message has to be mass, moles, or volume. Anything else, including an unknown spelling, is stored as `UnmeasuredAmount` (`CUSTOM`) with the original numeric value and unit in `details`. `gr` is left unmeasured: the UDM `unitMass` enumeration defines it as grain, while `units.py` uses that spelling for gram.

Omitted units follow the XSD defaults: `AMOUNT` → `mol` (`molType`), `SAMPLE_MASS` → `g` (`massType`), `VOLUME` → `L`. Combined value+unit strings (SURF) are split before resolution. Non-finite values (`inf`, `nan`, overflow) are silently dropped. When no amount element is present at all, the converter sets `UnmeasuredAmount` (`CUSTOM`, details `"amount not reported in UDM"`) so ORD validation passes.

### Molecule identifiers

Molecules with a `<MOLSTRUCTURE>` element get a `MOLBLOCK` identifier. The `<MOLSTRUCTURE>` text content is used; any `@format` or other attributes are ignored. Molecules without a structure get a `NAME` identifier containing the molecule name (`CUSTOM` is unusable here, since ORD requires a `details` string alongside it).

MolBlock text is kept verbatim by `etree_to_dict` rather than stripped, because MolBlock lines are column-sensitive and its first line is a title that is often blank. A newline after the opening tag and the per-line indentation an XML formatter adds are both indistinguishable from MolBlock content, so `_normalize_molblock()` tries the text as recorded, without its leading newline, dedented, and both — keeping the first form RDKit can read.

A `<MOLSTRUCTURE>` that RDKit cannot read at all is not recorded as a `MOLBLOCK`, since ORD validation rejects an unreadable one and would fail the whole dataset; the molecule name is recorded as a `NAME` identifier instead and a warning is logged.

---

## Conditions (`ReactionConditions` / `ReactionSetup`)

Mapped from `VARIATION/CONDITIONS/CONDITION_GROUP`, or from `VARIATION/CONDITIONS` directly when there is no `CONDITION_GROUP` (SURF).

### Multiple CONDITION_GROUPs

ORD models conditions as **one** `ReactionConditions` per `Reaction` (`reaction.conditions = 4` in the proto — not `repeated`). UDM may list several sibling `<CONDITION_GROUP>` elements (common in Reaxys literature exports).

| Step | Behavior |
| --- | --- |
| One group | Map supported fields structurally; preserve unsupported fields and one-sided bounds in `conditions.details` |
| Two or more groups | Set `conditions_are_dynamic = true`; write every group as a labeled stage in `conditions.details` |
| Structured setpoints for dynamic conditions | Leave unset; selecting or merging one stage would misrepresent the source |

UDM documents multiple groups as a dynamic multi-stage profile rather than independent reactions. ORD's `conditions_are_dynamic` and catch-all `details` fields represent that case without duplicating outcomes or fabricating a static composite.

SURF nests reactants/products/conditions under `VARIATION/SECTION`; the converter promotes those children to variation level before mapping (existing top-level keys win).

| UDM element | ORD field | Notes |
| --- | --- | --- |
| `TEMPERATURE/@unit(s)` + exact or min/max | `conditions.temperature.setpoint` | A complete range maps to midpoint + precision. Unit spellings come from `ord_schema.units` (`degC`, `C`, `K`, …). Missing unit → Celsius (XSD default `degC`) |
| `PRESSURE/@unit(s)` + exact or min/max | `conditions.pressure.setpoint` | A complete range maps to midpoint + precision. Unit spellings come from `ord_schema.units` (`torr`, `bar`, `atm`, `psi`, `Pa`, `kPa`, `mmHg`, …). **No unit → setpoint omitted**; raw value appended to `conditions.details`. The XSD default is torr, and it is not applied: Reaxys unitless values are not one scale |
| `PRESSURE/ATMOSPHERE` | `conditions.pressure.atmosphere.type` | `air`, `n2`/`nitrogen`, `ar`/`argon`, `o2`/`oxygen`, `h2`/`hydrogen`, `co`, `co2` |
| `STIRRING` (text) | `conditions.stirring.details` + `.type` + `.rate.rpm` | See [Stirring](#stirring) |
| `REFLUX` | `conditions.reflux` | True when value is `true`, `yes`, or `1` |
| `PH` exact or min/max; `<PH>7.0</PH>` | `conditions.ph` | Complete ranges use the midpoint; the range remains in `details` because ORD pH has no precision field |
| `PREPARATION` | `setup.environment.type`, or `notes.procedure_details` | Keyword match (`fume hood`, `bench top`, `glove box`, `glove bag`) sets the environment. Any other text, from `CONDITIONS` or `CONDITION_GROUP`, is procedure notes |
| `VESSEL/VESSEL_TYPE` | `setup.vessel.type` | `round bottom flask`/`rbf`, `vial`, `well plate`, `tube`, `microwave vial`, `nmr tube`, `pressure flask`, `pressure reactor` |
| `VESSEL/DETAILS` | `setup.vessel.details` | |

All numeric condition values guard against non-finite floats (`inf`, `nan`); value and units are assigned together. A lone `min` or `max` is retained as a labeled bound (`>=` / `<=`) in `details`, never fabricated as an exact value.

### Stirring

UDM records stirring as free text (e.g. `600 rpm magnetic stir bar`), but ORD requires a `StirringMethodType`, so the type is inferred from keywords in that text (first match wins) and the full text is kept in `details`:

| Keyword in `STIRRING` | ORD `stirring.type` |
| --- | --- |
| `stir bar`, `magnetic` | `STIR_BAR` |
| `overhead` | `OVERHEAD_MIXER` |
| `agitat` | `AGITATION` |
| `ball mill` | `BALL_MILLING` |
| `sonicat` | `SONICATION` |
| `unstirred`, `not stirred`, `none` | `NONE` |
| *(no match)* | `CUSTOM` (legal because `details` is populated) |

A `<N> rpm` substring additionally populates `conditions.stirring.rate.rpm`.

---

## Outcomes (`ReactionOutcome`)

An outcome is only created when the variation has a parseable `<DURATION>` (dict form with `exact` child) or at least one `<PRODUCT>`. Input-only variations produce zero outcomes.

| UDM element | ORD field | Notes |
| --- | --- | --- |
| `DURATION` / `CONDITIONS/.../TIME` exact or min/max | `outcome.reaction_time` | Complete range maps to midpoint + precision. Unit spellings come from `ord_schema.units` (`hr`, `min`, `s`, `d`, …). Missing unit → hour (XSD default `hr`) |
| `PRODUCT/MOLECULE/@MOL_ID` | `outcome.products[].identifiers` | Via molecule lookup |
| `REACTION/PRODUCT_ID` or `VARIATION/PRODUCT_ID` | `outcome.products[].identifiers` | Used when no `PRODUCT` blocks (Reaxys); resolved via `MOLECULES` |
| `PRODUCT/YIELD` exact or min/max; `<YIELD>85</YIELD>` | `outcome.products[].measurements[].percentage` | Complete range maps to midpoint + precision; lone bound goes to measurement `details`; non-finite values skipped |

---

## Notes and observations

| UDM element | ORD field | Notes |
| --- | --- | --- |
| `VARIATION/PROCEDURE` | `notes.procedure_details` | |
| `CONDITIONS/PREPARATION` or `CONDITION_GROUP/PREPARATION`, when not an environment keyword | `notes.procedure_details` | Environment keywords stay on `setup.environment` only |
| `VARIATION/COMMENT` | `observations[0].comment` | |

---

## Provenance (`ReactionProvenance`)

CLI precedence is asymmetric (see the user guide): `--name` / `--description` **override** UDM dataset metadata. When UDM has a `SCIENTIST`, that person is `experimenter` and `record_created.person`; CLI username, name, and ORCID are not merged onto them, and `--email` fills only a missing scientist email on `record_created` / `record_modified`. When UDM has no scientist, CLI depositor flags populate `record_created.person` and `experimenter` stays unset. `--created-date` still fills only a missing `CREATION_DATE`.

| UDM element | ORD field | Notes |
| --- | --- | --- |
| `LEGAL/PRODUCER` | `provenance.experimenter.organization` + `record_created.person.organization` | |
| `VARIATION/SCIENTIST` (bare string or `NAME` / `EMAIL`) | `experimenter`, and the same name/email on `record_created.person` | The experimenter is this scientist only. CLI username, name, and ORCID are not copied onto that person |
| *(no SCIENTIST)* | `record_created.person` from `--username` / `--person-name` / `--orcid` / `--email` | `experimenter` name, username, ORCID, and email stay unset. Literature exports must not record the depositor as the person who ran the reaction |
| `--email` when SCIENTIST has no email | `record_created.person.email` and `record_modified[].person.email` only | ORD requires that email. It is not written onto `experimenter` |
| `--username`, `--person-name`, `--orcid` | `record_created.person` and `record_modified` when UDM has no SCIENTIST | Depositor identity. `--person-name` is distinct from `--name` (dataset title) |
| `--created-date` | `provenance.record_created.time.value` when UDM has no `CREATION_DATE` | Gap-fill only; required by ORD validation (`RecordEvent.time`) |
| `LEGAL/DOI` | `provenance.doi` | Overridden by variation-level or reaction-level citation DOI if present |
| `VARIATION/CITATION/@CIT_ID` or `VARIATION/@CIT_ID` → `CITATIONS/CITATION/@ID/DOI` | `provenance.doi` | Per-variation citation lookup (SURF uses the attribute form). Wins over a reaction-level DOI |
| `REACTION/CITATIONS/CITATION/DOI` | `provenance.doi` | Reaction-level legacy path. Used only when the variation did not resolve a DOI; still overrides `LEGAL/DOI` |
| `REACTION/CITATIONS/CITATION/PATENT_NUMBER` | `provenance.patent` | |
| `VARIATION/CREATION_DATE` | `provenance.record_created.time.value` | Wins over `--created-date` |
| `VARIATION/MODIFICATION_DATE` | `provenance.record_modified[].time.value` | Plain string and list both supported |
| `REACTION/ORGANISATIONS[0]/ORGANISATION/ADDRESS` | `provenance.city` | |
| *(always)* | `provenance.is_mined = false` | |

### Scientist / email / created date

UDM v6 models `SCIENTIST` like `AUTHOR`: required `NAME`, optional `EMAIL`, `PHONE`, `ORGANISATION`. That person is the experimenter. ORD validation requires an email on `record_created.person` and every `record_modified` person, a `record_created.time`, and at least one of username/name/orcid on that person. When the export has a scientist, `record_created.person` is that scientist, and `--email` fills only a missing email. When the export has no scientist (common in SURF and literature exports), the CLI depositor is `record_created.person` and the experimenter is left unset. `--name` / `--description` still override dataset packaging fields. Bare-name strings (`<SCIENTIST>Alice</SCIENTIST>`) still map the name field only. Validation failures for these gaps print a Hint naming the flag to pass.

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

Text is stripped for every tag except those in `_RAW_TEXT_TAGS` (`MOLSTRUCTURE` and `RXNSTRUCTURE`), whose whitespace is significant; see [Molecule identifiers](#molecule-identifiers).

With `--include-udm-xml`, each ORD reaction stores its source `<REACTION>` element as `provenance.reaction_metadata["udm_reaction_xml"]` and shared document context (`UDM_VERSION`, `LEGAL`, `ORGANISATIONS`, `CITATIONS`, etc.) as `"udm_parent_xml"`. `REACTIONS` and the potentially large `MOLECULES` lookup are excluded from parent context. This is opt-in because source XML may have a different license.

---

## What is not converted

| UDM element / case | Reason |
| --- | --- |
| `<RXNSTRUCTURE format="cdxml">` value | CDX binary embedded in XML; no ORD SMILES equivalent |
| `<ANALYSIS>`, `<SPECTRUM>` | ORD `Analysis` proto exists but not yet wired up |
| `<SCALE>` | No direct ORD equivalent |
| `<PRESSURE>` value without a unit attribute | Not mapped to `pressure.setpoint`. The XSD default is torr, but Reaxys unitless values span more than one scale (about 760, and also ~2 and ~4.5e6), so the raw value stays in `conditions.details` |
| `<MODIFICATION_DATE>` nested structures | Plain strings and simple lists are handled; unusual nesting may vary |

---

## Converter policies for incomplete UDM (domain experts to review)

Policies applied when UDM is missing fields that ORD validation still requires. Intended for domain experts to review and challenge.

| Incomplete UDM pattern | ORD requirement | Converter policy |
| --- | --- | --- |
| Elsevier-style DOI with balanced `(…)` in the suffix; URL or `org/…` prefixes; wrappers like `(doi:10.…)` | `provenance.doi` must equal `parse_doi(doi)` | `parse_doi` keeps balanced parenthetical suffixes and trims an unmatched trailing `)` from the regex match; converter normalizes before storing |
| Products only as `REACTION/PRODUCT_ID` (or `VARIATION/PRODUCT_ID`), no `<PRODUCT>` block | ≥1 `ReactionOutcome` | Resolve IDs through `MOLECULES` into outcome products when no `PRODUCT` blocks exist |
| `<REACTANT>` / `<REAGENT>` / `<CATALYST>` / `<SOLVENT>` with no addition-order field | One `ReactionInput` | Components of one shared `combined` input; each block keeps its own role and amount |
| Reactants only as `REACTION/REACTANT_ID` (or `VARIATION/REACTANT_ID`), no role blocks | ≥1 reaction input | Resolve all IDs through `MOLECULES` as components of one shared `REACTANT_IDS` input |
| `REACTION` has no `VARIATION` | Preserve recoverable reaction-level data | Emit one ORD reaction using reaction-level identifiers (unattributed `<RXNSTRUCTURE>` → `rxn`) and `REACTANT_ID` / `PRODUCT_ID` fallbacks. A non-SMILES identifier-only record may require `--no-validate` |
| Free-text `PREPARATION` | Procedure text is not an environment | Write it to `notes.procedure_details`. Set `setup.environment` only when the text is a known environment keyword |
| Role compound without `AMOUNT` / `SAMPLE_MASS` / `VOLUME` | Every input component needs an `Amount` | `UnmeasuredAmount` with `type=CUSTOM` and details `"amount not reported in UDM"` |
| `PRESSURE/exact` without `@unit` / `@units` | If setpoint `value` is set, `units` is required | Omit setpoint. The XSD default (torr) is not applied, because unitless Reaxys values are not one scale. Record `UDM PRESSURE value=… (unit omitted; not mapped to setpoint)` in `conditions.details` |
| Empty `MOLECULE/NAME` and no readable `MOLSTRUCTURE` | Identifier `value` must be non-empty | Use molecule `@ID` as `NAME` |
| Several `CONDITION_GROUP`s | One `ReactionConditions` message | Set `conditions_are_dynamic`; summarize every group as a stage in `details` — see [Multiple CONDITION_GROUPs](#multiple-condition_groups) |

See also the user-facing summary in [`convert_udm_to_ord_guide.md`](convert_udm_to_ord_guide.md#converter-policies-for-incomplete-udm-domain-review).
