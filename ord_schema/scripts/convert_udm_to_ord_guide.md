# UDM → ORD Converter: User Guide

`convert_udm_to_ord.py` converts a UDM (Unified Data Model) v6.0.0 XML file exported from an ELN into an ORD `Dataset` protobuf file (`.pbtxt` or `.pb`).

Each UDM `<VARIATION>` element becomes one ORD `Reaction`.

---

## Requirements

```bash
uv sync --extra tests   # installs ord-schema with all dependencies
```

---

## Basic usage

```bash
uv run python ord_schema/scripts/convert_udm_to_ord.py \
    --input my_reactions.xml \
    --output my_dataset.pbtxt
```

If `--output` is omitted, the output filename is derived from the `<TITLE>` in the UDM file (unsafe characters are replaced with `_`). Falls back to `ord_dataset.pbtxt` when no title is present.

---

## Options

| Flag | Default | Description |
|---|---|---|
| `--input FILE` | _(required)_ | Path to the UDM v6.0.0 XML file |
| `--output FILE` | `<title>.pbtxt` | Output path; suffix determines format (`.pbtxt` = text proto, `.pb` = binary) |
| `--name TEXT` | UDM `LEGAL/TITLE` | **Override** dataset name (`Dataset.name`) |
| `--description TEXT` | DOI-derived text | **Override** dataset description (`Dataset.description`) |
| `--username TEXT` | _(none)_ | Depositor username → provenance `Person.username` (no UDM source) |
| `--person-name TEXT` | UDM `SCIENTIST/NAME` | **Fill gap** for depositor display name. Distinct from `--name` (dataset title) |
| `--orcid TEXT` | _(none)_ | Depositor ORCID iD → provenance `Person.orcid` (no UDM source) |
| `--email ADDRESS` | UDM `SCIENTIST/EMAIL` | **Fill gap** for depositor email; ORD requires email on `record_created` |
| `--created-date TEXT` | UDM `CREATION_DATE` | **Fill gap** for `record_created.time`; ORD requires a time |
| `--no-validate` | off | Skip ORD schema validation; useful for large batch jobs or partially complete data |

### CLI precedence (asymmetric on purpose)

Dataset packaging flags and reaction-provenance flags behave differently when both UDM and CLI supply a value:

| Kind | Flags | Both present? | Why |
|---|---|---|---|
| Dataset packaging | `--name`, `--description` | **CLI wins** (overrides UDM `TITLE` / DOI-derived text) | Naming the ORD Dataset is a deposit/packaging choice; same idea as other ORD scripts |
| Provenance gap-fill | `--email`, `--person-name`, `--created-date` | **UDM wins**; CLI used only when UDM omits the field | Scientist contact and creation time are recorded experiment facts — do not overwrite |
| Provenance only-on-CLI | `--username`, `--orcid` | CLI applied (no UDM equivalent in this converter) | |

`--name` is **not** a person identity flag. ORD still needs at least one of `--person-name`, `--username`, or `--orcid` (or UDM `SCIENTIST/NAME`) on `record_created.person`.

---

## Examples

**Convert and validate:**
```bash
uv run python ord_schema/scripts/convert_udm_to_ord.py \
    --input screen.xml \
    --output screen.pbtxt
```

**Fill depositor provenance when UDM has no SCIENTIST / CREATION_DATE (e.g. SURF exports):**
```bash
uv run python ord_schema/scripts/convert_udm_to_ord.py \
    --input surf_export.xml \
    --output surf_export.pbtxt \
    --name "SURF export" \
    --description "Literature Minisci screen" \
    --username ada \
    --person-name "Ada Lovelace" \
    --orcid 0000-0002-1825-0097 \
    --email ada@example.com \
    --created-date 2024-01-15
```

SURF files often nest chemistry under `VARIATION/SECTION` and write amounts as `0.3000 mmol`; the converter handles both. Empty `<LEGAL />` still needs `--name` / `--description` for ORD dataset validation.

**Convert without validation (faster, for drafts):**
```bash
uv run python ord_schema/scripts/convert_udm_to_ord.py \
    --input screen.xml \
    --output screen.pbtxt \
    --no-validate
```

**Override dataset name and write binary format:**
```bash
uv run python ord_schema/scripts/convert_udm_to_ord.py \
    --input screen.xml \
    --output screen.pb \
    --name "Suzuki Screen 2024-Q1" \
    --description "DOI: 10.1234/xyz"
```

**Batch convert a directory of UDM files (bash):**
```bash
for f in udm_exports/*.xml; do
    uv run python ord_schema/scripts/convert_udm_to_ord.py \
        --input "$f" \
        --output "ord_out/$(basename "${f%.xml}").pbtxt" \
        --no-validate
done
```

---

## Output format

The converter writes a standard ORD `Dataset` protobuf. Each reaction inside it has:

- A canonical `reaction_id` (auto-assigned as `ord-<sha256>`)
- Inputs keyed by `<MOL_ID>_<ROLE>` (e.g., `MOL-1_REACTANT`, `MOL-1_SOLVENT`)
- Conditions, outcomes, notes, and provenance populated where UDM data is present

To inspect the output:
```bash
# Print the text proto
cat my_dataset.pbtxt

# Load in Python
from ord_schema import datasets
ds = datasets.load_dataset("my_dataset.pbtxt", as_dataset=True)
print(len(ds.reactions), "reactions")
```

---

## Validation

By default the converter runs `validations.validate_datasets` before writing. If the dataset fails validation, the output file is **not** written and the converter exits with code 1:

```
ERROR - Validation failed (use --no-validate to write anyway):
  reaction[0]: ...
```

Common reasons for validation failure:
- `SCIENTIST/EMAIL` missing and `--email` not supplied (ORD requires an email on every provenance record)
- `CREATION_DATE` missing and `--created-date` not supplied (ORD requires `record_created.time`)
- No scientist identity and none of `--person-name` / `--username` / `--orcid` supplied
- Missing required fields (e.g., no amount on a component)
- Amount units unrecognised (logged as a warning during conversion)
- Reaction has no inputs or outcomes

On these provenance gaps the converter prints a short Hint naming the CLI flag to pass. Use `--no-validate` to write the file anyway and validate it separately later.

**Validate after the fact — CLI (single file):**
```bash
uv run python ord_schema/scripts/validate_dataset.py \
    --input_pattern "my_dataset.pbtxt"
```

**Validate multiple files at once:**
```bash
# All .pbtxt files in a directory
uv run python ord_schema/scripts/validate_dataset.py \
    --input_pattern "ord_out/*.pbtxt"

# With parallel workers for large batches
uv run python ord_schema/scripts/validate_dataset.py \
    --input_pattern "ord_out/*.pbtxt" \
    --n_jobs 4

# Also check that reaction_id / dataset_id are well-formed ord- hashes
uv run python ord_schema/scripts/validate_dataset.py \
    --input_pattern "ord_out/*.pbtxt" \
    --validate_ids
```

**Validate after the fact — Python:**
```python
from ord_schema import datasets, validations

ds = datasets.load_dataset("my_dataset.pbtxt", as_dataset=True)
validations.validate_datasets({"my_dataset": ds})  # raises ValidationError if invalid
```

---

## License notice

The ORD repository uses the **CC-BY-SA** license. Do not submit converted data to `ord-data` unless you hold the authority to relicense the source data under CC-BY-SA. This warning is logged at the start of every conversion run.

---

## Warnings and what they mean

Converter and RDKit may print warnings during conversion. Most do **not** stop the run; treat them as data-quality notes.

| Warning message | Source | Cause | Safe to ignore? |
|---|---|---|---|
| `could not find number of expected rings. Switching to an approximate ring finding algorithm.` | RDKit | Unusual / incomplete MolBlock ring perception | Usually yes, if the molecule still parses and ORD validation passes |
| `molecule is tagged as 2D, but at least one Z coordinate is not zero. Marking the mol as 3D.` | RDKit | MolBlock dimension flag disagrees with coordinates | Usually yes; RDKit retags and continues |
| `Element 'H+' not found` / `Post-condition Violation` … `Failed Expression: anum > -1` | RDKit | MolBlock uses a non-element symbol (often proton `H+`) RDKit cannot look up | Usually yes for conversion. Structure often fails parse; converter then falls back to `NAME` (see MOLSTRUCTURE row below). Common in SPRESI-style exports |
| `Cannot convert '…' to unsigned int on line 4` | RDKit | Malformed MolBlock counts line (bad atom/bond counts) | Same as above — often followed by NAME fallback |
| `Explicit valence for atom # … is greater than permitted` | RDKit | Drawn valence / charge inconsistent with RDKit rules | Usually yes if conversion finishes; structure may still parse or fall back to NAME |
| `atom … has specified valence (…) smaller than the drawn valence …` | RDKit | Valence annotation disagrees with bonds in the MolBlock | Usually yes; RDKit noise on legacy drawings |
| `Warning: ambiguous stereochemistry - overlapping neighbors … ignored` | RDKit | Stereo wedge/hash geometry is ambiguous | Usually yes; stereo may be dropped for that atom |
| `WARNING: not removing hydrogen atom without neighbors` (and similar H cleanup messages) | RDKit | Odd explicit hydrogens in the MolBlock | Usually yes |
| `Multiple CONDITION_GROUPs found; using the first.` | Converter | UDM has more than one `<CONDITION_GROUP>` under `CONDITIONS` | **Not silent data loss — review it.** Only the **first** group is mapped to ORD. Count how often this appears in the log (common in Reaxys). See [Dropped CONDITION_GROUPs](#dropped-condition_groups-domain-review) |
| `AMOUNT element has attributes but no text content; skipping.` | Converter | `<AMOUNT units="g"/>` with no value | No — fix the source UDM |
| `Non-finite AMOUNT value ...; skipping.` | Converter | Amount is `inf`, `nan`, or overflowing float | No — check the source value |
| `Molecule ... not found in MOLECULES lookup; skipping.` | Converter | `MOL_ID` has no matching `<MOLECULE>` | No — fix dangling references |
| `MOLSTRUCTURE for ... is not a readable MolBlock; recording NAME instead.` | Converter | RDKit cannot read the structure (often after `H+` / valence / counts-line errors above) | Often yes for conversion; **structure is lost**, molecule `@ID` or name kept as `NAME`. Same mol ID may warn many times if reused across reactions |
| `Reaction has no VARIATION elements; skipping.` | Converter | `<REACTION>` with no `<VARIATION>` children | Expected for templates; no ORD Reaction emitted |

To capture warnings (and errors) in a file, redirect **stderr** as well as stdout, e.g. `&> convert.log` or `2>&1 | tee convert.log`.

---

## Troubleshooting

**`sys.exit(1)` / "Input file not found"**
Check that `--input` points to an existing file.

**"Input file is not UDM format"**
The XML root element must be `<UDM>`. Check that you have a UDM v6.0.0 export, not a CambridgeSoft or other format.

**"Validation failed" with many errors after conversion**
Pass `--no-validate` to write the file, then load it and inspect the reactions manually. Consider filing a bug with a sanitised sample UDM file.

**Output filename contains `_` where you expected `/`**
Path separators in the UDM `<TITLE>` are replaced with `_` to avoid filesystem errors. Use `--output` to specify an exact path.

---

## Field mapping reference

See [`udm_to_ord_mapping.md`](udm_to_ord_mapping.md) for the full UDM → ORD field-by-field mapping table.

---

## Converter policies for incomplete UDM (domain experts to review)

Literature / ELN exports (especially Reaxys) often omit fields ORD validation requires. The converter applies the policies below so datasets can validate **without** inventing chemically wrong values. Domain experts: please comment if a policy should change.

| Situation | Converter policy | Rationale |
|---|---|---|
| DOI like `10.1016/S0022-328X(00)99569-X` (parentheses in suffix) | Keep the full DOI (`parse_doi` allows `(…)`) | Trimmed forms are regex artifacts and often do not resolve; the published DOI is kept |
| DOI prefixed with junk (`org/10.1016/…`, URL wrappers) | Normalize via `parse_doi` before writing `provenance.doi` | ORD requires the stored DOI to equal the parsed form |
| `VARIATION` has reagents but no `<PRODUCT>`; `REACTION` has `<PRODUCT_ID>` | Resolve `PRODUCT_ID` → `MOLECULES` and create an outcome product | ORD requires ≥1 outcome; Reaxys stores products as IDs at reaction level |
| `VARIATION` has no role compounds; `REACTION` has `<REACTANT_ID>` | Resolve `REACTANT_ID` → `MOLECULES` as `REACTANT` inputs (unmeasured amount) | ORD requires ≥1 input; Reaxys often lists reactants only as IDs |
| Free-text `<PREPARATION>` (not a known env keyword) | Set `setup.environment.type=CUSTOM` and put the text in `.details` | ORD requires `type` whenever environment fields are set |
| Input compound has no `<AMOUNT>` | Set `amount.unmeasured` with `type=CUSTOM`, `details="amount not reported in UDM"` | ORD requires an amount on every input; unmeasured is honest vs inventing a number |
| `<PRESSURE><exact>…</exact></PRESSURE>` with **no unit** | Do **not** set `pressure.setpoint`; append raw value to `conditions.details` | Reaxys mixes scales (≈760 torr vs large Pa-like numbers); guessing units is unreliable |
| `<MOLECULE>` has empty `<NAME/>` and no usable `MOLSTRUCTURE` | Use the molecule `@ID` string as a `NAME` identifier | Empty identifier values fail validation; ID preserves a stable handle |
| Multiple `<CONDITION_GROUP>` siblings | Use the **first only**; log `Multiple CONDITION_GROUPs found; using the first.` once per affected variation. 2nd+ groups are **not** written to ORD (not even into `conditions.details`) | ORD has a single `ReactionConditions` per Reaction — see [Dropped CONDITION_GROUPs](#dropped-condition_groups-domain-review) |

### Dropped CONDITION_GROUPs (domain experts to review)

UDM allows several `<CONDITION_GROUP>` children under one `<CONDITIONS>` (e.g. alternate pressure/temperature snapshots in Reaxys). ORD does **not**: each `Reaction` has one `ReactionConditions` message.

**What the converter does today**

1. Reads the first `<CONDITION_GROUP>` only (temperature, pressure, stirring, …).
2. Logs a warning for that variation: `Multiple CONDITION_GROUPs found; using the first.`
3. **Discards** every later group — contents are not copied into `conditions.details`, `conditions_are_dynamic`, outcomes, or notes.

**How to review a conversion log**

```bash
# How many variations lost 2nd+ groups?
rg -c 'Multiple CONDITION_GROUPs found' convert.log
```

Each hit ≈ one ORD reaction that may be missing alternate condition sets from the source UDM. Open the matching UDM `<VARIATION>` / `<CONDITIONS>` and check whether the dropped groups were duplicates, alternatives, or multi-stage steps you care about.

**Possible future improvements** (not implemented): merge groups; set `conditions_are_dynamic = true` and summarize extras in `conditions.details`; or emit one ORD reaction per group. Domain experts: comment which behavior you want.

Related: [`udm_to_ord_mapping.md` — Dropped CONDITION_GROUPs](udm_to_ord_mapping.md#dropped-condition_groups-domain-review) · [full policy table](udm_to_ord_mapping.md#converter-policies-for-incomplete-udm-domain-review).
