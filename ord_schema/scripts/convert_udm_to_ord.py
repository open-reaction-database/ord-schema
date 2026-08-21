# Copyright 2025 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Converts a UDM v6.0.0 XML file to an ORD Dataset (.pbtxt or .pb).

Each UDM VARIATION becomes a separate ORD Reaction.

Example usage:
    python convert_udm_to_ord.py --input my_dataset.xml --output my_dataset.pbtxt \
        --email me@example.com --person-name "Ada Lovelace" --created-date 2024-01-15
"""

import argparse
import contextlib
import math
import pathlib
import re
import sys
import textwrap
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import cast

from ord_schema import message_helpers, updates, validations
from ord_schema.logging import get_logger
from ord_schema.proto import dataset_pb2, reaction_pb2

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Unit normalisation maps
# ---------------------------------------------------------------------------

_MASS_UNITS: dict[str, reaction_pb2.Mass.MassUnit] = {
    "kg": reaction_pb2.Mass.KILOGRAM,
    "kilogram": reaction_pb2.Mass.KILOGRAM,
    "kilograms": reaction_pb2.Mass.KILOGRAM,
    "g": reaction_pb2.Mass.GRAM,
    "gram": reaction_pb2.Mass.GRAM,
    "grams": reaction_pb2.Mass.GRAM,
    "mg": reaction_pb2.Mass.MILLIGRAM,
    "milligram": reaction_pb2.Mass.MILLIGRAM,
    "milligrams": reaction_pb2.Mass.MILLIGRAM,
    "ug": reaction_pb2.Mass.MICROGRAM,
    "µg": reaction_pb2.Mass.MICROGRAM,
    "microgram": reaction_pb2.Mass.MICROGRAM,
    "micrograms": reaction_pb2.Mass.MICROGRAM,
}

_MOLES_UNITS: dict[str, reaction_pb2.Moles.MolesUnit] = {
    "mol": reaction_pb2.Moles.MOLE,
    "mole": reaction_pb2.Moles.MOLE,
    "moles": reaction_pb2.Moles.MOLE,
    "mmol": reaction_pb2.Moles.MILLIMOLE,
    "millimol": reaction_pb2.Moles.MILLIMOLE,
    "millimole": reaction_pb2.Moles.MILLIMOLE,
    "millimoles": reaction_pb2.Moles.MILLIMOLE,
    "umol": reaction_pb2.Moles.MICROMOLE,
    "µmol": reaction_pb2.Moles.MICROMOLE,
    "micromole": reaction_pb2.Moles.MICROMOLE,
    "micromoles": reaction_pb2.Moles.MICROMOLE,
    "nmol": reaction_pb2.Moles.NANOMOLE,
    "nanomole": reaction_pb2.Moles.NANOMOLE,
    "nanomoles": reaction_pb2.Moles.NANOMOLE,
}

_TIME_UNITS: dict[str, reaction_pb2.Time.TimeUnit] = {
    "h": reaction_pb2.Time.HOUR,
    "hr": reaction_pb2.Time.HOUR,
    "hour": reaction_pb2.Time.HOUR,
    "hours": reaction_pb2.Time.HOUR,
    "min": reaction_pb2.Time.MINUTE,
    "minute": reaction_pb2.Time.MINUTE,
    "minutes": reaction_pb2.Time.MINUTE,
    "s": reaction_pb2.Time.SECOND,
    "sec": reaction_pb2.Time.SECOND,
    "second": reaction_pb2.Time.SECOND,
    "seconds": reaction_pb2.Time.SECOND,
}

_VOLUME_UNITS: dict[str, reaction_pb2.Volume.VolumeUnit] = {
    "l": reaction_pb2.Volume.LITER,
    "liter": reaction_pb2.Volume.LITER,
    "liters": reaction_pb2.Volume.LITER,
    "litre": reaction_pb2.Volume.LITER,
    "litres": reaction_pb2.Volume.LITER,
    "ml": reaction_pb2.Volume.MILLILITER,
    "milliliter": reaction_pb2.Volume.MILLILITER,
    "milliliters": reaction_pb2.Volume.MILLILITER,
    "millilitre": reaction_pb2.Volume.MILLILITER,
    "millilitres": reaction_pb2.Volume.MILLILITER,
    "ul": reaction_pb2.Volume.MICROLITER,
    "µl": reaction_pb2.Volume.MICROLITER,
    "microliter": reaction_pb2.Volume.MICROLITER,
    "microliters": reaction_pb2.Volume.MICROLITER,
    "nl": reaction_pb2.Volume.NANOLITER,
    "nanoliter": reaction_pb2.Volume.NANOLITER,
    "nanoliters": reaction_pb2.Volume.NANOLITER,
}

_TEMP_UNITS: dict[str, reaction_pb2.Temperature.TemperatureUnit] = {
    "c": reaction_pb2.Temperature.CELSIUS,
    "celsius": reaction_pb2.Temperature.CELSIUS,
    "°c": reaction_pb2.Temperature.CELSIUS,
    "degc": reaction_pb2.Temperature.CELSIUS,
    "f": reaction_pb2.Temperature.FAHRENHEIT,
    "fahrenheit": reaction_pb2.Temperature.FAHRENHEIT,
    "°f": reaction_pb2.Temperature.FAHRENHEIT,
    "degf": reaction_pb2.Temperature.FAHRENHEIT,
    "k": reaction_pb2.Temperature.KELVIN,
    "kelvin": reaction_pb2.Temperature.KELVIN,
}

_PRESSURE_UNITS: dict[str, reaction_pb2.Pressure.PressureUnit] = {
    "bar": reaction_pb2.Pressure.BAR,
    "atm": reaction_pb2.Pressure.ATMOSPHERE,
    "atmosphere": reaction_pb2.Pressure.ATMOSPHERE,
    "atmospheres": reaction_pb2.Pressure.ATMOSPHERE,
    "psi": reaction_pb2.Pressure.PSI,
    "kpsi": reaction_pb2.Pressure.KPSI,
    "torr": reaction_pb2.Pressure.TORR,
}

_ATMOSPHERE_TYPES: dict[
    str, reaction_pb2.PressureConditions.Atmosphere.AtmosphereType
] = {
    "air": reaction_pb2.PressureConditions.Atmosphere.AIR,
    "n2": reaction_pb2.PressureConditions.Atmosphere.NITROGEN,
    "nitrogen": reaction_pb2.PressureConditions.Atmosphere.NITROGEN,
    "ar": reaction_pb2.PressureConditions.Atmosphere.ARGON,
    "argon": reaction_pb2.PressureConditions.Atmosphere.ARGON,
    "o2": reaction_pb2.PressureConditions.Atmosphere.OXYGEN,
    "oxygen": reaction_pb2.PressureConditions.Atmosphere.OXYGEN,
    "h2": reaction_pb2.PressureConditions.Atmosphere.HYDROGEN,
    "hydrogen": reaction_pb2.PressureConditions.Atmosphere.HYDROGEN,
    "co": reaction_pb2.PressureConditions.Atmosphere.CARBON_MONOXIDE,
    "co2": reaction_pb2.PressureConditions.Atmosphere.CARBON_DIOXIDE,
}

# Tags whose text must reach _normalize_molblock() unstripped; see etree_to_dict().
_RAW_TEXT_TAGS = frozenset({"MOLSTRUCTURE"})

# UDM STIRRING is free text ("600 rpm magnetic stir bar"), so the method type is
# inferred from keywords; order matters because the first match wins.
_STIRRING_TYPE_KEYWORDS: tuple[
    tuple[str, reaction_pb2.StirringConditions.StirringMethodType], ...
] = (
    ("stir bar", reaction_pb2.StirringConditions.STIR_BAR),
    ("magnetic", reaction_pb2.StirringConditions.STIR_BAR),
    ("overhead", reaction_pb2.StirringConditions.OVERHEAD_MIXER),
    ("agitat", reaction_pb2.StirringConditions.AGITATION),
    ("ball mill", reaction_pb2.StirringConditions.BALL_MILLING),
    ("sonicat", reaction_pb2.StirringConditions.SONICATION),
    ("unstirred", reaction_pb2.StirringConditions.NONE),
    ("not stirred", reaction_pb2.StirringConditions.NONE),
    ("none", reaction_pb2.StirringConditions.NONE),
)

_ENVIRONMENT_TYPES: dict[
    str, reaction_pb2.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
] = {
    "fume hood": reaction_pb2.ReactionSetup.ReactionEnvironment.FUME_HOOD,
    "fume_hood": reaction_pb2.ReactionSetup.ReactionEnvironment.FUME_HOOD,
    "bench": reaction_pb2.ReactionSetup.ReactionEnvironment.BENCH_TOP,
    "bench top": reaction_pb2.ReactionSetup.ReactionEnvironment.BENCH_TOP,
    "bench_top": reaction_pb2.ReactionSetup.ReactionEnvironment.BENCH_TOP,
    "glove box": reaction_pb2.ReactionSetup.ReactionEnvironment.GLOVE_BOX,
    "glove_box": reaction_pb2.ReactionSetup.ReactionEnvironment.GLOVE_BOX,
    "glovebox": reaction_pb2.ReactionSetup.ReactionEnvironment.GLOVE_BOX,
    "glove bag": reaction_pb2.ReactionSetup.ReactionEnvironment.GLOVE_BAG,
    "glove_bag": reaction_pb2.ReactionSetup.ReactionEnvironment.GLOVE_BAG,
    "glovebag": reaction_pb2.ReactionSetup.ReactionEnvironment.GLOVE_BAG,
}

_VESSEL_TYPES: dict[str, reaction_pb2.Vessel.VesselType] = {
    "round bottom flask": reaction_pb2.Vessel.ROUND_BOTTOM_FLASK,
    "round_bottom_flask": reaction_pb2.Vessel.ROUND_BOTTOM_FLASK,
    "rbf": reaction_pb2.Vessel.ROUND_BOTTOM_FLASK,
    "vial": reaction_pb2.Vessel.VIAL,
    "well plate": reaction_pb2.Vessel.WELL_PLATE,
    "well_plate": reaction_pb2.Vessel.WELL_PLATE,
    "microwave vial": reaction_pb2.Vessel.MICROWAVE_VIAL,
    "microwave_vial": reaction_pb2.Vessel.MICROWAVE_VIAL,
    "tube": reaction_pb2.Vessel.TUBE,
    "nmr tube": reaction_pb2.Vessel.NMR_TUBE,
    "nmr_tube": reaction_pb2.Vessel.NMR_TUBE,
    "pressure flask": reaction_pb2.Vessel.PRESSURE_FLASK,
    "pressure_flask": reaction_pb2.Vessel.PRESSURE_FLASK,
    "pressure reactor": reaction_pb2.Vessel.PRESSURE_REACTOR,
    "pressure_reactor": reaction_pb2.Vessel.PRESSURE_REACTOR,
}


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------


def _as_list(val: object) -> list:
    """Wraps a dict in a list; returns lists unchanged; returns [] for None."""
    if val is None:
        return []
    if isinstance(val, dict):
        return [val]
    if not isinstance(val, list):
        return []
    return val


def _text(val: object) -> str:
    """Extracts text from an etree_to_dict value (plain string or dict with '#text')."""
    if isinstance(val, dict):
        d = cast("dict[str, object]", val)
        return str(d.get("#text") or "")
    return str(val) if val is not None else ""


def _safe_filename(name: str, max_len: int = 200) -> str:
    """Converts a UDM title to a safe filesystem filename stem."""
    # Replace slash first so "Pd/C" → "Pd_C", not "C" (pathlib strips pre-slash part).
    safe = name.replace("/", "_").replace("\\", "_")
    safe = pathlib.Path(safe).name
    safe = re.sub(r"[^\w\s\-.]", "_", safe)
    return safe[:max_len] or "ord_dataset"


def _attr_unit(node: dict, *keys: str) -> str:
    """Returns the first non-empty unit attribute from a UDM dict node."""
    for key in keys:
        raw = node.get(key)
        if raw is not None and str(raw).strip():
            return str(raw).strip().lower()
    return ""


# SURF / literature exports often write "<AMOUNT>0.3000 mmol</AMOUNT>".
_COMBINED_AMOUNT_RE = re.compile(
    r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*([A-Za-zµμ°]+)\s*$"
)


def _split_combined_amount(text: str) -> tuple[str, str] | None:
    """Splits a combined amount string into (value, unit), or None if not matched."""
    match = _COMBINED_AMOUNT_RE.match(text)
    if not match:
        return None
    return match.group(1), match.group(2)


def _parse_amount(value_str: object, unit_str: object) -> reaction_pb2.Amount | None:
    """Converts a (value, unit) pair from UDM into an ORD Amount.

    Returns None if the value cannot be parsed as a finite float.
    """
    # Failure 1: <AMOUNT unit="g">1.5</AMOUNT> → {'@unit': 'g', '#text': '1.5'}
    # Legacy: <AMOUNT units="g">1.5</AMOUNT> or separate AMOUNT_UNIT child.
    # SURF: <AMOUNT>0.3000 mmol</AMOUNT> (value and unit in one text node).
    if isinstance(value_str, dict):
        d = cast("dict[str, object]", value_str)
        unit_str = unit_str or d.get("@unit") or d.get("@units")
        value_str = d.get("#text")
        if value_str is None:
            logger.warning(
                "AMOUNT element has attributes but no text content; skipping."
            )
            return None
    if value_str is None:
        return None
    text = str(value_str)
    try:
        value = float(text)
    except (ValueError, TypeError):
        split = _split_combined_amount(text)
        if split is None:
            return None
        value_text, combined_unit = split
        unit_str = unit_str or combined_unit
        try:
            value = float(value_text)
        except (ValueError, TypeError):
            return None
    if not math.isfinite(value):
        logger.warning("Non-finite AMOUNT value %r; skipping.", value_str)
        return None

    unit_key = str(unit_str).strip().lower() if unit_str else ""

    amount = reaction_pb2.Amount()
    if unit_key in _MASS_UNITS:
        amount.mass.value = value
        amount.mass.units = _MASS_UNITS[unit_key]
    elif unit_key in _MOLES_UNITS:
        amount.moles.value = value
        amount.moles.units = _MOLES_UNITS[unit_key]
    elif unit_key in _VOLUME_UNITS:
        amount.volume.value = value
        amount.volume.units = _VOLUME_UNITS[unit_key]
    else:
        # Unknown or absent unit: store as mass with value only so data is not lost.
        amount.mass.value = value
    return amount


def _compound_amount(compound_entry: dict) -> reaction_pb2.Amount | None:
    """Parses AMOUNT, SAMPLE_MASS, or VOLUME from a UDM compound role block."""
    for key in ("AMOUNT", "SAMPLE_MASS", "VOLUME"):
        raw = compound_entry.get(key)
        if raw is None:
            continue
        unit_raw = compound_entry.get("AMOUNT_UNIT") or compound_entry.get("UNIT")
        parsed = _parse_amount(raw, unit_raw)
        if parsed is not None:
            return parsed
    return None


def _molblock_parses(molblock: str) -> bool:
    """Tests whether RDKit can read a MolBlock, using ORD's own parse check."""
    return message_helpers.identifier_parses_unsanitized(
        reaction_pb2.CompoundIdentifier.MOLBLOCK, molblock
    )


def _normalize_molblock(text: str) -> str:
    """Recovers a MolBlock from XML text content.

    MolBlock lines are column-sensitive and its first line is a title that is often
    blank, so a newline after the opening tag and the indentation an XML formatter
    adds are both indistinguishable from the MolBlock's own content. Each repair is
    therefore only applied when it is what makes RDKit read the block; the text is
    returned as recorded when none of them helps.
    """
    as_recorded = text.rstrip() + "\n"
    dedented = textwrap.dedent(text)
    candidates = dict.fromkeys(  # Deduplicated so unindented text is parsed once.
        (
            as_recorded,
            text.removeprefix("\n").rstrip() + "\n",
            dedented.rstrip() + "\n",
            dedented.removeprefix("\n").rstrip() + "\n",
        )
    )
    return next(
        (candidate for candidate in candidates if _molblock_parses(candidate)),
        as_recorded,
    )


def _add_compound_identifier(
    component: reaction_pb2.Compound | reaction_pb2.ProductCompound,
    molval: dict,
    *,
    mol_id: str = "",
    fallback_name: str = "",
) -> None:
    """Populates identifiers on a Compound from the molecule lookup dict.

    A MOLSTRUCTURE that RDKit cannot read is not recorded as a MOLBLOCK, since ORD
    validation rejects an unreadable one; the molecule name is recorded instead.
    When name and structure are both missing (common in Reaxys stubs), the mol ID
    string is used as a NAME so validation does not see an empty identifier.
    """
    molblock = molval.get("molblock", "")
    if molblock and _molblock_parses(molblock):
        component.identifiers.add(
            type="MOLBLOCK",
            details="MOLECULE -> MOLSTRUCTURE from UDM",
            value=molblock,
        )
        return
    if molblock:
        logger.warning(
            "MOLSTRUCTURE for %r is not a readable MolBlock; recording NAME instead.",
            molval.get("name", "") or mol_id,
        )
    # NAME rather than CUSTOM: UDM MOLECULE/NAME is a common name, and a CUSTOM
    # identifier additionally requires a details string.
    name = (
        (molval.get("name") or "").strip()
        or (fallback_name or "").strip()
        or (mol_id or "").strip()
    )
    component.identifiers.add(type="NAME", value=name)


# ---------------------------------------------------------------------------
# XML → dict helper (StackOverflow: https://stackoverflow.com/q/7684333)
# ---------------------------------------------------------------------------


def etree_to_dict(t: ET.Element) -> dict:
    """Recursively converts an ElementTree element to a nested dict."""
    d: dict = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd: defaultdict = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(("@" + k, v) for k, v in t.attrib.items())
    if t.text:
        # MolBlock text is column- and line-sensitive, so it is kept verbatim here
        # and normalised in _normalize_molblock() instead.
        text = t.text if t.tag in _RAW_TEXT_TAGS else t.text.strip()
        if children or t.attrib:
            if text.strip():
                d[t.tag]["#text"] = text
        else:
            d[t.tag] = text
    return d


# ---------------------------------------------------------------------------
# Per-variation conversion
# ---------------------------------------------------------------------------


def _build_molecule_lookup(udm: dict) -> dict[str, dict]:
    """Builds a mol_id → {name, molblock?} lookup from UDM MOLECULES."""
    lookup: dict[str, dict] = {}
    molecules_root = udm.get("MOLECULES") or {}
    for molecule in _as_list(molecules_root.get("MOLECULE")):
        mol_id = molecule.get("@ID")
        if not mol_id:
            continue
        # SURF often has CAS and no NAME; CAS is still a usable display label.
        entry: dict = {
            "name": _text(molecule.get("NAME")) or _text(molecule.get("CAS"))
        }  # Failure 11
        molblock_raw = molecule.get("MOLSTRUCTURE")  # Failure 8
        if isinstance(molblock_raw, dict):
            molblock_raw = molblock_raw.get("#text", "")
        if molblock_raw:
            entry["molblock"] = _normalize_molblock(str(molblock_raw))
        lookup[mol_id] = entry
    return lookup


def _variation_with_section(variation: dict) -> dict:
    """Promotes VARIATION/SECTION children to variation level (SURF dialect).

    Strict UDM places REACTANT/PRODUCT/CONDITIONS directly under VARIATION. SURF
    nests them under a single SECTION extension element. Existing top-level
    variation keys win when both are present.
    """
    section = variation.get("SECTION")
    if section is None:
        return variation
    if isinstance(section, list):
        if len(section) > 1:
            logger.warning(
                "Multiple SECTION elements found under VARIATION; using the first."
            )
        section = section[0] if section else {}
    if not isinstance(section, dict):
        return variation

    merged = dict(variation)
    for key, value in section.items():
        existing = merged.get(key)
        if existing in (None, "", {}, []):
            merged[key] = value
    return merged


def _map_rxn_identifiers(reaction: dict, pb2_reaction: reaction_pb2.Reaction) -> None:
    """Maps UDM RXNSTRUCTURE entries to ORD ReactionIdentifiers."""
    for structure in _as_list(reaction.get("RXNSTRUCTURE")):
        if isinstance(structure, dict):
            # Schema: format attr + text content. Legacy: @value attribute.
            udmformat = structure.get("@format", "")
            value = _text(structure) or str(structure.get("@value") or "")
        else:
            udmformat, value = "", str(structure or "")
        if udmformat == "cdxml":
            ordtype, orddetails = 1, "cdxml"
        elif udmformat == "rinchi":
            ordtype, orddetails = 5, ""
        elif udmformat == "rsmiles":
            ordtype, orddetails = 2, ""
        elif udmformat == "rxn":
            ordtype, orddetails = 1, "rxn"
        else:
            ordtype, orddetails = 0, ""
        pb2_reaction.identifiers.add(type=ordtype, details=orddetails, value=value)


def _map_inputs(
    variation: dict,
    all_molecules: dict[str, dict],
    pb2_reaction: reaction_pb2.Reaction,
    *,
    reaction: dict | None = None,
) -> None:
    """Maps UDM REACTANT / REAGENT / CATALYST / SOLVENT to ORD ReactionInputs.

    When the variation has no role compound blocks (common in Reaxys), falls back
    to VARIATION/REACTANT_ID then REACTION/REACTANT_ID resolved through the
    molecule lookup as REACTANT inputs.
    """
    role_map = {
        "REACTANT": reaction_pb2.ReactionRole.REACTANT,
        "REAGENT": reaction_pb2.ReactionRole.REAGENT,
        "CATALYST": reaction_pb2.ReactionRole.CATALYST,
        "SOLVENT": reaction_pb2.ReactionRole.SOLVENT,
    }
    mapped_any = False
    for udm_key, ord_role in role_map.items():
        for compound_entry in _as_list(variation.get(udm_key)):
            mol_ref = compound_entry.get("MOLECULE")
            if not isinstance(mol_ref, dict):
                continue
            mol_id = mol_ref.get("@MOL_ID", "")
            molval = all_molecules.get(mol_id)
            if molval is None:
                logger.warning(
                    "Molecule %r not found in MOLECULES lookup; skipping.", mol_id
                )
                continue

            # Failure 6: role-qualified key prevents cross-role slot contamination.
            input_key = f"{mol_id}_{udm_key}"
            molinput = pb2_reaction.inputs[input_key]
            molcomponent = molinput.components.add()
            _add_compound_identifier(
                molcomponent,
                molval,
                mol_id=str(mol_id),
                fallback_name=_text(compound_entry.get("NAME")),
            )
            molcomponent.reaction_role = ord_role

            parsed = _compound_amount(compound_entry)
            if parsed is not None:
                molcomponent.amount.CopyFrom(parsed)
            else:
                # ORD requires an amount on every input component; Reaxys often
                # omits quantities. Record an explicit unmeasured placeholder.
                molcomponent.amount.unmeasured.type = (
                    reaction_pb2.UnmeasuredAmount.CUSTOM
                )
                molcomponent.amount.unmeasured.details = "amount not reported in UDM"
            mapped_any = True

    if mapped_any:
        return

    reactant_ids = _mol_ids_from(variation, "REACTANT_ID")
    if not reactant_ids and reaction is not None:
        reactant_ids = _mol_ids_from(reaction, "REACTANT_ID")
    for mol_id in reactant_ids:
        molval = all_molecules.get(mol_id)
        input_key = f"{mol_id}_REACTANT"
        molinput = pb2_reaction.inputs[input_key]
        molcomponent = molinput.components.add()
        if molval is not None:
            _add_compound_identifier(molcomponent, molval, mol_id=mol_id)
        else:
            logger.warning(
                "REACTANT_ID %r not found in MOLECULES lookup; recording ID as NAME.",
                mol_id,
            )
            molcomponent.identifiers.add(type="NAME", value=mol_id)
        molcomponent.reaction_role = reaction_pb2.ReactionRole.REACTANT
        molcomponent.amount.unmeasured.type = reaction_pb2.UnmeasuredAmount.CUSTOM
        molcomponent.amount.unmeasured.details = "amount not reported in UDM"


_CONDITION_FIELDS = frozenset(
    {
        "TEMPERATURE",
        "PRESSURE",
        "TIME",
        "STIRRING",
        "PH",
        "REFLUX",
        "ATMOSPHERE",
        "PREPARATION",
        "VESSEL",
    }
)


def _condition_group(conditions: dict) -> dict:
    """Returns the condition-field dict for a UDM CONDITIONS element.

    Prefers CONDITION_GROUP when present. SURF places TEMPERATURE/TIME directly
    under CONDITIONS with no CONDITION_GROUP wrapper.
    """
    cg_raw = conditions.get("CONDITION_GROUP")
    if isinstance(cg_raw, list):
        if len(cg_raw) > 1:
            logger.warning("Multiple CONDITION_GROUPs found; using the first.")
        cg_raw = cg_raw[0] if cg_raw else {}
    if isinstance(cg_raw, dict) and cg_raw:
        return cg_raw
    if _CONDITION_FIELDS & conditions.keys():
        return conditions
    return {}


def _map_conditions(
    variation: dict,
    pb2_reaction: reaction_pb2.Reaction,
) -> None:
    """Maps UDM CONDITIONS to ORD ReactionConditions and ReactionSetup."""
    conditions = variation.get("CONDITIONS") or {}
    if not isinstance(conditions, dict):
        return
    # Failure 2: normalise multiple CONDITION_GROUPs to the first entry.
    cg = _condition_group(conditions)
    if not cg:
        return

    # Temperature — Failure 7: guard against non-finite values.
    temp = cg.get("TEMPERATURE")
    if isinstance(temp, dict):
        exact = temp.get("exact")
        if exact is not None:
            with contextlib.suppress(ValueError, TypeError):
                val = float(exact)
                if math.isfinite(val):
                    pb2_reaction.conditions.temperature.setpoint.value = val
        unit_key = _attr_unit(temp, "@unit", "@units")
        if unit_key in _TEMP_UNITS:
            pb2_reaction.conditions.temperature.setpoint.units = _TEMP_UNITS[unit_key]
        elif pb2_reaction.conditions.temperature.setpoint.HasField("value"):
            # SURF omits unit; Celsius is the lab default for these exports.
            pb2_reaction.conditions.temperature.setpoint.units = (
                reaction_pb2.Temperature.CELSIUS
            )

    # Pressure — Failure 7: guard against non-finite values.
    # Reaxys often omits units and mixes scales (torr vs Pa vs atm); guessing is
    # unreliable, so setpoint is only written when a recognised unit is present.
    # The raw number is preserved in conditions.details for domain review.
    pressure = cg.get("PRESSURE")
    if isinstance(pressure, dict):
        exact = pressure.get("exact")
        unit_key = _attr_unit(pressure, "@unit", "@units")
        if exact is not None:
            with contextlib.suppress(ValueError, TypeError):
                val = float(exact)
                if math.isfinite(val):
                    if unit_key in _PRESSURE_UNITS:
                        pb2_reaction.conditions.pressure.setpoint.value = val
                        pb2_reaction.conditions.pressure.setpoint.units = (
                            _PRESSURE_UNITS[unit_key]
                        )
                    else:
                        note = (
                            f"UDM PRESSURE exact={val} "
                            "(unit omitted; not mapped to setpoint)"
                        )
                        if pb2_reaction.conditions.details:
                            pb2_reaction.conditions.details += f"; {note}"
                        else:
                            pb2_reaction.conditions.details = note

    # ATMOSPHERE is a CONDITION_GROUP sibling in UDM v6 (legacy: under PRESSURE).
    atm_raw = (
        str(
            cg.get("ATMOSPHERE")
            or (pressure.get("ATMOSPHERE") if isinstance(pressure, dict) else "")
            or ""
        )
        .strip()
        .lower()
    )
    if atm_raw in _ATMOSPHERE_TYPES:
        pb2_reaction.conditions.pressure.atmosphere.type = _ATMOSPHERE_TYPES[atm_raw]

    # Stirring — schema uses stirringRange; older files may use free text.
    stirring = cg.get("STIRRING")
    if isinstance(stirring, dict):
        exact = stirring.get("exact")
        unit_key = _attr_unit(stirring, "@unit", "@units") or "rpm"
        if exact is not None:
            with contextlib.suppress(ValueError, TypeError):
                rpm_val = float(exact)
                if math.isfinite(rpm_val):
                    pb2_reaction.conditions.stirring.rate.rpm = int(rpm_val)
                    pb2_reaction.conditions.stirring.details = (
                        f"{int(rpm_val)} {unit_key}"
                    )
        pb2_reaction.conditions.stirring.type = (
            reaction_pb2.StirringConditions.CUSTOM
            if pb2_reaction.conditions.stirring.details
            else reaction_pb2.StirringConditions.UNSPECIFIED
        )
    elif stirring is not None:
        text = str(stirring)
        pb2_reaction.conditions.stirring.details = text
        lowered = text.lower()
        stirring_type = next(
            (value for keyword, value in _STIRRING_TYPE_KEYWORDS if keyword in lowered),
            reaction_pb2.StirringConditions.CUSTOM,
        )
        pb2_reaction.conditions.stirring.type = stirring_type
        rpm_match = re.search(r"(\d+)\s*rpm", lowered)
        if rpm_match:
            pb2_reaction.conditions.stirring.rate.rpm = int(rpm_match.group(1))

    # Reflux
    reflux_raw = cg.get("REFLUX")
    if reflux_raw is not None:
        pb2_reaction.conditions.reflux = str(reflux_raw).strip().lower() in (
            "true",
            "yes",
            "1",
        )

    # pH — Failure 12: accept plain string; Failure 7: guard non-finite.
    ph = cg.get("PH")
    if isinstance(ph, dict):
        exact = ph.get("exact")
    elif ph is not None:
        exact = ph  # plain string like "7.0"
    else:
        exact = None
    if exact is not None:
        with contextlib.suppress(ValueError, TypeError):
            val = float(exact)
            if math.isfinite(val):
                pb2_reaction.conditions.ph = val

    # Setup environment from CONDITION_GROUP PREPARATION when it matches a known env.
    prep_raw = cg.get("PREPARATION")
    preparations = [prep_raw] if isinstance(prep_raw, str) else _as_list(prep_raw)
    for preparation in preparations:
        env_key = str(preparation).strip().lower()
        env_type = _ENVIRONMENT_TYPES.get(env_key)
        if env_type is not None:
            pb2_reaction.setup.environment.type = env_type
            break
        if preparation and not pb2_reaction.setup.environment.details:
            # Free-text prep is not a known environment keyword; CUSTOM + details
            # satisfies ORD type/details validation.
            pb2_reaction.setup.environment.type = (
                reaction_pb2.ReactionSetup.ReactionEnvironment.CUSTOM
            )
            pb2_reaction.setup.environment.details = str(preparation)

    # Vessel
    vessel_raw = cg.get("VESSEL") or variation.get("VESSEL")
    if isinstance(vessel_raw, dict):
        vessel_type_raw = str(vessel_raw.get("VESSEL_TYPE", "")).strip().lower()
        vessel_type = _VESSEL_TYPES.get(vessel_type_raw)
        if vessel_type is not None:
            pb2_reaction.setup.vessel.type = vessel_type
        details = vessel_raw.get("DETAILS")
        if details:
            pb2_reaction.setup.vessel.details = str(details)
    elif vessel_raw:
        vessel_type = _VESSEL_TYPES.get(str(vessel_raw).strip().lower())
        if vessel_type is not None:
            pb2_reaction.setup.vessel.type = vessel_type


def _map_notes(variation: dict, pb2_reaction: reaction_pb2.Reaction) -> None:
    """Maps UDM procedure text to ORD notes.procedure_details.

    Prefers legacy VARIATION/PROCEDURE, else CONDITIONS/PREPARATION (UDM v6).
    """
    procedure = variation.get("PROCEDURE")
    if not procedure:
        conditions = variation.get("CONDITIONS") or {}
        if isinstance(conditions, dict):
            prep_raw = conditions.get("PREPARATION")
            preps = [prep_raw] if isinstance(prep_raw, str) else _as_list(prep_raw)
            # Prefer a PREPARATION that is not an environment keyword.
            procedure = next(
                (p for p in preps if str(p).strip().lower() not in _ENVIRONMENT_TYPES),
                preps[0] if preps else None,
            )
    if procedure:
        pb2_reaction.notes.procedure_details = _text(procedure) or str(procedure)


def _map_observations(variation: dict, pb2_reaction: reaction_pb2.Reaction) -> None:
    """Maps UDM COMMENT to ORD observations."""
    comment_raw = variation.get("COMMENT")
    comments = [comment_raw] if isinstance(comment_raw, str) else _as_list(comment_raw)
    for comment in comments:
        text = _text(comment) if isinstance(comment, dict) else str(comment)
        if text:
            obs = pb2_reaction.observations.add()
            obs.comment = text


def _condition_time(variation: dict) -> dict | None:
    """Returns a TIME/DURATION dict from VARIATION or CONDITIONS."""
    duration = variation.get("DURATION")
    if isinstance(duration, dict):
        return duration
    conditions = variation.get("CONDITIONS") or {}
    if not isinstance(conditions, dict):
        return None
    cg = _condition_group(conditions)
    time_node = cg.get("TIME")
    if isinstance(time_node, dict):
        return time_node
    return None


def _mol_ids_from(node: dict, key: str) -> list[str]:
    """Returns ID strings for a UDM ID field such as PRODUCT_ID or REACTANT_ID."""
    ids: list[str] = []
    raw = node.get(key)
    # Plain string is dropped by _as_list (same Failure 4 pattern as MODIFICATION_DATE).
    entries = [raw] if isinstance(raw, str) else _as_list(raw)
    for entry in entries:
        text = _text(entry) if isinstance(entry, dict) else str(entry or "")
        text = text.strip()
        if text:
            ids.append(text)
    return ids


def _map_outcomes(
    variation: dict,
    all_molecules: dict[str, dict],
    pb2_reaction: reaction_pb2.Reaction,
    *,
    reaction: dict | None = None,
) -> None:
    """Maps UDM PRODUCT entries and reaction time to ORD ReactionOutcomes.

    Prefers VARIATION/PRODUCT blocks. When none are present (common in Reaxys),
    falls back to VARIATION/PRODUCT_ID then REACTION/PRODUCT_ID resolved through
    the molecule lookup.
    """
    # Failure 5: don't create an empty outcome for input-only variations.
    duration = _condition_time(variation)
    product_entries = _as_list(variation.get("PRODUCT"))
    product_ids: list[str] = []
    if not product_entries:
        product_ids = _mol_ids_from(variation, "PRODUCT_ID")
        if not product_ids and reaction is not None:
            product_ids = _mol_ids_from(reaction, "PRODUCT_ID")
    if not (duration or product_entries or product_ids):
        return

    outcome = pb2_reaction.outcomes.add()

    # Reaction time — Failure 7: guard non-finite. UDM v6 uses CONDITIONS/TIME.
    if isinstance(duration, dict):
        dur_val = duration.get("exact") or duration.get("value")
        unit_key = _attr_unit(duration, "@unit", "@units")
        if dur_val is not None:
            with contextlib.suppress(ValueError, TypeError):
                val = float(dur_val)
                if math.isfinite(val):
                    outcome.reaction_time.value = val
            if unit_key in _TIME_UNITS:
                outcome.reaction_time.units = _TIME_UNITS[unit_key]
            elif outcome.reaction_time.HasField("value"):
                # SURF omits unit; hour is the usual duration scale in these exports.
                outcome.reaction_time.units = reaction_pb2.Time.HOUR

    # Products from VARIATION/PRODUCT blocks.
    for udm_product in product_entries:
        mol_ref = udm_product.get("MOLECULE")
        if not isinstance(mol_ref, dict):
            continue
        mol_id = mol_ref.get("@MOL_ID", "")
        molval = all_molecules.get(mol_id)
        product = outcome.products.add()
        if molval is not None:
            _add_compound_identifier(
                product,
                molval,
                mol_id=str(mol_id),
                fallback_name=_text(udm_product.get("NAME")),
            )
        else:
            logger.warning("Product molecule %r not found in MOLECULES lookup.", mol_id)
            # SURF PRODUCT blocks often carry a CAS-like NAME when MOLECULES lookup misses.
            product_name = _text(udm_product.get("NAME")) or str(mol_id)
            if product_name:
                product.identifiers.add(type="NAME", value=product_name)

        yield_data = udm_product.get("YIELD")
        if yield_data is not None:
            yield_val = (
                yield_data.get("exact") if isinstance(yield_data, dict) else yield_data
            )
            if yield_val is not None:
                # Failure 7: add() only after confirming a finite value.
                y: float | None = None
                with contextlib.suppress(ValueError, TypeError):
                    raw_y = float(yield_val)
                    if math.isfinite(raw_y):
                        y = raw_y
                if y is not None:
                    measurement = product.measurements.add()
                    measurement.type = (
                        reaction_pb2.ProductMeasurement.ProductMeasurementType.YIELD
                    )
                    measurement.percentage.value = y

    # Reaxys-style PRODUCT_ID fallback when no PRODUCT blocks were present.
    for mol_id in product_ids:
        product = outcome.products.add()
        molval = all_molecules.get(mol_id)
        if molval is not None:
            _add_compound_identifier(product, molval, mol_id=mol_id)
        else:
            logger.warning(
                "PRODUCT_ID %r not found in MOLECULES lookup; recording ID as NAME.",
                mol_id,
            )
            product.identifiers.add(type="NAME", value=mol_id)


def _scientist_fields(scientist: object) -> tuple[str, str]:
    """Returns (name, email) from a UDM SCIENTIST element.

    UDM v6 allows either a bare name string or an AUTHOR-shaped block with
    ``NAME`` / ``EMAIL`` children; both forms are accepted.
    """
    if isinstance(scientist, dict):
        d = cast("dict[str, object]", scientist)
        return _text(d.get("NAME")), _text(d.get("EMAIL"))
    if scientist is not None:
        return str(scientist), ""
    return "", ""


def _fill_person(
    person: reaction_pb2.Person,
    *,
    username: str = "",
    name: str = "",
    orcid: str = "",
    email: str = "",
) -> None:
    """Writes non-empty fields onto a Person, leaving already-set values alone."""
    if username and not person.username:
        person.username = username
    if name and not person.name:
        person.name = name
    if orcid and not person.orcid:
        person.orcid = orcid
    if email and not person.email:
        person.email = email


def _normalize_doi(raw: object) -> str:
    """Returns an ORD-valid DOI string, or empty if none can be parsed."""
    text = _text(raw) if isinstance(raw, dict) else str(raw or "").strip()
    if not text:
        return ""
    with contextlib.suppress(ValueError):
        return message_helpers.parse_doi(text)
    return text


def _map_provenance(
    reaction: dict,
    variation: dict,
    udm: dict,
    pb2_reaction: reaction_pb2.Reaction,
    *,
    username: str = "",
    person_name: str = "",
    orcid: str = "",
    email: str = "",
    created_date: str = "",
) -> None:
    """Maps UDM provenance fields to ORD ReactionProvenance.

    CLI depositor fields fill gaps when UDM has no SCIENTIST name/email or
    CREATION_DATE (common in SURF exports). UDM wins when both are present.
    Username and ORCID have no UDM equivalent in this converter. Dataset
    ``--name`` / ``--description`` are separate packaging overrides (CLI wins).
    """
    legal = udm.get("LEGAL") or {}

    # Failure 11: PRODUCER may carry XML attributes, producing a dict.
    producer = _text(legal.get("PRODUCER"))
    if producer:
        pb2_reaction.provenance.experimenter.organization = producer
        pb2_reaction.provenance.record_created.person.organization = producer

    scientist_name, scientist_email = _scientist_fields(variation.get("SCIENTIST"))
    # UDM first, then CLI fallbacks for missing person fields.
    resolved_name = scientist_name or person_name
    resolved_email = scientist_email or email
    if resolved_name:
        pb2_reaction.provenance.experimenter.name = resolved_name
        pb2_reaction.provenance.record_created.person.name = resolved_name
    if resolved_email:
        pb2_reaction.provenance.experimenter.email = resolved_email
        pb2_reaction.provenance.record_created.person.email = resolved_email
    _fill_person(
        pb2_reaction.provenance.experimenter,
        username=username,
        orcid=orcid,
    )
    _fill_person(
        pb2_reaction.provenance.record_created.person,
        username=username,
        orcid=orcid,
    )

    orgs = _as_list(reaction.get("ORGANISATIONS"))
    if orgs:
        address = (orgs[0].get("ORGANISATION") or {}).get("ADDRESS")
        if address:
            pb2_reaction.provenance.city = str(address)

    # DOI resolution: per-variation citation overrides global DOI.
    # Failure 11: DOI element may carry XML attributes.
    # SURF: VARIATION/@CIT_ID; strict UDM: VARIATION/CITATION/@CIT_ID.
    global_doi = _normalize_doi(legal.get("DOI"))
    variation_doi = ""
    cit_id = ""
    var_citation = variation.get("CITATION")
    if var_citation is not None:
        cit_id = (
            var_citation if isinstance(var_citation, dict) else var_citation[0]
        ).get("@CIT_ID", "")
    cit_id = cit_id or str(variation.get("@CIT_ID") or "")
    if cit_id:
        for citation in _as_list((udm.get("CITATIONS") or {}).get("CITATION")):
            if citation.get("@ID") == cit_id and "DOI" in citation:
                variation_doi = _normalize_doi(citation.get("DOI"))
                break

    pb2_reaction.provenance.doi = variation_doi or global_doi

    # Patent from reaction-level citations (legacy path)
    reaction_citations = _as_list(reaction.get("CITATIONS"))
    if reaction_citations:
        first_cit = reaction_citations[0].get("CITATION") or {}
        if "DOI" in first_cit:
            pb2_reaction.provenance.doi = _normalize_doi(first_cit["DOI"])
        if "PATENT_NUMBER" in first_cit:
            pb2_reaction.provenance.patent = first_cit["PATENT_NUMBER"]

    # UDM CREATION_DATE wins; --created-date fills when absent.
    creation_date = variation.get("CREATION_DATE") or created_date
    if creation_date:
        pb2_reaction.provenance.record_created.time.value = str(creation_date)

    # Failure 4: plain-string MODIFICATION_DATE is silently dropped by _as_list.
    mod_raw = variation.get("MODIFICATION_DATE")
    mod_dates = [mod_raw] if isinstance(mod_raw, str) else _as_list(mod_raw)
    for mod_date in mod_dates:
        event = pb2_reaction.provenance.record_modified.add()
        event.time.value = str(mod_date)
        _fill_person(
            event.person,
            username=username,
            name=resolved_name,
            orcid=orcid,
            email=resolved_email,
        )

    pb2_reaction.provenance.is_mined = False


def _validation_flag_hints(error_text: str) -> str:
    """Returns CLI flag hints for common ORD validation gaps UDM often omits."""
    text = error_text.lower()
    hints: list[str] = []
    if "email" in text and ("required" in text or "must have" in text):
        hints.append(
            "Pass --email when the UDM file has no SCIENTIST/EMAIL "
            "(ORD requires record_created.person.email)."
        )
    if "time" in text and ("recordevent" in text or "must have" in text):
        hints.append(
            "Pass --created-date when the UDM file has no CREATION_DATE "
            "(ORD requires record_created.time)."
        )
    if "username" in text or "orcid" in text:
        hints.append(
            "Pass --person-name, --username, or --orcid so record_created.person "
            "is identifiable (ORD requires at least one)."
        )
    if not hints:
        return ""
    return "Hint:\n" + "\n".join(f"  • {h}" for h in hints)


# ---------------------------------------------------------------------------
# Main conversion logic
# ---------------------------------------------------------------------------


def convert(
    input_path: pathlib.Path,
    *,
    name: str = "",
    description: str = "",
    username: str = "",
    person_name: str = "",
    orcid: str = "",
    email: str = "",
    created_date: str = "",
) -> dataset_pb2.Dataset:
    """Parses a UDM XML file and returns an ORD Dataset.

    Args:
        input_path: Path to the UDM v6.0.0 XML file.
        name: Dataset name. When set, overrides UDM LEGAL/TITLE (packaging
            override, unlike provenance gap-fill flags below).
        description: Dataset description. When set, overrides the DOI-derived
            default (packaging override).
        username: Depositor username for provenance (no UDM equivalent).
        person_name: Depositor display name; fills only when UDM SCIENTIST has
            no name (UDM wins if both present). Distinct from ``name``.
        orcid: Depositor ORCID iD for provenance (no UDM equivalent).
        email: Depositor email; fills only when UDM SCIENTIST has no email (UDM
            wins if both present). ORD validation requires email on
            record_created.
        created_date: Depositor record-created timestamp; fills only when UDM
            has no CREATION_DATE (UDM wins if both present). ORD validation
            requires record_created.time.

    Returns:
        A populated dataset_pb2.Dataset.

    Raises:
        SystemExit: On unrecoverable input errors (missing file, bad XML, wrong root).
    """
    logger.info(
        "Starting conversion from UDM v6.0.0 to ORD. "
        "*** The ORD repository uses the CC-BY-SA license. Do not push converted data "
        "to ord-data unless you hold the authority to relicense it. ***"
    )

    if not input_path.is_file():
        logger.error("Input file not found: %s", input_path)
        sys.exit(1)

    try:
        root = ET.parse(input_path).getroot()  # noqa: S314  (parses user-supplied UDM files; defusedxml not a dependency)
    except ET.ParseError:
        logger.exception("Failed to parse XML: %s", input_path)
        sys.exit(1)

    raw = etree_to_dict(root)
    if "UDM" not in raw:
        logger.error(
            "Input file is not UDM format — <UDM> must be the root element: %s",
            input_path,
        )
        sys.exit(1)

    udm = raw["UDM"]

    # Dataset-level metadata — Failure 11: TITLE/DOI may carry XML attributes.
    legal = udm.get("LEGAL") or {}
    dataset_name = name or _text(legal.get("TITLE"))
    global_doi = _text(legal.get("DOI"))
    dataset_description = description or (
        f"UDM dataset DOI: {global_doi}" if global_doi else ""
    )

    if "REACTIONS" not in udm:
        logger.error("<REACTIONS> element not found in %s", input_path)
        sys.exit(1)

    all_molecules = _build_molecule_lookup(udm)
    pb2_reactions: list[reaction_pb2.Reaction] = []

    # Failure 10: <REACTIONS/> (self-closing) produces None, not {}.
    for reaction in _as_list((udm["REACTIONS"] or {}).get("REACTION")):
        _map_rxn_identifiers(reaction, _scratch := reaction_pb2.Reaction())
        rxn_identifiers = _scratch.identifiers[:]  # carry forward to each variation

        for variation in _as_list(reaction.get("VARIATION")):
            variation = _variation_with_section(variation)
            pb2_reaction = reaction_pb2.Reaction()

            # Copy reaction-level identifiers into each variation's Reaction.
            for ident in rxn_identifiers:
                new_id = pb2_reaction.identifiers.add()
                new_id.CopyFrom(ident)

            _map_inputs(variation, all_molecules, pb2_reaction, reaction=reaction)
            _map_conditions(variation, pb2_reaction)
            _map_notes(variation, pb2_reaction)
            _map_observations(variation, pb2_reaction)
            _map_outcomes(
                variation, all_molecules, pb2_reaction, reaction=reaction
            )
            _map_provenance(
                reaction,
                variation,
                udm,
                pb2_reaction,
                username=username,
                person_name=person_name,
                orcid=orcid,
                email=email,
                created_date=created_date,
            )

            pb2_reactions.append(pb2_reaction)

        if not _as_list(reaction.get("VARIATION")):
            logger.warning("Reaction has no VARIATION elements; skipping.")

    dataset = dataset_pb2.Dataset(
        name=dataset_name,
        description=dataset_description,
        reactions=pb2_reactions,
    )
    updates.update_dataset(dataset)
    return dataset


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert a UDM v6.0.0 XML file to an ORD Dataset (.pbtxt or .pb)."
    )
    parser.add_argument("--input", required=True, help="Path to the UDM XML file.")
    parser.add_argument(
        "--output",
        default=None,
        help="Output path for the ORD Dataset (*.pbtxt or *.pb). "
        "Defaults to <UDM TITLE>.pbtxt or ord_dataset.pbtxt.",
    )
    parser.add_argument(
        "--name",
        default="",
        help="Dataset name. Overrides UDM LEGAL/TITLE when both are set "
        "(packaging override; not a person identity field).",
    )
    parser.add_argument(
        "--description",
        default="",
        help="Dataset description. Overrides the DOI-derived default when set "
        "(packaging override).",
    )
    parser.add_argument(
        "--username",
        default="",
        help="Depositor username written into provenance Person.username "
        "(no UDM equivalent).",
    )
    parser.add_argument(
        "--person-name",
        default="",
        help="Depositor display name for provenance. Fills only when UDM has no "
        "SCIENTIST/NAME (UDM wins if both set). Distinct from --name.",
    )
    parser.add_argument(
        "--orcid",
        default="",
        help="Depositor ORCID iD written into provenance Person.orcid "
        "(no UDM equivalent).",
    )
    parser.add_argument(
        "--email",
        default="",
        help="Depositor email for provenance. Fills only when UDM has no "
        "SCIENTIST/EMAIL (UDM wins if both set). ORD requires email on "
        "record_created.",
    )
    parser.add_argument(
        "--created-date",
        default="",
        help="record_created.time (e.g. 2024-01-15). Fills only when UDM has no "
        "CREATION_DATE (UDM wins if both set). ORD requires a time on "
        "record_created.",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip ORD schema validation of the converted dataset.",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Entry point for the UDM → ORD converter."""
    input_path = pathlib.Path(args.input)
    dataset = convert(
        input_path,
        name=args.name,
        description=args.description,
        username=args.username,
        person_name=args.person_name,
        orcid=args.orcid,
        email=args.email,
        created_date=args.created_date,
    )

    # Failure 9: catch ValidationError and exit cleanly instead of showing a traceback.
    if not args.no_validate:
        try:
            validations.validate_datasets({"_COMBINED": dataset})
        except validations.ValidationError as exc:
            hint = _validation_flag_hints(str(exc))
            message = f"Validation failed (use --no-validate to write anyway):\n{exc}"
            if hint:
                message = f"{message}\n{hint}"
            logger.error(  # noqa: TRY400 — traceback is noise for user-facing validation errors
                "%s", message
            )
            sys.exit(1)

    # Failure 3: sanitise dataset name before using as a filename.
    if args.output:
        output_path = pathlib.Path(args.output)
    elif dataset.name:
        output_path = pathlib.Path(_safe_filename(dataset.name) + ".pbtxt")
    else:
        output_path = pathlib.Path("ord_dataset.pbtxt")

    message_helpers.save_message(dataset, output_path)
    logger.info("Wrote %d reaction(s) to %s", len(dataset.reactions), output_path)


if __name__ == "__main__":
    main(parse_args())
