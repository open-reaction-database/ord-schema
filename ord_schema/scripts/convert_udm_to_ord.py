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

UDM (Unified Data Model) is a Pistoia Alliance format (MIT license); see
https://github.com/PistoiaAlliance/UDM for the schema, documentation, and
example data.

This converter's handling of the SURF dialect of UDM is based on
https://github.com/alexarnimueller/surf (MIT license). Please cite SURF as:
Nippa, Mueller, Atz, Konrad, Grether, Martin & Schneider (2023), "Simple
User-Friendly Reaction Format," ChemRxiv, https://doi.org/10.26434/chemrxiv-2023-nfq7h

Example usage::

    python convert_udm_to_ord.py \\
        --input my_dataset.xml \\
        --output my_dataset.pbtxt \\
        --email me@example.com --person-name "Ada Lovelace" \\
        --created-date 2024-01-15
"""

import argparse
import contextlib
import copy
import math
import pathlib
import re
import sys
import textwrap
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import cast

from ord_schema import message_helpers, units, validations
from ord_schema.logging import get_logger
from ord_schema.proto import dataset_pb2, reaction_pb2

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Unit resolution
# ---------------------------------------------------------------------------

# Spellings live in ord_schema.units so this converter cannot drift from them.
_UNIT_RESOLVER = units.UnitResolver()

# UDM unitMass defines "gr" as grain. units.py maps that spelling to gram.
_UDM_GRAIN_UNIT = "gr"

_AMOUNT_DEFAULT_UNITS = {
    "AMOUNT": "mol",  # molType
    "SAMPLE_MASS": "g",  # massType
    "VOLUME": "L",  # VOLUME / volumeRange
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

# Tags whose text must stay unstripped; see etree_to_dict(). MolBlock and RXN
# counts lines are column-sensitive.
_RAW_TEXT_TAGS = frozenset({"MOLSTRUCTURE", "RXNSTRUCTURE"})

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


def _finite_float(value: object) -> float | None:
    """Returns a finite float parsed from an XML value, or None."""
    with contextlib.suppress(ValueError, TypeError):
        parsed = float(_text(value))
        if math.isfinite(parsed):
            return parsed
    return None


def _parse_range(value: object) -> tuple[float, float | None] | None:
    """Returns (midpoint, precision) for a UDM exact/min/max range.

    Lone bounds cannot be represented structurally in ORD and return None.
    Inverted complete ranges are normalized before calculating precision.
    """
    if not isinstance(value, dict):
        parsed = _finite_float(value)
        return (parsed, None) if parsed is not None else None
    mapping = cast("dict[str, object]", value)
    if "exact" in mapping:
        parsed = _finite_float(mapping["exact"])
        return (parsed, None) if parsed is not None else None
    minimum = _finite_float(mapping.get("min"))
    maximum = _finite_float(mapping.get("max"))
    if minimum is None or maximum is None:
        return None
    low, high = sorted((minimum, maximum))
    precision = (high - low) / 2
    return (low + precision, precision or None)


_CONDITION_DEFAULT_UNITS = {
    "TEMPERATURE": "degC",
    "TIME": "hr",
    "STIRRING": "rpm",
    "REACTION_MOLARITY": "mol/L",
    "BUFFER_CONCENTRATION": "mol/L",
    "TOTAL_VOLUME": "L",
}


def _format_number(value: float) -> str:
    """Formats a numeric condition value compactly."""
    return f"{value:g}"


def _format_range_text(value: dict, default_unit: str = "") -> str | None:
    """Formats a UDM range, including lone bounds, for condition details."""
    unit = _attr_unit(value, "@unit", "@units") or default_unit
    parsed = _parse_range(value)
    if parsed is not None:
        midpoint, precision = parsed
        text = _format_number(midpoint)
        if precision is not None:
            text += f"±{_format_number(precision)}"
    else:
        minimum = _finite_float(value.get("min"))
        maximum = _finite_float(value.get("max"))
        if minimum is not None:
            text = f">={_format_number(minimum)}"
        elif maximum is not None:
            text = f"<={_format_number(maximum)}"
        else:
            return None
    return f"{text} {unit}".strip()


def _format_detail_value(value: object, *, default_unit: str = "") -> str:
    """Formats an XML-derived scalar, list, range, or nested dictionary."""
    if isinstance(value, list):
        return ", ".join(
            text
            for item in value
            if (text := _format_detail_value(item, default_unit=default_unit))
        )
    if isinstance(value, dict):
        mapping = cast("dict[str, object]", value)
        range_text = _format_range_text(mapping, default_unit)
        if range_text is not None:
            increment = _format_detail_value(mapping.get("incr"))
            return f"{range_text}; ramp={increment}" if increment else range_text
        text = _text(mapping)
        if text:
            unit = _attr_unit(mapping, "@unit", "@units")
            return f"{text} {unit}".strip()
        return ", ".join(
            f"{key.lstrip('@')}={formatted}"
            for key, item in mapping.items()
            if key != "SECTION" and (formatted := _format_detail_value(item))
        )
    return _text(value).strip()


def _condition_group_details(
    group: dict,
    *,
    exclude: frozenset[str] = frozenset(),
) -> str:
    """Formats CONDITION_GROUP fields not represented structurally."""
    parts = []
    for key, value in group.items():
        if key in exclude or key == "SECTION" or key.startswith("@"):
            continue
        formatted = _format_detail_value(
            value,
            default_unit=_CONDITION_DEFAULT_UNITS.get(key, ""),
        )
        if formatted:
            parts.append(f"{key.lower().replace('_', ' ')} {formatted}")
    return "; ".join(parts)


def _append_condition_details(
    pb2_reaction: reaction_pb2.Reaction,
    details: str,
) -> None:
    """Appends non-empty text to ReactionConditions.details."""
    if not details:
        return
    if pb2_reaction.conditions.details:
        pb2_reaction.conditions.details += f"; {details}"
    else:
        pb2_reaction.conditions.details = details


# SURF / literature exports often write "<AMOUNT>0.3000 mmol</AMOUNT>".
_COMBINED_AMOUNT_RE = re.compile(
    r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*([A-Za-zµμ°]+)\s*$"
)


def _resolve_unit(unit_key: str) -> tuple[type, int] | None:
    """Resolves a UDM unit spelling to an ORD message class and enum value.

    Returns None when the spelling is unknown, forbidden, or ``gr``. UDM
    ``unitMass`` defines ``gr`` as grain; ``units.py`` uses that spelling for gram.
    """
    if not unit_key or unit_key == _UDM_GRAIN_UNIT:
        return None
    try:
        message_cls, enum_value = _UNIT_RESOLVER.resolve_unit(unit_key)
    except KeyError:
        return None
    return message_cls, enum_value


def _unit_enum(unit_key: str, message_cls: type) -> int | None:
    """Returns the enum value when ``unit_key`` resolves to ``message_cls``."""
    resolved = _resolve_unit(unit_key)
    if resolved is not None and resolved[0] is message_cls:
        return resolved[1]
    return None


def _split_combined_amount(text: str) -> tuple[str, str] | None:
    """Splits a combined amount string into (value, unit), or None if not matched."""
    match = _COMBINED_AMOUNT_RE.match(text)
    if not match:
        return None
    return match.group(1), match.group(2)


def _parse_amount(
    value_str: object,
    unit_str: object,
    *,
    default_unit: str = "",
) -> reaction_pb2.Amount | None:
    """Converts a (value, unit) pair from UDM into an ORD Amount.

    Unknown units are preserved in a CUSTOM unmeasured amount instead of being
    assigned an incorrect physical quantity type.

    Returns None if the value cannot be parsed as a finite float.
    """
    # <AMOUNT unit="g">1.5</AMOUNT> → {'@unit': 'g', '#text': '1.5'} after XML parse.
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

    unit_text = _text(unit_str).strip() if unit_str else ""
    if not unit_text:
        unit_text = default_unit
    unit_key = unit_text.lower()

    amount = reaction_pb2.Amount()
    resolved = _resolve_unit(unit_key)
    if resolved is not None and resolved[0] is reaction_pb2.Mass:
        amount.mass.value = value
        amount.mass.units = resolved[1]
    elif resolved is not None and resolved[0] is reaction_pb2.Moles:
        amount.moles.value = value
        amount.moles.units = resolved[1]
    elif resolved is not None and resolved[0] is reaction_pb2.Volume:
        amount.volume.value = value
        amount.volume.units = resolved[1]
    else:
        amount.unmeasured.type = reaction_pb2.UnmeasuredAmount.CUSTOM
        if unit_text:
            amount.unmeasured.details = (
                f"UDM amount {value:g} has unsupported unit {unit_text!r}"
            )
        else:
            amount.unmeasured.details = f"UDM amount {value:g} has no unit"
    return amount


def _compound_amount(compound_entry: dict) -> reaction_pb2.Amount | None:
    """Parses AMOUNT, SAMPLE_MASS, or VOLUME from a UDM compound role block."""
    for key in ("AMOUNT", "SAMPLE_MASS", "VOLUME"):
        raw = compound_entry.get(key)
        if raw is None:
            continue
        unit_raw = compound_entry.get("AMOUNT_UNIT") or compound_entry.get("UNIT")
        # XSD defaults: molType → mol, massType → g, VOLUME → L.
        parsed = _parse_amount(
            raw,
            unit_raw,
            default_unit=_AMOUNT_DEFAULT_UNITS[key],
        )
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
        }
        molblock_raw = molecule.get("MOLSTRUCTURE")
        if isinstance(molblock_raw, dict):
            # Attributes such as format="mol" leave the MolBlock in #text.
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
    """Maps UDM RXNSTRUCTURE entries to ORD ReactionIdentifiers.

    ``format`` is optional and defaults to ``rxn``. With no attributes,
    ``etree_to_dict`` returns a plain string, which ``_as_list`` would drop.
    """
    raw = reaction.get("RXNSTRUCTURE")
    structures = [raw] if isinstance(raw, str) else _as_list(raw)
    identifier = reaction_pb2.ReactionIdentifier
    for structure in structures:
        if isinstance(structure, dict):
            # Schema: format attr + text content. Legacy: @value attribute.
            # Missing format is the XSD default, rxn, not an unknown type.
            udmformat = str(structure.get("@format") or "rxn")
            value = _text(structure) or str(structure.get("@value") or "")
        else:
            udmformat, value = "rxn", str(structure or "")
        if udmformat == "cdxml":
            ordtype, orddetails = identifier.CUSTOM, "cdxml"
        elif udmformat == "rinchi":
            ordtype, orddetails = identifier.RINCHI, ""
        elif udmformat == "rsmiles":
            ordtype, orddetails = identifier.REACTION_SMILES, ""
        elif udmformat == "rxn":
            ordtype, orddetails = identifier.CUSTOM, "rxn"
        else:
            ordtype, orddetails = identifier.UNSPECIFIED, ""
        pb2_reaction.identifiers.add(type=ordtype, details=orddetails, value=value)


def _map_inputs(
    variation: dict,
    all_molecules: dict[str, dict],
    pb2_reaction: reaction_pb2.Reaction,
    *,
    reaction: dict | None = None,
) -> None:
    """Maps UDM REACTANT / REAGENT / CATALYST / SOLVENT to ORD ReactionInputs.

    UDM role blocks do not record addition order or grouping, so every compound
    in a variation shares one ReactionInput (key ``combined``). Each block stays
    its own component, with its own reaction_role and amount. A separate
    ReactionInput would claim a separate addition event.

    When the variation has no role compound blocks (common in Reaxys), falls back
    to VARIATION/REACTANT_ID then REACTION/REACTANT_ID resolved through the
    molecule lookup as components of one ``REACTANT_IDS`` input.
    """
    role_map = {
        "REACTANT": reaction_pb2.ReactionRole.REACTANT,
        "REAGENT": reaction_pb2.ReactionRole.REAGENT,
        "CATALYST": reaction_pb2.ReactionRole.CATALYST,
        "SOLVENT": reaction_pb2.ReactionRole.SOLVENT,
    }
    molinput = None
    for udm_key, ord_role in role_map.items():
        for compound_entry in _as_list(variation.get(udm_key)):
            mol_ref = compound_entry.get("MOLECULE")
            if not isinstance(mol_ref, dict):
                continue
            mol_id = str(mol_ref.get("@MOL_ID", ""))
            molval = all_molecules.get(mol_id)
            local_name = _text(compound_entry.get("NAME")) or mol_id
            if molval is None and not local_name:
                logger.warning(
                    "Molecule %r not found in MOLECULES lookup; skipping.", mol_id
                )
                continue

            if molinput is None:
                molinput = pb2_reaction.inputs["combined"]
            molcomponent = molinput.components.add()
            if molval is not None:
                _add_compound_identifier(
                    molcomponent,
                    molval,
                    mol_id=mol_id,
                    fallback_name=local_name,
                )
            else:
                logger.warning(
                    "Molecule %r not found in MOLECULES lookup; recording ID as NAME.",
                    mol_id,
                )
                molcomponent.identifiers.add(type="NAME", value=local_name)
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

    if molinput is not None:
        return

    reactant_ids = _mol_ids_from(variation, "REACTANT_ID")
    if not reactant_ids and reaction is not None:
        reactant_ids = _mol_ids_from(reaction, "REACTANT_ID")
    if not reactant_ids:
        return
    molinput = pb2_reaction.inputs["REACTANT_IDS"]
    for mol_id in reactant_ids:
        molval = all_molecules.get(mol_id)
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
        "PROCESS",
        "REACTION_MOLARITY",
        "BUFFER_TYPE",
        "BUFFER_CONCENTRATION",
        "TOTAL_VOLUME",
        "REACTANT_ID",
        "REAGENT_ID",
        "CATALYST_ID",
        "SOLVENT_ID",
    }
)


def _condition_groups(conditions: dict) -> list[dict]:
    """Returns condition groups, including SURF's unwrapped representation.

    Strict UDM uses one or more CONDITION_GROUP children. SURF places condition
    fields directly under CONDITIONS without that wrapper.
    """
    groups = [
        group
        for group in _as_list(conditions.get("CONDITION_GROUP"))
        if isinstance(group, dict) and group
    ]
    if groups:
        return groups
    if _CONDITION_FIELDS & conditions.keys():
        return [conditions]
    return []


def _map_conditions(
    variation: dict,
    pb2_reaction: reaction_pb2.Reaction,
) -> None:
    """Maps UDM CONDITIONS to ORD ReactionConditions and ReactionSetup."""
    conditions = variation.get("CONDITIONS") or {}
    if not isinstance(conditions, dict):
        return
    groups = _condition_groups(conditions)
    if not groups:
        return
    if len(groups) > 1:
        pb2_reaction.conditions.conditions_are_dynamic = True
        for index, group in enumerate(groups, start=1):
            details = _condition_group_details(group)
            _append_condition_details(
                pb2_reaction,
                f"Stage {index}: {details or '(no data)'}",
            )
        return
    cg = groups[0]
    captured: set[str] = set()

    # Temperature: value and units are assigned atomically.
    temp = cg.get("TEMPERATURE")
    if isinstance(temp, dict):
        parsed = _parse_range(temp)
        unit_key = _attr_unit(temp, "@unit", "@units")
        # Unitless TEMPERATURE defaults to Celsius (UDM XSD temperatureRange → degC).
        units = (
            _unit_enum(unit_key, reaction_pb2.Temperature)
            if unit_key
            else reaction_pb2.Temperature.CELSIUS
        )
        if parsed is not None and units is not None:
            value, precision = parsed
            setpoint = pb2_reaction.conditions.temperature.setpoint
            setpoint.value = value
            setpoint.units = units
            if precision is not None:
                setpoint.precision = precision
            if "incr" not in temp:
                captured.add("TEMPERATURE")
        elif parsed is not None:
            logger.warning(
                "Unsupported TEMPERATURE unit %r; skipping setpoint.",
                unit_key,
            )

    # The XSD default for pressureRange is torr, but Reaxys omits the unit on
    # values that are not one scale (about 760, and also ~2 and ~4.5e6). Applying
    # torr would mislabel that file, so a setpoint is written only when a
    # recognised unit is present. The raw number stays in conditions.details.
    pressure = cg.get("PRESSURE")
    if isinstance(pressure, dict):
        parsed = _parse_range(pressure)
        unit_key = _attr_unit(pressure, "@unit", "@units")
        units = _unit_enum(unit_key, reaction_pb2.Pressure) if unit_key else None
        if parsed is not None and units is not None:
            value, precision = parsed
            setpoint = pb2_reaction.conditions.pressure.setpoint
            setpoint.value = value
            setpoint.units = units
            if precision is not None:
                setpoint.precision = precision
            captured.add("PRESSURE")
        elif parsed is not None and not unit_key:
            value, _ = parsed
            _append_condition_details(
                pb2_reaction,
                f"UDM PRESSURE value={value:g} (unit omitted; not mapped to setpoint)",
            )
            captured.add("PRESSURE")

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
        captured.add("ATMOSPHERE")

    # Stirring — schema uses stirringRange; older files may use free text.
    stirring = cg.get("STIRRING")
    if isinstance(stirring, dict):
        parsed = _parse_range(stirring)
        unit_key = _attr_unit(stirring, "@unit", "@units") or "rpm"
        if parsed is not None and unit_key == "rpm":
            rpm_value, precision = parsed
            rpm = round(rpm_value)
            if 0 <= rpm <= 2**31 - 1:
                pb2_reaction.conditions.stirring.type = (
                    reaction_pb2.StirringConditions.CUSTOM
                )
                pb2_reaction.conditions.stirring.rate.rpm = rpm
                pb2_reaction.conditions.stirring.details = f"{rpm} rpm"
                if precision is None:
                    captured.add("STIRRING")
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
        captured.add("STIRRING")

    # Reflux
    reflux_raw = cg.get("REFLUX")
    if reflux_raw is not None:
        pb2_reaction.conditions.reflux = str(reflux_raw).strip().lower() in (
            "true",
            "yes",
            "1",
        )
        captured.add("REFLUX")

    # pH accepts either a range dictionary or a plain string.
    ph = cg.get("PH")
    parsed_ph = _parse_range(ph) if ph is not None else None
    if parsed_ph is not None:
        value, precision = parsed_ph
        pb2_reaction.conditions.ph = value
        if precision is None:
            captured.add("PH")

    time_value = cg.get("TIME")
    if isinstance(time_value, dict):
        parsed_time = _parse_range(time_value)
        unit_key = _attr_unit(time_value, "@unit", "@units")
        if parsed_time is not None and (
            not unit_key or _unit_enum(unit_key, reaction_pb2.Time) is not None
        ):
            captured.add("TIME")

    # Environment only when PREPARATION is a known keyword. Procedure text goes
    # to notes.procedure_details (_map_notes), including this group's text.
    prep_raw = cg.get("PREPARATION")
    preparations = [prep_raw] if isinstance(prep_raw, str) else _as_list(prep_raw)
    for preparation in preparations:
        env_type = _ENVIRONMENT_TYPES.get(str(preparation).strip().lower())
        if env_type is not None:
            pb2_reaction.setup.environment.type = env_type
            break
    if preparations:
        captured.add("PREPARATION")

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
        if vessel_type is not None or details:
            captured.add("VESSEL")
    elif vessel_raw:
        vessel_type = _VESSEL_TYPES.get(str(vessel_raw).strip().lower())
        if vessel_type is not None:
            pb2_reaction.setup.vessel.type = vessel_type
            captured.add("VESSEL")

    _append_condition_details(
        pb2_reaction,
        _condition_group_details(cg, exclude=frozenset(captured)),
    )


def _preparation_texts(node: dict) -> list[str]:
    """Returns PREPARATION strings from a CONDITIONS or CONDITION_GROUP node."""
    raw = node.get("PREPARATION")
    preps = [raw] if isinstance(raw, str) else _as_list(raw)
    texts = []
    for prep in preps:
        text = _text(prep).strip()
        if text:
            texts.append(text)
    return texts


def _map_notes(variation: dict, pb2_reaction: reaction_pb2.Reaction) -> None:
    """Maps UDM procedure text to ORD notes.procedure_details.

    Prefers legacy VARIATION/PROCEDURE. Otherwise uses PREPARATION text from
    CONDITIONS or CONDITION_GROUP when it is not an environment keyword.
    """
    procedure = variation.get("PROCEDURE")
    if not procedure:
        conditions = variation.get("CONDITIONS") or {}
        preps: list[str] = []
        if isinstance(conditions, dict):
            preps.extend(_preparation_texts(conditions))
            for group in _as_list(conditions.get("CONDITION_GROUP")):
                if isinstance(group, dict):
                    preps.extend(_preparation_texts(group))
        procedure = next(
            (text for text in preps if text.lower() not in _ENVIRONMENT_TYPES),
            None,
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
    groups = _condition_groups(conditions)
    if len(groups) != 1:
        return None
    time_node = groups[0].get("TIME")
    if isinstance(time_node, dict):
        return time_node
    return None


def _mol_ids_from(node: dict, key: str) -> list[str]:
    """Returns ID strings for a UDM ID field such as PRODUCT_ID or REACTANT_ID."""
    ids: list[str] = []
    raw = node.get(key)
    # Plain-string ID fields are handled before _as_list (which would drop them).
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
    # Do not create an empty outcome for input-only variations.
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

    # Reaction time from CONDITIONS/TIME; non-finite values are dropped in _parse_range.
    if isinstance(duration, dict):
        parsed = _parse_range(duration)
        if parsed is None and "value" in duration:
            parsed = _parse_range(duration["value"])
        if parsed is not None:
            unit_key = _attr_unit(duration, "@unit", "@units")
            # Unitless TIME defaults to hour (UDM XSD timeRange → hr).
            units = (
                _unit_enum(unit_key, reaction_pb2.Time)
                if unit_key
                else reaction_pb2.Time.HOUR
            )
            if units is not None:
                value, precision = parsed
                outcome.reaction_time.value = value
                outcome.reaction_time.units = units
                if precision is not None:
                    outcome.reaction_time.precision = precision
            else:
                logger.warning(
                    "Unsupported reaction TIME unit %r; skipping reaction time.",
                    unit_key,
                )

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
            # SURF PRODUCT blocks often carry a CAS-like NAME
            # when the MOLECULES lookup misses.
            product_name = _text(udm_product.get("NAME")) or str(mol_id)
            if product_name:
                product.identifiers.add(type="NAME", value=product_name)

        yield_data = udm_product.get("YIELD")
        if yield_data is not None:
            parsed = _parse_range(yield_data)
            if parsed is not None:
                value, precision = parsed
                measurement = product.measurements.add()
                measurement.type = (
                    reaction_pb2.ProductMeasurement.ProductMeasurementType.YIELD
                )
                measurement.percentage.value = value
                if precision is not None:
                    measurement.percentage.precision = precision
            elif isinstance(yield_data, dict):
                details = _format_range_text(yield_data, "percent")
                if details:
                    measurement = product.measurements.add()
                    measurement.type = (
                        reaction_pb2.ProductMeasurement.ProductMeasurementType.YIELD
                    )
                    measurement.details = f"yield {details}"

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

    ``experimenter`` is the UDM SCIENTIST only, and stays unset when the export
    has none. ``record_created.person`` is that scientist, or the CLI depositor
    when there is no scientist — never a mix of the two identities. ``--email``
    still fills ``record_created.person.email`` when the scientist has no email,
    because ORD validation requires it. Dataset ``--name`` / ``--description``
    are separate packaging overrides (CLI wins).
    """
    legal = udm.get("LEGAL") or {}

    # PRODUCER may carry XML attributes; _text reads the #text node.
    producer = _text(legal.get("PRODUCER"))
    if producer:
        pb2_reaction.provenance.experimenter.organization = producer
        pb2_reaction.provenance.record_created.person.organization = producer

    scientist_name, scientist_email = _scientist_fields(variation.get("SCIENTIST"))
    if scientist_name or scientist_email:
        if scientist_name:
            pb2_reaction.provenance.experimenter.name = scientist_name
            pb2_reaction.provenance.record_created.person.name = scientist_name
        if scientist_email:
            pb2_reaction.provenance.experimenter.email = scientist_email
            pb2_reaction.provenance.record_created.person.email = scientist_email
        elif email:
            pb2_reaction.provenance.record_created.person.email = email
        record_username = ""
        record_name = scientist_name
        record_orcid = ""
        record_email = scientist_email or email
    else:
        _fill_person(
            pb2_reaction.provenance.record_created.person,
            username=username,
            name=person_name,
            orcid=orcid,
            email=email,
        )
        record_username = username
        record_name = person_name
        record_orcid = orcid
        record_email = email

    orgs = _as_list(reaction.get("ORGANISATIONS"))
    if orgs:
        address = (orgs[0].get("ORGANISATION") or {}).get("ADDRESS")
        if address:
            pb2_reaction.provenance.city = str(address)

    # DOI resolution: per-variation citation overrides global DOI.
    # DOI may carry XML attributes; _normalize_doi reads via _text.
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
        # Variation citation wins. Reaction-level DOI fills only when that
        # lookup did not resolve one; it still overrides LEGAL/DOI.
        if "DOI" in first_cit and not variation_doi:
            pb2_reaction.provenance.doi = _normalize_doi(first_cit["DOI"])
        if "PATENT_NUMBER" in first_cit:
            pb2_reaction.provenance.patent = first_cit["PATENT_NUMBER"]

    # UDM CREATION_DATE wins; --created-date fills when absent.
    creation_date = variation.get("CREATION_DATE") or created_date
    if creation_date:
        pb2_reaction.provenance.record_created.time.value = str(creation_date)

    # Plain-string MODIFICATION_DATE is handled without _as_list.
    mod_raw = variation.get("MODIFICATION_DATE")
    mod_dates = [mod_raw] if isinstance(mod_raw, str) else _as_list(mod_raw)
    for mod_date in mod_dates:
        event = pb2_reaction.provenance.record_modified.add()
        event.time.value = str(mod_date)
        _fill_person(
            event.person,
            username=record_username,
            name=record_name,
            orcid=record_orcid,
            email=record_email,
        )

    pb2_reaction.provenance.is_mined = False


def _validation_flag_hints(error_text: str) -> str:
    """Returns CLI flag hints for common ORD validation gaps UDM often omits."""
    text = error_text.lower()
    hints: list[str] = []
    missing_provenance_email = (
        "user email is required for record_created" in text
        or "user email is required for record_modified" in text
    )
    if missing_provenance_email:
        hints.append(
            "Pass --email when the UDM file has no SCIENTIST/EMAIL "
            "(ORD requires record_created.person.email)."
        )
    if "recordevent" in text and "time" in text and "must have" in text:
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


def _document_context_xml(root: ET.Element) -> str:
    """Serializes shared UDM context without mutating the parsed source tree."""
    context = ET.Element(root.tag, root.attrib)
    for child in root:
        if child.tag not in ("REACTIONS", "MOLECULES"):
            context.append(copy.deepcopy(child))
    return ET.tostring(context, encoding="unicode")


def _set_xml_metadata(
    pb2_reaction: reaction_pb2.Reaction,
    *,
    reaction_xml: str,
    parent_xml: str,
) -> None:
    """Stores opt-in UDM source XML on reaction provenance."""
    reaction_data = pb2_reaction.provenance.reaction_metadata["udm_reaction_xml"]
    reaction_data.string_value = reaction_xml
    reaction_data.format = "xml"
    reaction_data.description = (
        "Raw UDM REACTION element from which this ORD Reaction was converted."
    )
    parent_data = pb2_reaction.provenance.reaction_metadata["udm_parent_xml"]
    parent_data.string_value = parent_xml
    parent_data.format = "xml"
    parent_data.description = (
        "Shared UDM document context excluding REACTIONS and MOLECULES."
    )


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
    include_udm_xml: bool = False,
) -> dataset_pb2.Dataset:
    """Parses a UDM XML file and returns an ORD Dataset.

    Args:
        input_path: Path to the UDM v6.0.0 XML file.
        name: Dataset name. When set, overrides UDM LEGAL/TITLE (packaging
            override, unlike provenance gap-fill flags below).
        description: Dataset description. When set, overrides the DOI-derived
            default (packaging override).
        username: Depositor username. Written to record_created.person only when
            UDM has no SCIENTIST. Distinct from the experimenter.
        person_name: Depositor display name. Written to record_created.person
            only when UDM has no SCIENTIST. Distinct from ``name``.
        orcid: Depositor ORCID iD. Written to record_created.person only when
            UDM has no SCIENTIST.
        email: Depositor email. Fills record_created.person.email when UDM has
            no SCIENTIST/EMAIL. Not copied onto experimenter. ORD validation
            requires email on record_created.
        created_date: Depositor record-created timestamp; fills only when UDM
            has no CREATION_DATE (UDM wins if both present). ORD validation
            requires record_created.time.
        include_udm_xml: Whether to embed reaction and shared document XML in
            provenance.reaction_metadata.

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

    # Dataset-level metadata; TITLE/DOI may carry XML attributes (_text reads #text).
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
    reactions_element = root.find("REACTIONS")
    if reactions_element is None:
        logger.error("<REACTIONS> element not found in %s", input_path)
        sys.exit(1)
    parent_xml = _document_context_xml(root) if include_udm_xml else ""

    # Empty <REACTIONS/> is valid: findall returns no REACTION children.
    for reaction_element in reactions_element.findall("REACTION"):
        reaction = etree_to_dict(reaction_element)["REACTION"]
        reaction_xml = (
            ET.tostring(reaction_element, encoding="unicode") if include_udm_xml else ""
        )
        _map_rxn_identifiers(reaction, _scratch := reaction_pb2.Reaction())
        rxn_identifiers = _scratch.identifiers[:]  # carry forward to each variation

        variations = _as_list(reaction.get("VARIATION")) or [{}]
        for raw_variation in variations:
            variation = _variation_with_section(raw_variation)
            pb2_reaction = reaction_pb2.Reaction()

            # Copy reaction-level identifiers into each variation's Reaction.
            for ident in rxn_identifiers:
                new_id = pb2_reaction.identifiers.add()
                new_id.CopyFrom(ident)

            _map_inputs(variation, all_molecules, pb2_reaction, reaction=reaction)
            _map_conditions(variation, pb2_reaction)
            _map_notes(variation, pb2_reaction)
            _map_observations(variation, pb2_reaction)
            _map_outcomes(variation, all_molecules, pb2_reaction, reaction=reaction)
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
            if include_udm_xml:
                _set_xml_metadata(
                    pb2_reaction,
                    reaction_xml=reaction_xml,
                    parent_xml=parent_xml,
                )

            pb2_reactions.append(pb2_reaction)

    return dataset_pb2.Dataset(
        name=dataset_name,
        description=dataset_description,
        reactions=pb2_reactions,
    )


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
        help="Depositor username for record_created.person. Applied only when "
        "UDM has no SCIENTIST (not copied onto experimenter).",
    )
    parser.add_argument(
        "--person-name",
        default="",
        help="Depositor display name for record_created.person. Applied only "
        "when UDM has no SCIENTIST. Distinct from --name.",
    )
    parser.add_argument(
        "--orcid",
        default="",
        help="Depositor ORCID iD for record_created.person. Applied only when "
        "UDM has no SCIENTIST (not copied onto experimenter).",
    )
    parser.add_argument(
        "--email",
        default="",
        help="Depositor email for record_created.person. Fills when UDM has no "
        "SCIENTIST/EMAIL. Not copied onto experimenter. ORD requires email on "
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
        "--include-udm-xml",
        action="store_true",
        help=(
            "Embed source REACTION XML and shared document context in "
            "provenance.reaction_metadata. Off by default because source XML "
            "may have a different license."
        ),
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
        include_udm_xml=args.include_udm_xml,
    )

    # Catch ValidationError and exit cleanly with optional CLI hints.
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

    # Sanitise dataset name before using as a default filename.
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
