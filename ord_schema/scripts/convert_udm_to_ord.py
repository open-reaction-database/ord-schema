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
"""Builds a Dataset from a set of (UDM) Reactions.

Given a UDM file, converts it to a Dataset .pbtxt for ORD. A UDM <REACTION>
may contain multiple <VARIATION> elements; each variation is emitted as its
own ORD Reaction, since UDM variations represent distinct experiments run
under the same reaction identifiers.

UDM (Unified Data Model) is a Pistoia Alliance format; see
https://github.com/PistoiaAlliance/UDM for the schema and example data.
"""

import argparse
import getpass
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Any

from ord_schema import message_helpers, updates, validations
from ord_schema.logging import get_logger
from ord_schema.proto import dataset_pb2, reaction_pb2

logger = get_logger(__name__)

UdmDict = dict[str, Any]

_MOLBLOCK_DETAILS = "MOLECULE -> MOLSTRUCTURE from UDM"

# UDM roles that map directly onto an ORD ReactionRole.
_ROLES = {
    "REACTANT": reaction_pb2.ReactionRole.REACTANT,
    "REAGENT": reaction_pb2.ReactionRole.REAGENT,
    "CATALYST": reaction_pb2.ReactionRole.CATALYST,
    "SOLVENT": reaction_pb2.ReactionRole.SOLVENT,
}

# UDM RXNSTRUCTURE @format values that map onto an ORD ReactionIdentifierType.
# Anything not listed here (including UDM's default, unlabeled format) falls
# back to UNSPECIFIED; that reaction identifier is then dropped rather than
# emitted with an invalid type (see _build_base_reaction).
_RXNSTRUCTURE_FORMATS: dict[str, tuple[int, str]] = {
    "cdxml": (reaction_pb2.ReactionIdentifier.CUSTOM, "cdxml"),
    "rinchi": (reaction_pb2.ReactionIdentifier.RINCHI, ""),
    "rsmiles": (reaction_pb2.ReactionIdentifier.REACTION_SMILES, ""),
    "rxn": (reaction_pb2.ReactionIdentifier.CUSTOM, "rxn"),
}

# Unit mappings below are taken from the UDM v6.0.0 XSD schema
# (udm_6_0_0_units.xsd, udm_6_0_0.xsd), not guessed:
# https://github.com/PistoiaAlliance/UDM/blob/master/udm_6_0_0_units.xsd

# <AMOUNT> is UDM's molType: a decimal with an optional @unit attribute
# defaulting to "mol". Units not listed here fall back to UNSPECIFIED.
_MOLES_UNITS: dict[str, int] = {
    "mol": reaction_pb2.Moles.MOLE,
    "mmol": reaction_pb2.Moles.MILLIMOLE,
    "umol": reaction_pb2.Moles.MICROMOLE,
    "μmol": reaction_pb2.Moles.MICROMOLE,  # UDM allows the literal mu character.
    "nmol": reaction_pb2.Moles.NANOMOLE,
}

# temperatureRange's @unit attribute defaults to "degC".
_TEMPERATURE_UNITS: dict[str, int] = {
    "degC": reaction_pb2.Temperature.CELSIUS,
    "degF": reaction_pb2.Temperature.FAHRENHEIT,
    "K": reaction_pb2.Temperature.KELVIN,
}

# pressureRange's @unit attribute defaults to "torr", not "atm".
_PRESSURE_UNITS: dict[str, int] = {
    "atm": reaction_pb2.Pressure.ATMOSPHERE,
    "bar": reaction_pb2.Pressure.BAR,
    "psi": reaction_pb2.Pressure.PSI,
    "KPa": reaction_pb2.Pressure.KILOPASCAL,
    "Pa": reaction_pb2.Pressure.PASCAL,
    "torr": reaction_pb2.Pressure.TORR,
    "mmHg": reaction_pb2.Pressure.MM_HG,
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Convert a UDM dataset to ORD")
    parser.add_argument("--input", required=True, help="XML filename in UDM format")
    parser.add_argument("--output", help="Output Dataset filename (*.pbtxt)")
    parser.add_argument("--name", help="Name for this dataset")
    parser.add_argument("--description", help="Description for this dataset")
    parser.add_argument(
        "--email",
        help=(
            "Email address of the person running this conversion. Recorded as "
            "the record_created/record_modified provenance, since that person "
            "-- not the original UDM scientist, who is recorded separately as "
            "the experimenter -- is who is creating this ORD record."
        ),
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="If set, do not run validations on reactions",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> None:
    """Converts a UDM XML dataset to an ORD Dataset and writes it to disk."""
    logger.info("Starting conversion from UDM v6.0.0 to ORD.")
    logger.info(
        "*** Important message - The Open Reaction Database repository, ord-data, "
        "uses the CC-BY-SA license for all data. Please do not push your converted "
        "dataset to ord-data unless you have the authority to change the license "
        "on the data to CC-BY-SA. ***"
    )

    if not Path(args.input).is_file():
        logger.error("Conversion failed: --input file %s does not exist.", args.input)
        raise SystemExit(1)

    udm_tree = ET.parse(args.input)  # noqa: S314 -- UDM inputs are trusted local files, not untrusted network input.
    udm_root = etree_to_dict(udm_tree.getroot())
    if "UDM" not in udm_root:
        logger.error("Input file is not UDM format: expected a <UDM> root element.")
        raise SystemExit(1)
    udm_reactions: UdmDict = udm_root["UDM"]

    if "REACTIONS" not in udm_reactions:
        logger.error("No <REACTIONS> element found in input .xml file.")
        raise SystemExit(1)
    if "MOLECULES" not in udm_reactions:
        logger.error("No <MOLECULES> element found in input .xml file.")
        raise SystemExit(1)

    legal: UdmDict = udm_reactions.get("LEGAL", {})
    dataset_name = args.name or legal.get("TITLE", "")
    dataset_description = args.description or (
        f"UDM dataset DOI: {legal['DOI']}" if "DOI" in legal else ""
    )
    outputfile = args.output or (
        f"{legal['TITLE']}.pbtxt" if "TITLE" in legal else "ord_dataset.pbtxt"
    )

    all_molecules: dict[str, UdmDict] = {}
    for molecule in _as_list(udm_reactions["MOLECULES"]["MOLECULE"]):
        mol_data: UdmDict = {"name": molecule.get("NAME")}
        if "MOLSTRUCTURE" in molecule:
            mol_data["molblock"] = molecule["MOLSTRUCTURE"]
        all_molecules[molecule["@ID"]] = mol_data

    # The person running this conversion is the one creating this ORD record,
    # distinct from the original UDM scientist who performed the reaction
    # (recorded separately as the experimenter).
    submitter_username = getpass.getuser()

    pb2_reactions = []
    for reaction in _as_list(udm_reactions["REACTIONS"]["REACTION"]):
        base_reaction = _build_base_reaction(reaction)
        # Each variation becomes its own Reaction, seeded with the
        # reaction-level identifiers built above. A reaction with no
        # <VARIATION> elements still gets a single Reaction, since it has no
        # variation-specific data to add.
        for variation in _as_list(reaction.get("VARIATION")) or [{}]:
            pb2_reaction = reaction_pb2.Reaction()
            pb2_reaction.CopyFrom(base_reaction)
            _add_variation(
                pb2_reaction,
                reaction,
                variation,
                udm_reactions,
                all_molecules,
                submitter_username,
                args.email,
            )
            pb2_reactions.append(pb2_reaction)

    dataset = dataset_pb2.Dataset(
        name=dataset_name, description=dataset_description, reactions=pb2_reactions
    )
    # Assigns canonical dataset_id/reaction_ids and provenance updates.
    updates.update_dataset(dataset)

    if not args.no_validate:
        validations.validate_datasets({"_COMBINED": dataset})

    message_helpers.save_message(dataset, outputfile)
    logger.info("Conversion completed successfully: wrote %s.", outputfile)


def _as_list(value: Any) -> list[Any]:
    """Normalizes an etree_to_dict value that may be absent, a single item, or a list.

    etree_to_dict collapses a single occurrence of a repeatable element to a
    bare dict or string rather than a one-element list, so a naive `for x in
    value` silently iterates a string's characters instead of treating it as
    one item. Only an already-list value is left as-is.
    """
    if not value:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _text_and_unit(value: Any) -> tuple[str | None, str | None]:
    """Splits an etree_to_dict value into (text, unit).

    Several UDM elements (AMOUNT, PREPARATION) carry an optional @unit or
    @format attribute alongside their text. etree_to_dict represents such an
    element as plain text when no attribute is present, or as a dict with
    "#text" plus "@..." keys when one is. This normalizes both shapes.
    """
    if isinstance(value, dict):
        return value.get("#text"), value.get("@unit") or value.get("@format")
    return value, None


def _build_base_reaction(reaction: UdmDict) -> reaction_pb2.Reaction:
    """Builds the reaction-level data shared by every variation of `reaction`.

    Covers the UDM fields that live directly on <REACTION> rather than on a
    specific <VARIATION>: REACTANT_ID/PRODUCT_ID placeholders and RXNSTRUCTURE
    identifiers.
    """
    base_reaction = reaction_pb2.Reaction()

    if "REACTANT_ID" in reaction:
        molinput = base_reaction.inputs["REACTANT_IDS"]
        for reactant_id in _as_list(reaction["REACTANT_ID"]):
            identifier = molinput.components.add().identifiers.add(
                type="CUSTOM", details="REACTANT_ID from UDM"
            )
            identifier.value = reactant_id

    for structure in _as_list(reaction.get("RXNSTRUCTURE")):
        udm_format = structure.get("format", 1)
        ordtype, orddetails = _RXNSTRUCTURE_FORMATS.get(
            udm_format, (reaction_pb2.ReactionIdentifier.UNSPECIFIED, "")
        )
        if (
            "value" in structure
            and ordtype != reaction_pb2.ReactionIdentifier.UNSPECIFIED
        ):
            base_reaction.identifiers.add(
                type=ordtype, details=orddetails, value=structure["value"]
            )

    if "PRODUCT_ID" in reaction:
        outcome = base_reaction.outcomes.add()
        for product_id in _as_list(reaction["PRODUCT_ID"]):
            identifier = outcome.products.add().identifiers.add(
                type="CUSTOM", details="PRODUCT_ID from UDM"
            )
            identifier.value = product_id

    return base_reaction


def _add_variation(
    pb2_reaction: reaction_pb2.Reaction,
    reaction: UdmDict,
    variation: UdmDict,
    udm_reactions: UdmDict,
    all_molecules: dict[str, UdmDict],
    submitter_username: str,
    submitter_email: str | None,
) -> None:
    """Populates a Reaction already seeded with base data from one VARIATION."""
    for role_name, role in _ROLES.items():
        for entry in _as_list(variation.get(role_name)):
            _add_component(pb2_reaction, all_molecules, entry, role, role_name)

    conditions: UdmDict = variation.get("CONDITIONS", {})
    # <PREPARATION> is a direct child of <CONDITIONS>, not <CONDITION_GROUP>.
    preparations = []
    for prep in _as_list(conditions.get("PREPARATION")):
        text, _ = _text_and_unit(prep)
        if text:
            preparations.append(text)
    if preparations:
        pb2_reaction.setup.environment.type = (
            reaction_pb2.ReactionSetup.ReactionEnvironment.CUSTOM
        )
        pb2_reaction.setup.environment.details = "; ".join(preparations)
    # <CONDITIONS> may repeat <CONDITION_GROUP>; ORD has only one
    # ReactionConditions per Reaction, so later groups overwrite fields set
    # by earlier ones rather than being merged or dropped.
    for condition_group in _as_list(conditions.get("CONDITION_GROUP")):
        _add_conditions(pb2_reaction, condition_group)

    comment = variation.get("COMMENT")
    if comment:
        pb2_reaction.observations.add().comment = comment

    for udm_product in _as_list(variation.get("PRODUCT")):
        _add_product(pb2_reaction, all_molecules, udm_product)

    _add_provenance(
        pb2_reaction,
        reaction,
        variation,
        udm_reactions,
        submitter_username,
        submitter_email,
    )
    pb2_reaction.provenance.is_mined = False


def _add_component(
    pb2_reaction: reaction_pb2.Reaction,
    all_molecules: dict[str, UdmDict],
    entry: UdmDict,
    role: int,
    role_label: str,
) -> None:
    """Adds one UDM REACTANT/REAGENT/CATALYST/SOLVENT entry as a ReactionInput."""
    if "MOLECULE" not in entry:
        return
    mol_id = entry["MOLECULE"]["@MOL_ID"]
    component = pb2_reaction.inputs[mol_id].components.add()
    _set_molecule_identifier(component, all_molecules, mol_id, role_label)
    component.reaction_role = role
    if "AMOUNT" in entry:
        # UDM's <AMOUNT> (molType) is molar, defaulting to "mol"; it is not a
        # mass. See _MOLES_UNITS.
        text, unit = _text_and_unit(entry["AMOUNT"])
        if text:
            component.amount.moles.value = float(text)
            component.amount.moles.units = _MOLES_UNITS.get(
                unit or "mol", reaction_pb2.Moles.UNSPECIFIED
            )


def _set_molecule_identifier(
    component: reaction_pb2.Compound | reaction_pb2.ProductCompound,
    all_molecules: dict[str, UdmDict],
    mol_id: str,
    role_label: str,
) -> None:
    """Sets a CompoundIdentifier on `component`, falling back if `mol_id` is unknown."""
    molval = all_molecules.get(mol_id)
    if molval and molval.get("molblock"):
        identifier = component.identifiers.add(
            type="MOLBLOCK", details=_MOLBLOCK_DETAILS
        )
        identifier.value = molval["molblock"]
    elif molval and molval.get("name"):
        identifier = component.identifiers.add(
            type="CUSTOM", details=f"{role_label} MOLECULE NAME from UDM"
        )
        identifier.value = molval["name"]
    else:
        identifier = component.identifiers.add(
            type="CUSTOM",
            details=f"{role_label} MOL_ID from UDM (not found in <MOLECULES>)",
        )
        identifier.value = mol_id


def _add_conditions(
    pb2_reaction: reaction_pb2.Reaction, condition_group: UdmDict
) -> None:
    """Populates ReactionConditions from a UDM CONDITION_GROUP."""
    temperature = condition_group.get("TEMPERATURE", {})
    if "exact" in temperature:
        pb2_reaction.conditions.temperature.setpoint.value = float(temperature["exact"])
        pb2_reaction.conditions.temperature.setpoint.units = _TEMPERATURE_UNITS.get(
            temperature.get("@unit", "degC"), reaction_pb2.Temperature.UNSPECIFIED
        )

    pressure = condition_group.get("PRESSURE", {})
    if "exact" in pressure:
        pb2_reaction.conditions.pressure.setpoint.value = float(pressure["exact"])
        pb2_reaction.conditions.pressure.setpoint.units = _PRESSURE_UNITS.get(
            pressure.get("@unit", "torr"), reaction_pb2.Pressure.UNSPECIFIED
        )

    # UDM's STIRRING is a numeric rate (stirringRange, @unit defaulting to
    # "rpm"), not free text; ORD's StirringConditions.rate.rpm is the closest
    # match. Units other than rpm can't be represented and are dropped.
    stirring = condition_group.get("STIRRING", {})
    if "exact" in stirring and stirring.get("@unit", "rpm") == "rpm":
        rpm = round(float(stirring["exact"]))
        pb2_reaction.conditions.stirring.type = reaction_pb2.StirringConditions.CUSTOM
        pb2_reaction.conditions.stirring.details = f"{rpm} rpm"
        pb2_reaction.conditions.stirring.rate.rpm = rpm

    if "REFLUX" in condition_group:
        pb2_reaction.conditions.reflux = _parse_bool(condition_group["REFLUX"])

    ph = condition_group.get("PH", {})
    if "exact" in ph:
        pb2_reaction.conditions.ph = float(ph["exact"])


def _parse_bool(value: str) -> bool:
    """Parses a UDM boolean-ish string (e.g. "true"/"false"/"1"/"0")."""
    return str(value).strip().lower() in ("true", "1", "yes")


def _add_product(
    pb2_reaction: reaction_pb2.Reaction,
    all_molecules: dict[str, UdmDict],
    udm_product: UdmDict,
) -> None:
    """Adds one UDM VARIATION/PRODUCT entry as a ReactionOutcome product."""
    molecule = udm_product.get("MOLECULE")
    if not molecule:
        return
    component = pb2_reaction.outcomes.add().products.add()
    _set_molecule_identifier(
        component, all_molecules, molecule.get("@MOL_ID"), "PRODUCT"
    )
    yield_ = udm_product.get("YIELD", {})
    if "exact" in yield_:
        measurement = component.measurements.add()
        measurement.type = reaction_pb2.ProductMeasurement.ProductMeasurementType.YIELD
        measurement.float_value.value = float(yield_["exact"])


def _add_provenance(
    pb2_reaction: reaction_pb2.Reaction,
    reaction: UdmDict,
    variation: UdmDict,
    udm_reactions: UdmDict,
    submitter_username: str,
    submitter_email: str | None,
) -> None:
    """Populates ReactionProvenance from reaction- and variation-level UDM fields."""
    # record_created/record_modified describe who is creating this ORD
    # record -- i.e. whoever is running this conversion -- not the original
    # UDM scientist, who is recorded separately below as the experimenter.
    pb2_reaction.provenance.record_created.person.username = submitter_username
    if submitter_email:
        pb2_reaction.provenance.record_created.person.email = submitter_email

    legal: UdmDict = udm_reactions.get("LEGAL", {})
    producer = legal.get("PRODUCER")
    if producer:
        # Fallback organization; overridden below if the scientist has their
        # own affiliation, since PRODUCER is the data vendor (e.g. "REAXYS"),
        # not necessarily the scientist's research institution.
        pb2_reaction.provenance.experimenter.organization = producer

    # SCIENTIST is UDM's authorDetails type (NAME/EMAIL/PHONE/ORGANISATION),
    # not a plain string, and may repeat; the first is taken as the
    # experimenter of record.
    scientists = _as_list(variation.get("SCIENTIST"))
    scientist = scientists[0] if scientists else None
    if isinstance(scientist, dict):
        if scientist.get("NAME"):
            pb2_reaction.provenance.experimenter.name = scientist["NAME"]
        if scientist.get("EMAIL"):
            pb2_reaction.provenance.experimenter.email = scientist["EMAIL"]
        organisation = scientist.get("ORGANISATION")
        if isinstance(organisation, dict):
            if organisation.get("NAME"):
                pb2_reaction.provenance.experimenter.organization = organisation["NAME"]
            if organisation.get("ADDRESS"):
                pb2_reaction.provenance.city = organisation["ADDRESS"]
    elif scientist:
        # Defensive fallback in case a non-conformant UDM file has a bare
        # string SCIENTIST rather than the schema's authorDetails structure.
        pb2_reaction.provenance.experimenter.name = scientist

    doi = legal.get("DOI")
    if doi and "CITATIONS" not in reaction:
        pb2_reaction.provenance.doi = doi
    elif "CITATIONS" in udm_reactions and "CITATION" in variation:
        variation_citation = variation["CITATION"]
        if isinstance(variation_citation, list):
            variation_citation = variation_citation[0]
        for citation in _as_list(udm_reactions["CITATIONS"]["CITATION"]):
            if (
                citation.get("@ID") == variation_citation.get("@CIT_ID")
                and "DOI" in citation
            ):
                pb2_reaction.provenance.doi = citation["DOI"]
                break

    if "CITATIONS" in reaction:
        reaction_citation = _as_list(reaction["CITATIONS"])[0]["CITATION"]
        if "DOI" in reaction_citation:
            pb2_reaction.provenance.doi = reaction_citation["DOI"]
        if "PATENT_NUMBER" in reaction_citation:
            pb2_reaction.provenance.patent = reaction_citation["PATENT_NUMBER"]

    creation_date = variation.get("CREATION_DATE")
    if creation_date:
        pb2_reaction.provenance.record_created.time.value = creation_date

    modification_date = variation.get("MODIFICATION_DATE")
    if modification_date:
        event = pb2_reaction.provenance.record_modified.add()
        event.time.value = modification_date
        event.person.username = submitter_username
        if submitter_email:
            event.person.email = submitter_email


# Converts XML tree into Python dictionary.
# https://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree
def etree_to_dict(t: ET.Element) -> UdmDict:
    """Recursively converts an ElementTree Element into a nested dict."""
    d: UdmDict = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(("@" + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]["#text"] = text
        else:
            d[t.tag] = text
    return d


if __name__ == "__main__":
    main(parse_args())
