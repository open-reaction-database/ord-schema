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
"""

import argparse
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


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Convert a UDM dataset to ORD")
    parser.add_argument("--input", required=True, help="XML filename in UDM format")
    parser.add_argument("--output", help="Output Dataset filename (*.pbtxt)")
    parser.add_argument("--name", help="Name for this dataset")
    parser.add_argument("--description", help="Description for this dataset")
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
                pb2_reaction, reaction, variation, udm_reactions, all_molecules
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


def _as_list(value: Any) -> list[UdmDict]:
    """Normalizes an etree_to_dict value that may be absent, a dict, or a list."""
    if not value:
        return []
    if isinstance(value, dict):
        return [value]
    return value


def _build_base_reaction(reaction: UdmDict) -> reaction_pb2.Reaction:
    """Builds the reaction-level data shared by every variation of `reaction`.

    Covers the UDM fields that live directly on <REACTION> rather than on a
    specific <VARIATION>: REACTANT_ID/PRODUCT_ID placeholders and RXNSTRUCTURE
    identifiers.
    """
    base_reaction = reaction_pb2.Reaction()

    if "REACTANT_ID" in reaction:
        molinput = base_reaction.inputs["REACTANT_IDS"]
        for reactant_id in reaction["REACTANT_ID"]:
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
        for product_id in reaction["PRODUCT_ID"]:
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
) -> None:
    """Populates a Reaction already seeded with base data from one VARIATION."""
    for role_name, role in _ROLES.items():
        for entry in _as_list(variation.get(role_name)):
            _add_component(pb2_reaction, all_molecules, entry, role, role_name)

    condition_group: UdmDict = variation.get("CONDITIONS", {}).get(
        "CONDITION_GROUP", {}
    )
    if condition_group:
        _add_conditions(pb2_reaction, condition_group)

    comment = variation.get("COMMENT")
    if comment:
        pb2_reaction.observations.add().comment = comment

    for udm_product in _as_list(variation.get("PRODUCT")):
        _add_product(pb2_reaction, all_molecules, udm_product)

    _add_provenance(pb2_reaction, reaction, variation, udm_reactions)
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
        # UDM's parsed <AMOUNT> is a bare number with no unit attribute visible
        # in the source data available at the time of writing; grams are
        # assumed. Revisit once real UDM sample data confirms the unit.
        component.amount.mass.value = float(entry["AMOUNT"])
        component.amount.mass.units = reaction_pb2.Mass.GRAM


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
    """Populates ReactionConditions/ReactionSetup from a UDM CONDITION_GROUP."""
    temperature = condition_group.get("TEMPERATURE", {})
    if "exact" in temperature:
        # No unit is present in the source data available at the time of
        # writing, so Celsius is assumed.
        pb2_reaction.conditions.temperature.setpoint.value = float(temperature["exact"])
        pb2_reaction.conditions.temperature.setpoint.units = (
            reaction_pb2.Temperature.CELSIUS
        )

    pressure = condition_group.get("PRESSURE", {})
    if "exact" in pressure:
        # Same caveat as temperature: unit assumed, not present in the source data.
        pb2_reaction.conditions.pressure.setpoint.value = float(pressure["exact"])
        pb2_reaction.conditions.pressure.setpoint.units = (
            reaction_pb2.Pressure.ATMOSPHERE
        )

    if "STIRRING" in condition_group:
        pb2_reaction.conditions.stirring.type = reaction_pb2.StirringConditions.CUSTOM
        pb2_reaction.conditions.stirring.details = condition_group["STIRRING"]

    if "REFLUX" in condition_group:
        pb2_reaction.conditions.reflux = _parse_bool(condition_group["REFLUX"])

    ph = condition_group.get("PH", {})
    if "exact" in ph:
        pb2_reaction.conditions.ph = float(ph["exact"])

    if "PREPARATION" in condition_group:
        pb2_reaction.setup.environment.type = (
            reaction_pb2.ReactionSetup.ReactionEnvironment.CUSTOM
        )
        pb2_reaction.setup.environment.details = condition_group["PREPARATION"]


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
) -> None:
    """Populates ReactionProvenance from reaction- and variation-level UDM fields."""
    legal: UdmDict = udm_reactions.get("LEGAL", {})
    producer = legal.get("PRODUCER")
    if producer:
        pb2_reaction.provenance.experimenter.organization = producer
        pb2_reaction.provenance.record_created.person.organization = producer

    scientist = variation.get("SCIENTIST")
    if scientist:
        pb2_reaction.provenance.experimenter.name = scientist
        pb2_reaction.provenance.record_created.person.name = scientist

    if "ORGANISATIONS" in reaction:
        pb2_reaction.provenance.city = reaction["ORGANISATIONS"][0]["ORGANISATION"][
            "ADDRESS"
        ]

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
        if scientist:
            event.person.name = scientist
        if producer:
            event.person.organization = producer


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
