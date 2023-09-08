# Copyright 2020 Open Reaction Database Project Authors
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
"""Parses CML from the NRD.

Data is at:
https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873

See also:
https://depth-first.com/articles/2019/01/28/the-nextmove-patent-reaction-dataset/

Usage:
    parse_uspto.py --input_pattern=<str> --name=<str> --output=<str> [--n_jobs=<int>]

Options:
    --input_pattern=<str>       Input pattern for CML files
    --name=<str>                Dataset name
    --output=<str>              Output Dataset filename
    --n_jobs=<str>              Number of parallel workers [default: 1]
"""
# pylint: disable=import-error
# pylint: disable=too-many-branches

import datetime
import glob
import os
import re
from typing import Union
from xml.etree import ElementTree

import docopt
import joblib
from rdkit import RDLogger

import ord_schema
from ord_schema.logging import get_logger
from ord_schema import message_helpers
from ord_schema import units
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

logger = get_logger(__name__)
RDLogger.DisableLog("rdApp.*")  # Disable RDKit logging.

# XML namespaces.
NAMESPACES = {
    "cml": "http://www.xml-cml.org/schema",
    "cmlDict": "http://www.xml-cml.org/dictionary/cml/",
    "dl": "http://bitbucket.org/dan2097",
    "nameDict": "http://www.xml-cml.org/dictionary/cml/name/",
    "unit": "http://www.xml-cml.org/unit/",
}

UNIT_RESOLVER = units.UnitResolver()

REACTION_ROLES = {
    "catalyst": reaction_pb2.ReactionRole.CATALYST,
    "product": reaction_pb2.ReactionRole.PRODUCT,
    "reactant": reaction_pb2.ReactionRole.REACTANT,
    "solvent": reaction_pb2.ReactionRole.SOLVENT,
}

PRODUCT_STATES = {
    "crystal": reaction_pb2.Texture.CRYSTAL,
    "crystalline-solid": reaction_pb2.Texture.CRYSTAL,
    "crystals": reaction_pb2.Texture.CRYSTAL,
    "foam": reaction_pb2.Texture.FOAM,
    "foam-like": reaction_pb2.Texture.FOAM,
    "foam/very": reaction_pb2.Texture.FOAM,
    "needle": reaction_pb2.Texture.CRYSTAL,
    "needle-like": reaction_pb2.Texture.CRYSTAL,
    "needles": reaction_pb2.Texture.CRYSTAL,
    "oil": reaction_pb2.Texture.OIL,
    "oil-like": reaction_pb2.Texture.OIL,
    "oils": reaction_pb2.Texture.OIL,
    "powder": reaction_pb2.Texture.POWDER,
    "powders": reaction_pb2.Texture.POWDER,
    "prisms": reaction_pb2.Texture.CRYSTAL,
    "semisolid": reaction_pb2.Texture.SEMI_SOLID,
    "semi-solid": reaction_pb2.Texture.SEMI_SOLID,
}

WORKUP_TYPES = {
    "Add": reaction_pb2.ReactionWorkup.ADDITION,
    "Cool": reaction_pb2.ReactionWorkup.TEMPERATURE,
    "Concentrate": reaction_pb2.ReactionWorkup.CONCENTRATION,
    "Dissolve": reaction_pb2.ReactionWorkup.DISSOLUTION,
    "Distill": reaction_pb2.ReactionWorkup.DISTILLATION,
    "Dry with material": reaction_pb2.ReactionWorkup.DRY_WITH_MATERIAL,
    "Extract": reaction_pb2.ReactionWorkup.EXTRACTION,
    "Filter": reaction_pb2.ReactionWorkup.FILTRATION,
    "Heat": reaction_pb2.ReactionWorkup.TEMPERATURE,
    "Sample": reaction_pb2.ReactionWorkup.ALIQUOT,
    "Stir": reaction_pb2.ReactionWorkup.STIRRING,
    "Wait": reaction_pb2.ReactionWorkup.WAIT,
    "Wash": reaction_pb2.ReactionWorkup.WASH,
}

UNIT_REPLACEMENTS = {
    "˜": " ",
    "≈": " ",
    "about ": " ",
    "approx ": " ",
    "approx. ": " ",
    "approximately ": " ",
    "around ": " ",
    "ca. ": " ",
    "yield ": " ",
    "yield: ": " ",
    "Yield ": " ",
    "Yield: ": " ",
    "1/2 ": " 0.5 ",
    "½ ": " 0.5 ",
    "half ": " 0.5 ",
    " to ": "-",
    "×10": "e",
    "−": "-",
    # Digits.
    "zero ": "0 ",
    "one ": "1 ",
    "two ": "2 ",
    "three ": "3 ",
    "four ": "4 ",
    "five ": "5 ",
    "six ": "6 ",
    "seven ": "7 ",
    "eight ": "8 ",
    "nine ": "9 ",
}


def get_tag(element: ElementTree.Element) -> str:
    """Fetches the namespace-translated element tag."""
    for key, value in NAMESPACES.items():
        if value in element.tag:
            return element.tag.replace(f"{{{value}}}", f"{key}:")
    return element.tag


def get_molecule_id(root: ElementTree.Element) -> str:
    """Fetches the ID of the molecule under this element."""
    molecules = root.findall("cml:molecule", namespaces=NAMESPACES)
    if len(molecules) != 1:
        raise NotImplementedError(len(molecules))
    return molecules[0].attrib["id"]


def resolve_units(value: str) -> ord_schema.UnitMessage:
    """Resolves a value/unit string."""
    for key, replacement in UNIT_REPLACEMENTS.items():
        value = value.replace(key, replacement)
    return UNIT_RESOLVER.resolve(value, allow_range=True)


def get_component_map(root: ElementTree.Element) -> dict[str, str]:
    """Builds a mapping of components to inputs."""
    reaction_inputs = {}
    # Start with a separate input for each component.
    for component in root.find("cml:reactantList", namespaces=NAMESPACES):
        molecule_id = get_molecule_id(component)
        reaction_inputs[molecule_id] = molecule_id
    for component in root.find("cml:spectatorList", namespaces=NAMESPACES):
        molecule_id = get_molecule_id(component)
        reaction_inputs[molecule_id] = molecule_id
    for action in root.find("dl:reactionActionList", namespaces=NAMESPACES):
        components = []
        for component in action.findall("dl:chemical", namespaces=NAMESPACES):
            if "ref" in component.attrib:
                components.append(component.attrib["ref"])
        if components:
            key = "_".join(components)
            for component in components:
                reaction_inputs[component] = key
    return reaction_inputs


def parse_reaction(root: ElementTree.Element) -> reaction_pb2.Reaction:
    """Parses reaction CML into a Reaction message."""
    reaction = reaction_pb2.Reaction()
    component_map = get_component_map(root)
    for child in root:
        tag = get_tag(child)
        if tag == "dl:source":
            parse_source(child, reaction)
        elif tag == "dl:reactionSmiles":
            reaction.identifiers.add(type="REACTION_CXSMILES", value=child.text, is_mapped=True)
        elif tag == "cml:productList":
            outcome = reaction.outcomes.add()  # Add a single outcome.
            for product in child:
                parse_product(product, outcome.products.add())
        elif tag in ["cml:reactantList", "cml:spectatorList"]:
            for compound in child:
                molecule_id = get_molecule_id(compound)
                reaction_input = reaction.inputs[component_map[molecule_id]]
                parse_reactant(compound, reaction_input.components.add())
        elif tag == "dl:reactionActionList":
            for action in child:
                parse_conditions(action, reaction)
                parse_workup(action, reaction)
        else:
            raise NotImplementedError(child)
    # We surely didn't capture everything.
    reaction.conditions.conditions_are_dynamic = True
    reaction.conditions.details = "See reaction.notes.procedure_details."
    return reaction


def parse_source(root: ElementTree.Element, reaction: reaction_pb2.Reaction):
    """Adds provenance information to a Reaction."""
    for child in root:
        tag = get_tag(child)
        if tag == "dl:documentId":
            reaction.provenance.patent = child.text
        elif tag == "dl:paragraphText":
            reaction.notes.procedure_details = child.text
        elif tag in ["dl:headingText", "dl:paragraphNum"]:
            continue  # Ignored.
        else:
            raise NotImplementedError(child)


def parse_product(root: ElementTree.Element, product_compound: reaction_pb2.ProductCompound):
    """Adds product information to a ProductCompound."""
    role = root.attrib.get("role")
    if role:
        product_compound.reaction_role = REACTION_ROLES[role]
    for child in root:
        tag = get_tag(child)
        if tag == "dl:entityType":
            pass
        elif tag == "cml:molecule":
            parse_molecule(child, product_compound)
        elif tag == "cml:amount":
            parse_product_amount(child, product_compound)
        elif tag == "cml:identifier":
            parse_identifier(child, product_compound)
        elif tag == "dl:state":
            texture = PRODUCT_STATES.get(child.text.lower(), reaction_pb2.ProductCompound.Texture.CUSTOM)
            product_compound.texture.type = texture
            product_compound.texture.details = child.text
        elif tag == "dl:appearance":
            product_compound.isolated_color = child.text
        else:
            raise NotImplementedError(child)


def parse_molecule(
    root: ElementTree.Element,
    compound: Union[reaction_pb2.Compound, reaction_pb2.ProductCompound],
):
    """Adds NAME identifiers to a Compound."""
    for child in root:
        tag = get_tag(child)
        if tag in ["cml:name", "dl:nameResolved"]:
            compound.identifiers.add(type="NAME", value=child.text)
        else:
            raise NotImplementedError(child)


def parse_product_amount(root: ElementTree.Element, product_compound: reaction_pb2.ProductCompound):
    """Adds amount information to a ProductCompound."""
    property_type = root.attrib[f'{{{NAMESPACES["dl"]}}}propertyType']
    if "PERCENTYIELD" in property_type:
        match = re.search(r"(\d+\.?\d*)", root.text)
        if not match:
            logger.debug(f"AMOUNT: {root.text}")
            return
        measurement = product_compound.measurements.add()
        measurement.type = measurement.YIELD
        measurement.percentage.value = float(match.group(1))
    else:
        measurement = product_compound.measurements.add()
        measurement.type = measurement.AMOUNT
        parse_amount(root, measurement)
    measurement.details = property_type


def parse_identifier(
    root: ElementTree.Element,
    compound: Union[reaction_pb2.Compound, reaction_pb2.ProductCompound],
):
    """Adds a SMILES or INCHI identifier to a Compound."""
    kind = root.attrib["dictRef"]
    value = root.attrib["value"]
    if kind == "cml:smiles":
        compound.identifiers.add(type="SMILES", value=value)
    elif kind == "cml:inchi":
        compound.identifiers.add(type="INCHI", value=value)
    else:
        raise NotImplementedError(kind)


def parse_reactant(root: ElementTree.Element, compound: reaction_pb2.Compound):
    """Populates an input Compound."""
    role = root.attrib.get("role")
    if role:
        compound.reaction_role = REACTION_ROLES[role]
    for child in root:
        tag = get_tag(child)
        if tag == "dl:entityType":
            if child.text not in ["exact", "chemicalClass", "definiteReference"]:
                raise NotImplementedError(child.text)
        elif tag == "cml:molecule":
            parse_molecule(child, compound)
        elif tag == "cml:amount":
            parse_amount(child, compound)
        elif tag == "cml:identifier":
            parse_identifier(child, compound)
        elif tag in ["dl:state", "dl:appearance"]:
            continue  # Not supported by Compound.
        else:
            raise NotImplementedError(child)


def parse_amount(
    root: ElementTree.Element,
    compound: Union[reaction_pb2.Compound, reaction_pb2.ProductMeasurement],
):
    """Parses an amount."""
    property_type = root.attrib[f'{{{NAMESPACES["dl"]}}}propertyType']
    if property_type in ["MOLARITY", "PH"]:
        return
    if compound.amount.WhichOneof("kind") not in ["mass", "volume"]:
        # Don't overwrite existing Mass or Volume with Moles.
        try:
            amount = resolve_units(root.text)
        except (KeyError, ValueError) as error:
            logger.debug(f'UNITS: {error} ("{root.text}")')
            return
        if isinstance(amount, reaction_pb2.Mass):
            compound.amount.mass.CopyFrom(amount)
        elif isinstance(amount, reaction_pb2.Moles):
            compound.amount.moles.CopyFrom(amount)
        elif isinstance(amount, reaction_pb2.Volume):
            compound.amount.volume.CopyFrom(amount)
        else:
            raise NotImplementedError(amount)


def parse_conditions(root: ElementTree.Element, reaction: reaction_pb2.Reaction):
    """Parses reaction conditions."""
    del reaction  # Unused.
    if not root.findall("cml:chemical", namespaces=NAMESPACES):
        return  # Refers to an added component; is a workup.
    action = root.attrib["action"]
    del action  # Unused.
    return  # TODO(kearnes): Implement this?


def parse_workup(root: ElementTree.Element, reaction: reaction_pb2.Reaction):
    """Parses a workup step."""
    if root.findall("dl:chemical", namespaces=NAMESPACES):
        return  # Refers to an input component; not a workup.
    action = root.attrib["action"]
    components = []
    for component in root.findall("cml:chemical", namespaces=NAMESPACES):
        entity_type = component.find("dl:entityType", namespaces=NAMESPACES)
        assert entity_type is not None  # Type hint.
        if entity_type.text in ["exact", "chemicalClass", "definiteReference"]:
            components.append(component)
    if action == "Dry" and components:
        action = "Dry with material"
    details = root.find("dl:phraseText", namespaces=NAMESPACES)
    assert details is not None  # Type hint.
    details = details.text
    if action in ["Purify", "Recover"] and details:
        # Make some actions more specific.
        # TODO(kearnes): This could be expanded.
        if "distill" in details.lower():
            action = "Distill"
        elif "filtration" in details.lower():
            action = "Filter"
    workup_type = WORKUP_TYPES.get(action, reaction_pb2.ReactionWorkup.CUSTOM)
    if workup_type == reaction_pb2.ReactionWorkup.CUSTOM:
        logger.debug(f'CUSTOM workup "{action}": {details}')
    workup = reaction.workups.add(type=workup_type)
    for component in components:
        compound = workup.input.components.add()
        parse_reactant(component, compound)
        compound.reaction_role = reaction_pb2.ReactionRole.WORKUP
    for child in root:
        tag = get_tag(child)
        if tag == "cml:chemical":
            pass  # Components are handled all together above.
        elif tag == "dl:phraseText":
            workup.details = child.text
        elif tag == "dl:parameter":
            try:
                parse_parameter(child, workup)
            except ValueError as error:
                logger.debug(f"PARAMETER: {ElementTree.tostring(child)}: {error}")
        elif tag in ["dl:atmosphere"]:
            pass  # Not supported by ReactionWorkup.
        else:
            raise NotImplementedError(child)


def parse_parameter(root: ElementTree.Element, workup: reaction_pb2.ReactionWorkup):
    """Parses a workup value."""
    kind = root.attrib["propertyType"]
    if kind == "Time":
        if root.text in ["overnight", "Overnight"]:
            workup.duration.value = 8
            workup.duration.precision = 8
            workup.duration.units = reaction_pb2.Time.HOUR
        else:
            try:
                workup.duration.CopyFrom(resolve_units(root.text))
            except (KeyError, ValueError) as error:
                logger.debug(f'UNITS: {error} ("{root.text}")')
    elif kind == "Temperature":
        if root.text.lower() in [
            "room temperature",
            "room temp",
            "rt",
            "r.t",
            "r.t.",
            "ambient temperature",
            "ambient temp",
        ]:
            workup.temperature.control.type = reaction_pb2.TemperatureConditions.TemperatureControl.AMBIENT
        else:
            value = root.text.rstrip(".").replace("° ", "°")
            try:
                temperature = resolve_units(value)
                if temperature.units == temperature.CELSIUS and temperature.value < -274:
                    raise ValueError("bad temperature")
                if temperature.precision < 0:
                    raise ValueError("bad temperature")
                workup.temperature.setpoint.CopyFrom(temperature)
            except (KeyError, ValueError) as error:
                logger.debug(f'UNITS: {error} ("{root.text}")')
    elif kind in ["Frequency", "Length", "Pressure"]:
        return  # Not supported in ReactionWorkup.
    else:
        raise NotImplementedError(kind)


def clean_reaction(reaction: reaction_pb2.Reaction):
    """Cleans a reaction so it will pass validations."""
    # Add a placeholder amount to components with no amount information.
    empty_amount = reaction_pb2.Moles(value=0, precision=1, units="MOLE")
    components = message_helpers.find_submessages(reaction, reaction_pb2.Compound)
    for component in components:
        if not component.amount.WhichOneof("kind"):
            component.amount.moles.CopyFrom(empty_amount)
    # Adjust identifier types as needed.
    identifiers = message_helpers.find_submessages(reaction, reaction_pb2.CompoundIdentifier)
    for identifier in identifiers:
        output = validations.validate_message(identifier, raise_on_error=False)
        if output.errors:
            old_type = reaction_pb2.CompoundIdentifier.IdentifierType.Name(identifier.type)
            identifier.details = f"Originally defined as {old_type}"
            identifier.type = reaction_pb2.CompoundIdentifier.CUSTOM
    # Adjust workup types as needed.
    for workup in reaction.workups:
        output = validations.validate_message(workup, raise_on_error=False)
        if output.errors:
            workup.type = reaction_pb2.ReactionWorkup.CUSTOM
    # Move some workups into ReactionConditions/ReactionOutcome.
    workups = [workup for workup in reaction.workups]  # pylint: disable=unnecessary-comprehension
    del reaction.workups[:]
    temperature_conditions = False
    stirring_conditions = False
    reaction_time = False
    for workup in workups:
        if workup.type == workup.STIRRING and not stirring_conditions:
            reaction.conditions.stirring.CopyFrom(workup.stirring)
            reaction.conditions.stirring.type = reaction_pb2.StirringConditions.StirringMethodType.CUSTOM
            reaction.conditions.stirring.details = workup.details
            stirring_conditions = True
            if workup.HasField("temperature") and not temperature_conditions:
                reaction.conditions.temperature.CopyFrom(workup.temperature)
                temperature_conditions = True
            if workup.HasField("duration") and not reaction_time:
                reaction.outcomes[0].reaction_time.CopyFrom(workup.duration)
                reaction_time = True
        elif workup.type == workup.WAIT and workup.HasField("duration") and not reaction_time:
            reaction.outcomes[0].reaction_time.CopyFrom(workup.duration)
            reaction_time = True
        elif workup.type == workup.TEMPERATURE and workup.HasField("temperature") and not temperature_conditions:
            reaction.conditions.temperature.CopyFrom(workup.temperature)
            temperature_conditions = True
        else:
            reaction.workups.add().CopyFrom(workup)


def run(filename: str) -> tuple[list[reaction_pb2.Reaction], list[reaction_pb2.Reaction]]:
    """Parses reactions from a single CML file."""
    RDLogger.DisableLog("rdApp.*")  # Disable RDKit logging.
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    reactions = []
    failures = []
    for reaction_cml in root.iterfind("cml:reaction", namespaces=NAMESPACES):
        try:
            reaction = parse_reaction(reaction_cml)
        except (KeyError, NotImplementedError) as error:
            raise ValueError(ElementTree.dump(reaction_cml)) from error
        event = reaction_pb2.RecordEvent(
            time={"value": str(datetime.datetime.now())},
            person={
                "username": "skearnes",
                "name": "Steven Kearnes",
                "orcid": "0000-0003-4579-4388",
                "organization": "Relay Therapeutics",
                "email": "skearnes@relaytx.com",
            },
        )
        reaction.provenance.record_created.CopyFrom(event)
        reaction.provenance.doi = "10.6084/m9.figshare.5104873.v1"
        clean_reaction(reaction)
        # Check reaction SMILES.
        try:
            message_helpers.validate_reaction_smiles(reaction.identifiers[0].value.split()[0])
        except ValueError:
            failures.append(reaction)
            continue
        output = validations.validate_message(reaction, raise_on_error=False)
        if output.errors:
            message = "\n".join(output.errors)
            raise validations.ValidationError(f"{reaction}\n{message}")
        reactions.append(reaction)
    return reactions, failures


def main(kwargs):
    filenames = sorted(glob.glob(kwargs["--input_pattern"]))
    all_reactions = joblib.Parallel(n_jobs=int(kwargs["--n_jobs"]), verbose=True)(
        joblib.delayed(run)(filename) for filename in filenames
    )
    reactions = []
    failures = []
    for file_reactions, file_failures in all_reactions:
        reactions.extend(file_reactions)
        failures.extend(file_failures)
    dataset = dataset_pb2.Dataset(reactions=reactions, name=kwargs["--name"])
    basenames = [os.path.basename(filename) for filename in filenames]
    dataset.description = f'CML filenames: {",".join(basenames)}'
    if kwargs["--output"] and reactions:
        logger.info(f"Writing {len(reactions)} reactions to {kwargs['--output']}")
        message_helpers.write_message(dataset, kwargs["--output"])
        if failures:
            failure_dataset = dataset_pb2.Dataset(reactions=failures, name=kwargs["--name"])
            message_helpers.write_message(failure_dataset, kwargs["--output"] + ".failures.pb")


if __name__ == "__main__":
    main(docopt.docopt(__doc__))
