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

Given a UDM file, convert to .pbtxt for ORD

Usage:
    convert_udm_to_ord.py --input=<str> [--output=<str>] [--name=<str>] [--description=<str>] [--no-validate]

Options:
    --input=<str>           XML filename in UDM format
    --output=<str>          Output Dataset filename (*.pbtxt)
    --name=<str>            Name for this dataset
    --description=<str>     Description for this dataset
    --no-validate           If set, do run validations on reactions
"""
# import glob

import json
import os.path
import xml.etree.ElementTree as ET
from collections import defaultdict

import docopt
from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.logging import get_logger
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema import updates

import random

logger = get_logger(__name__)


def main(kwargs):
    logger.info("Starting conversion from UDM v6.0.0 to ORD.")
    logger.info("*** Important message - The Open Reaction Database repository, ord-data, uses the CC-BY-SA license"
                " for all data. "
                "Please do not push your converted dataset to ord-data unless you have the authority to "
                "change the license on the data to CC-BY-SA. ***")

    # Default input and output file names
    inputfile = "udm_dataset.xml"
    outputfile = "ord_dataset.pbtxt"

    # Pull input filename from argument parameters
    if kwargs["--input"]:
        inputfile = kwargs["--input"]

    # If the input file is not found, exit with an error message
    if not os.path.isfile(inputfile):
        logger.info("Error - Conversion failed. Please check that the --input file exists at the specified location.")
        exit(1)

    # Load UDM file and parse XML into a dictionary
    udm_tree = ET.parse(inputfile)
    root = udm_tree.getroot()
    # for child in root:
        # print(child.tag, child.attrib)
    udm_reactions = etree_to_dict(root)

    # <UDM> is expected to be the root element of the UDM file - if it isn't, we exit with an error.
    if "UDM" not in udm_reactions:
        logger.info("Error - Input file is not UDM format - please ensure the <UDM> tag is used as the root of your "
                    ".xml file.")
        exit(1)

    # Access the data in the .xml file.
    udm_reactions = udm_reactions["UDM"]

    # Determine dataset name and description - if given, the <TITLE> will be used as the name.
    # If a global DOI is given, this will be included in the dataset description.
    dataset_name = ""
    dataset_description = ""
    if "LEGAL" in udm_reactions:
        if "TITLE" in udm_reactions["LEGAL"]:
            outputfile = udm_reactions["LEGAL"]["TITLE"] + ".pbtxt"
            dataset_name = udm_reactions["LEGAL"]["TITLE"]
        if "DOI" in udm_reactions["LEGAL"]:
            dataset_description = "UDM dataset DOI: " + udm_reactions["LEGAL"]["DOI"]

    # Set the dataset name if explicitly provided (overrides what's in the UDM file)
    if kwargs["--name"]:
        dataset_name = kwargs["--name"]

    # Set the dataset description if explicitly provided (overrides what's in the UDM file)
    if kwargs["--description"]:
        dataset_description = kwargs["--description"]

    # Set the output file name if explicitly provided (overrides what's in the UDM file)
    if kwargs["--output"]:
        outputfile = kwargs["--output"]

    # Inspect UDM reactions
    # print(str(udm_reactions))

    # List of ORD reactions-to-be
    ord_reactions = []
    pb2_reactions = []

    # This is unlikely to be hit, but just in case
    if "REACTIONS" not in udm_reactions:
        logger.info("Error - <REACTIONS> element not found in input .xml file.")
        exit(1)

    # If there is just one reaction in the entire file, we have to change our approach.
    # Create a list of every reaction in the file.
    udm_reactions_list = []
    if isinstance(udm_reactions["REACTIONS"]["REACTION"], dict):
        udm_reactions_list.append(udm_reactions["REACTIONS"]["REACTION"])
    else:
        for reaction in udm_reactions["REACTIONS"]["REACTION"]:
            udm_reactions_list.append(reaction)

    all_molecules = {}
    for molecule in udm_reactions["MOLECULES"]["MOLECULE"]:
        mol_data = dict()
        mol_data["name"] = molecule["NAME"]
        if "MOLSTRUCTURE" in molecule:
            mol_data["molblock"] = molecule["MOLSTRUCTURE"]
        all_molecules[molecule["@ID"]] = mol_data

    # Loop over each UDM reaction. UDM reactions may have one or more Variations.
    # We will treat each Variation as a separate ORD Reaction.
    # May be 1 or multiple, with null checks
    # Return errors if UDM format is wrong
    # Return warnings if data has to be omitted if not supported in ORD
    for reaction in udm_reactions_list:
        # Create a blank dictionary for each reaction.
        pb2_reaction = reaction_pb2.Reaction()

        # Step 1 of 9: Identifiers
        # Used to identify molecules

        pb2_inputs = []

        if "REACTANT_ID" in reaction:
            molinput = pb2_reaction.inputs[""]
            reactant_list = reaction["REACTANT_ID"]
            for reactant_id in reactant_list:
                molcomponent = molinput.components.add()
                molcomponent.identifiers.add(type="CUSTOM", details="REACTANT_ID from UDM")
                molcomponent.identifiers[0].value = reactant_id

        # May be multiple RXNSTRUCTs
        if "RXNSTRUCTURE" in reaction:
            rxnstructures = []
            if isinstance(reaction["RXNSTRUCTURE"], dict):
                rxnstructures.append(reaction["RXNSTRUCTURE"])
            else:
                for structure in reaction["RXNSTRUCTURE"]:
                    rxnstructures.append(structure)

            for structure in rxnstructures:
                udmformat = 1
                if "format" in structure:
                    udmformat = structure["format"]

                # Default reaction identifier type will be unspecified.
                ordtype = 0
                orddetails = ""

                # For types not supported in ORD, the type will be CUSTOM
                # and the details will contain the name of the type.
                if udmformat == "cdxml":
                    ordtype = 1
                    orddetails = "cdxml"
                elif udmformat == "rinchi":
                    ordtype = 5
                elif udmformat == "rsmiles":
                    ordtype = 2
                elif udmformat == "rxn":
                    ordtype = 1
                    orddetails = "rxn"

                udmvalue = None
                if "value" in structure:
                    udmvalue = structure["value"]

                pb2_reaction.identifiers.add(type=ordtype, details=orddetails, value=udmvalue)

        # Step 2 of 9: Inputs
        # Used for reaction inputs, i.e. reactants, catalysts, solvents etc.

        # Reactions in UDM are more VARIATIONS on a reaction conducted during EXPERIMENTS
        # Roles for Compounds and Products together in one enum for ORD.
        # UNSPECIFIED, REACTANT, REAGENT, SOLVENT, CATALYST, WORKUP, INTERNAL_STANDARD, AUTHENTIC_STANDARD
        # PRODUCT, BYPRODUCT, SIDE_PRODUCT
        # Of these, REACTANT, REAGENT, SOLVENT, CATALYST, and PRODUCT are represented in UDM as separate data structures.
        # They allow MULTIPLE of each one in each reaction with MULTIPLE compounds.
        ord_inputs = []

        # If there is one variation, it is a dict
        variations = []
        if "VARIATION" in reaction:
            if isinstance(reaction["VARIATION"], dict):
                variations.append(reaction["VARIATION"])
            else:
                for variation in reaction["VARIATION"]:
                    variations.append(variation)

        # If there are multiple variations, each will be its own reaction.
        for variation in variations:
            pb2_reaction = reaction_pb2.Reaction()

            reactants = []
            if "REACTANT" in variation:
                if isinstance(variation["REACTANT"], dict):
                    reactants.append(variation["REACTANT"])
                else:
                    for reactant in variation["REACTANT"]:
                        reactants.append(reactant)

            reagents = []
            if "REAGENT" in variation:
                if isinstance(variation["REAGENT"], dict):
                    reagents.append(variation["REAGENT"])
                else:
                    for reagent in variation["REAGENT"]:
                        reagents.append(reagent)

            catalysts = []
            if "CATALYST" in variation:
                if isinstance(variation["CATALYST"], dict):
                    catalysts.append(variation["CATALYST"])
                else:
                    for catalyst in variation["CATALYST"]:
                        catalysts.append(catalyst)

            solvents = []
            if "SOLVENT" in variation:
                if isinstance(variation["SOLVENT"], dict):
                    solvents.append(variation["SOLVENT"])
                else:
                    for solvent in variation["SOLVENT"]:
                        solvents.append(solvent)

            for reactant in reactants:
                pb2_compound = reaction_pb2.Compound()
                compound = dict()
                # TODO: Check the following translation from MOLECULE to identifiers.
                compound["identifiers"] = reactant["MOLECULE"]
                # compound["amount"] = reactant["AMOUNT"]
                compound["reaction_role"] = []
                compound["is_limiting"] = []
                compound["preparations"] = []
                compound["source"] = []
                compound["features"] = []
                compound["analyses"] = []
                compound["texture"] = []

                if "MOLECULE" in reactant:
                    molval = all_molecules.get(reactant["MOLECULE"]["@MOL_ID"])
                    # pb2_reaction.identifiers.add(type=0, details='', value=molval)
                    molinput = pb2_reaction.inputs[reactant["MOLECULE"]["@MOL_ID"]]
                    molcomponent = molinput.components.add()
                    if "molblock" in molval:
                        molcomponent.identifiers.add(type="MOLBLOCK", details="MOLECULE -> MOLSTRUCTURE from UDM")
                        molcomponent.identifiers[0].value = molval["molblock"]
                    else:
                        molcomponent.identifiers.add(type="CUSTOM")
                        if molval is not None:
                            molcomponent.identifiers[0].value = molval["name"]
                if "MOLECULE" in reactant and "NAME" in reactant["MOLECULE"]:
                    pb2_compound.identifiers.add(value=reactant["MOLECULE"]["NAME"])
                if "AMOUNT" in reactant:
                    convertedAmount = float(reactant["AMOUNT"])
                    pb2_compound.amount.mass.value = convertedAmount

                pb2_components = []
                pb2_components.append(pb2_compound)

                pb2_input = reaction_pb2.ReactionInput()
                pb2_input.components.append(pb2_compound)
                pb2_inputs.append(pb2_input)

                # pb2_reaction.inputs.update(key=pb2_input, value=pb2_input)
                # pb2_input.components.add(pb2_compound)
                # pb2_reaction.inputs.update(key=reactant["MOLECULE"]["@MOL_ID"], value=pb2_input)
                # updates.update_reaction(pb2_reaction)

                # Step 4 of 9: Conditions
                # Describes the reaction conditions.

            for reagent in reagents:
                if "MOLECULE" in reagent:
                    molval = all_molecules.get(reagent["MOLECULE"]["@MOL_ID"])
                    # print(molval)
                    # pb2_reaction.identifiers.add(type=0, details='', value=molval)
                    molinput = pb2_reaction.inputs[reagent["MOLECULE"]["@MOL_ID"]]
                    molcomponent = molinput.components.add()
                    if "molblock" in molval:
                        molcomponent.identifiers.add(type="MOLBLOCK", details="MOLECULE -> MOLSTRUCTURE from UDM")
                        molcomponent.identifiers[0].value = molval["molblock"]
                    else:
                        molcomponent.identifiers.add(type="CUSTOM")
                        if molval["name"] is not None:
                            molcomponent.identifiers[0].value = molval["name"]
                    molcomponent.reaction_role = reaction_pb2.ReactionRole.REAGENT

            for catalyst in catalysts:
                if "MOLECULE" in catalyst:
                    molval = all_molecules.get(catalyst["MOLECULE"]["@MOL_ID"])
                    # print(molval)
                    # pb2_reaction.identifiers.add(type=0, details='', value=molval)
                    molinput = pb2_reaction.inputs[catalyst["MOLECULE"]["@MOL_ID"]]
                    molcomponent = molinput.components.add()
                    if "molblock" in molval:
                        molcomponent.identifiers.add(type="MOLBLOCK", details="MOLECULE -> MOLSTRUCTURE from UDM")
                        molcomponent.identifiers[0].value = molval["molblock"]
                    else:
                        molcomponent.identifiers.add(type="CUSTOM")
                        if molval["name"] is not None:
                            molcomponent.identifiers[0].value = molval["name"]
                    molcomponent.reaction_role = reaction_pb2.ReactionRole.CATALYST

            for solvent in solvents:
                if "MOLECULE" in solvent:
                    molval = all_molecules.get(solvent["MOLECULE"]["@MOL_ID"])
                    # print(molval)
                    # pb2_reaction.identifiers.add(type=0, details='', value=molval)
                    molinput = pb2_reaction.inputs[solvent["MOLECULE"]["@MOL_ID"]]
                    molcomponent = molinput.components.add()
                    if "molblock" in molval:
                        molcomponent.identifiers.add(type="MOLBLOCK", details="MOLECULE -> MOLSTRUCTURE from UDM")
                        molcomponent.identifiers[0].value = molval["molblock"]
                    else:
                        molcomponent.identifiers.add(type="CUSTOM")
                        if molval["name"] is not None:
                            molcomponent.identifiers[0].value = molval["name"]
                    molcomponent.reaction_role = reaction_pb2.ReactionRole.SOLVENT


            if "CONDITIONS" in variation and "CONDITION_GROUP" in variation["CONDITIONS"]:
                pb2_conditions = reaction_pb2.ReactionConditions()
                if "TEMPERATURE" in variation["CONDITIONS"]["CONDITION_GROUP"]:
                    if "exact" in variation["CONDITIONS"]["CONDITION_GROUP"]["TEMPERATURE"]:
                        pb2_reaction.conditions.temperature.setpoint.value = float(variation["CONDITIONS"]["CONDITION_GROUP"]["TEMPERATURE"]["exact"])
                if "PRESSURE" in variation["CONDITIONS"]["CONDITION_GROUP"]:
                    if "exact" in variation["CONDITIONS"]["CONDITION_GROUP"]["PRESSURE"]:
                        pb2_reaction.conditions.pressure.setpoint.value = float(variation["CONDITIONS"]["CONDITION_GROUP"]["PRESSURE"]["exact"])
                if "STIRRING" in variation["CONDITIONS"]["CONDITION_GROUP"]:
                    pb2_reaction.conditions.stirring.details = variation["CONDITIONS"]["CONDITION_GROUP"]["STIRRING"]
                if "REFLUX" in variation["CONDITIONS"]["CONDITION_GROUP"]:
                    pb2_reaction.conditions.reflux = variation["CONDITIONS"]["CONDITION_GROUP"]["REFLUX"]
                if "PH" in variation["CONDITIONS"]["CONDITION_GROUP"]:
                    if "exact" in variation["CONDITIONS"]["CONDITION_GROUP"]["PH"]:
                        pb2_reaction.conditions.ph = float(variation["CONDITIONS"]["CONDITION_GROUP"]["PH"]["exact"])

                pb2_setup = reaction_pb2.ReactionSetup()
                if "PREPARATION" in variation["CONDITIONS"]["CONDITION_GROUP"]:
                    pb2_reaction.setup.environment = variation["CONDITIONS"]["CONDITION_GROUP"]["PREPARATION"]

                # pb2_reaction.conditions = pb2_conditions

        # Step 3 of 9: Setup
        # Describes the reaction setup.

        pb2_setup = reaction_pb2.ReactionSetup()
        # pb2_reaction.setup = pb2_setup

        # Step 5 of 9: Notes
        # Extra notes about the reaction.

        pb2_notes = reaction_pb2.ReactionNotes()
        # pb2_notes.safety_notes
        # pb2_notes.procedure_details
        pb2_reaction.notes.procedure_details = ""

        # Step 6 of 9: Observations
        # Notes about observations during the reaction.

        pb2_observations = reaction_pb2.ReactionObservation()
        # pb2_observations.time
        if "VARIATION" in reaction and "COMMENT" in reaction["VARIATION"]:
            pb2_observations.comment = reaction["VARIATION"]["COMMENT"]
        pb2_reaction.observations.append(pb2_observations)

        # Step 7 of 9: Workups
        # Reaction workups

        pb2_workups = reaction_pb2.ReactionWorkup()
        # type, details, duration, input, amount, temperature, stirring, target_ph
        pb2_reaction.workups.append(pb2_workups)

        # Step 8 of 9: Outcomes
        # The outcomes of the reaction - time taken, the products etc.

        pb2_outcomes = reaction_pb2.ReactionOutcome()
        # reaction_time, conversion, analyses

        products = []

        if "PRODUCT_ID" in reaction:
            molinput = pb2_reaction.inputs[""]
            product_list = reaction["PRODUCT_ID"]
            outcome = pb2_reaction.outcomes.add()
            for product_id in product_list:
                molcomponent = outcome.products.add()
                molcomponent.identifiers.add(type="CUSTOM", details="PRODUCT_ID from UDM")
                molcomponent.identifiers[0].value = product_id

        if "VARIATION" in reaction and "PRODUCT" in reaction["VARIATION"]:
            if isinstance(reaction["VARIATION"]["PRODUCT"], dict):
                products.append(reaction["VARIATION"]["PRODUCT"])
            else:
                for product in reaction["VARIATION"]["PRODUCT"]:
                    products.append(product)

            for udm_product in products:
                molecule = udm_product["MOLECULE"]
                outcome = pb2_reaction.outcomes.add()
                product = outcome.products.add()
                if all_molecules.get(molecule["@MOL_ID"]) is not None:
                    molval = all_molecules.get(molecule["@MOL_ID"])
                    if "molblock" in molval:
                        product.identifiers.add(type="MOLBLOCK", details="MOLECULE -> MOLSTRUCTURE from UDM")
                        product.identifiers[0].value = molval["molblock"]
                    else:
                        product.identifiers.add(type="CUSTOM")
                        if molval["name"] is not None:
                            product.identifiers[0].value = molval["name"]
                if "YIELD" in udm_product:
                    product_measurement = product.measurements.add()
                    product_measurement.type = reaction_pb2.ProductMeasurement.ProductMeasurementType.YIELD
                    product_measurement.float_value.value = float(udm_product["YIELD"]["exact"])

        # Step 9 of 9: Provenance
        # Publication and patent details, attribution, other metadata

        ord_provenance = dict()
        pb2_provenance = reaction_pb2.ReactionProvenance()
        ord_provenance_experimenter = dict()
        if "LEGAL" in udm_reactions:
            pb2_reaction.provenance.experimenter.organization = udm_reactions["LEGAL"]["PRODUCER"]
        if "VARIATION" in reaction and "SCIENTIST" in reaction["VARIATION"]:
            pb2_reaction.provenance.experimenter.name = reaction["VARIATION"]["SCIENTIST"]
        if "LEGAL" in udm_reactions:
            pb2_reaction.provenance.record_created.person.organization = udm_reactions["LEGAL"]["PRODUCER"]
        if "VARIATION" in reaction and "SCIENTIST" in reaction["VARIATION"]:
            pb2_reaction.provenance.record_created.person.name = reaction["VARIATION"]["SCIENTIST"]
        if "ORGANISATIONS" in reaction:
            pb2_reaction.provenance.city = reaction["ORGANISATIONS"][0]["ORGANISATION"]["ADDRESS"]

        # experiment start

        if "LEGAL" in udm_reactions and "DOI" in udm_reactions["LEGAL"] and "CITATIONS" not in reaction:
            pb2_reaction.provenance.doi = udm_reactions["LEGAL"]["DOI"]
        elif "CITATIONS" in udm_reactions and "VARIATION" in reaction and "CITATION" in reaction["VARIATION"]:
            variation_citation = reaction["VARIATION"]["CITATION"]
            variation_doi = ""
            for citation in udm_reactions["CITATIONS"]["CITATION"]:
                if isinstance(variation_citation, list):
                    variation_citation = variation_citation[0]
                if citation["@ID"] == variation_citation["@CIT_ID"] and "DOI" in citation:
                    variation_doi = citation["DOI"]
            if variation_doi != "":
                pb2_reaction.provenance.doi = variation_doi

        if "CITATIONS" in reaction:
            pb2_reaction.provenance.doi = reaction["CITATIONS"][0]["CITATION"]["DOI"]
            pb2_reaction.provenance.patent = reaction["CITATIONS"][0]["CITATION"]["PATENT_NUMBER"]

        # publication url is not provided

        if "VARIATION" in reaction and "CREATION_DATE" in reaction["VARIATION"]:
            pb2_reaction.provenance.record_created.time = reaction["VARIATION"]["CREATION_DATE"]

        if "VARIATION" in reaction and "MODIFICATION_DATE" in reaction["VARIATION"]:
            pb2_reaction.provenance.record_modified.time = reaction["VARIATION"]["MODIFICATION_DATE"]

        # ord_provenance["reaction_metadata"] = ""
        pb2_reaction.provenance.is_mined = False

        # Generate REACTION ID
        pb2_reaction.reaction_id = generate_reaction_id()

        # Finally, update the reaction and add it to the list of reactions
        updates.update_reaction(pb2_reaction)
        pb2_reactions.append(pb2_reaction)

    # Create the ORD dataset
    dataset = dataset_pb2.Dataset(name=dataset_name, description=dataset_description, reactions=pb2_reactions)

    # Validate the dataset
    if not kwargs["--no-validate"]:
        validations.validate_datasets({"_COMBINED": dataset})

    # Create the .pbtxt file
    message_helpers.write_message(dataset, outputfile)

    # Write JSON file for debug
    debugfile = "ORD_UDM_CONVERSION_" + inputfile + ".debug.json"
    f = open(debugfile, "w")
    f.write(json.dumps({}, indent=4, separators=(', ', ': ')))
    f.close()

    # Successfully converted UDM .xml file into ORD dataset!
    logger.info("Conversion completed successfully! Check ORD_UDM_CONVERSION_"
    + inputfile + ".debug.json file for errors.")

    exit(0)


# Converts XML tree into Python dictionary.
# https://stackoverflow.com/questions/7684333/converting-xml-to-dictionary-using-elementtree
def etree_to_dict(t):
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v
                     for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v)
                        for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d

def generate_reaction_id():
    random_id_parts = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f']
    new_reaction_id = "ord-"
    for i in range(0, 31):
        new_reaction_id += random.choice(random_id_parts)
    return new_reaction_id

# Shows the Usage information if parameters not provided (self documenting)
if __name__ == "__main__":
    main(docopt.docopt(__doc__))

