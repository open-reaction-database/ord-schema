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
    convert_udm_to_ord.py --input=<str> --output=<str>

Options:
    --input=<str>           XML filename in UDM format
    --output=<str>          Output Dataset filename (*.pbtxt)
    --debugfile=<str>       Output Debug filename (*.json)
    --name=<str>            Name for this dataset
    --description=<str>     Description for this dataset
    --no-validate           If set, do run validations on reactions
"""
# import glob

import docopt

import xml.etree.ElementTree as ET

from ord_schema.logging import get_logger

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

from collections import defaultdict
import json
import os.path

import message_helpers

logger = get_logger(__name__)


def main(kwargs):
    # Add message helpers and reaction pb2 from build_dataset script
    # Write to pbtxt file
    # Use test UDM files from real source
    # Test with 1 reaction and multiple, same with other structures
    # Tolerate or warn if values are null or missing
    # Fill data into custom columns where no field is appropriate
    # Finish schema matching

    # reactions.append(message_helpers.load_message(filename, reaction_pb2.Reaction))

    logger.info("Starting conversion from UDM v6.0.0 to ORD.")

    logger.info("Important message - The Open Reaction Database repository, ord-data, uses the CC-BY-SA license"
                " for all data. "
                "Please do not push your converted dataset to ord-data unless you have the authority to "
                "change the license on the data to CC-BY-SA.")

    # Default input and output file names
    inputfile = "udm_dataset.xml"
    # Output - protobuf
    outputfile = "ord_dataset.json"

    # Pull input and output file names from parameters
    if kwargs["--input"]:
        inputfile = kwargs["--input"]

    if not os.path.isfile(inputfile):
        logger.info("Error - Conversion failed. Please check that the --input file exists at the specified location.")
        exit(1)

    # Load UDM file and parse XML into a dictionary
    udm_tree = ET.parse(inputfile)
    root = udm_tree.getroot()
    # for child in root:
        # print(child.tag, child.attrib)
    udm_reactions = etree_to_dict(root)

    if "UDM" not in udm_reactions:
        logger.info("Error - Input file is not UDM format - please ensure the <UDM> tag is used as the root of your "
                    ".xml file.")
        exit(1)

    udm_reactions = udm_reactions["UDM"]

    dataset_name = ""
    dataset_description = ""
    if "LEGAL" in udm_reactions:
        if "TITLE" in udm_reactions["LEGAL"]:
            outputfile = udm_reactions["LEGAL"]["TITLE"] + ".json"
            dataset_name = udm_reactions["LEGAL"]["TITLE"]
        if "DOI" in udm_reactions["LEGAL"]:
            dataset_description = udm_reactions["LEGAL"]["DOI"]

    if kwargs["--name"]:
        dataset_name = kwargs["--name"]

    if kwargs["--description"]:
        dataset_description = kwargs["--description"]

    # If converter user specified a different output filename:
    if kwargs["--output"]:
        outputfile = kwargs["--output"]

    # Inspect UDM reactions
    # print(str(udm_reactions))

    # Load UDM file and parse XML
    # Convert fields
    # Return errors if UDM format is wrong
    # Return warnings if data has to be omitted if not supported in ORD
    # Todo: What should we do if we need to change data in UDM vs ORD?

    # List of ORD reactions.
    ord_reactions = []
    pb2_reactions = []

    # Loop over each UDM reaction.
    # Loop over VARIATIONs.
    for reaction in udm_reactions["REACTIONS"]["REACTION"]:
    # Inner loop over udm_reactions["REACTIONS"]["REACTION"]["VARIATION"]
        # Create a blank dictionary for each reaction.
        ord_reaction = dict()
        pb2_reaction = message_helpers.create_message("Reaction")

        # Initialize.
        ord_reaction["identifiers"] = []
        ord_reaction["inputs"] = []
        ord_reaction["setup"] = []
        ord_reaction["conditions"] = []
        ord_reaction["notes"] = []
        ord_reaction["observations"] = []
        ord_reaction["workups"] = []
        ord_reaction["outcomes"] = []
        ord_reaction["provenance"] = []

        # Step 1 of 9: Identifiers
        # Used to identify molecules

        # for each identifier, add to list
        ord_identifiers = []

        # May be multiple RXNSTRUCTs
        if "RXNSTRUCTURE" in reaction:
            rxnstructure = reaction["RXNSTRUCTURE"]
            udmformat = rxnstructure["format"]
            # Default reaction identifier type will be unspecified.
            ordtype = 0
            orddetails = ""

            # For types not supported in ORD, the type will be CUSTOM and the details will contain the name of the type.
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

            ord_identifier = dict()
            ord_identifier["type"] = ordtype
            ord_identifier["details"] = orddetails
            # The actual identifier for the reaction
            ord_identifier["value"] = rxnstructure["value"]
            ord_identifier["is_mapped"] = False
            ord_identifiers.append(ord_identifier)

            ord_reaction["identifiers"] = ord_identifiers

        # Step 2 of 9: Inputs
        # Used for reaction inputs, i.e. reactants, catalysts, solvents etc.

        # Reactions in UDM are more VARIATIONS on a reaction conducted during EXPERIMENTS
        # Roles for Compounds and Products together in one enum for ORD.
        # UNSPECIFIED, REACTANT, REAGENT, SOLVENT, CATALYST, WORKUP, INTERNAL_STANDARD, AUTHENTIC_STANDARD
        # PRODUCT, BYPRODUCT, SIDE_PRODUCT
        # Of these, REACTANT, REAGENT, SOLVENT, CATALYST, and PRODUCT are represented in UDM as separate data structures.
        # The thing is, they allow MULTIPLE of each one in each reaction with MULTIPLE compounds.
        ord_inputs = []

        # may be multiple variations on one reaction
        # variation = reaction["VARIATION"]
        # may be multiple reactants on one variation
        # reactant = variation["REACTANT"]
        # print(reaction["VARIATION"])
        for variation in reaction["VARIATION"]:
            reactant = None
            reagent = None
            catalyst = None
            solvent = None

            if "REACTANT" in variation:
                reactant = variation["REACTANT"]
            if "REAGENT" in variation:
                reagent = variation["REAGENT"]
            if "CATALYST" in variation:
                catalyst = variation["CATALYST"]
            if "SOLVENT" in variation:
                solvent = variation["SOLVENT"]

            if reactant:
                # compound = reaction_pb2.Compound()
                compound = dict()
                # TODO: Check the following translation from MOLECULE to identifiers.
                compound["identifiers"] = reactant["MOLECULE"]
                compound["amount"] = reactant["AMOUNT"]
                compound["reaction_role"] = []
                compound["is_limiting"] = []
                compound["preparations"] = []
                compound["source"] = []
                compound["features"] = []
                compound["analyses"] = []
                compound["texture"] = []

                # Missing from UDM - addition details such as time, speed, duration, and device.
                ord_input = dict()
                ord_input["components"] = []
                ord_input["crude_components"] = []
                ord_input["addition_order"] = 0
                ord_input["addition_time"] = ""
                ord_input["addition_speed"] = 0
                ord_input["addition_duration"] = ""
                ord_input["flow_rate"] = 0
                ord_input["addition_device"] = 0
                ord_input["addition_temperature"] = ""
                ord_input["texture"] = ""

                ord_inputs.append(ord_input)

                ord_reaction["inputs"] = ord_inputs

        # Step 3 of 9: Setup
        # Describes the reaction setup.

        ord_setup = dict()
        ord_setup["vessel"] = ""
        ord_setup["is_automated"] = True
        ord_setup["automation_platform"] = ""
        ord_setup["automation_code"] = ""
        ord_setup["environment"] = ""
        ord_reaction["setup"] = ord_setup

        # Step 4 of 9: Conditions
        # Describes the reaction conditions.

        # Add extra data from UDM concatenate - into free text custom field.
        ord_conditions = dict()
        if "CONDITIONS" in reaction:
            ord_conditions["temperature"] = reaction["CONDITIONS"][0]["CONDITION_GROUP"]["TEMPERATURE"]
            ord_conditions["pressure"] = reaction["CONDITIONS"][0]["CONDITION_GROUP"]["PRESSURE"]
            ord_conditions["stirring"] = reaction["CONDITIONS"][0]["CONDITION_GROUP"]["STIRRING"]
            ord_conditions["illumination"] = ""
            ord_conditions["electrochemistry"] = ""
            ord_conditions["flow"] = ""
            ord_conditions["reflux"] = reaction["CONDITIONS"][0]["CONDITION_GROUP"]["REFLUX"]
            ord_conditions["ph"] = reaction["CONDITIONS"][0]["CONDITION_GROUP"]["PH"]
            ord_conditions["conditions_are_dynamic"] = True
            ord_conditions["details"] = ""
            ord_reaction["conditions"] = ord_conditions

        # Step 5 of 9: Notes
        # Extra notes about the reaction.

        # TODO: Remove undefined properties if not determined by original data
        ord_notes = dict()
        ord_notes["is_heterogeneous"] = True
        ord_notes["forms_precipitate"] = False
        ord_notes["is_exothermic"] = False
        ord_notes["offgasses"] = False
        ord_notes["is_sensitive_to_moisture"] = False
        ord_notes["is_sensitive_to_oxygen"] = False
        ord_notes["is_sensitive_to_light"] = False
        ord_notes["safety_notes"] = ""
        ord_notes["procedure_details"] = ""
        ord_reaction["notes"] = ord_notes

        # Step 6 of 9: Observations
        # Notes about observations during the reaction.

        ord_observations = dict()
        ord_observations["time"] = ""
        if "COMMENT" in reaction["VARIATION"]:
            ord_observations["comment"] = reaction["VARIATION"]["COMMENT"]
        ord_observations["image"] = ""
        ord_reaction["observations"] = ord_observations

        # Step 7 of 9: Workups
        # Reaction workups

        ord_workups = dict()
        ord_workups["type"] = 0
        ord_workups["details"] = ""
        ord_workups["duration"] = ""
        ord_workups["input"] = ""
        ord_workups["amount"] = ""
        ord_workups["temperature"] = ""
        ord_workups["keep_phase"] = ""
        ord_workups["stirring"] = ""
        ord_workups["target_ph"] = 0.0
        ord_workups["is_automated"] = False
        ord_reaction["workups"] = ord_workups

        # Step 8 of 9: Outcomes
        # The outcomes of the reaction - time taken, the products etc.
        ord_outcomes = dict()
        ord_outcomes["reaction_time"] = ""
        ord_outcomes["conversion"] = ""
        if "PRODUCT" in reaction["VARIATION"]:
            ord_outcomes["products"] = reaction["VARIATION"]["PRODUCT"]
        ord_outcomes["analyses"] = []
        ord_reaction["outcomes"] = ord_outcomes

        # Step 9 of 9: Provenance
        # Publication and patent details, attribution, other metadata

        ord_provenance = dict()
        ord_provenance_experimenter = dict()
        if "LEGAL" in udm_reactions:
            ord_provenance_experimenter["organization"] = udm_reactions["LEGAL"]["PRODUCER"]

        if "SCIENTIST" in reaction["VARIATION"]:
            ord_provenance_experimenter["name"] = reaction["VARIATION"]["SCIENTIST"]

        ord_provenance["experimenter"] = ord_provenance_experimenter
        ord_record_created = dict()
        ord_record_created["person"] = ord_provenance_experimenter
        ord_provenance["record_created"] = ord_record_created

        if "ORGANISATIONS" in reaction:
            ord_provenance["city"] = reaction["ORGANISATIONS"][0]["ORGANISATION"]["ADDRESS"]

        ord_provenance["experiment_start"] = ""

        if "LEGAL" in udm_reactions and "DOI" in udm_reactions["LEGAL"] and "CITATIONS" not in reaction:
            ord_provenance["doi"] = udm_reactions["LEGAL"]["DOI"]

        if "CITATIONS" in reaction:
            ord_provenance["doi"] = reaction["CITATIONS"][0]["CITATION"]["DOI"]
            ord_provenance["patent"] = reaction["CITATIONS"][0]["CITATION"]["PATENT_NUMBER"]

        ord_provenance["publication_url"] = ""

        if "CREATION_DATE" in reaction["VARIATION"]:
            ord_provenance["record_created"] = reaction["VARIATION"]["CREATION_DATE"]

        if "MODIFICATION_DATE" in reaction["VARIATION"]:
            ord_provenance["record_modified"] = reaction["VARIATION"]["MODIFICATION_DATE"]

        ord_provenance["reaction_metadata"] = ""
        ord_provenance["is_mined"] = False
        ord_reaction["provenance"] = ord_provenance

        ord_reactions.append(ord_reaction)
        pb2_reactions.append(pb2_reaction)

    dataset = dataset_pb2.Dataset(name=dataset_name, description=dataset_description, reactions=pb2_reactions)
    if not kwargs["--no-validate"]:
        validations.validate_datasets({"_COMBINED": dataset})
    message_helpers.write_message(dataset, kwargs["--output"])

    f = open(outputfile, "w")
    f.write(json.dumps(ord_reactions, indent=4, separators=(', ', ': ')))
    f.close()

    # Write JSON file for debug
    debugfile = "ORD_UDM_CONVERSION_" + inputfile + ".debug.json"
    f = open(debugfile, "w")
    f.write(json.dumps({}, indent=4, separators=(', ', ': ')))
    f.close()

    logger.info("Conversion completed successfully! Check ord_data.json file for errors.")

    return


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


# Shows the Usage information if parameters not provided (self documenting)
if __name__ == "__main__":
    main(docopt.docopt(__doc__))

