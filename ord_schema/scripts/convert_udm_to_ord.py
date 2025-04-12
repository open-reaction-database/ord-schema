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
import glob

import docopt

import xml.etree.ElementTree as ET

from ord_schema.logging import get_logger
# from ord_schema import message_helpers
# from ord_schema import validations
# from ord_schema.proto import dataset_pb2
# from ord_schema.proto import reaction_pb2

from collections import defaultdict
import json

logger = get_logger(__name__)


def main(kwargs):
    # Add dataset name and description parameters
    # Add message helpers and reaction pb2 from build_dataset script
    # Write to pbtxt file
    # Load pbtxt files into ORD database and display in interface
    # Use test UDM files from real source - Github or Pistoia
    # TODO Test with 1 reaction and multiple, same with other structures
    # Tolerate or warn if values are null or missing
    # Fill data into custom columns where no field is appropriate
    # Finish schema matching

    logger.info("Starting conversion from UDM to ORD.")

    # Default input and output file names
    inputfile = "udm_dataset.xml"
    # Output - protobuf
    outputfile = "ord_dataset.json"

    # Pull input and output file names from parameters
    if kwargs["--input"]:
        inputfile = kwargs["--input"]
    if kwargs["--output"]:
        outputfile = kwargs["--output"]

    # Load UDM file and parse XML into a dictionary
    udm_tree = ET.parse(inputfile)
    root = udm_tree.getroot()
    # for child in root:
        # print(child.tag, child.attrib)
    udm_reactions = etree_to_dict(root)

    # Inspect UDM reactions
    # print(str(udm_reactions))

    # Load UDM file and parse XML
    # Convert fields
    # Return errors if UDM format is wrong
    # Return warnings if data has to be omitted if not supported in ORD
    # Todo: What should we do if we need to change data in UDM vs ORD?

    # List of ORD reactions.
    ord_reactions = []

    # Loop over each UDM reaction.
    for reaction in udm_reactions["REACTIONS"]["REACTION"]:
        # Create a blank dictionary for each reaction.
        ord_reaction = dict()

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
        reactant = reaction["VARIATION"]["REACTANT"]

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
        ord_outcomes["products"] = reaction["VARIATION"]["PRODUCT"]
        ord_outcomes["analyses"] = []
        ord_reaction["outcomes"] = ord_outcomes

        # Step 9 of 9: Provenance
        # Publication and patent details, attribution, other metadata

        ord_provenance = dict()
        ord_provenance["experimenter"] = reaction["ORGANISATIONS"][0]["ORGANISATION"]["NAME"]
        ord_provenance["city"] = reaction["ORGANISATIONS"][0]["ORGANISATION"]["ADDRESS"]
        ord_provenance["experiment_start"] = ""
        ord_provenance["doi"] = reaction["CITATIONS"][0]["CITATION"]["DOI"]
        ord_provenance["patent"] = reaction["CITATIONS"][0]["CITATION"]["PATENT_NUMBER"]
        ord_provenance["publication_url"] = ""
        ord_provenance["record_created"] = ""
        ord_provenance["record_modified"] = ""
        ord_provenance["reaction_metadata"] = ""
        ord_provenance["is_mined"] = False
        ord_reaction["provenance"] = ord_provenance

        ord_reactions.append(ord_reaction)

    # Write JSON file for debug
    f = open(outputfile, "w")
    f.write(json.dumps(ord_reactions, indent=4, separators=(', ', ': ')))
    f.close()

    logger.info("Conversion completed successfully! Check ord_data.json file for errors.")

    # TODO - Run script with xml file from UDM repo passed in and see errors
    # TODO - Allow name and description parameters to create an ORD dataset

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

