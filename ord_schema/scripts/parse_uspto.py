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

Data is at https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873.

See also https://depth-first.com/articles/2019/01/28/the-nextmove-patent-reaction-dataset/.
"""

import glob
import re
from typing import Dict
from xml.etree import ElementTree

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema import units
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input_pattern', None, 'Input pattern for CML files.')
flags.DEFINE_string('output', None, 'Output Dataset filename.')

# XML namespaces.
NAMESPACES = {
    'cmlDict': 'http://www.xml-cml.org/dictionary/cml/',
    'nameDict': 'http://www.xml-cml.org/dictionary/cml/name/',
    'unit': 'http://www.xml-cml.org/unit/',
    'cml': 'http://www.xml-cml.org/schema',
    'dl': 'http://bitbucket.org/dan2097',
}

UNIT_RESOLVER = units.UnitResolver()

WORKUP_TYPES = {
    'Extract': reaction_pb2.ReactionWorkup.EXTRACTION,
    'Stir': reaction_pb2.ReactionWorkup.STIRRING,
    'Wash': reaction_pb2.ReactionWorkup.WASH,
    'Dry with material': reaction_pb2.ReactionWorkup.DRY_WITH_MATERIAL,
    'Remove': reaction_pb2.ReactionWorkup.CUSTOM,
}


def get_tag(element):
    for key, value in NAMESPACES.items():
        if value in element.tag:
            return element.tag.replace(f'{{{value}}}', f'{key}:')
    return element.tag


def get_molecule_id(root):
    molecules = root.findall('cml:molecule', namespaces=NAMESPACES)
    if len(molecules) != 1:
        raise NotImplementedError(len(molecules))
    return molecules[0].attrib['id']


def get_component_map(root: ElementTree.ElementTree) -> Dict[str, str]:
    """Builds a map of components to inputs."""
    reaction_inputs = {}
    # Start with a separate input for each component.
    for component in root.find('cml:reactantList', namespaces=NAMESPACES):
        molecule_id = get_molecule_id(component)
        reaction_inputs[molecule_id] = molecule_id
    for action in root.find('dl:reactionActionList', namespaces=NAMESPACES):
        components = []
        for component in action.findall('dl:chemical', namespaces=NAMESPACES):
            if 'ref' in component.attrib:
                components.append(component.attrib['ref'])
        if components:
            key = '_'.join(components)
            for component in components:
                reaction_inputs[component] = key
    return reaction_inputs


def parse_reaction(root):
    """Parses reaction CML into a Reaction message.

    Input components are grouped by reading the reactionActionList.
    """
    reaction = reaction_pb2.Reaction()
    component_map = get_component_map(root)
    for child in root:
        tag = get_tag(child)
        if tag == 'dl:source':
            parse_source(child, reaction)
        elif tag == 'dl:reactionSmiles':
            reaction.identifiers.add(type='REACTION_SMILES', value=child.text)
        elif tag == 'cml:productList':
            outcome = reaction.outcomes.add()  # Add a single outcome.
            for product in child:
                parse_product(product, outcome.products.add())
        elif tag in ['cml:reactantList', 'cml:spectatorList']:
            for compound in child:
                molecule_id = get_molecule_id(compound)
                reaction_input = reaction.inputs[component_map[molecule_id]]
                parse_reactant(compound, reaction_input.components.add())
        elif tag == 'dl:reactionActionList':
            for action in child:
                parse_action(action, reaction)
        else:
            raise NotImplementedError(child)
    return reaction


def parse_source(root, reaction):
    for child in root:
        tag = get_tag(child)
        if tag == 'dl:documentId':
            reaction.provenance.patent = child.text
        elif tag == 'dl:paragraphText':
            reaction.notes.procedure_details = child.text
        else:
            raise NotImplementedError(child)


def parse_product(root, product_compound):
    role = root.attrib.get('role')
    if role == 'product':
        product_compound.reaction_role = reaction_pb2.ReactionRole.PRODUCT
    elif role:
        raise NotImplementedError(role)
    for child in root:
        tag = get_tag(child)
        if tag == 'dl:entityType':
            if child.text not in ['exact', 'chemicalClass']:
                raise NotImplementedError(child.text)
            continue
        elif tag == 'cml:molecule':
            parse_molecule(child, product_compound)
        elif tag == 'cml:amount':
            parse_product_amount(child, product_compound)
        elif tag == 'cml:identifier':
            parse_identifier(child, product_compound)
        elif tag == 'dl:state':
            parse_product_state(child, product_compound)
        else:
            raise NotImplementedError(child)


def parse_molecule(root, compound):
    for child in root:
        tag = get_tag(child)
        if tag == 'cml:name':
            compound.identifiers.add(type='NAME', value=child.text)
        else:
            raise NotImplementedError(child)


def parse_product_amount(root, product_compound):
    measurement = product_compound.measurements.add()
    match = re.fullmatch(r'([\d.]+)% yield', root.text)
    if match:
        measurement.type = measurement.YIELD
        measurement.percentage.value = float(match.group(1))
    else:
        parse_amount(root, measurement)
        measurement.type = measurement.AMOUNT


def parse_identifier(root, compound):
    kind = root.attrib['dictRef']
    value = root.attrib['value']
    if kind == 'cml:smiles':
        compound.identifiers.add(type='SMILES', value=value)
    elif kind == 'cml:inchi':
        compound.identifiers.add(type='INCHI', value=value)
    else:
        raise NotImplementedError(kind)


def parse_product_state(root, product_compound):
    if root.text == 'oil':
        product_compound.texture.type = product_compound.texture.OIL
    else:
        raise NotImplementedError(root.text)


def parse_reactant(root, compound):
    role = root.attrib.get('role')
    if role == 'reactant':
        compound.reaction_role = reaction_pb2.ReactionRole.REACTANT
    elif role == 'solvent':
        compound.reaction_role = reaction_pb2.ReactionRole.SOLVENT
    elif role:
        raise NotImplementedError(role)
    for child in root:
        tag = get_tag(child)
        if tag == 'dl:entityType':
            if child.text not in ['exact', 'chemicalClass']:
                raise NotImplementedError(child.text)
            continue
        elif tag == 'cml:molecule':
            parse_molecule(child, compound)
        elif tag == 'cml:amount':
            parse_amount(child, compound)
        elif tag == 'cml:identifier':
            parse_identifier(child, compound)
        else:
            raise NotImplementedError(child)


def parse_amount(root, compound):
    if compound.amount.WhichOneof('kind') not in ['mass', 'volume']:
        # Don't overwrite existing Mass or Volume with Moles.
        try:
            amount = UNIT_RESOLVER.resolve(root.text)
        except (KeyError, ValueError) as error:
            logging.info(f'{root.text}: {error}')
            return
        if isinstance(amount, reaction_pb2.Mass):
            compound.amount.mass.CopyFrom(amount)
        elif isinstance(amount, reaction_pb2.Moles):
            compound.amount.moles.CopyFrom(amount)
        elif isinstance(amount, reaction_pb2.Volume):
            compound.amount.volume.CopyFrom(amount)
        else:
            raise NotImplementedError(amount)


def parse_action(root, reaction):
    if root.findall('dl:chemical', namespaces=NAMESPACES):
        return  # Refers to an input component; not a workup.
    action = root.attrib['action']
    if action == 'Dry':
        if root.findall('cml:chemical', namespaces=NAMESPACES):
            action = 'Dry with material'
        else:
            raise NotImplementedError(action)
    if not WORKUP_TYPES[action]:
        return
    workup = reaction.workups.add(type=WORKUP_TYPES[action])
    for child in root:
        tag = get_tag(child)
        if tag == 'cml:chemical':
            compound = workup.input.components.add()
            parse_reactant(child, compound)
            compound.reaction_role = reaction_pb2.ReactionRole.WORKUP
        elif tag == 'dl:phraseText':
            workup.details = child.text
        elif tag == 'dl:parameter':
            parse_parameter(child, workup)
        else:
            raise NotImplementedError(child)


def parse_parameter(root: ElementTree.Element, workup: reaction_pb2.ReactionWorkup):
    kind = root.attrib['propertyType']
    if kind == 'Time':
        workup.duration.CopyFrom(UNIT_RESOLVER.resolve(root.text))
    elif kind == 'Temperature':
        if root.text == 'room temperature':
            workup.temperature.control.type = (
                reaction_pb2.TemperatureConditions.TemperatureControl.AMBIENT)
        else:
            raise NotImplementedError(root.text)
    else:
        raise NotImplementedError(kind)


def main(argv):
    del argv  # Only used by app.run().
    reactions = []
    for filename in glob.glob(FLAGS.input_pattern):
        logging.info(filename)
        tree = ElementTree.parse(filename)
        root = tree.getroot()
        for reaction_cml in root.iterfind('cml:reaction',
                                          namespaces=NAMESPACES):
            reaction = parse_reaction(reaction_cml)
            reactions.append(reaction)
            assert False, reaction
    dataset = dataset_pb2.Dataset(reactions=reactions)
    if FLAGS.output:
        message_helpers.write_message(dataset, FLAGS.output)


if __name__ == '__main__':
    flags.mark_flag_as_required('input_pattern')
    app.run(main)
