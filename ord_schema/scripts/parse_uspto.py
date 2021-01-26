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
    'cml': 'http://www.xml-cml.org/schema',
    'cmlDict': 'http://www.xml-cml.org/dictionary/cml/',
    'dl': 'http://bitbucket.org/dan2097',
    'nameDict': 'http://www.xml-cml.org/dictionary/cml/name/',
    'unit': 'http://www.xml-cml.org/unit/',
}

UNIT_RESOLVER = units.UnitResolver()

REACTION_ROLES = {
    'catalyst': reaction_pb2.ReactionRole.CATALYST,
    'product': reaction_pb2.ReactionRole.PRODUCT,
    'reactant': reaction_pb2.ReactionRole.REACTANT,
    'solvent': reaction_pb2.ReactionRole.SOLVENT,
}

PRODUCT_STATES = {
    'crystal': reaction_pb2.ProductCompound.Texture.CRYSTAL,
    'crystalline-solid': reaction_pb2.ProductCompound.Texture.CRYSTAL,
    'crystals': reaction_pb2.ProductCompound.Texture.CRYSTAL,
    'foam': reaction_pb2.ProductCompound.Texture.FOAM,
    'foam-like': reaction_pb2.ProductCompound.Texture.FOAM,
    'foam/very': reaction_pb2.ProductCompound.Texture.FOAM,
    'needle': reaction_pb2.ProductCompound.Texture.CRYSTAL,
    'needle-like': reaction_pb2.ProductCompound.Texture.CRYSTAL,
    'needles': reaction_pb2.ProductCompound.Texture.CRYSTAL,
    'oil': reaction_pb2.ProductCompound.Texture.OIL,
    'oil-like': reaction_pb2.ProductCompound.Texture.OIL,
    'oils': reaction_pb2.ProductCompound.Texture.OIL,
    'powder': reaction_pb2.ProductCompound.Texture.POWDER,
    'powders': reaction_pb2.ProductCompound.Texture.POWDER,
    'prisms': reaction_pb2.ProductCompound.Texture.CRYSTAL,
    'semisolid': reaction_pb2.ProductCompound.Texture.SEMI_SOLID,
    'semi-solid': reaction_pb2.ProductCompound.Texture.SEMI_SOLID,
}

WORKUP_TYPES = {
    'Add': reaction_pb2.ReactionWorkup.ADDITION,
    'Cool': reaction_pb2.ReactionWorkup.TEMPERATURE,
    'Concentrate': reaction_pb2.ReactionWorkup.CONCENTRATION,
    'Dissolve': reaction_pb2.ReactionWorkup.DISSOLUTION,
    'Distill': reaction_pb2.ReactionWorkup.DISTILLATION,
    'Dry with material': reaction_pb2.ReactionWorkup.DRY_WITH_MATERIAL,
    'Extract': reaction_pb2.ReactionWorkup.EXTRACTION,
    'Filter': reaction_pb2.ReactionWorkup.FILTRATION,
    'Heat': reaction_pb2.ReactionWorkup.TEMPERATURE,
    'Sample': reaction_pb2.ReactionWorkup.ALIQUOT,
    'Stir': reaction_pb2.ReactionWorkup.STIRRING,
    'Wait': reaction_pb2.ReactionWorkup.WAIT,
    'Wash': reaction_pb2.ReactionWorkup.WASH,
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


def resolve_units(value):
    value = value.replace('˜', ' ')
    value = value.replace('yield ', ' ')
    return UNIT_RESOLVER.resolve(value, allow_range=True)


def get_component_map(root: ElementTree.ElementTree) -> Dict[str, str]:
    """Builds a map of components to inputs."""
    reaction_inputs = {}
    # Start with a separate input for each component.
    for component in root.find('cml:reactantList', namespaces=NAMESPACES):
        molecule_id = get_molecule_id(component)
        reaction_inputs[molecule_id] = molecule_id
    for component in root.find('cml:spectatorList', namespaces=NAMESPACES):
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
                parse_conditions(action, reaction)
                parse_workup(action, reaction)
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
        elif tag in ['dl:headingText', 'dl:paragraphNum']:
            continue  # Ignored.
        else:
            raise NotImplementedError(child)


def parse_product(root, product_compound):
    role = root.attrib.get('role')
    if role:
        product_compound.reaction_role = REACTION_ROLES[role]
    for child in root:
        tag = get_tag(child)
        if tag == 'dl:entityType':
            pass
        elif tag == 'cml:molecule':
            parse_molecule(child, product_compound)
        elif tag == 'cml:amount':
            parse_product_amount(child, product_compound)
        elif tag == 'cml:identifier':
            parse_identifier(child, product_compound)
        elif tag == 'dl:state':
            texture = PRODUCT_STATES.get(
                child.text.lower(), reaction_pb2.ProductCompound.Texture.CUSTOM)
            product_compound.texture.type = texture
            product_compound.texture.details = child.text
        elif tag == 'dl:appearance':
            product_compound.isolated_color = child.text
        else:
            raise NotImplementedError(child)


def parse_molecule(root, compound):
    for child in root:
        tag = get_tag(child)
        if tag in ['cml:name', 'dl:nameResolved']:
            compound.identifiers.add(type='NAME', value=child.text)
        else:
            raise NotImplementedError(child)


def parse_product_amount(root, product_compound):
    property_type = root.attrib[f'{{{NAMESPACES["dl"]}}}propertyType']
    if 'PERCENTYIELD' in property_type:
        match = re.search(r'(\d+\.?\d*)', root.text)
        if not match:
            logging.warning(f'AMOUNT: {root.text}')
            return
        measurement = product_compound.measurements.add()
        measurement.type = measurement.YIELD
        measurement.percentage.value = float(match.group(1))
    else:
        measurement = product_compound.measurements.add()
        measurement.type = measurement.AMOUNT
        parse_amount(root, measurement)


def parse_identifier(root, compound):
    kind = root.attrib['dictRef']
    value = root.attrib['value']
    if kind == 'cml:smiles':
        compound.identifiers.add(type='SMILES', value=value)
    elif kind == 'cml:inchi':
        compound.identifiers.add(type='INCHI', value=value)
    else:
        raise NotImplementedError(kind)


def parse_reactant(root, compound):
    role = root.attrib.get('role')
    if role:
        compound.reaction_role = REACTION_ROLES[role]
    for child in root:
        tag = get_tag(child)
        if tag == 'dl:entityType':
            if child.text not in [
                    'exact', 'chemicalClass', 'definiteReference'
            ]:
                raise NotImplementedError(child.text)
            continue
        elif tag == 'cml:molecule':
            parse_molecule(child, compound)
        elif tag == 'cml:amount':
            parse_amount(child, compound)
        elif tag == 'cml:identifier':
            parse_identifier(child, compound)
        elif tag in ['dl:state', 'dl:appearance']:
            continue  # Not supported by Compound.
        else:
            raise NotImplementedError(child)


def parse_amount(root, compound):
    property_type = root.attrib[f'{{{NAMESPACES["dl"]}}}propertyType']
    if property_type in ['MOLARITY', 'PH']:
        return
    if compound.amount.WhichOneof('kind') not in ['mass', 'volume']:
        # Don't overwrite existing Mass or Volume with Moles.
        try:
            amount = resolve_units(root.text)
        except (KeyError, ValueError) as error:
            logging.warning(f'UNITS: {error} ("{root.text}")')
            return
        if isinstance(amount, reaction_pb2.Mass):
            compound.amount.mass.CopyFrom(amount)
        elif isinstance(amount, reaction_pb2.Moles):
            compound.amount.moles.CopyFrom(amount)
        elif isinstance(amount, reaction_pb2.Volume):
            compound.amount.volume.CopyFrom(amount)
        else:
            raise NotImplementedError(amount)


def parse_conditions(root, reaction):
    del reaction  # Unused.
    if not root.findall('cml:chemical', namespaces=NAMESPACES):
        return  # Refers to an added component; is a workup.
    action = root.attrib['action']
    return action  # TODO(kearnes): Implement this?


def parse_workup(root, reaction):
    if root.findall('dl:chemical', namespaces=NAMESPACES):
        return  # Refers to an input component; not a workup.
    action = root.attrib['action']
    components = []
    for component in root.findall('cml:chemical', namespaces=NAMESPACES):
        entity_type = component.find('dl:entityType', namespaces=NAMESPACES)
        if entity_type.text in ['exact', 'chemicalClass', 'definiteReference']:
            components.append(component)
    if action == 'Dry' and components:
        action = 'Dry with material'
    details = root.find('dl:phraseText', namespaces=NAMESPACES)
    if action in ['Purify', 'Recover'] and details:
        # Make some actions more specific.
        # TODO(kearnes): This could be expanded.
        if 'distill' in details.lower():
            action = 'Distill'
        elif 'filtration' in details.lower():
            action = 'Filter'
        else:
            logging.info(f'{action}: {details}')
    workup_type = WORKUP_TYPES.get(action, reaction_pb2.ReactionWorkup.CUSTOM)
    workup = reaction.workups.add(type=workup_type)
    for component in components:
        compound = workup.input.components.add()
        parse_reactant(component, compound)
        compound.reaction_role = reaction_pb2.ReactionRole.WORKUP
    for child in root:
        tag = get_tag(child)
        if tag == 'cml:chemical':
            continue  # Components are handled all together above.
        elif tag == 'dl:phraseText':
            workup.details = child.text
        elif tag == 'dl:parameter':
            try:
                parse_parameter(child, workup)
            except ValueError as error:
                logging.warning(
                    f'PARAMETER: {ElementTree.tostring(child)}: {error}')
        elif tag in ['dl:atmosphere']:
            continue  # Not supported by ReactionWorkup.
        else:
            raise NotImplementedError(child)


def parse_parameter(root: ElementTree.Element,
                    workup: reaction_pb2.ReactionWorkup):
    kind = root.attrib['propertyType']
    if kind == 'Time':
        if root.text in ['overnight']:
            workup.duration.value = 8
            workup.duration.precision = 8
            workup.duration.units = reaction_pb2.Time.HOUR
        else:
            try:
                workup.duration.CopyFrom(resolve_units(root.text))
            except (KeyError, ValueError) as error:
                logging.warning(f'UNITS: {error} ("{root.text}")')
    elif kind == 'Temperature':
        if root.text.lower() in [
                'room temperature',
                'room temp',
                'rt',
                'r.t',
                'r.t.',
                'ambient temperature',
        ]:
            workup.temperature.control.type = (
                reaction_pb2.TemperatureConditions.TemperatureControl.AMBIENT)
        else:
            value = root.text.rstrip('.').replace('° ', '°')
            try:
                workup.temperature.setpoint.CopyFrom(resolve_units(value))
            except (KeyError, ValueError) as error:
                if str(error) not in ['unrecognized units: eq']:
                    logging.warning(f'UNITS: {error} ("{root.text}")')
    elif kind in ['Frequency', 'Length', 'Pressure']:
        return  # Not supported in ReactionWorkup.
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
            try:
                reaction = parse_reaction(reaction_cml)
            except (KeyError, NotImplementedError) as error:
                raise ValueError(ElementTree.dump(reaction_cml)) from error
            reactions.append(reaction)
    dataset = dataset_pb2.Dataset(reactions=reactions)
    if FLAGS.output:
        message_helpers.write_message(dataset, FLAGS.output)


if __name__ == '__main__':
    flags.mark_flag_as_required('input_pattern')
    app.run(main)
