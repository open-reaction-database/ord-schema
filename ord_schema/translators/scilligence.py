# Copyright 2022 Open Reaction Database Project Authors
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
"""Helper functions for constructing Protocol Buffer messages."""

import os
import warnings
import datetime
from typing import List
from collections import defaultdict
import xmltodict

from ord_schema.proto import reaction_pb2
from ord_schema.units import UnitResolver

_resolver = UnitResolver()

__SUPPORTED_VERSIONS = ['3.0']
__ANALYSIS_DICTIONARY = {'HPLC': 'LC', 'MS': 'MS', 'LCMS': 'LCMS',
            'IR': 'IR', 'C-NMR': 'NMR_13C', 'H-NMR': 'NMR_1H'}
__TEXTURE_DICTIONARY = {'solid': 'SOLID', 'liquid': 'LIQUID', 'oil': 'OIL'}

def convert_reaction(xml_string: str) -> reaction_pb2.Reaction:
    # pylint: disable=too-many-locals,too-many-branches,too-many-statements

    """Builds a Reaction message from an xml string originating from a
    Scilligence notebook export.

    Args:
        xml_string: XML-formatted string from a Scilligence ELN export.

    Returns:
        Reaction message.
    """

    reaction = reaction_pb2.Reaction()
    try:
        dic = xmltodict.parse(xml_string)
    except xmltodict.expat.ExpatError as exception:
        raise ValueError('Invalid XML string') from exception

    # Meta-data about the conversion
    try:
        assert dic['eln']['@createdby'] == 'eln'
        eln_version = dic['eln']['@version']
        assert eln_version in __SUPPORTED_VERSIONS # supported versions
    except Exception as exception:
        raise ValueError('Parsed XML is not a supported Scilligence ELN export'
            f'; supported versions: {__SUPPORTED_VERSIONS}') from exception

    # Reformat parsed xml into a more convenient dictionary
    reaction_dict = defaultdict(lambda: None)
    for field in dic['eln']['field']:
        field_id = field['@id']
        if 'jssdf' in field:
            if 'r' not in field['jssdf']:
                warnings.warn(f'Field {field_id} is empty, skipping...')
                continue
            field_data = field['jssdf']['r']

            if isinstance(field_data, dict):
                reaction_dict[field_id] = xml_datafield_to_dict(field_data['i'])
            elif isinstance(field_data, list):
                reaction_dict[field_id] = [
                    xml_datafield_to_dict(datum['i']) for datum in field_data
                ]
        elif field_id == 'procedure':
            reaction_dict[field_id] = field.get('#text', None)

    # Expected fields in reaction_dict are as follows
    # "request" : dict
        # TODO: figure out what this means
    # "summary" : dict
        # Expecting temperature, time, pressure, yield, succeed,
        # rating
    # "summary2" : dict
        # Expecting customerid, compoundnovelty, color, state,
        # stability, status (of the reported product)
    # "rxn" : dict
        # cannot capture JSDraw information currently!!
    # "reactants" : list
        # Expecting name, compid, structure (JSDraw), saltratio,
        # salt, mf, mw, eq, casno, mass, moles, hazardflags,
        # safety, stockcon, solvent, volume, density, molar,
        # limiting, notes
        # Note: certain fields, like 'moles', may also have an
        # "@num" key using a standard unit like moles. However,
        # the "#text" key's value seems easy to parse.
        # TODO -- learn to parse JSDraw structures
    # "reagents": list
        # Expecting name,
        # mf, mw, eq, casno, mass, moles, source, hazardflags,
        # safety, catalyst, stockcon, solvent, volume, density,
        # purity, molar, notes
    # "solvents": list
        # Expecting name, volume, ratio, source, purity,
        # hazardflags, safety, notes
    # "products": dict
        # Expecting name, compid, regno, structure (JSDraw),
        # compoundtype, saltratio, salt, eq, mf, mw, casno,
        # mass, theorymass, purity, yield, chiralpurity,
        # moles, state, volume, density, notes
        # Note: "state" may contain notes about yield quantification
        # and can be one of Experiment Discarded,Purified Product,
        # Crude Product,Product Not Isolated,Yield From HPLC,
        # Yield from Conversion Rate
    # "containers" : dict
        # TODO
        # Expecting select, compid, amount, barcode
    # "analyticalsamples" : dict
        # TODO
        # Expecting select, sampleid, amount, concentration,
        # analyticaltype, analyticalmethod, solvent, mw, requestid,
        # status, result_mw, result_purity, result_summary,
        # result_report, comment
    # "analyticalresults" : dict (?)
        # Expecting sampleid, filetype, mw, purity, summary, file
    # "procedure" : string
    # "operations" : dict (?)
        # Expecting time, operation, observation, reffiles
        # Note: this appears to be interventions/observations during
        # the course of the reaction
    # "reference" : dict
        # Expecting reffiles, authors, title, journal, year, volume, pages

    # Populate the reaction message
    reaction.notes.procedure_details = reaction_dict['procedure']

    # Populate all input species
    def add_input(name, compound_dict, role='UNSPECIFIED'):
        # TODO: figure out how/if stock solutions are defined in Scilligence;
        # this assumes that all inputs are single component, which I think
        # is all we can discern from the ELN export
        _input = reaction.inputs[name]
        compound = _input.components.add(reaction_role=role)
        if compound_dict.get('limiting', None) == '1':
            compound.is_limiting = True
        if compound_dict.get('source', None):
            compound.source.vendor = compound_dict['source']
        # TODO: consider adding purity field to compounds in the schema

        # Identifiers
        if compound_dict.get('name', None):
            compound.identifiers.add(type='NAME',
                value=compound_dict['name'])
        if compound_dict.get('casno', None):
            compound.identifiers.add(type='CAS_NUMBER',
                value=compound_dict['casno'])
        if compound_dict.get('structure', None):
            # TODO: parse JSDraw
            pass
        if compound_dict.get('mw', None):
            compound.features['mw'].float_value = float(compound_dict['mw'])
        if compound_dict.get('salt', None):
            compound.preparations.add(type='CUSTOM',
                details=(f'Used as the {compound_dict["salt"]} salt with '
                    f'a salt ratio of {compound_dict["saltratio"]}'))

        # For amounts, prefer mass > moles > volume
        if compound_dict.get('mass', None):
            compound.amount.mass.CopyFrom(_resolver.resolve(
                compound_dict['mass']))
        elif compound_dict.get('moles', None):
            compound.amount.moles.CopyFrom(_resolver.resolve(
                compound_dict['moles']))
        elif compound_dict.get('volume', None):
            compound.amount.volume.CopyFrom(_resolver.resolve(
                compound_dict['volume']))

    for i,reactant_dict in enumerate(reaction_dict['reactants']):
        add_input(f'reactant{i+1}', reactant_dict, 'REACTANT')
    for i,reagent_dict in enumerate(reaction_dict['reagents']):
        add_input(f'reagent{i+1}', reagent_dict, 'REAGENT')
    for i,solvent_dict in enumerate(reaction_dict['solvents']):
        add_input(f'solvent{i+1}', solvent_dict, 'SOLVENT')

    # Conditions
    if reaction_dict['summary']['temperature']:
        reaction.conditions.temperature.setpoint.CopyFrom(
            _resolver.resolve(reaction_dict['summary']['temperature']))

    if reaction_dict['summary']['pressure']:
        reaction.conditions.pressure.setpoint.CopyFrom(
            _resolver.resolve(reaction_dict['summary']['pressure']))

    # Outcome (only single outcome possible)
    outcome = reaction.outcomes.add()
    if reaction_dict['summary']['time']:
        outcome.reaction_time.CopyFrom(
            _resolver.resolve(reaction_dict['summary']['time']))

    # Products (only single product possible (?) TODO: check
    product_dict = reaction_dict['products']
    product = outcome.products.add()
    if product_dict['name']:
        product.identifiers.add(type='NAME',
            value=product_dict['name'])
    if product_dict['casno']:
        product.identifiers.add(type='CAS_NUMBER',
            value=product_dict['casno'])
    if product_dict['structure']:
        # TODO: parse JSDraw
        pass

    if product_dict['yield']:
        assert product_dict['yield'][-1] == '%'
        yield_pct = float(product_dict['yield'][:-1])
        measurement = product.measurements.add(type='YIELD',
            percentage=dict(value=yield_pct))
    if product_dict['purity']:
        assert product_dict['purity'][-1] == '%'
        purity_pct = float(product_dict['purity'][:-1])
        measurement = product.measurements.add(type='PURITY',
            percentage=dict(value=purity_pct))
    if product_dict['moles']:
        measurement = product.measurements.add(type='AMOUNT')
        measurement.amount.moles.CopyFrom(_resolver.resolve(
            product_dict['moles']))
    # TODO: use state, e.g., "Yield from HPLC"
    # TODO: use chiralpurity when we have an example

    if reaction_dict['summary2']['color']:
        product.isolated_color = reaction_dict['summary2']['color']
    if reaction_dict['summary2']['state']:
        texture_string = reaction_dict['summary2']['state']
        try:
            product.texture.CopyFrom(reaction_pb2.ProductCompound.Texture(
                type=__TEXTURE_DICTIONARY[texture_string.lower()]))
        except KeyError:
            product.texture.CopyFrom(reaction_pb2.ProductCompound.Texture(
                type='CUSTOM', details=texture_string))

    # TODO: figure out if analyticalsamples is sometimes a list, or if there is
    # always only a single entry
    for i,analysis_dict in enumerate([reaction_dict['analyticalsamples']]):
        analysis = outcome.analyses[f'analysis{i}']
        try:
            analysis.CopyFrom(reaction_pb2.Analysis(
                type=__ANALYSIS_DICTIONARY[analysis_dict['analyticalmethod']]))
        except KeyError:
            analysis.CopyFrom(reaction_pb2.Analysis(type='CUSTOM',
                details=analysis_dict['analyticalmethod']))

    # TODO: figure out if operatoins is sometimes a list, or if there is
    # always only a single entry
    for operation_dict in [reaction_dict['operations']]:
        if operation_dict['operation']: # there was an intervention
            reaction.conditions.conditions_are_dynamic = True
            reaction.conditions.details += (f'@ {operation_dict["time"]}, '
                f' {operation_dict["operation"]}. ')

    # Safety notes are copied from starting materials only
    safety_flags = set()
    for compound_dict in (reaction_dict['reactants'] +
            reaction_dict['reagents'] + reaction_dict['solvents']):
        if compound_dict.get('safety', None):
            safety_flags |= set(compound_dict['safety'].split(','))
        if compound_dict.get('hazardflags', None):
            safety_flags |= set(compound_dict['hazardflags'].split(','))
    if safety_flags:
        reaction.notes.safety_notes = ('Starting material handling notes: ' +
            '; '.join(list(safety_flags)))

    # TODO: figure out how to map provenance onto our schema
    event = reaction_pb2.RecordEvent(
        time=dict(value=str(datetime.datetime.now())),
        details=('Automatically translated from Scilligence ELN ' +
            f'version {eln_version}'),
    )
    reaction.provenance.record_created.CopyFrom(event)

    return reaction


def xml_datafield_to_dict(subfield_list: List[dict]) -> dict:
    """Reformats a list of data fields into a dictionary.

    Args:
        subfield_list: List of dictionaries with keys "@n" and "#text"

    Returns:
        Dictionary equivalent of the datafield list.
    """
    dic = {}
    for subfield in subfield_list:
        dic[subfield['@n']] = subfield.get('#text', None)
    return dic


if __name__ == '__main__':
    with open(os.path.join(os.path.dirname(__file__),
            'scilligence_example.xml')) as fid:
        print(convert_reaction(fid.read()))
