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
"""Template filters to use with Jinja to format reaction messages into text."""

import collections

from ord_schema import units
from ord_schema import message_helpers
from ord_schema.visualization import drawing


def _is_true(boolean):
    return bool(boolean)


def _count_addition_order(inputs):
    """Returns the number of inputs for each addition_order value."""
    counts = collections.defaultdict(int)
    for value in inputs.values():
        counts[value.addition_order] += len(value.components)
    for order in sorted(counts):
        yield order, counts[order]


def _sort_addition_order(inputs):
    """Sort inputs by addition order, sorting again within stages/steps.

    Args:
        inputs: Map of ReactionInput messages.

    Returns:
        key: Text key into `inputs`.
        value: ReactionInput message.
    """
    orders = collections.defaultdict(list)
    for key in sorted(inputs):
        value = inputs[key]
        orders[value.addition_order].append((key, value))
    for order in sorted(orders):
        for i, (key, value) in enumerate(orders[order]):
            yield key, value


def _get_input_borders(components):
    for i, component in enumerate(components):
        if i == 0 and i + 1 == len(components):
            border = 'both'
        elif i == 0:
            border = 'left'
        elif i + 1 == len(components):
            border = 'right'
        else:
            border = 'clean'
        yield component, border


def _stirring_conditions(stirring):
    if stirring.method.type == stirring.method.NONE:
        return 'No stirring was used.'
    txt = ''
    txt += {
        stirring.rate.UNSPECIFIED: '',
        stirring.rate.HIGH: 'at a high rate',
        stirring.rate.MEDIUM: 'at a medium rate',
        stirring.rate.LOW: 'at a low rate',
    }[stirring.rate.type]
    if stirring.rate.rpm:
        txt += f' ({stirring.rate.rpm} rpm)'
    if stirring.method.type != stirring.method.UNSPECIFIED:
        txt += ' using '
        txt += {
            stirring.method.CUSTOM: 'a custom setup',
            stirring.method.STIR_BAR: 'a stir bar',
            stirring.method.OVERHEAD_MIXER: 'an overhead mixer',
            stirring.method.AGITATION: 'external agitation',
        }[stirring.method.type]
        txt += f' {_parenthetical_if_def(stirring.method.details)}'
    if txt.strip():
        txt = 'The reaction mixture was stirred ' + txt + '.'
    return txt


def _stirring_conditions_html(stirring):
    if stirring.method.type == stirring.method.NONE:
        return ''
    txt = ''
    if stirring.method.type != stirring.method.UNSPECIFIED:
        txt += {
            stirring.method.CUSTOM: stirring.method.details,
            stirring.method.STIR_BAR: 'stir bar',
            stirring.method.OVERHEAD_MIXER: 'overhead mixer',
            stirring.method.AGITATION: 'agitation',
        }[stirring.method.type]
    if stirring.rate.rpm:
        txt += f' ({stirring.rate.rpm} rpm)'
    return txt


def _pressure_conditions(pressure):
    txt = ''
    if pressure.atmosphere.type != pressure.atmosphere.UNSPECIFIED:
        txt += {
            pressure.atmosphere.CUSTOM: 'under a custom atmosphere',
            pressure.atmosphere.AIR: 'under air',
            pressure.atmosphere.NITROGEN: 'under nitrogen',
            pressure.atmosphere.ARGON: 'under argon',
            pressure.atmosphere.OXYGEN: 'under oxygen',
            pressure.atmosphere.HYDROGEN: 'under hydrogen',
        }[pressure.atmosphere.type]
        txt += f' {_parenthetical_if_def(pressure.atmosphere.details)}'
    if pressure.control.type != pressure.control.UNSPECIFIED:
        txt += ' '
        txt += {
            pressure.control.CUSTOM: 'using a custom pressure controller',
            pressure.control.AMBIENT: 'using ambient pressure',
            pressure.control.SEALED: 'after fully sealing the reaction vessel',
            pressure.control.PRESSURIZED: 'using pressurization',
        }[pressure.control.type]
        txt += f' {_parenthetical_if_def(pressure.control.details)}'
        setpoint = units.format_message(pressure.setpoint)
        if setpoint:
            txt += f' with a setpoint of {setpoint}'
    if txt:
        txt = 'The reaction was run ' + txt + '.'
    return txt


def _pressure_conditions_html(pressure):
    txt = ''
    if pressure.atmosphere.type != pressure.atmosphere.UNSPECIFIED:
        txt += {
            pressure.atmosphere.CUSTOM: pressure.atmosphere.details,
            pressure.atmosphere.AIR: 'in air',
            pressure.atmosphere.NITROGEN: 'under nitrogen',
            pressure.atmosphere.ARGON: 'under argon',
            pressure.atmosphere.OXYGEN: 'under oxygen',
            pressure.atmosphere.HYDROGEN: 'under hydrogen',
        }[pressure.atmosphere.type]
    if pressure.atmosphere.type != pressure.atmosphere.UNSPECIFIED:
        setpoint = units.format_message(pressure.setpoint)
        if setpoint:
            txt += f' ({setpoint})'
    return txt


def _temperature_conditions(temperature):
    txt = ''
    if temperature.control.type != temperature.control.UNSPECIFIED:
        txt += 'The reaction was run '
        txt += {
            temperature.control.CUSTOM:
                'under custom temperature conditions',
            temperature.control.AMBIENT:
                'under ambient temperature conditions',
            temperature.control.OIL_BATH:
                'in an oil bath',
            temperature.control.WATER_BATH:
                'in a water bath',
            temperature.control.SAND_BATH:
                'in a sand bath',
            temperature.control.ICE_BATH:
                'in an ice bath',
            temperature.control.DRY_ALUMINUM_PLATE:
                'using an aluminum heating block',
            temperature.control.MICROWAVE:
                'in a microwave reactor',
            temperature.control.DRY_ICE_BATH:
                'in a dry ice bath',
            temperature.control.AIR_FAN:
                'using a fan for temperautre control',
            temperature.control.LIQUID_NITROGEN:
                'using liquid nitrogen for temperature control',
        }[temperature.control.type]
        txt += f' {_parenthetical_if_def(temperature.control.details)}'
        setpoint = units.format_message(temperature.setpoint)
        if setpoint:
            txt += f' with a setpoint of {setpoint}'
    return txt + '.'


def _temperature_conditions_html(temperature):
    txt = ''
    if (temperature.control.type == temperature.control.UNSPECIFIED
            or temperature.control.type == temperature.control.AMBIENT):
        return 'ambient temperature'
    setpoint = units.format_message(temperature.setpoint)
    if setpoint:
        txt += f'{setpoint}'
    return txt


def _product_color_texture(product):
    txt = ''
    txt += f'{product.isolated_color} '
    txt += {
        product.Texture.UNSPECIFIED:
            '',
        product.Texture.CUSTOM:
            product.texture.details,
        product.Texture.POWDER:
            f'powder {_parenthetical_if_def(product.texture.details)}',
        product.Texture.CRYSTAL:
            f'set of crystals {_parenthetical_if_def(product.texture.details)}',
        product.Texture.OIL:
            f'oil {_parenthetical_if_def(product.texture.details)}',
    }[product.texture.type]
    if not txt.strip():
        return ''
    return f'It appeared as a {txt}.'


def _selectivity_type(selectivity):
    return {
        selectivity.CUSTOM: selectivity.details,
        selectivity.EE: 'e.e.',
        selectivity.ER: 'e.r.',
        selectivity.DE: 'd.e.',
    }[selectivity.type]


def _analysis_format(analysis):
    # TODO(ccoley) include data?
    return {
        analysis.UNSPECIFIED: '<UNK_ANALYSIS>',
        analysis.CUSTOM: 'a custom analysis',
        analysis.LC: 'liquid chromatography',
        analysis.GC: 'gas chromatography',
        analysis.IR: 'IR spectroscopy',
        analysis.NMR_1H: '1H NMR',
        analysis.NMR_13C: '13C NMR',
        analysis.NMR_OTHER: 'NMR (other)',
        analysis.MP: 'melting point characterization',
        analysis.UV: 'UV spectroscopy',
        analysis.TLC: 'thin-layer chromatography',
        analysis.MS: 'mass spectrometry',
        analysis.HRMS: 'high-resolution mass spectrometry',
        analysis.MSMS: 'two-dimensional MS',
        analysis.WEIGHT: 'isolated weight',
        analysis.LCMS: 'LCMS',
        analysis.GCMS: 'GCMS',
        analysis.ELSD: 'ELSD',
        analysis.CD: 'circular dichroism',
        analysis.SFC: 'supercritical fluid chromatography',
    }[analysis.type]


def _compound_svg(compound):
    mol = message_helpers.get_compound_mol(compound)
    if mol:
        return drawing.mol_to_svg(mol)
    return 'no RDKIT_BINARY'


def _compound_png(compound):
    mol = message_helpers.get_compound_mol(compound)
    if mol:
        return drawing.mol_to_png(mol)
    return 'no RDKIT_BINARY'


def _compound_amount(compound):
    amount = compound.WhichOneof('amount')
    if not amount:
        return ''
    return units.format_message(getattr(compound, amount))


def _compound_name(compound):
    for identifier in compound.identifiers:
        if identifier.type == identifier.NAME:
            return identifier.value
    return ''


def _compound_smiles(compound):
    for identifier in compound.identifiers:
        if identifier.type == identifier.SMILES:
            return identifier.value
    return ''


def _compound_role(compound):
    limiting_if_true = {
        True: 'limiting',
        False: '',
        None: '',
    }
    return {
        compound.ReactionRole.UNSPECIFIED:
            '',
        compound.ReactionRole.REACTANT:
            f'{limiting_if_true[compound.is_limiting]} reactant',
        compound.ReactionRole.REAGENT:
            'reagent',
        compound.ReactionRole.SOLVENT:
            'solvent',
        compound.ReactionRole.CATALYST:
            'catalyst',
        compound.ReactionRole.INTERNAL_STANDARD:
            'internal standard',
        compound.ReactionRole.WORKUP:
            '',
        compound.ReactionRole.PRODUCT:
            'product',
    }[compound.reaction_role]


def _compound_source_prep(compound):
    txt = []
    if compound.vendor_source:
        txt.append(f'purchased from {compound.vendor_source}')
    if compound.vendor_id:
        txt.append(f'catalog #{compound.vendor_id}')
    if compound.vendor_lot:
        txt.append(f'lot #{compound.vendor_lot}')
    for preparation in compound.preparations:
        txt.append({
            preparation.UNSPECIFIED: '',
            preparation.CUSTOM: '',
            preparation.NONE: '',
            preparation.REPURIFIED: 'repurified',
            preparation.SPARGED: 'sparged',
            preparation.DRIED: 'dried',
            preparation.SYNTHESIZED: 'synthesized in-house'
        }[preparation.type])
        txt.append(preparation.details)
    if any(elem for elem in txt):
        return '(' + ', '.join([elem for elem in txt if elem]) + ')'
    return ''


def _parenthetical_if_def(string):
    if not string:
        return ''
    return f'({string})'


def _vessel_prep(vessel):
    preparation_strings = []
    for preparation in vessel.preparations:
        preparation_strings.append({
            preparation.UNSPECIFIED: '',
            preparation.CUSTOM: 'prepared',
            preparation.NONE: '',
            preparation.OVEN_DRIED: 'oven-dried',
        }[preparation.type])
    return ', '.join(preparation_strings)


def _vessel_size(vessel):
    if vessel.volume.value:
        return f'{units.format_message(vessel.volume)}'
    return ''


def _vessel_material(vessel):
    return {
        vessel.material.UNSPECIFIED: '',
        vessel.material.CUSTOM: f'{vessel.material.details}',
        vessel.material.GLASS: 'glass',
        vessel.material.POLYPROPYLENE: 'polypropylene',
        vessel.material.PLASTIC: 'plastic',
    }[vessel.material.type]


def _vessel_type(vessel):
    return {
        vessel.type.UNSPECIFIED:
            'vessel',
        vessel.type.CUSTOM:
            'vessel',
        vessel.type.ROUND_BOTTOM_FLASK:
            'round bottom flask',
        vessel.type.VIAL:
            'vial',
        vessel.type.WELL_PLATE:
            'well-plate',
        vessel.type.MICROWAVE_VIAL:
            'microwave vial',
        vessel.type.TUBE:
            'tube',
        vessel.type.CONTINUOUS_STIRRED_TANK_REACTOR:
            'continuous stirred-tank reactor',
        vessel.type.PACKED_BED_REACTOR:
            'packed bed reactor',
    }[vessel.type.type]


def _input_addition(reaction_input):
    txt = []
    if reaction_input.addition_time.value:
        txt.append(
            f'after {units.format_message(reaction_input.addition_time)}')
    txt.append({
        reaction_input.addition_speed.UNSPECIFIED: '',
        reaction_input.addition_speed.ALL_AT_ONCE: 'all at once',
        reaction_input.addition_speed.FAST: 'quickly',
        reaction_input.addition_speed.SLOW: 'slowly',
        reaction_input.addition_speed.DROPWISE: 'dropwise',
    }[reaction_input.addition_speed.type])
    if reaction_input.addition_duration.value:
        txt.append(
            f'over {units.format_message(reaction_input.addition_duration)}')
    if any(elem for elem in txt):
        return '(' + ', '.join([elem for elem in txt if elem]) + ')'
    return ''


def _uses_addition_order(reaction):
    return any(value.addition_order for value in reaction.inputs.values())


def _round(value, places=2):
    fstring = '{:.%gg}' % places
    return fstring.format(value)


def _datetimeformat(value, format_string='%H:%M / %d-%m-%Y'):
    return value.strftime(format_string)


def _width(values):
    return 100 / len(values)


TEMPLATE_FILTERS = {
    'round': _round,
    'is_true': _is_true,
    'datetimeformat': _datetimeformat,
    'uses_addition_order': _uses_addition_order,
    'input_addition': _input_addition,
    'compound_svg': _compound_svg,
    'compound_png': _compound_png,
    'compound_amount': _compound_amount,
    'compound_name': _compound_name,
    'compound_smiles': _compound_smiles,
    'compound_role': _compound_role,
    'compound_source_prep': _compound_source_prep,
    'vessel_prep': _vessel_prep,
    'vessel_type': _vessel_type,
    'vessel_material': _vessel_material,
    'vessel_size': _vessel_size,
    'unit_format': units.format_message,
    'parenthetical_if_def': _parenthetical_if_def,
    'analysis_format': _analysis_format,
    'selectivity_type': _selectivity_type,
    'product_color_texture': _product_color_texture,
    'temperature_conditions': _temperature_conditions,
    'temperature_conditions_html': _temperature_conditions_html,
    'pressure_conditions': _pressure_conditions,
    'pressure_conditions_html': _pressure_conditions_html,
    'stirring_conditions': _stirring_conditions,
    'stirring_conditions_html': _stirring_conditions_html,
    'count_addition_order': _count_addition_order,
    'sort_addition_order': _sort_addition_order,
    'get_input_borders': _get_input_borders,
    'width': _width,
}
