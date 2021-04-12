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
"""Template filters to use with Jinja to format reaction messages into text.

NOTE(kearnes): Flask escapes HTML by default, so the HTML-producing filters
in this module do not include any HTML tags, only their contents.
"""

import collections
from typing import Any, Iterable, List, Mapping, Optional, Tuple

from dateutil import parser
from google.protobuf import text_format  # pytype: disable=import-error

from ord_schema import units
from ord_schema import message_helpers
from ord_schema.visualization import drawing
from ord_schema.proto import reaction_pb2


def _is_true(boolean: Any) -> bool:
    """Returns whether a value is True."""
    return bool(boolean)


def _count_addition_order(
    inputs: Mapping[str,
                    reaction_pb2.ReactionInput]) -> Iterable[Tuple[int, int]]:
    """Returns the number of inputs for each addition_order value.

    Args:
        inputs: Map of ReactionInput messages.

    Yields:
        order: Integer addition order.
        count: Integer number of Compounds with that addition order.
    """
    counts = collections.defaultdict(int)
    for value in inputs.values():
        counts[value.addition_order] += len(value.components)
    for order in sorted(counts):
        yield order, counts[order]


def _sort_addition_order(
    inputs: Mapping[str, reaction_pb2.ReactionInput]
) -> Iterable[Tuple[str, reaction_pb2.ReactionInput]]:
    """Sorts inputs by addition order, sorting again within stages/steps.

    Args:
        inputs: Map of ReactionInput messages.

    Yields:
        key: Text key into `inputs`.
        value: ReactionInput message.
    """
    orders = collections.defaultdict(list)
    for key in sorted(inputs):
        value = inputs[key]
        orders[value.addition_order].append((key, value))
    for order in sorted(orders):
        for key, value in orders[order]:
            yield key, value


def _get_input_borders(
    components: List[reaction_pb2.Compound]
) -> Iterable[Tuple[reaction_pb2.Compound, str]]:
    """Returns the CSS class for a Compound cell.

    The HTML representation of a Reaction groups Compounds by their parent
    ReactionInput. This function assigns a CSS class to each component for use
    in defining borders for table cells.

    Args:
        components: List of Compound messages.

    Yields:
        component: Compound message.
        border: String CSS class describing the border for this cell.
    """
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


def _stirring_conditions(stirring: reaction_pb2.StirringConditions) -> str:
    """Generates a text description of stirring conditions.

    Args:
        stirring: StirringConditions message.

    Returns:
        String description of the stirring conditions.
    """
    if stirring.type == stirring.NONE:
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
    if stirring.type != stirring.UNSPECIFIED:
        txt += ' using '
        txt += {
            stirring.CUSTOM: 'a custom setup',
            stirring.STIR_BAR: 'a stir bar',
            stirring.OVERHEAD_MIXER: 'an overhead mixer',
            stirring.AGITATION: 'external agitation',
        }[stirring.type]
        txt += f' {_parenthetical_if_def(stirring.details)}'
    if txt.strip():
        txt = 'The reaction mixture was stirred ' + txt + '.'
    return txt


def _stirring_conditions_html(stirring: reaction_pb2.StirringConditions) -> str:
    """Generates an HTML-ready description of stirring conditions.

    Args:
        stirring: StirringConditions message.

    Returns:
        String description of the stirring conditions.
    """
    if stirring.type == stirring.NONE:
        return ''
    txt = ''
    if stirring.type != stirring.UNSPECIFIED:
        txt += {
            stirring.CUSTOM: stirring.details,
            stirring.STIR_BAR: 'stir bar',
            stirring.OVERHEAD_MIXER: 'overhead mixer',
            stirring.AGITATION: 'agitation',
        }[stirring.type]
    if stirring.rate.rpm:
        txt += f' ({stirring.rate.rpm} rpm)'
    return txt


def _pressure_conditions(pressure: reaction_pb2.PressureConditions) -> str:
    """Generates a text description of pressure conditions.

    Args:
        pressure: PressureConditions message.

    Returns:
        String description of the pressure conditions.
    """
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


def _pressure_conditions_html(pressure: reaction_pb2.PressureConditions) -> str:
    """Generates an HTML-ready description of pressure conditions.

    Args:
        pressure: PressureConditions message.

    Returns:
        String description of the pressure conditions.
    """
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


def _temperature_conditions(
        temperature: reaction_pb2.TemperatureConditions) -> str:
    """Generates a text description of temperature conditions.

    Args:
        temperature: TemperatureConditions message.

    Returns:
        String description of the temperature conditions.
    """
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


def _temperature_conditions_html(
        temperature: reaction_pb2.TemperatureConditions) -> str:
    """Generates an HTML-ready description of temperature conditions.

    Args:
        temperature: TemperatureConditions message.

    Returns:
        String description of the temperature conditions.
    """
    txt = ''
    if temperature.control.type == temperature.control.UNSPECIFIED:
        return ''
    if temperature.control.type == temperature.control.AMBIENT:
        return 'ambient temperature'
    setpoint = units.format_message(temperature.setpoint)
    if setpoint:
        txt += f'{setpoint}'
    return txt


def _product_color_texture(product: reaction_pb2.ProductCompound) -> str:
    """Generates a text description of the color and texture of a product.

    Args:
        product: ProductCompound message.

    Returns:
        String description of the product.
    """
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


def _selectivity_type(
    selectivity: reaction_pb2.ProductMeasurement.Selectivity.SelectivityType
) -> str:
    """Returns a string version of the selectivity type."""
    return {
        selectivity.CUSTOM: selectivity.details,
        selectivity.EE: 'e.e.',
        selectivity.ER: 'e.r.',
        selectivity.DE: 'd.e.',
    }[selectivity.type]


def _analysis_format(analysis: reaction_pb2.Analysis.AnalysisType) -> str:
    """Returns a string version of the analysis type."""
    # TODO(ccoley) include data?
    return {
        analysis.UNSPECIFIED: 'an UNSPECIFIED analysis',
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


def _compound_svg(compound: reaction_pb2.Compound,
                  bond_length: int = 25) -> str:
    """Returns an SVG string for the given compound.

    If the compound does not have a structural identifier, a sentinel value
    is returned instead.

    Args:
        compound: Compound message.
        bond_length: Bond length in pixels.

    Returns:
        String SVG or sentinel value.
    """
    try:
        mol = message_helpers.mol_from_compound(compound)
        if mol:
            svg = drawing.mol_to_svg(mol, bond_length=bond_length)
            if svg is None:
                return (message_helpers.get_compound_smiles(compound) or
                        '[Compound]')
            return svg
    except ValueError:
        pass
    return message_helpers.get_compound_smiles(compound) or '[Compound]'


def _compound_png(compound: reaction_pb2.Compound) -> str:
    """Returns a PNG string for the given compound.

    If the compound does not have a structural identifier, a sentinel value
    is returned instead.

    Args:
        compound: Compound message.

    Returns:
        String PNG or sentinel value.
    """
    try:
        mol = message_helpers.mol_from_compound(compound)
        if mol:
            return drawing.mol_to_png(mol)
    except ValueError:
        pass
    return message_helpers.get_compound_smiles(compound) or '[Compound]'


def _amount(amount: reaction_pb2.Amount) -> Optional[str]:
    """Returns a string representing an Amount."""
    kind = amount.WhichOneof('kind')
    if not kind:
        return ''
    return units.format_message(getattr(amount, kind))


def _compound_name(compound: reaction_pb2.Compound) -> str:
    """Returns the compound name, if defined."""
    for identifier in compound.identifiers:
        if identifier.type == identifier.NAME:
            return identifier.value
    return ''


def _compound_smiles(compound: reaction_pb2.Compound) -> str:
    """Returns the compound SMILES, if defined."""
    for identifier in compound.identifiers:
        if identifier.type == identifier.SMILES:
            return identifier.value
    return ''


def _compound_role(compound: reaction_pb2.Compound, text: bool = False) -> str:
    """Returns a description of the compound role.

    Args:
        compound: Compound message.
        text: If True, return a text description. If False, return an
            HTML-ready description.

    Returns:
        String compound role description.
    """
    limiting_if_true = {
        True: 'limiting',
        False: '',
        None: '',
    }
    try:
        limiting = limiting_if_true[compound.is_limiting]
    except AttributeError:
        limiting = ''
    if text:
        options = {
            reaction_pb2.ReactionRole.UNSPECIFIED:
                '',
            reaction_pb2.ReactionRole.REACTANT:
                f'as a {limiting} reactant',
            reaction_pb2.ReactionRole.REAGENT:
                'as a reagent',
            reaction_pb2.ReactionRole.SOLVENT:
                'as a solvent',
            reaction_pb2.ReactionRole.CATALYST:
                'as a catalyst',
            reaction_pb2.ReactionRole.INTERNAL_STANDARD:
                'as an internal standard',
            reaction_pb2.ReactionRole.WORKUP:
                '',
            reaction_pb2.ReactionRole.AUTHENTIC_STANDARD:
                'as an authentic standard',
            reaction_pb2.ReactionRole.PRODUCT:
                'as a product',
        }
    else:
        options = {
            reaction_pb2.ReactionRole.UNSPECIFIED: '',
            reaction_pb2.ReactionRole.REACTANT: f'{limiting} reactant',
            reaction_pb2.ReactionRole.REAGENT: 'reagent',
            reaction_pb2.ReactionRole.SOLVENT: 'solvent',
            reaction_pb2.ReactionRole.CATALYST: 'catalyst',
            reaction_pb2.ReactionRole.INTERNAL_STANDARD: 'internal standard',
            reaction_pb2.ReactionRole.WORKUP: '',
            reaction_pb2.ReactionRole.AUTHENTIC_STANDARD: 'authentic standard',
            reaction_pb2.ReactionRole.PRODUCT: 'product',
        }
    return options[compound.reaction_role]


def _compound_source_prep(compound: reaction_pb2.Compound) -> str:
    """Returns a string describing the compound source and preparation.

    Args:
        compound: Compound message.

    Returns:
        String description of the source and preparation.
    """
    txt = []
    if compound.source.vendor:
        txt.append(f'purchased from {compound.source.vendor}')
    if compound.source.id:
        txt.append(f'catalog #{compound.source.id}')
    if compound.source.lot:
        txt.append(f'lot #{compound.source.lot}')
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


def _product_yield(product: reaction_pb2.ProductCompound) -> str:
    """Returns a string describing how a product yield was calculated."""
    for measurement in product.measurements:
        if measurement.type == measurement.YIELD:
            if measurement.percentage.HasField('value'):
                string = f'{measurement.percentage.value:0.3f}%'
                if measurement.percentage.HasField('precision'):
                    string += f' (± {measurement.percentage.precision:0.3f}%)'
                return string
    return ''


def _parenthetical_if_def(string: str) -> str:
    """Returns a parenthesized version of a string, if defined."""
    if not string:
        return ''
    return f'({string})'


def _vessel_prep(vessel: reaction_pb2.Vessel) -> str:
    """Returns a description of the vessel preparation."""
    preparation_strings = []
    for preparation in vessel.preparations:
        preparation_strings.append({
            preparation.UNSPECIFIED: '',
            preparation.CUSTOM: 'prepared',
            preparation.NONE: '',
            preparation.OVEN_DRIED: 'oven-dried',
        }[preparation.type])
    return ', '.join(preparation_strings)


def _vessel_size(vessel: reaction_pb2.Vessel) -> str:
    """Returns a description of the vessel volume."""
    if vessel.volume.value:
        return f'{units.format_message(vessel.volume)}'
    return ''


def _vessel_material(vessel: reaction_pb2.Vessel) -> str:
    """Returns a description of the vessel material."""
    return {
        vessel.material.UNSPECIFIED: '',
        vessel.material.CUSTOM: f'{vessel.material.details}',
        vessel.material.GLASS: 'glass',
        vessel.material.POLYPROPYLENE: 'polypropylene',
        vessel.material.PLASTIC: 'plastic',
    }[vessel.material.type]


def _vessel_type(vessel: reaction_pb2.Vessel) -> str:
    """Returns a description of the vessel type.

    Args:
        vessel: Vessel message.

    Returns:
        String description of the vessel type.
    """
    return {
        vessel.UNSPECIFIED:
            'vessel',
        vessel.CUSTOM:
            'vessel',
        vessel.ROUND_BOTTOM_FLASK:
            'round bottom flask',
        vessel.VIAL:
            'vial',
        vessel.WELL_PLATE:
            'well-plate',
        vessel.MICROWAVE_VIAL:
            'microwave vial',
        vessel.TUBE:
            'tube',
        vessel.CONTINUOUS_STIRRED_TANK_REACTOR:
            'continuous stirred-tank reactor',
        vessel.PACKED_BED_REACTOR:
            'packed bed reactor',
    }[vessel.type]


def _input_addition(reaction_input: reaction_pb2.ReactionInput) -> str:
    """Returns a description of the addition of a ReactionInput.

    Args:
        reaction_input: ReactionInput message.

    Returns:
        String description of the addition of this input.
    """
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


def _uses_addition_order(reaction: reaction_pb2.Reaction) -> bool:
    """Returns whether any ReactionInput has a non-zero addition_order."""
    return any(value.addition_order for value in reaction.inputs.values())


def _round(value: float, places=2) -> str:
    """Rounds a value to the given number of decimal places."""
    fstring = '{:.%gg}' % places
    return fstring.format(value)


def _datetimeformat(message: reaction_pb2.DateTime,
                    format_string: str = '%Y-%m-%d / %H:%M') -> str:
    """Formats a date/time string."""
    value = parser.parse(message.value)
    return value.strftime(format_string)


def _get_compact_components(
    inputs: Mapping[str, reaction_pb2.ReactionInput]
) -> Iterable[Tuple[reaction_pb2.Compound, bool]]:
    """Returns a list of input components for 'compact' visualization.

    Args:
        inputs: Reaction inputs map.

    Yields:
        component: Compound message.
        is_last: Whether this is the last compound that will be displayed.
    """
    roles_to_keep = [
        reaction_pb2.ReactionRole.REACTANT,
        reaction_pb2.ReactionRole.UNSPECIFIED,
    ]
    compounds = []
    for _, value in _sort_addition_order(inputs):
        for component in value.components:
            if component.reaction_role in roles_to_keep:
                compounds.append(component)
    for i, compound in enumerate(compounds):
        if i + 1 == len(compounds):
            yield compound, True
        else:
            yield compound, False


def _get_compact_products(
    products: Iterable[reaction_pb2.ProductCompound]
) -> List[reaction_pb2.ProductCompound]:
    """Returns a list of product compounds for 'compact' visualization."""
    roles_to_keep = [
        reaction_pb2.ReactionRole.PRODUCT,
        reaction_pb2.ReactionRole.UNSPECIFIED,
    ]
    return [
        compound for compound in products
        if compound.reaction_role in roles_to_keep
    ]


def _value_and_precision(message) -> str:
    """Returns value +/- precision."""
    txt = f'{message.value:.7g}'
    if message.precision:
        txt += f' ± {message.precision:.7g}'
    return txt


def _product_measurement_value(message) -> str:
    """Returns the value for a product measurement."""
    if isinstance(message, reaction_pb2.Percentage):
        return _value_and_precision(message) + '%'
    if isinstance(message, reaction_pb2.FloatValue):
        return _value_and_precision(message)
    if isinstance(message, str):
        return message
    if isinstance(message, reaction_pb2.Amount):
        return _amount(message) or ''
    return ''


def _pbtxt(reaction: reaction_pb2.Reaction) -> str:
    """Converts a message to text format."""
    return text_format.MessageToString(reaction)


def _product_pbtxt(product: reaction_pb2.ProductCompound) -> str:
    """Converts a ProductCompound to text format without measurements."""
    trimmed = reaction_pb2.ProductCompound()
    trimmed.CopyFrom(product)
    del trimmed.measurements[:]
    return _pbtxt(trimmed)


def _oneof(message, name='kind'):
    """Retrieves the proper oneof value."""
    return getattr(message, message.WhichOneof(name))


def _defined(message):
    """Returns whether the message is defined (not empty)."""
    return message != type(message)()


def _type_and_details(message):
    """Returns type (details)."""
    assert len(message.DESCRIPTOR.enum_types) == 1
    value = message.DESCRIPTOR.enum_types[0].values_by_number[message.type].name
    if message.details:
        value += f' ({message.details})'
    return value


def _events(
    message: reaction_pb2.ReactionProvenance
) -> Iterable[reaction_pb2.RecordEvent]:
    """Returns a generator of RecordEvent messages."""
    events = [message.record_created]
    events.extend(message.record_modified)
    for event in events:
        yield event


TEMPLATE_FILTERS = {
    'round': _round,
    'is_true': _is_true,
    'datetimeformat': _datetimeformat,
    'uses_addition_order': _uses_addition_order,
    'get_compact_components': _get_compact_components,
    'get_compact_products': _get_compact_products,
    'input_addition': _input_addition,
    'compound_svg': _compound_svg,
    'compound_png': _compound_png,
    'amount': _amount,
    'compound_name': _compound_name,
    'compound_smiles': _compound_smiles,
    'compound_role': _compound_role,
    'compound_source_prep': _compound_source_prep,
    'product_yield': _product_yield,
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
    'pbtxt': _pbtxt,
    'product_pbtxt': _product_pbtxt,
    'oneof': _oneof,
    'defined': _defined,
    'type_and_details': _type_and_details,
    'product_measurement_value': _product_measurement_value,
    'events': _events,
}
