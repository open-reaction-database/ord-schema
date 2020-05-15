"""Template filters to use with Jinja to format reaction messages into text.
"""

from ord_schema import units
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import drawing


def _is_true(boolean):
    return boolean == reaction_pb2.Boolean.TRUE


def _sort_addition_order(inputs):
    for k in sorted(inputs.keys(), key=lambda k: inputs[k].addition_order):
        yield k, inputs[k]


def _stirring_conditions(stirring):
    if stirring.type == stirring.StirringMethod.NONE:
        return 'No stirring was used.'
    txt = ''
    txt += {
        stirring.StirringRate.UNSPECIFIED: '',
        stirring.StirringRate.HIGH: 'at a high rate',
        stirring.StirringRate.MEDIUM: 'at a medium rate',
        stirring.StirringRate.LOW: 'at a low rate',
    }[stirring.rate]
    if stirring.rpm:
        txt += f' ({stirring.rpm} rpm)'
    if stirring.type != stirring.StirringMethod.UNSPECIFIED:
        txt += ' using '
        txt += {
            stirring.StirringMethod.CUSTOM: 'a custom setup',
            stirring.StirringMethod.STIR_BAR: 'a stir bar',
            stirring.StirringMethod.OVERHEAD_MIXER: 'an overhead mixer',
            stirring.StirringMethod.AGITATION: 'external agitation',
        }[stirring.type]
        txt += f' {_parenthetical_if_def(stirring.details)}'
    if txt.strip():
        txt = 'The reaction mixture was stirred ' + txt + '.'
    return txt


def _stirring_conditions_html(stirring):
    if stirring.type == stirring.StirringMethod.NONE:
        return ''
    txt = ''
    if stirring.type != stirring.StirringMethod.UNSPECIFIED:
        txt += {
            stirring.StirringMethod.CUSTOM: stirring.details,
            stirring.StirringMethod.STIR_BAR: 'stir bar',
            stirring.StirringMethod.OVERHEAD_MIXER: 'overhead mixer',
            stirring.StirringMethod.AGITATION: 'agitation',
        }[stirring.type]
    if stirring.rpm:
        txt += f' ({stirring.rpm} rpm)'
    txt += '<br>'
    return txt


def _pressure_conditions(pressure):
    txt = ''
    if pressure.atmosphere != pressure.Atmosphere.UNSPECIFIED:
        txt += {
            pressure.Atmosphere.CUSTOM: 'under a custom atmosphere',
            pressure.Atmosphere.AIR: 'under air',
            pressure.Atmosphere.NITROGEN: 'under nitrogen',
            pressure.Atmosphere.ARGON: 'under argon',
            pressure.Atmosphere.OXYGEN: 'under oxygen',
            pressure.Atmosphere.HYDROGEN: 'under hydrogen',
        }[pressure.atmosphere]
        txt += f' {_parenthetical_if_def(pressure.atmosphere_details)}'
    if pressure.type != pressure.PressureControl.UNSPECIFIED:
        txt += ' '
        txt += {
            pressure.PressureControl.CUSTOM:
                'using a custom pressure controller',
            pressure.PressureControl.AMBIENT: 'using ambient pressure',
            pressure.PressureControl.BALLOON:
                'using a balloon for pressure control',
            pressure.PressureControl.SEALED:
                'after fully sealing the reaction vessel',
            pressure.PressureControl.SEPTUM_WITH_NEEDLE:
                'using a needle to pierce the vessel septum',
            pressure.PressureControl.RELEASEVALVE:
                'using a pressure release valve',
            pressure.PressureControl.BPR: 'using a backpressure regulator',
        }[pressure.type]
        txt += f' {_parenthetical_if_def(pressure.details)}'
        setpoint = units.format_message(pressure.setpoint)
        if setpoint:
            txt += f' with a setpoint of {setpoint}'
    if txt:
        txt = 'The reaction was run ' + txt + '.'
    return txt


def _pressure_conditions_html(pressure):
    txt = ''
    if pressure.atmosphere != pressure.Atmosphere.UNSPECIFIED:
        txt += {
            pressure.Atmosphere.CUSTOM: pressure.atmosphere_details,
            pressure.Atmosphere.AIR: 'in air',
            pressure.Atmosphere.NITROGEN: 'under nitrogen',
            pressure.Atmosphere.ARGON: 'under argon',
            pressure.Atmosphere.OXYGEN: 'under oxygen',
            pressure.Atmosphere.HYDROGEN: 'under hydrogen',
        }[pressure.atmosphere]
    if pressure.type != pressure.PressureControl.UNSPECIFIED:
        setpoint = units.format_message(pressure.setpoint)
        if setpoint:
            txt += f' ({setpoint})'
    if txt:
        txt += '<br>'
    return txt


def _temperature_conditions(temperature):
    txt = ''
    if temperature.type != temperature.TemperatureControl.UNSPECIFIED:
        txt += 'The reaction was run '
        txt += {
            temperature.TemperatureControl.CUSTOM:
                'under custom temperature conditions',
            temperature.TemperatureControl.AMBIENT:
                'under ambient temperature conditions',
            temperature.TemperatureControl.OIL_BATH: 'in an oil bath',
            temperature.TemperatureControl.WATER_BATH: 'in a water bath',
            temperature.TemperatureControl.SAND_BATH: 'in a sand bath',
            temperature.TemperatureControl.ICE_BATH: 'in an ice bath',
            temperature.TemperatureControl.DRY_ALUMINUM_PLATE:
                'using an aluminum heating block',
            temperature.TemperatureControl.MICROWAVE:
                'in a microwave reactor',
            temperature.TemperatureControl.DRY_ICE_BATH: 'in a dry ice bath',
            temperature.TemperatureControl.AIR_FAN:
                'using a fan for temperautre control',
            temperature.TemperatureControl.LIQUID_NITROGEN:
                'using liquid nitrogen for temperature control',
        }[temperature.type]
        txt += f' {_parenthetical_if_def(temperature.details)}'
        setpoint = units.format_message(temperature.setpoint)
        if setpoint:
            txt += f' with a setpoint of {setpoint}'
    return txt + '.'


def _temperature_conditions_html(temperature):
    txt = ''
    if (temperature.type == temperature.TemperatureControl.UNSPECIFIED or
            temperature.type == temperature.TemperatureControl.AMBIENT):
        return 'ambient temperature<br>'
    setpoint = units.format_message(temperature.setpoint)
    if setpoint:
        txt += f'{setpoint}'
    if txt:
        txt += '<br>'
    return txt


def _product_color_texture(product):
    txt = ''
    txt += f'{product.isolated_color} '
    txt += {
        product.Texture.UNSPECIFIED: '',
        product.Texture.CUSTOM: product.texture_details,
        product.Texture.POWDER:
            f'powder {_parenthetical_if_def(product.texture_details)}',
        product.Texture.CRYSTAL:
            f'set of crystals {_parenthetical_if_def(product.texture_details)}',
        product.Texture.OIL:
            f'oil {_parenthetical_if_def(product.texture_details)}',
    }[product.texture]
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
        analysis.NMR: 'NMR',
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


def _compound_name(compound, use_br=False):
    txt = ''
    for identifier in compound.identifiers:
        if identifier.type == identifier.NAME:
            txt += f'{identifier.value}'
    for identifier in compound.identifiers:
        if identifier.type == identifier.SMILES:
            if use_br:
                txt += '<br>'
            txt += f' "{identifier.value}"'
    if not txt:
        return '<UNK_COMPOUND>'
    return txt


def _compound_role(compound):
    limiting_if_true = {
        reaction_pb2.Boolean.UNSPECIFIED: '',
        reaction_pb2.Boolean.TRUE: 'limiting',
        reaction_pb2.Boolean.FALSE: '',
    }
    return {
        compound.ReactionRole.UNSPECIFIED: '',
        compound.ReactionRole.REACTANT:
            f'as a {limiting_if_true[compound.is_limiting]} reactant',
        compound.ReactionRole.REAGENT: 'as a reagent',
        compound.ReactionRole.SOLVENT: 'as a solvent',
        compound.ReactionRole.CATALYST: 'as a catalyst',
        compound.ReactionRole.INTERNAL_STANDARD: 'as an internal standard',
        compound.ReactionRole.WORKUP: '',
        compound.ReactionRole.PRODUCT: 'as a product',
    }[compound.reaction_role]


def _compound_source_prep(compound):
    txt = []
    if compound.vendor_source:
        txt.append(f'purchased from {compound.vendor_source}')
    if compound.vendor_id:
        txt.append(f'catalog #{compound.vendor_id}')
    if compound.vendor_lot:
        txt.append(f'lot #{compound.vendor_lot}')
    txt.append({
        compound.preparation.UNSPECIFIED: '',
        compound.preparation.CUSTOM: '',
        compound.preparation.NONE: '',
        compound.preparation.REPURIFIED: 'repurified',
        compound.preparation.SPARGED: 'sparged',
        compound.preparation.DRIED: 'dried',
        compound.preparation.SYNTHESIZED: 'synthesized in-house'
    }[compound.preparation.type])
    txt.append(compound.preparation.details)
    if any(elem for elem in txt):
        return '(' + ', '.join([elem for elem in txt if elem]) + ')'
    return ''


def _parenthetical_if_def(string):
    if not string:
        return ''
    return f'({string})'


def _vessel_prep(vessel):
    return {
        vessel.VesselPreparation.UNSPECIFIED: '',
        vessel.VesselPreparation.CUSTOM: 'prepared',
        vessel.VesselPreparation.NONE: '',
        vessel.VesselPreparation.OVEN_DRIED: 'oven-dried',
    }[vessel.preparation]


def _vessel_size(vessel):
    if vessel.volume.value:
        return f'{units.format_message(vessel.volume)}'
    return ''


def _vessel_material(vessel):
    return {
        vessel.VesselMaterial.UNSPECIFIED: '',
        vessel.VesselMaterial.CUSTOM: f'{vessel.material_details}',
        vessel.VesselMaterial.GLASS: f'glass',
        vessel.VesselMaterial.POLYPROPYLENE: f'polypropylene',
        vessel.VesselMaterial.PLASTIC: f'plastic',
    }[vessel.material]


def _vessel_type(vessel):
    return {
        vessel.VesselType.UNSPECIFIED: 'vessel',
        vessel.VesselType.CUSTOM: f'vessel',
        vessel.VesselType.ROUND_BOTTOM_FLASK: f'round bottom flask',
        vessel.VesselType.VIAL: f'vial',
        vessel.VesselType.WELL_PLATE: f'well-plate',
        vessel.VesselType.MICROWAVE_VIAL: f'microwave vial',
        vessel.VesselType.TUBE: f'tube',
        vessel.VesselType.CONTINUOUS_STIRRED_TANK_REACTOR:
            f'continuous stirred-tank reactor',
        vessel.VesselType.PACKED_BED_REACTOR: f'packed bed reactor',
    }[vessel.type]


def _input_addition(reaction_input):
    txt = []
    if reaction_input.addition_time.value:
        txt.append(
            f'after {units.format_message(reaction_input.addition_time)}')
    txt.append({
        reaction_input.AdditionSpeed.UNSPECIFIED: '',
        reaction_input.AdditionSpeed.ALL_AT_ONCE: 'all at once',
        reaction_input.AdditionSpeed.FAST: 'quickly',
        reaction_input.AdditionSpeed.SLOW: 'slowly',
        reaction_input.AdditionSpeed.DROPWISE: 'dropwise',
    }[reaction_input.addition_speed])
    if reaction_input.addition_duration.value:
        txt.append(
            f'over {units.format_message(reaction_input.addition_duration)}')
    if any(elem for elem in txt):
        return '(' + ', '.join([elem for elem in txt if elem]) + ')'
    return ''


def _uses_addition_order(reaction):
    return any(input.addition_order for input in reaction.inputs.values())


def _round(value, places=2):
    fstring = '{:.%gg}' % places
    return fstring.format(value)


def _datetimeformat(value, format_string='%H:%M / %d-%m-%Y'):
    return value.strftime(format_string)


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
    'sort_addition_order': _sort_addition_order,
}
