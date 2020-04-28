"""Text generation script for Reaction messages.

This script is meant to convert a Reaction message into a full-fledged text
paragraph that would be appropriate for inclusion in a Supplemental Information
document for a publication.

Example usage:
* For normal operation from a pbtxt
  $ python generate_text.py --input_reaction=reaction.pb --input_format=pbtxt

# TODO(ccoley) figure out what ideal output format is
"""

import os
import jinja2
import re

from absl import app
from absl import flags
from absl import logging

from ord_schema import units
from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input_file', None,
                    'File containing a Reaction message.')
flags.DEFINE_enum('input_format', 'binary',
                  [f.value for f in message_helpers.MessageFormats],
                  'Input message format.')
flags.DEFINE_string('output', None, 'Filename for output Dataset.')


_TEMPLATE = '''
To a 

{# VESSEL #}
{{ reaction.setup.vessel|vessel_prep }} 
{{ reaction.setup.vessel.preparation_details|parenthetical_if_def }} 
{{ reaction.setup.vessel|vessel_size }} 
{{ reaction.setup.vessel|vessel_material }} 
{{ reaction.setup.vessel.material_details|parenthetical_if_def }} 
{{ reaction.setup.vessel|vessel_type }} 
{{ reaction.setup.vessel.details|parenthetical_if_def }} 

was 

{# AUTOMATION #}
{% if reaction.setup.is_automated %} automatically {% endif %}
{{ reaction.setup.automation_platform|parenthetical_if_def }}

added 

{# REACTION INPUTS - SORT BY ADDITION ORDER #}
{% for group in reaction.inputs.values()|groupby('addition_order') %}
    {% if loop.last %} and {% endif %}

    {% if reaction|uses_addition_order %}
        ({{ group.list[0].addition_order }})
    {% endif %}

    {% for input in group.list %}
        {% for compound in input.components %}

            {# COMPOUND DETAILS #}
            {{ compound|compound_amount }}
            {{ compound|compound_name }} 
            {{ compound|compound_role }} 
            {{ compound|compound_source_prep }} 

            {% if not loop.last %} + {% endif %}
        {% endfor %}
        {{ input|input_addition }}{% if not loop.last %}, {% endif %}
    {% endfor %}

    {% if not loop.last %}; {% endif %}
{% endfor %}. 

{# TEMPERATURE #}
{{ reaction.conditions.temperature|temperature_conditions }} 

{# PRESSURE #}
{{ reaction.conditions.pressure|pressure_conditions }} 

{# STIRRING #}
{{ reaction.conditions.stirring|stirring_conditions }} 

{# TODO(ccoley) ILLUMINATION #}

{# TODO(ccoley) ELECTROCHEMISTRY #}

{# TODO(ccoley) FLOW #}

{# OBSERVATIONS - images not handled #}
{% if reaction.observations %}
    During the reaction, the following observations were noted: 
    {% for observation in reaction.observations %}
        {% if loop.last  and loop.index > 1 %} and {% endif %}
        {% if observation.comment %}
            after 
            {% if observation.time.value %}
                {{ observation.time|unit_format }}
            {% else %}
                an unspecified amount of time
            {% endif %}, 
            {{ observation.comment }}
            {% if not loop.last %}; {% endif %}
        {% endif %}
    {% endfor %}.
{% endif %}

{# TODO(ccoley) WORKUP #}
{% if reaction.workup %}
    The workup procedure consisted of 
    {% for workup in reaction.workup %}
        {% if loop.last  and loop.index > 1 %} and {% endif %}
        ({{ loop.index }}) 
            <WORKUP_STEP>
        {% if not loop.last %}; {% endif %}
    {% endfor %}
{% endif %}

{# OUTCOME #}
{% for outcome in reaction.outcomes %}
    The reaction was {% if not loop.first %} also {% endif %} analyzed after 
    {% if outcome.reaction_time.value %}
        {{ outcome.reaction_time|unit_format }}
    {% else %}
        an unspecified amount of time
    {% endif %} 
    by  
    {% for analysis in outcome.analyses.values() %}
        {% if loop.last and loop.index > 1 %} and {% endif %}
        {{ analysis|analysis_format }}
        {{ analysis.details|parenthetical_if_def }}
        {% if not loop.last %}, {% endif %}
    {% endfor %}.

    {# PRODUCTS #}
    {% for product in outcome.products %}
        {{ product.compound|compound_name }}
        {% if product.compound.is_desired_product %}
            , the desired product,
        {% endif %}

        {# PRODUCT IDENTITY #}
        {% if product.analysis_identity %}
            was identified through
            {% for key in product.analysis_identity %}
                {% if loop.last and loop.index > 1 %} and {% endif %}
                {{ outcome.analyses[key]|analysis_format }}
                {% if not loop.last %}, {% endif %}
            {% endfor %}
        {% else %}
            was observed
        {% endif %}
        
        {# PRODUCT YIELD #}
        {% if product.compound_yield.value %}
            with a yield of {{ product.compound_yield.value|round(3) }}%
            {% if product.compound_yield.precision %}
                (p/m {{ product.compound_yield.precision|round(3) }}%)
            {% endif %}
            {% if product.analysis_yield %}
                (assayed through
                {% for key in product.analysis_yield %}
                    {% if loop.last and loop.index > 1 %} and {% endif %}
                    {{ outcome.analyses[key]|analysis_format }}
                    {% if not loop.last %}, {% endif %}
                {% endfor %})
            {% endif %}
        {% endif %}

        {# PRODUCT YIELD #}
        {% if product.selectivity.value %}
            with a selectivity of {{ product.selectivity.value|round(3) }}
            {{ product.selectivity|selectivity_type }}
            {% if product.selectivity.precision %}
                (p/m {{ product.selectivity.precision|round(3) }})
            {% endif %}
            {% if product.analysis_selectivity %}
                (assayed through
                {% for key in product.analysis_selectivity %}
                    {% if loop.last and loop.index > 1 %} and {% endif %}
                    {{ outcome.analyses[key]|analysis_format }}
                    {% if not loop.last %}, {% endif %}
                {% endfor %})
            {% endif %}
        {% endif %}

        {# PRODUCT PURITY #}
        {% if product.purity.value %}
            with a purity of {{ product.purity.value|round(3) }}%
            {% if product.purity.precision %}
                (p/m {{ product.purity.precision|round(3) }}%)
            {% endif %}
            {% if product.analysis_purity %}
                (assayed through
                {% for key in product.analysis_purity %}
                    {% if loop.last and loop.index > 1 %} and {% endif %}
                    {{ outcome.analyses[key]|analysis_format }}
                    {% if not loop.last %}, {% endif %}
                {% endfor %})
            {% endif %}
        {% endif %}.

        {# PRODUCT COLOR/TEXTURE #}
        {{ product|product_color_texture }}.


    {% endfor %}

{% endfor %}

{# DICLAIMER ABOUT SCHEMA LIMITATIONS #}
{% if reaction.conditions.conditions_are_dynamic %}
    Note: the conditions for this reaction may not precisely fit the schema, so
    additional care should be taken when reading this description.
{% endif %}
{% if reaction.conditions.details %}
    Condition details: {{ reaction.conditions.details }}.
{% endif %}
{% if reaction.notes.procedure_details %}
    Procedure details: {{ reaction.notes.procedure_details }}.
{% endif %}

'''


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
    if not txt.strip():
        return ''
    return 'The reaction mixture was stirred ' + txt + '.'


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
        return 'The reaction was run ' + txt + '.'


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


def _reaction_has_notes(reaction):
    return (reaction.notes)


def _compound_amount(compound):
    amount = compound.WhichOneof('amount')
    if not amount:
        return ''
    return units.format_message(getattr(compound, amount))


def _compound_name(compound):
    txt = ''
    for identifier in compound.identifiers:
        if identifier.type == identifier.NAME:
            txt += f'{identifier.value}'
    for identifier in compound.identifiers:
        if identifier.type == identifier.SMILES:
            txt += f' "{identifier.value}"'
    if not txt:
        return '<UNK_COMPOUND>'
    return txt


def _compound_role(compound):
    limiting_if_true = {True: 'limiting', False: ''}
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
        txt.append(f'lot #{vendor_lot}')
    txt.append({
        compound.preparation.UNSPECIFIED: '',
        compound.preparation.CUSTOM: '',
        compound.preparation.NONE: 'used as received',
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


def _input_addition(input):
    txt = []
    if input.addition_time.value:
        txt.append(f'after {units.format_message(input.addition_time)}')
    txt.append({
        input.AdditionSpeed.UNSPECIFIED: '',
        input.AdditionSpeed.ALL_AT_ONCE: 'all at once',
        input.AdditionSpeed.FAST: 'quickly',
        input.AdditionSpeed.SLOW: 'slowly',
        input.AdditionSpeed.DROPWISE: 'dropwise',
    }[input.addition_speed])
    if input.addition_duration.value:
        txt.append(f'over {units.format_message(input.addition_duration)}')
    if any(elem for elem in txt):
        return '(' + ', '.join([elem for elem in txt if elem]) + ')'
    return ''


def _uses_addition_order(reaction):
    return any(input.addition_order for input in reaction.inputs.values())


def _round(value, places=2):
    fstring = '{:.%gg}' % places
    return fstring.format(value)


def _datetimeformat(value, format='%H:%M / %d-%m-%Y'):
    return value.strftime(format)


def generate_text(reaction, template_string=_TEMPLATE):
    """Generates a reaction description.

    Args:
        reaction: a Reaction message.

    Raises:
        #TODO

    Returns:
        String description
    """
    env = jinja2.Environment(loader=jinja2.BaseLoader())
    env.filters['round'] = _round
    env.filters['datetimeformat'] = _datetimeformat
    env.filters['uses_addition_order'] = _uses_addition_order
    env.filters['input_addition'] = _input_addition
    env.filters['compound_amount'] = _compound_amount
    env.filters['compound_name'] = _compound_name
    env.filters['compound_role'] = _compound_role
    env.filters['compound_source_prep'] = _compound_source_prep
    env.filters['vessel_prep'] = _vessel_prep
    env.filters['vessel_type'] = _vessel_type
    env.filters['vessel_material'] = _vessel_material
    env.filters['vessel_size'] = _vessel_size
    env.filters['unit_format'] = units.format_message
    env.filters['parenthetical_if_def'] = _parenthetical_if_def
    env.filters['analysis_format'] = _analysis_format
    env.filters['selectivity_type'] = _selectivity_type
    env.filters['product_color_texture'] = _product_color_texture
    env.filters['temperature_conditions'] = _temperature_conditions
    env.filters['pressure_conditions'] = _pressure_conditions
    env.filters['stirring_conditions'] = _stirring_conditions
    template = env.from_string(template_string)
    text = template.render(reaction=reaction)

    # Fix line breaks, extra spaces, "a" versus "an"
    text = ''.join(text.strip().splitlines())
    text = re.sub(r'[ ]{2,}', ' ', text)
    text = re.sub(r' a ([aeiouAEIOU])', r' an \1', text)
    text = re.sub(r'[ ]\.', '.', text)
    text = re.sub(r'[ ]\;', ';', text)
    text = re.sub(r'[ ]\)', ')', text)
    text = re.sub(r'[ ]\,', ',', text)
    text = re.sub(r'\.\.', '.', text)
    return text


def main(argv):
    del argv  # Only used by app.run()
    reaction = message_helpers.load_message(
        FLAGS.input_file, reaction_pb2.Reaction, FLAGS.input_format)
    text = generate_text(reaction)
    if FLAGS.output:
        with open(FLAGS.output, 'w') as fid:
            fid.write(text)
    else:
        print(text)


if __name__ == '__main__':
    flags.mark_flag_as_required('input_file')
    app.run(main)
