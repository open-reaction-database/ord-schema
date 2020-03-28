"""Helpers validating specific Message types."""

from ord_schema.proto import ord_schema_pb2 as schema

import re
import math
import dateutil.parser


def validate_message(message, recurse=True):
    """Template function for validating custom messages in the schema.

    Messages are not validated to check enum values, since these are enforced
    by the schema. Instead, we only check for validity of items that cannot be
    enforced in the schema (e.g., non-negativity of certain measurements,
    consistency of cross-referenced keys).

    Args:
        message: A message to validate.
        recurse: A boolean that controls whether submessages of message (i.e.,
            fields that are messages) should also be validated. Defaults to
            True.

    Returns:
        The input message, with any unambiguous changes made as
        needed to ensure validity.

    Raises:
        ValueError: If any fields are invalid.
    """

    # Recurse through submessages
    if recurse:
        for field, value in message.ListFields():
            if field.type == field.TYPE_MESSAGE:  # need to recurse
                if field.label == field.LABEL_REPEATED:
                    if field.message_type.GetOptions().map_entry:  # map
                        # value is message
                        if field.message_type.fields_by_name['value'].type == \
                                field.TYPE_MESSAGE:
                            for key, submessage in value.items():
                                submessage.CopyFrom(
                                    validate_message(submessage)
                                )
                        else:  # value is a primitive
                            pass
                    else:  # Just a repeated message
                        for submessage in value:
                            submessage.CopyFrom(
                                validate_message(submessage)
                            )
                else:  # no recursion needed
                    submessage = value
                    submessage.CopyFrom(
                        validate_message(submessage)
                    )

    # Message-specific validation
    try:
        return _VALIDATOR_SWITCH[type(message)](message)
    except KeyError:
        # NOTE(ccoley): I made the conscious decision to raise an error here,
        # rather than assume that the message is valid. If a message does not
        # require any message-level checks (not uncommon), then it should still
        # be listed in the dictionary switch above withpass. This will force
        # us to think about what is necessary if/when new messages are added.
        raise NotImplementedError(f"Don't know how to validate {type(message)}")


class ValidationWarning(Warning):
    pass


def ensure_float_nonnegative(message, field):
    if getattr(message, field) < 0:
        raise ValueError(f'Field {field} of message '
                         f'{type(message).DESCRIPTOR.name} must be'
                         ' non-negative')


def ensure_float_range(message, field, min=-math.inf, max=math.inf):
    if getattr(message, field) < min or getattr(message, field) > max:
        raise ValueError(f'Field {field} of message '
                         f'{type(message).DESCRIPTOR.name} must be between'
                         f' {min} and {max}')


def ensure_units_specified_if_value_defined(message):
    if message.value and message.units == message.UNSPECIFIED:
        raise ValueError(f'Unspecified units for {type(message)} with '
                         f'value defined ({message.value})')


def ensure_details_specified_if_type_custom(message):
    if message.type == message.CUSTOM and not message.details:
        raise ValueError(f'Custom type defined for {type(message)}, '
                         'but details field is empty')
    return message


def validate_reaction(message):
    if len(message.inputs) == 0:
        raise ValueError('Reactions should have at least 1 reaction input')
    # TODO(ccoley) Should outcomes also have a minimum length?
    return message


def validate_reaction_identifier(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_reaction_input(message):
    if len(message.components) == 0:
        raise ValueError('Reaction inputs must have at least one component')
    return message


def validate_compound(message):
    if len(message.identifiers) == 0:
        raise ValueError('Compounds must have at least one identifier')
    return message


def validate_compound_feature(message):
    return message


def validate_compound_preparation(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_compound_identifier(message):
    ensure_details_specified_if_type_custom(message)
    # TODO(ccoley): Add identifier-specific validation, e.g., by using
    # RDKit to try to parse SMILES, looking up NAMEs using online resolvers
    return message


def validate_vessel(message):
    if message.type == message.VesselType.CUSTOM and not message.details:
        raise ValueError('VesselType custom, but no details provided')
    if message.material == message.VesselMaterial.CUSTOM and \
            not message.material_details:
        raise ValueError('VesselMaterial custom, but no details provided')
    if message.preparation == message.VesselPreparation.CUSTOM and \
            not message.preparation_details:
        raise ValueError('VesselPreparation custom, but no details provided')
    return message


def validate_reaction_setup(message):
    return message


def validate_reaction_conditions(message):
    if message.conditions_are_dynamic and not message.details:
        raise ValueError('Reaction conditions are dynamic, but no details'
                         ' provided to explain how procedure deviates from'
                         ' normal single-step reaction conditions.')
    return message


def validate_temperature_conditions(message):
    if message.type == message.TemperatureControl.CUSTOM and \
            not message.details:
        raise ValueError('Temperature control custom, but no details provided')
    return message


def validate_temperature_measurement(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_pressure_conditions(message):
    if message.type == message.PressureControl.CUSTOM and not message.details:
        raise ValueError('Pressure control custom, but no details provided')
    if message.atmosphere == message.Atmosphere.CUSTOM and \
            not message.atmosphere_details:
        raise ValueError(
            'Atmosphere custom, but no atmosphere_details provided')
    return message


def validate_pressure_measurement(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_stirring_conditions(message):
    ensure_float_nonnegative(message, 'rpm')
    if message.type == message.StirringMethod.CUSTOM and not message.details:
        raise ValueError('Stirring method custom, but no details provided')
    return message


def validate_illumination_conditions(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_electrochemistry_conditions(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_electrochemistry_measurement(message):
    return message


def validate_flow_conditions(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_tubing(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_reaction_notes(message):
    return message


def validate_reaction_observation(message):
    return message


def validate_reaction_workup(message):
    ensure_details_specified_if_type_custom(message)
    return message


def validate_reaction_outcome(message):
    # Can only have one desired product
    if sum(product.is_desired_product for product in message.products) > 1:
        raise ValueError('Cannot have more than one desired product!')
    # Check key values for product analyses
    # NOTE(ccoley): Could use any(), but using expanded loops for clarity
    analysis_keys = list(message.analyses.keys())
    for product in message.products:
        for field in ['analysis_identity', 'analysis_yield', 'analysis_purity',
                      'analysis_selectivity']:
            for key in getattr(product, field):
                if key not in analysis_keys:
                    raise ValueError(f'Undefined analysis key {key} '
                                     'in ReactionProduct')
    return message


def validate_reaction_product(message):
    if message.texture == message.Texture.CUSTOM and \
            not message.texture_details:
        raise ValueError(f'Custom texture defined for {type(message)}, '
                         'but texture_details field is empty')
    return message


def validate_selectivity(message):
    ensure_float_nonnegative(message, 'precision')
    if message.type == message.EE:
        ensure_float_range(message, 'value', 0, 100)
        if 0 < message.value < 1:
            raise ValidationWarning('EE selectivity values are 0-100, '
                                    f'not fractions ({message.value} used)')
    ensure_details_specified_if_type_custom(message)
    return message


def validate_dateTime(message):
    if message.value:
        try:
            message.value = dateutil.parser.parse(message.value).ctime()
        except dateutil.parser.ParserError:
            raise ValueError(f'Could not parse DateTime string {message.value}')
    return message


def validate_reaction_analysis(message):
    # TODO(ccoley): Will be lots to expand here if we add structured data.
    ensure_details_specified_if_type_custom(message)
    return message


def validate_reaction_provenance(message):
    # Prepare datetimes
    if message.experiment_start.value:
        experiment_start = dateutil.parser.parse(
            message.experiment_start.value)
    if message.record_created.value:
        record_created = dateutil.parser.parse(message.record_created.value)
    if message.record_modified.value:
        record_modified = dateutil.parser.parse(message.record_created.value)
    # Check if record_created undefined
    if message.record_modified.value and not message.record_created.value:
        raise ValidationWarning('record_created not defined, but '
                                'record_modified is')
    # Check signs of time differences
    if message.experiment_start.value and message.record_created.value:
        if (record_created - experiment_start).total_seconds() < 0:
            raise ValueError('Record creation time should be after experiment')
    if message.record_modified.value and message.record_created.value:
        if (record_modified - record_created).total_seconds() < 0:
            raise ValueError('Record modified time should be after creation')
    # TODO(ccoley) could check if publication_url is valid, etc.
    return message


def validate_person(message):
    # NOTE(ccoley): final character is checksum, but ignoring that for now
    if message.orcid:
        if not re.match('[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]',
                        message.orcid):
            raise ValueError('Invalid ORCID: Enter as 0000-0000-0000-0000')
    return message


def validate_time(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_mass(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_moles(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_volume(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_concentration(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_pressure(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_temperature(message):
    if message.units == message.CELSIUS:
        ensure_float_range(message, 'value', min=-273.15)
    elif message.units == message.FAHRENHEIT:
        ensure_float_range(message, 'value', min=-459)
    elif message.units == message.KELVIN:
        ensure_float_range(message, 'value', min=0)
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_current(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_voltage(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_length(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_wavelength(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_flowRate(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)
    return message


def validate_percentage(message):
    if 0 < message.value < 1:
        raise ValidationWarning('Percentage values are 0-100, '
                                f'not fractions ({message.value} used)')
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_float_range(message, 'value', 0, 105)  # generous upper bound
    return message

_VALIDATOR_SWITCH = {
    schema.Reaction: validate_reaction,
    # Basics
    schema.ReactionIdentifier: validate_reaction_identifier,
    schema.ReactionInput: validate_reaction_input,
    # Compounds
    schema.Compound: validate_compound,
    schema.Compound.Feature: validate_compound_feature,
    schema.CompoundPreparation: validate_compound_preparation,
    schema.CompoundIdentifier: validate_compound_identifier,
    # Setup
    schema.Vessel: validate_vessel,
    schema.ReactionSetup: validate_reaction_setup,
    # Conditions
    schema.ReactionConditions: validate_reaction_conditions,
    schema.TemperatureConditions: validate_temperature_conditions,
    schema.TemperatureConditions.Measurement: validate_temperature_measurement,
    schema.PressureConditions: validate_pressure_conditions,
    schema.PressureConditions.Measurement: validate_pressure_measurement,
    schema.StirringConditions: validate_stirring_conditions,
    schema.IlluminationConditions: validate_illumination_conditions,
    schema.ElectrochemistryConditions: validate_electrochemistry_conditions,
    schema.ElectrochemistryConditions.Measurement: 
        validate_electrochemistry_measurement,
    schema.FlowConditions: validate_flow_conditions,
    schema.FlowConditions.Tubing: validate_tubing,
    # Annotations
    schema.ReactionNotes: validate_reaction_notes,
    schema.ReactionObservation: validate_reaction_observation,
    # Outcome
    schema.ReactionWorkup: validate_reaction_workup,
    schema.ReactionOutcome: validate_reaction_outcome,
    schema.ReactionProduct: validate_reaction_product,
    schema.Selectivity: validate_selectivity,
    schema.DateTime: validate_dateTime,
    schema.ReactionAnalysis: validate_reaction_analysis,
    # Metadata
    schema.ReactionProvenance: validate_reaction_provenance,
    schema.Person: validate_person,
    # Units
    schema.Time: validate_time,
    schema.Mass: validate_mass,
    schema.Moles: validate_moles,
    schema.Volume: validate_volume,
    schema.Concentration: validate_concentration,
    schema.Pressure: validate_pressure,
    schema.Temperature: validate_temperature,
    schema.Current: validate_current,
    schema.Voltage: validate_voltage,
    schema.Length: validate_length,
    schema.Wavelength: validate_wavelength,
    schema.FlowRate: validate_flowRate,
    schema.Percentage: validate_percentage,
}
