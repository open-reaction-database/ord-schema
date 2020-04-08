"""Helpers validating specific Message types."""

import math
import re
from dateutil import parser

from ord_schema.proto import reaction_pb2


def validate_message(message, recurse=True):
    """Template function for validating custom messages in the reaction_pb2.

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
    # pylint: disable=too-many-nested-blocks
    if recurse:
        for field, value in message.ListFields():
            if field.type == field.TYPE_MESSAGE:  # need to recurse
                if field.label == field.LABEL_REPEATED:
                    if field.message_type.GetOptions().map_entry:  # map
                        # value is message
                        if field.message_type.fields_by_name['value'].type == \
                                field.TYPE_MESSAGE:
                            for submessage in value.values():
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
    # pylint: enable=too-many-nested-blocks

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


# pylint: disable=missing-function-docstring
def ensure_float_nonnegative(message, field):
    if getattr(message, field) < 0:
        raise ValueError(f'Field {field} of message '
                         f'{type(message).DESCRIPTOR.name} must be'
                         ' non-negative')


def ensure_float_range(message, field, min_value=-math.inf, max_value=math.inf):
    if (getattr(message, field) < min_value or
            getattr(message, field) > max_value):
        raise ValueError(f'Field {field} of message '
                         f'{type(message).DESCRIPTOR.name} must be between'
                         f' {min_value} and {max_value}')


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
    if not message.value and not message.bytes_value:
        raise ValueError('{bytes_}value must be set')
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
    if not message.value and not message.bytes_value:
        raise ValueError('{bytes_}value must be set')
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


def validate_date_time(message):
    if message.value:
        try:
            message.value = parser.parse(message.value).ctime()
        except parser.ParserError:
            raise ValueError(f'Could not parse DateTime string {message.value}')
    return message


def validate_reaction_analysis(message):
    # TODO(ccoley): Will be lots to expand here if we add structured data.
    ensure_details_specified_if_type_custom(message)
    return message


def validate_reaction_provenance(message):
    # Prepare datetimes
    # TODO(kearnes): Require these to be set?
    experiment_start = None
    record_created = None
    record_modified = None
    if message.experiment_start.value:
        experiment_start = parser.parse(message.experiment_start.value)
    if message.record_created.time.value:
        record_created = parser.parse(message.record_created.time.value)
    for record in message.record_modified:
        # Use the last record as the most recent modification time.
        record_modified = parser.parse(record.time.value)
    # Check if record_created undefined
    if record_modified and not record_created:
        raise ValidationWarning('record_created not defined, but '
                                'record_modified is')
    # Check signs of time differences
    if experiment_start and record_created:
        if (record_created - experiment_start).total_seconds() < 0:
            raise ValueError('Record creation time should be after experiment')
    if record_modified and record_created:
        if (record_modified - record_created).total_seconds() < 0:
            raise ValueError('Record modified time should be after creation')
    # TODO(ccoley) could check if publication_url is valid, etc.
    return message


def validate_record_event(message):
    if not message.time.value:
        raise ValueError('RecordEvent must have `time` specified')
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
        ensure_float_range(message, 'value', min_value=-273.15)
    elif message.units == message.FAHRENHEIT:
        ensure_float_range(message, 'value', min_value=-459)
    elif message.units == message.KELVIN:
        ensure_float_range(message, 'value', min_value=0)
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


def validate_flow_rate(message):
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


def validate_binary_data(message):
    if not message.value:
        raise ValueError('value is required for BinaryData')
    if not message.format:
        raise ValidationWarning('No format specified for BinaryData')


# pylint: enable=missing-function-docstring

_VALIDATOR_SWITCH = {
    reaction_pb2.Reaction: validate_reaction,
    # Basics
    reaction_pb2.ReactionIdentifier: validate_reaction_identifier,
    reaction_pb2.ReactionInput: validate_reaction_input,
    # Compounds
    reaction_pb2.Compound: validate_compound,
    reaction_pb2.Compound.Feature: validate_compound_feature,
    reaction_pb2.CompoundPreparation: validate_compound_preparation,
    reaction_pb2.CompoundIdentifier: validate_compound_identifier,
    # Setup
    reaction_pb2.Vessel: validate_vessel,
    reaction_pb2.ReactionSetup: validate_reaction_setup,
    # Conditions
    reaction_pb2.ReactionConditions: validate_reaction_conditions,
    reaction_pb2.TemperatureConditions: validate_temperature_conditions,
    reaction_pb2.TemperatureConditions.Measurement: (
        validate_temperature_measurement),
    reaction_pb2.PressureConditions: validate_pressure_conditions,
    reaction_pb2.PressureConditions.Measurement: validate_pressure_measurement,
    reaction_pb2.StirringConditions: validate_stirring_conditions,
    reaction_pb2.IlluminationConditions: validate_illumination_conditions,
    reaction_pb2.ElectrochemistryConditions: (
        validate_electrochemistry_conditions),
    reaction_pb2.ElectrochemistryConditions.Measurement:
        validate_electrochemistry_measurement,
    reaction_pb2.FlowConditions: validate_flow_conditions,
    reaction_pb2.FlowConditions.Tubing: validate_tubing,
    # Annotations
    reaction_pb2.ReactionNotes: validate_reaction_notes,
    reaction_pb2.ReactionObservation: validate_reaction_observation,
    # Outcome
    reaction_pb2.ReactionWorkup: validate_reaction_workup,
    reaction_pb2.ReactionOutcome: validate_reaction_outcome,
    reaction_pb2.ReactionProduct: validate_reaction_product,
    reaction_pb2.Selectivity: validate_selectivity,
    reaction_pb2.DateTime: validate_date_time,
    reaction_pb2.ReactionAnalysis: validate_reaction_analysis,
    # Metadata
    reaction_pb2.ReactionProvenance: validate_reaction_provenance,
    reaction_pb2.ReactionProvenance.RecordEvent: validate_record_event,
    reaction_pb2.Person: validate_person,
    # Units
    reaction_pb2.Time: validate_time,
    reaction_pb2.Mass: validate_mass,
    reaction_pb2.Moles: validate_moles,
    reaction_pb2.Volume: validate_volume,
    reaction_pb2.Concentration: validate_concentration,
    reaction_pb2.Pressure: validate_pressure,
    reaction_pb2.Temperature: validate_temperature,
    reaction_pb2.Current: validate_current,
    reaction_pb2.Voltage: validate_voltage,
    reaction_pb2.Length: validate_length,
    reaction_pb2.Wavelength: validate_wavelength,
    reaction_pb2.FlowRate: validate_flow_rate,
    reaction_pb2.Percentage: validate_percentage,
    reaction_pb2.BinaryData: validate_binary_data,
}
