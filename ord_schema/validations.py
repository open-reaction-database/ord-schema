"""Helpers validating specific Message types."""

from ord_schema.proto import ord_schema_pb2 as schema

import re
from dateutil import parser as dateparser

def ValidateMessage(message, recurse=True):
    """Template function for validating custom messages in the schema.

        Messages are not validated to check enum values, since these
        are enforced by the schema. Instead, we only check for validity
        of items that cannot be enforced in the schema (e.g., non-negativity
        of certain measurements, consistency of cross-referenced keys)
        
        Args:
            message: a message to validate.

        Returns:
            The input message, with any unambiguous changes made as
            needed to ensure validity.

        Raises:
            ValueError: If any fields are invalid."""

    # Recurse through submessages
    if recurse:
        for field in message.DESCRIPTOR.fields:
            if field.type == field.TYPE_MESSAGE: # recurse
                if field.label == field.LABEL_REPEATED:
                    if field.message_type.GetOptions().map_entry: # map
                        for key,submessage in getattr(message, field.name).items():
                            submessage.CopyFrom(
                                ValidateMessage(submessage)
                            )
                    else: # Just a repeated message
                        for submessage in getattr(message, field.name):
                            submessage.CopyFrom(
                                ValidateMessage(submessage)
                            )
                else: # no recursion needed
                    submessage = getattr(message, field.name)
                    submessage.CopyFrom(
                        ValidateMessage(submessage)
                    )

    # Message-specific validation
    try:
        return {
            schema.Reaction: ValidateReaction,
            # Basics
            schema.ReactionIdentifier: ValidateReactionIdentifier,
            schema.ReactionInput: ValidateReactionInput,
            # Compounds
            schema.Compound: ValidateCompound,
            schema.Compound.Feature: ValidateCompoundFeature,
            schema.CompoundPreparation: ValidateCompoundPreparation,
            schema.CompoundIdentifier: ValidateCompoundIdentifier,
            # Setup
            schema.Vessel: ValidateVessel,
            schema.ReactionSetup: ValidateReactionSetup,
            # Conditions
            schema.ReactionConditions: ValidateReactionConditions,
            schema.TemperatureConditions: ValidateTemperatureConditions,
            schema.TemperatureConditions.Measurement: ValidateTemperatureMeasurement,
            schema.PressureConditions: ValidatePressureConditions,
            schema.PressureConditions.Measurement: ValidatePressureMeasurement,
            schema.StirringConditions: ValidateStirringConditions,
            schema.IlluminationConditions: ValidateIlluminationConditions,
            schema.ElectrochemistryConditions: ValidateElectrochemistryConditions,
            schema.ElectrochemistryConditions.Measurement: ValidateElectrochemistryMeasurement,
            schema.FlowConditions: ValidateFlowConditions,
            schema.FlowConditions.Tubing: ValidateTubing,
            # Annotations
            schema.ReactionNotes: ValidateReactionNotes,
            schema.ReactionObservation: ValidateReactionObservation,
            # Outcome
            schema.ReactionWorkup: ValidateReactionWorkup,
            schema.ReactionOutcome: ValidateReactionOutcome,
            schema.ReactionProduct: ValidateReactionProduct,
            schema.Selectivity: ValidateSelectivity,
            schema.DateTime: ValidateDateTime,
            schema.ReactionAnalysis: ValidateReactionAnalysis,
            # Metadata
            schema.ReactionProvenance: ValidateReactionProvenance,
            schema.Person: ValidatePerson,
            # Units
            schema.Time: ValidateTime,
            schema.Mass: ValidateMass,
            schema.Moles: ValidateMoles,
            schema.Volume: ValidateVolume,
            schema.Concentration: ValidateConcentration,
            schema.Pressure: ValidatePressure,
            schema.Temperature: ValidateTemperature,
            schema.Current: ValidateCurrent,
            schema.Voltage: ValidateVoltage,
            schema.Length: ValidateLength,
            schema.Wavelength: ValidateWavelength,
            schema.FlowRate: ValidateFlowRate,
            schema.Percentage: ValidatePercentage,
        }[type(message)](message)
    except KeyError:
        # NOTE(ccoley): I made the conscious decision to raise an error here,
        # rather than assume that the message is valid. If a message does not
        # require any message-level checks (not uncommon), then it should still
        # be listed in the dictionary switch above withpass. This will force
        # us to think about what is necessary if/when new messages are added.
        raise NotImplementedError(f"Don't know how to validate {type(message)}")


class ValidationWarning(Warning):
    pass

def return_message_if_valid(func):
    def wrapper(message):
        func(message)
        return message
    return wrapper

def ensure_float_nonnegative(message, field):
    if getattr(message, field) < 0:
        raise ValueError(f'Field {field} of message '\
            f'{type(message).DESCRIPTOR.name} must be non-negative')

def ensure_float_range(message, field, min, max):
    if getattr(message, field) < min or getattr(message, field) > max:
        raise ValueError(f'Field {field} of message '\
            f'{type(message).DESCRIPTOR.name} must be between {min} and {max}')

def ensure_units_specified_if_value_defined(message):
    if message.value: # can't distinguish 0 and unspecified anyway
        if message.units == message.UNSPECIFIED:
            raise ValueError(f'Unspecified units for {type(message)} with ' \
                f'value defined ({message.value})')

def ensure_details_specified_if_type_custom(message):
    if message.type == message.CUSTOM:
        if not message.details:
            raise ValueError(f'Custom type defined for {type(message)}, ' \
                'but details field is empty')

@return_message_if_valid
def ValidateReaction(message):
    if len(message.inputs) == 0:
        raise ValueError('Reactions should have at least 1 reaction input')
    # TODO(ccoley) Should outcomes also have a minimum length?

@return_message_if_valid
def ValidateReactionIdentifier(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateReactionInput(message):
    if len(message.components) == 0:
        raise ValueError('Reaction inputs must have at least one component')

@return_message_if_valid
def ValidateCompound(message):
    if len(message.identifiers) == 0:
        raise ValueError('Compounds must have at least one identifier')

@return_message_if_valid
def ValidateCompoundFeature(message):
    pass

@return_message_if_valid
def ValidateCompoundPreparation(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateCompoundIdentifier(message):
    ensure_details_specified_if_type_custom(message)
    # TODO(ccoley): Add identifier-specific validation, e.g., by using
    # RDKit to try to parse SMILES, looking up NAMEs using online resolvers

@return_message_if_valid
def ValidateVessel(message):
    if message.type == message.VesselType.CUSTOM:
        if not message.details:
            raise ValueError('VesselType custom, but no details provided')
    if message.material == message.VesselMaterial.CUSTOM:
        if not message.material_details:
            raise ValueError('VesselMaterial custom, but no details provided')
    if message.preparation == message.VesselPreparation.CUSTOM:
        if not message.preparation_details:
            raise ValueError('VesselPreparation custom, but no details provided')

@return_message_if_valid
def ValidateReactionSetup(message):
    pass

@return_message_if_valid
def ValidateReactionConditions(message):
    if message.conditions_are_dynamic:
        if not message.details:
            raise ValueError('Reaction conditions are dynamic, but no details' \
                ' provided to explain how procedure deviates from normal ' \
                'single-step reaction conditions.')

@return_message_if_valid
def ValidateTemperatureConditions(message):
    if message.type == message.TemperatureControl.CUSTOM:
        if not message.details:
            raise ValueError('Temperature control custom, but no details provided')

@return_message_if_valid
def ValidateTemperatureMeasurement(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidatePressureConditions(message):
    if message.type == message.PressureControl.CUSTOM:
        if not message.details:
            raise ValueError('Pressure control custom, but no details provided')
    if message.atmosphere == message.Atmosphere.CUSTOM:
        if not message.atmosphere_details:
            raise ValueError('Atmosphere custom, but no atmosphere_details provided')

@return_message_if_valid
def ValidatePressureMeasurement(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateStirringConditions(message):
    ensure_float_nonnegative(message, 'rpm')
    if message.type == message.StirringMethod.CUSTOM:
        if not message.details:
            raise ValueError('Stirring method custom, but no details provided')

@return_message_if_valid
def ValidateIlluminationConditions(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateElectrochemistryConditions(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateElectrochemistryMeasurement(message):
    pass

@return_message_if_valid
def ValidateFlowConditions(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateTubing(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateReactionNotes(message):
    pass

@return_message_if_valid
def ValidateReactionObservation(message):
    pass

@return_message_if_valid
def ValidateReactionWorkup(message):
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateReactionOutcome(message):
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
                    raise ValueError(f'Undefined analysis key {key} '\
                        'in ReactionProduct')

@return_message_if_valid
def ValidateReactionProduct(message):
    if message.texture == message.Texture.CUSTOM:
        if not message.texture_details:
            raise ValueError(f'Custom texture defined for {type(message)}, ' \
                'but texture_details field is empty')

@return_message_if_valid
def ValidateSelectivity(message):
    ensure_float_nonnegative(message, 'precision')
    if message.type == message.EE:
        ensure_float_range(message, 'value', 0, 100)
        if message.value > 0 and message.value < 1:
            raise ValidationWarning('EE selectivity values are 0-100, ' \
                f'not fractions ({message.value} used)')
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateDateTime(message):
    if message.value:
        try:
            message.value = dateparser.parse(message.value).ctime()
        except dateparser.ParserError:
            raise ValueError(f'Could not parse DateTime string {message.value}')

@return_message_if_valid
def ValidateReactionAnalysis(message):
    # TODO(ccoley): Will be lots to expand here if we add structured data.
    ensure_details_specified_if_type_custom(message)

@return_message_if_valid
def ValidateReactionProvenance(message):
    # Prepare datetimes
    if message.experiment_start.value:
        experiment_start = dateparser.parse(message.experiment_start.value)
    if message.record_created.value:
        record_created = dateparser.parse(message.record_created.value)
    if message.record_modified.value:
        record_modified = dateparser.parse(message.record_created.value)
    # Check if record_created undefined
    if message.record_modified.value and not message.record_created.value:
        raise ValidationWarning('record_created not defined, but ' \
            'record_modified is')
    # Check signs of time differences
    if message.experiment_start.value and message.record_created.value:
        if (record_created - experiment_start) < 0:
            raise ValueError('Record creation time should be after experiment')
    if message.record_modified.value and message.record_created.value:
        if (record_modified - record_created) < 0:
            raise ValueError('Record modified time should be after creation')
    # TODO(ccoley) could check if publication_url is valid, etc.

@return_message_if_valid
def ValidatePerson(message):
    # NOTE(ccoley): final character is checksum, but ignoring that for now
    if message.orcid:
        if not re.match('[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]', 
                message.orcid):
            raise ValueError('Invalid ORCID: Enter as 0000-0000-0000-0000')

@return_message_if_valid
def ValidateTime(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateMass(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateMoles(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateVolume(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateConcentration(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidatePressure(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateTemperature(message):
    if message.units == message.CELSIUS:
        ensure_float_nonnegative(message, 'value')
    elif message.units == message.FAHRENHEIT:
        ensure_float_nonnegative(message + 459.67, 'value')
    elif message.units == message.KELVIN:
        ensure_float_nonnegative(message + 273.15, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateCurrent(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateVoltage(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateLength(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateWavelength(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidateFlowRate(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)

@return_message_if_valid
def ValidatePercentage(message):
    if message.value > 0 and message.value < 1:
        raise ValidationWarning('Percentage values are 0-100, ' \
            f'not fractions ({message.value} used)')
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_float_range(message, 'value', 0, 105)  # generous upper bound


