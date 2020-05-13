"""Helpers validating specific Message types."""

import math
import os
import re
import warnings

from absl import logging
from dateutil import parser
from rdkit import Chem
from rdkit import __version__ as RDKIT_VERSION

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


def validate_datasets(datasets, write_errors=False):
    """Runs validation for a set of datasets.

    Args:
        datasets: Dict mapping text filenames to Dataset protos.
        write_errors: If True, errors are written to disk.

    Raises:
        ValidationError: if any Dataset does not pass validation.
    """
    all_errors = []
    for filename, dataset in datasets.items():
        errors = _validate_dataset(filename, dataset)
        if errors:
            for error in errors:
                all_errors.append(f'{filename}: {error}')
            if write_errors:
                with open(f'{filename}.error', 'w') as f:
                    for error in errors:
                        f.write(f'{error}\n')
    # NOTE(kearnes): We run validation for all datasets before exiting if there
    # are errors.
    if all_errors:
        error_string = '\n'.join(all_errors)
        raise ValidationError(f'validation encountered errors:\n{error_string}')


def _validate_dataset(filename, dataset):
    """Validates Reaction messages in a Dataset.

    Note that validation may change the message. For example, NAME
    identifiers will be resolved to structures.

    Args:
        filename: Text filename; the dataset source.
        dataset: dataset_pb2.Dataset message.

    Returns:
        List of validation error messages.
    """
    basename = os.path.basename(filename)
    errors = []
    num_bad_reactions = 0
    for i, reaction in enumerate(dataset.reactions):
        reaction_errors = validate_message(reaction, raise_on_error=False)
        if reaction_errors:
            num_bad_reactions += 1
        for error in reaction_errors:
            errors.append(error)
            logging.warning('Validation error for %s[%d]: %s',
                            basename, i, error)
    logging.info('Validation summary for %s: %d/%d successful (%d failures)',
                 basename,
                 len(dataset.reactions) - num_bad_reactions,
                 len(dataset.reactions),
                 num_bad_reactions)
    return errors


# pylint: disable=too-many-branches
# pylint: disable=too-many-nested-blocks
def validate_message(message, recurse=True, raise_on_error=True):
    """Template function for validating custom messages in the reaction_pb2.

    Messages are not validated to check enum values, since these are enforced
    by the schema. Instead, we only check for validity of items that cannot be
    enforced in the schema (e.g., non-negativity of certain measurements,
    consistency of cross-referenced keys).

    Note that the message may be modified in-place with any unambiguous changes
    needed to ensure validity.

    Args:
        message: A message to validate.
        recurse: A boolean that controls whether submessages of message (i.e.,
            fields that are messages) should also be validated. Defaults to
            True.
        raise_on_error: If True, raises a ValidationError exception when errors
            are encountered. If False, the user must manually check the return
            value to identify validation errors.

    Returns:
        List of text validation errors, if any.

    Raises:
        ValidationError: If any fields are invalid.
    """
    errors = []
    # Recurse through submessages
    if recurse:
        for field, value in message.ListFields():
            if field.type == field.TYPE_MESSAGE:  # need to recurse
                if field.label == field.LABEL_REPEATED:
                    if field.message_type.GetOptions().map_entry:  # map
                        # value is message
                        if field.message_type.fields_by_name['value'].type == \
                                field.TYPE_MESSAGE:
                            for submessage in value.values():
                                errors.extend(validate_message(
                                    submessage, raise_on_error=raise_on_error))
                        else:  # value is a primitive
                            pass
                    else:  # Just a repeated message
                        for submessage in value:
                            errors.extend(validate_message(
                                submessage, raise_on_error=raise_on_error))
                else:  # no recursion needed
                    errors.extend(
                        validate_message(value, raise_on_error=raise_on_error))

    # Message-specific validation
    if not isinstance(message, tuple(_VALIDATOR_SWITCH.keys())):
        # NOTE(ccoley): I made the conscious decision to raise an error here,
        # rather than assume that the message is valid. If a message does not
        # require any message-level checks (not uncommon), then it should still
        # be listed in the dictionary switch above withpass. This will force
        # us to think about what is necessary if/when new messages are added.
        raise NotImplementedError(f"Don't know how to validate {type(message)}")

    with warnings.catch_warnings(record=True) as tape:
        _VALIDATOR_SWITCH[type(message)](message)
    for warning in tape:
        if issubclass(warning.category, ValidationError):
            if raise_on_error:
                raise warning.message
            errors.append(str(warning.message))
        else:
            warnings.warn(warning.message)
    return errors


# pylint: enable=too-many-branches
# pylint: enable=too-many-nested-blocks


class ValidationError(Warning):
    pass


class ValidationWarning(Warning):
    pass


# pylint: disable=missing-function-docstring
def ensure_float_nonnegative(message, field):
    if getattr(message, field) < 0:
        warnings.warn(f'Field {field} of message '
                      f'{type(message).DESCRIPTOR.name} must be'
                      ' non-negative', ValidationError)


def ensure_float_range(message, field, min_value=-math.inf, max_value=math.inf):
    if (getattr(message, field) < min_value or
            getattr(message, field) > max_value):
        warnings.warn(f'Field {field} of message '
                      f'{type(message).DESCRIPTOR.name} must be between'
                      f' {min_value} and {max_value}', ValidationError)


def ensure_units_specified_if_value_defined(message):
    if message.value and message.units == message.UNSPECIFIED:
        warnings.warn(f'Unspecified units for {type(message)} with '
                      f'value defined ({message.value})', ValidationError)


def ensure_details_specified_if_type_custom(message):
    if message.type == message.CUSTOM and not message.details:
        warnings.warn(f'Custom type defined for {type(message)}, '
                      'but details field is empty', ValidationError)


def reaction_has_internal_standard(message):
    """Whether any reaction component uses the internal standard role."""
    for reaction_input in message.inputs.values():
        for compound in reaction_input.components:
            if (compound.reaction_role ==
                    compound.ReactionRole.INTERNAL_STANDARD):
                return True
    for workup in message.workup:
        for compound in workup.components:
            if (compound.reaction_role ==
                    compound.ReactionRole.INTERNAL_STANDARD):
                return True
    return False


def reaction_has_limiting_component(message):
    """Whether any reaction input compound is limiting."""
    for reaction_input in message.inputs.values():
        for compound in reaction_input.components:
            if compound.is_limiting:
                return True
    return False


def reaction_needs_internal_standard(message):
    """Whether any analysis uses an internal standard."""
    for outcome in message.outcomes:
        for analysis in outcome.analyses.values():
            if analysis.uses_internal_standard:
                return True
    return False


def validate_dataset(message):
    if message.dataset_id:
        # The dataset_id is a 32-character uuid4 hex string.
        if not re.fullmatch('^ord_dataset-[0-9a-f]{32}$', message.dataset_id):
            warnings.warn('Dataset ID is malformed', ValidationError)


def validate_dataset_example(message):
    if not message.description:
        warnings.warn('DatasetExample.description is required', ValidationError)
    if not message.url:
        warnings.warn('DatasetExample.url is required', ValidationError)
    if not message.HasField('created'):
        warnings.warn('DatasetExample.created is required', ValidationError)


def validate_reaction(message):
    if len(message.inputs) == 0:
        warnings.warn('Reactions should have at least 1 reaction input',
                      ValidationError)
    if len(message.outcomes) == 0:
        warnings.warn('Reactions should have at least 1 reaction outcome',
                      ValidationError)
    if (reaction_needs_internal_standard(message)
            and not reaction_has_internal_standard(message)):
        warnings.warn('Reaction analysis uses an internal standard, but no '
                      'component (as reaction input or workup) uses the '
                      'reaction role INTERNAL_STANDARD', ValidationError)
    if (any(outcome.HasField('conversion') for outcome in message.outcomes)
            and not reaction_has_limiting_component(message)):
        warnings.warn('If reaction conversion is specified, at least one '
                      'reaction input component must be labeled is_limiting',
                      ValidationError)


def validate_reaction_identifier(message):
    ensure_details_specified_if_type_custom(message)
    if not message.value and not message.bytes_value:
        warnings.warn('{bytes_}value must be set', ValidationError)


def validate_reaction_input(message):
    if len(message.components) == 0:
        warnings.warn('Reaction inputs must have at least one component',
                      ValidationError)
    for component in message.components:
        if not component.WhichOneof('amount'):
            warnings.warn('Reaction input\'s components require an amount',
                          ValidationError)


def validate_compound(message):
    if len(message.identifiers) == 0:
        warnings.warn('Compounds must have at least one identifier',
                      ValidationError)
    if all(identifier.type == identifier.NAME
           for identifier in message.identifiers):
        warnings.warn('Compounds should have more specific identifiers than '
                      'NAME whenever possible', ValidationWarning)
    for identifier in message.identifiers:
        if Chem and identifier.type == identifier.SMILES:
            mol = Chem.MolFromSmiles(identifier.value)
            if mol is None:
                warnings.warn(f'RDKit {RDKIT_VERSION} could not validate'
                              f' SMILES identifier {identifier.value}',
                              ValidationError)
        elif identifier.type == identifier.INCHI:
            mol = Chem.MolFromInchi(identifier.value)
            if mol is None:
                warnings.warn(f'RDKit {RDKIT_VERSION} could not validate'
                              f' InChI identifier {identifier.value}',
                              ValidationError)
        elif identifier.type == identifier.MOLBLOCK:
            mol = Chem.MolFromMolBlock(identifier.value)
            if mol is None:
                warnings.warn(f'RDKit {RDKIT_VERSION} could not validate'
                              ' MolBlock identifier', ValidationError)
        elif identifier.type == identifier.RDKIT_BINARY:
            mol = Chem.Mol(identifier.bytes_value)
            if mol is None:
                warnings.warn(f'RDKit {RDKIT_VERSION} could not validate'
                              ' RDKit Binary identifier',
                              ValidationError)


def validate_compound_feature(message):
    if not message.name:
        warnings.warn('Compound features must have names', ValidationError)


def validate_compound_preparation(message):
    ensure_details_specified_if_type_custom(message)


def validate_compound_identifier(message):
    ensure_details_specified_if_type_custom(message)
    if not message.value and not message.bytes_value:
        warnings.warn('{bytes_}value must be set', ValidationError)


def validate_vessel(message):
    if message.type == message.VesselType.CUSTOM and not message.details:
        warnings.warn('VesselType custom, but no details provided',
                      ValidationError)
    if message.material == message.VesselMaterial.CUSTOM and \
            not message.material_details:
        warnings.warn('VesselMaterial custom, but no details provided',
                      ValidationError)
    if message.preparation == message.VesselPreparation.CUSTOM and \
            not message.preparation_details:
        warnings.warn('VesselPreparation custom, but no details provided',
                      ValidationError)


def validate_reaction_setup(message):
    del message  # Unused.


def validate_reaction_conditions(message):
    if message.conditions_are_dynamic and not message.details:
        warnings.warn('Reaction conditions are dynamic, but no details'
                      ' provided to explain how procedure deviates from'
                      ' normal single-step reaction conditions.',
                      ValidationError)
    if message.details and not message.conditions_are_dynamic:
        warnings.warn('Reaction condition details provided but field '
                      'conditions_are_dynamic is False. If the conditions '
                      'cannot be fully captured by the schema, set to True.',
                      ValidationWarning)


def validate_temperature_conditions(message):
    if message.type == message.TemperatureControl.CUSTOM and \
            not message.details:
        warnings.warn('Temperature control custom, but no details provided',
                      ValidationError)
    if not message.setpoint.value:
        warnings.warn('Temperature setpoints should be specified; even if '
                      'using ambient conditions, estimate room temperature and '
                      'the precision of your estimate.', ValidationWarning)


def validate_temperature_measurement(message):
    ensure_details_specified_if_type_custom(message)


def validate_pressure_conditions(message):
    if message.type == message.PressureControl.CUSTOM and not message.details:
        warnings.warn('Pressure control custom, but no details provided',
                      ValidationError)
    if message.atmosphere == message.Atmosphere.CUSTOM and \
            not message.atmosphere_details:
        warnings.warn(
            'Atmosphere custom, but no atmosphere_details provided',
            ValidationError)


def validate_pressure_measurement(message):
    ensure_details_specified_if_type_custom(message)


def validate_stirring_conditions(message):
    ensure_float_nonnegative(message, 'rpm')
    if message.type == message.StirringMethod.CUSTOM and not message.details:
        warnings.warn('Stirring method custom, but no details provided',
                      ValidationError)


def validate_illumination_conditions(message):
    ensure_details_specified_if_type_custom(message)


def validate_electrochemistry_conditions(message):
    ensure_details_specified_if_type_custom(message)


def validate_electrochemistry_measurement(message):
    del message  # Unused.


def validate_flow_conditions(message):
    ensure_details_specified_if_type_custom(message)


def validate_tubing(message):
    ensure_details_specified_if_type_custom(message)


def validate_reaction_notes(message):
    del message  # Unused.


def validate_reaction_observation(message):
    del message  # Unused.


def validate_reaction_workup(message):
    ensure_details_specified_if_type_custom(message)
    if (message.type == reaction_pb2.ReactionWorkup.WAIT and
            not message.duration.value):
        warnings.warn('"WAIT" workup steps require a defined duration',
                      ValidationError)
    if (message.type == reaction_pb2.ReactionWorkup.TEMPERATURE and
            not message.HasField('temperature')):
        warnings.warn('"TEMPERATURE" workup steps require defined '
                      'temperature conditions', ValidationError)
    if (message.type in (reaction_pb2.ReactionWorkup.EXTRACTION,
                         reaction_pb2.ReactionWorkup.FILTRATION) and
            not message.keep_phase):
        warnings.warn('Workup step EXTRACTION or FILTRATION missing '
                      'required field keep_phase', ValidationError)
    if (message.type in (reaction_pb2.ReactionWorkup.ADDITION,
                         reaction_pb2.ReactionWorkup.WASH,
                         reaction_pb2.ReactionWorkup.DRY_WITH_MATERIAL,
                         reaction_pb2.ReactionWorkup.SCAVENGING,
                         reaction_pb2.ReactionWorkup.DISSOLUTION,
                         reaction_pb2.ReactionWorkup.PH_ADJUST) and
            not message.components):
        warnings.warn('Workup step missing required components definition',
                      ValidationError)
    if (message.type == reaction_pb2.ReactionWorkup.STIRRING and
            not message.stirring):
        warnings.warn('Stirring workup step missing stirring definition',
                      ValidationError)
    if (message.type == reaction_pb2.ReactionWorkup.PH_ADJUST and
            not message.target_ph):
        warnings.warn('pH adjustment workup missing target pH',
                      ValidationError)


def validate_reaction_outcome(message):
    # Can only have one desired product
    if sum(product.is_desired_product for product in message.products) > 1:
        warnings.warn('Cannot have more than one desired product!',
                      ValidationError)
    # Check key values for product analyses
    # NOTE(ccoley): Could use any(), but using expanded loops for clarity
    analysis_keys = list(message.analyses.keys())
    for product in message.products:
        for field in ['analysis_identity', 'analysis_yield', 'analysis_purity',
                      'analysis_selectivity']:
            for key in getattr(product, field):
                if key not in analysis_keys:
                    warnings.warn(f'Undefined analysis key {key} '
                                  'in ReactionProduct', ValidationError)
    if not message.products and not message.HasField('conversion'):
        warnings.warn('No products or conversion are specified for reaction;'
                      ' at least one must be specified', ValidationError)


def validate_reaction_product(message):
    if message.texture == message.Texture.CUSTOM and \
            not message.texture_details:
        warnings.warn(f'Custom texture defined for {type(message)}, '
                      'but texture_details field is empty', ValidationError)


def validate_selectivity(message):
    ensure_float_nonnegative(message, 'precision')
    if message.type == message.EE:
        ensure_float_range(message, 'value', 0, 100)
        if 0 < message.value < 1:
            warnings.warn('EE selectivity values are 0-100, not fractions '
                          f'({message.value} used)', ValidationWarning)
    ensure_details_specified_if_type_custom(message)


def validate_date_time(message):
    if message.value:
        try:
            parser.parse(message.value).ctime()
        except parser.ParserError:
            warnings.warn(f'Could not parse DateTime string {message.value}',
                          ValidationError)


def validate_reaction_analysis(message):
    # TODO(ccoley): Will be lots to expand here if we add structured data.
    ensure_details_specified_if_type_custom(message)


def validate_reaction_provenance(message):
    # Prepare datetimes
    if not message.HasField('record_created'):
        warnings.warn('Reactions must have record_created defined.',
                      ValidationError)
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
    # Check signs of time differences
    if experiment_start and record_created:
        if (record_created - experiment_start).total_seconds() < 0:
            warnings.warn('Record creation time should be after experiment',
                          ValidationError)
    if record_modified and record_created:
        if (record_modified - record_created).total_seconds() < 0:
            warnings.warn('Record modified time should be after creation',
                          ValidationError)
    if message.record_id:
        # The record_id suffix is a 32-character uuid4 hex string.
        if not re.fullmatch('^ord-[0-9a-f]{32}$', message.record_id):
            warnings.warn('Record ID is malformed', ValidationError)
    # TODO(ccoley) could check if publication_url is valid, etc.


def validate_record_event(message):
    if not message.time.value:
        warnings.warn('RecordEvent must have `time` specified',
                      ValidationError)


def validate_person(message):
    # NOTE(ccoley): final character is checksum, but ignoring that for now
    if message.orcid:
        if not re.match('[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]',
                        message.orcid):
            warnings.warn('Invalid ORCID: Enter as 0000-0000-0000-0000',
                          ValidationError)


def validate_time(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_mass(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_moles(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_volume(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_concentration(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_pressure(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_temperature(message):
    if message.units == message.CELSIUS:
        ensure_float_range(message, 'value', min_value=-273.15)
    elif message.units == message.FAHRENHEIT:
        ensure_float_range(message, 'value', min_value=-459)
    elif message.units == message.KELVIN:
        ensure_float_range(message, 'value', min_value=0)
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_current(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_voltage(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_length(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_wavelength(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_flow_rate(message):
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_units_specified_if_value_defined(message)


def validate_percentage(message):
    if 0 < message.value < 1:
        warnings.warn('Percentage values are 0-100, not fractions '
                      f'({message.value} used)', ValidationWarning)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')
    ensure_float_range(message, 'value', 0, 105)  # generous upper bound


def validate_data(message):
    # TODO(kearnes): Validate/ping URLs?
    if not message.WhichOneof('kind'):
        warnings.warn('Data requires one of {value, bytes_value, url}',
                      ValidationError)
    if message.bytes_value and not message.format:
        warnings.warn('Data format is required for bytes_data',
                      ValidationError)


# pylint: enable=missing-function-docstring

_VALIDATOR_SWITCH = {
    dataset_pb2.Dataset: validate_dataset,
    dataset_pb2.DatasetExample: validate_dataset_example,
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
    reaction_pb2.RecordEvent: validate_record_event,
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
    reaction_pb2.Data: validate_data,
}
