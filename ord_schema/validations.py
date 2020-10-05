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
"""Helpers validating specific Message types."""

import dataclasses
import math
import os
import re
from typing import List
import warnings

from absl import logging
from dateutil import parser
from rdkit import Chem
from rdkit import __version__ as RDKIT_VERSION

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


@dataclasses.dataclass
class ValidationOptions:
    """Options for message validation."""
    # Check that Dataset and Reaction IDs are well-formed.
    validate_ids: bool = False
    # Require ReactionProvenance for Reactions.
    require_provenance: bool = False


@dataclasses.dataclass
class ValidationOutput:
    """Validation output: errors and warnings."""
    errors: List[str] = dataclasses.field(default_factory=list)
    warnings: List[str] = dataclasses.field(default_factory=list)

    def extend(self, other):
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)


def validate_datasets(datasets, write_errors=False, options=None):
    """Runs validation for a set of datasets.

    Args:
        datasets: Dict mapping text filenames to Dataset protos.
        write_errors: If True, errors are written to disk.
        options: ValidationOptions.

    Raises:
        ValidationError: if any Dataset does not pass validation.
    """
    all_errors = []
    for filename, dataset in datasets.items():
        basename = os.path.basename(filename)
        errors = _validate_datasets(dataset, label=basename, options=options)
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


def _validate_datasets(dataset, label='dataset', options=None):
    """Validates Reaction messages and cross-references in a Dataset.

    Args:
        dataset: dataset_pb2.Dataset message.
        label: string label for logging purposes only.
        options: ValidationOptions.

    Returns:
        List of validation error messages.
    """
    errors = []
    # Reaction-level validation.
    num_bad_reactions = 0
    for i, reaction in enumerate(dataset.reactions):
        reaction_output = validate_message(reaction,
                                           raise_on_error=False,
                                           options=options)
        if reaction_output.errors:
            num_bad_reactions += 1
        for error in reaction_output.errors:
            errors.append(error)
            logging.warning('Validation error for %s[%d]: %s', label, i, error)
    logging.info('Validation summary for %s: %d/%d successful (%d failures)',
                 label,
                 len(dataset.reactions) - num_bad_reactions,
                 len(dataset.reactions), num_bad_reactions)
    # Dataset-level validation of cross-references.
    dataset_output = validate_message(dataset,
                                      raise_on_error=False,
                                      recurse=False,
                                      options=options)
    for error in dataset_output.errors:
        errors.append(error)
        logging.warning('Validation error for %s: %s', label, error)

    return errors


def validate_message(message,
                     recurse=True,
                     raise_on_error=True,
                     options=None,
                     trace=None):
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
        options: ValidationOptions.
        trace: Tuple containing a string "stack trace" to track the position of
            the current message relative to the recursion root.

    Returns:
        ValidationOutput.

    Raises:
        ValidationError: If any fields are invalid.
    """
    if trace is None:
        trace = (message.DESCRIPTOR.name,)
    output = ValidationOutput()
    # Recurse through submessages
    if recurse:
        for field, value in message.ListFields():
            if field.type == field.TYPE_MESSAGE:  # need to recurse
                _validate_message(field=field,
                                  value=value,
                                  output=output,
                                  raise_on_error=raise_on_error,
                                  options=options,
                                  trace=trace)

    # Message-specific validation
    if not isinstance(message, tuple(_VALIDATOR_SWITCH.keys())):
        # NOTE(ccoley): I made the conscious decision to raise an error here,
        # rather than assume that the message is valid. If a message does not
        # require any message-level checks (not uncommon), then it should still
        # be listed in the dictionary switch above withpass. This will force
        # us to think about what is necessary if/when new messages are added.
        raise NotImplementedError(f"Don't know how to validate {type(message)}")

    with warnings.catch_warnings(record=True) as tape:
        if isinstance(message, (reaction_pb2.Reaction, dataset_pb2.Dataset)):
            _VALIDATOR_SWITCH[type(message)](message, options=options)
        else:
            _VALIDATOR_SWITCH[type(message)](message)
    stack = '.'.join(trace)
    for warning in tape:
        message = f'{stack}: {warning.message}'
        if issubclass(warning.category, ValidationError):
            if raise_on_error:
                raise ValidationError(message)
            output.errors.append(message)
        else:
            output.warnings.append(message)
    return output


def _validate_message(field, value, output, raise_on_error, options, trace):
    """Validates a single message field and its children.

    Args:
        field: FieldDescriptor instance.
        value: The value of the current message field.
        output: ValidationOutput.
        raise_on_error: If True, raises a ValidationError exception when errors
            are encountered. If False, the user must manually check the return
            value to identify validation errors.
        options: ValidationOptions.
        trace: Tuple containing a string "stack trace" to track the position of
            the current message relative to the recursion root.
    """
    if field.label == field.LABEL_REPEATED:
        if field.message_type.GetOptions().map_entry:  # map
            # value is message
            if field.message_type.fields_by_name['value'].type == \
                    field.TYPE_MESSAGE:
                for key, submessage in value.items():
                    this_trace = trace + (f'{field.name}["{key}"]',)
                    this_output = validate_message(
                        submessage,
                        raise_on_error=raise_on_error,
                        options=options,
                        trace=this_trace)
                    output.extend(this_output)
            else:  # value is a primitive
                pass
        else:  # Just a repeated message
            for index, submessage in enumerate(value):
                this_trace = trace + (f'{field.name}[{index}]',)
                this_output = validate_message(submessage,
                                               raise_on_error=raise_on_error,
                                               options=options,
                                               trace=this_trace)
                output.extend(this_output)
    else:  # no recursion needed
        this_trace = trace + (field.name,)
        this_output = validate_message(value,
                                       raise_on_error=raise_on_error,
                                       options=options,
                                       trace=this_trace)
        output.extend(this_output)


class ValidationError(Warning):
    pass


class ValidationWarning(Warning):
    pass


# pylint: disable=missing-function-docstring
def ensure_float_nonnegative(message, field):
    if getattr(message, field) < 0:
        warnings.warn(
            f'Field {field} of message '
            f'{type(message).DESCRIPTOR.name} must be'
            ' non-negative', ValidationError)


def ensure_float_range(message, field, min_value=-math.inf, max_value=math.inf):
    if (getattr(message, field) < min_value or
            getattr(message, field) > max_value):
        warnings.warn(
            f'Field {field} of message '
            f'{type(message).DESCRIPTOR.name} must be between'
            f' {min_value} and {max_value}', ValidationError)


def check_value_and_units(message):
    """Checks that value/units messages are complete."""
    if not message.HasField('value'):
        warnings.warn(f'{type(message)} requires `value` to be set',
                      ValidationError)
    if message.units == message.UNSPECIFIED:
        warnings.warn(f'{type(message)} requires `units` to be set',
                      ValidationError)


def check_type_and_details(message):
    """Checks that type/details messages are complete."""
    if message.type == message.UNSPECIFIED:
        warnings.warn(f'{type(message)} requires `type` to be set',
                      ValidationError)
    if message.type == message.CUSTOM and not message.details:
        warnings.warn(
            f'{type(message)} has type CUSTOM but details field is empty',
            ValidationError)


def reaction_has_internal_standard(message):
    """Whether any reaction component uses the internal standard role."""
    for reaction_input in message.inputs.values():
        for compound in reaction_input.components:
            if (compound.reaction_role ==
                    compound.ReactionRole.INTERNAL_STANDARD):
                return True
    for workup in message.workup:
        if workup.input:
            for compound in workup.input.components:
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


def get_referenced_reaction_ids(message):
    """Return the set of reaction IDs that are referenced in a Reaction."""
    referenced_ids = set()
    for reaction_input in message.inputs.values():
        for component in reaction_input.components:
            for preparation in component.preparations:
                if preparation.reaction_id:
                    referenced_ids.add(preparation.reaction_id)
        for crude_component in reaction_input.crude_components:
            referenced_ids.add(crude_component.reaction_id)
    return referenced_ids


def validate_dataset(message, options=None):
    # pylint: disable=too-many-branches,too-many-nested-blocks
    if options is None:
        options = ValidationOptions()
    if not message.reactions and not message.reaction_ids:
        warnings.warn('Dataset requires reactions or reaction_ids',
                      ValidationError)
    elif message.reactions and message.reaction_ids:
        warnings.warn('Dataset requires reactions or reaction_ids, not both',
                      ValidationError)
    if message.reaction_ids:
        for reaction_id in message.reaction_ids:
            if not re.fullmatch('^ord-[0-9a-f]{32}$', reaction_id):
                warnings.warn('Reaction ID is malformed', ValidationError)
    if options.validate_ids:
        # The dataset_id is a 32-character uuid4 hex string.
        if not re.fullmatch('^ord_dataset-[0-9a-f]{32}$', message.dataset_id):
            warnings.warn('Dataset ID is malformed', ValidationError)
    # Check cross-references
    dataset_referenced_ids = set()
    dataset_defined_ids = set()
    for reaction in message.reactions:
        if reaction.reaction_id:
            if reaction.reaction_id in dataset_defined_ids:
                warnings.warn(
                    'Multiple Reactions should never have the same '
                    'IDs', ValidationError)
            dataset_defined_ids.add(reaction.reaction_id)
        referenced_ids = get_referenced_reaction_ids(reaction)
        if any(_id == reaction.reaction_id for _id in referenced_ids):
            warnings.warn('A Reaction should not reference its own ID',
                          ValidationError)
        dataset_referenced_ids |= referenced_ids
    if len(dataset_referenced_ids - dataset_defined_ids) > 0:
        warnings.warn(
            'Reactions in the Dataset refer to undefined '
            f'reaction_ids {dataset_referenced_ids - dataset_defined_ids}',
            ValidationError)


def validate_dataset_example(message):
    if not message.description:
        warnings.warn('DatasetExample.description is required', ValidationError)
    if not message.url:
        warnings.warn('DatasetExample.url is required', ValidationError)
    if not message.HasField('created'):
        warnings.warn('DatasetExample.created is required', ValidationError)


def validate_reaction(message, options=None):
    if options is None:
        options = ValidationOptions()
    if len(message.inputs) == 0:
        warnings.warn('Reactions should have at least 1 reaction input',
                      ValidationError)
    if len(message.outcomes) == 0:
        warnings.warn('Reactions should have at least 1 reaction outcome',
                      ValidationError)
    if (reaction_needs_internal_standard(message) and
            not reaction_has_internal_standard(message)):
        warnings.warn(
            'Reaction analysis uses an internal standard, but no '
            'component (as reaction input or workup) uses the '
            'reaction role INTERNAL_STANDARD', ValidationError)
    if (any(outcome.HasField('conversion') for outcome in message.outcomes) and
            not reaction_has_limiting_component(message)):
        warnings.warn(
            'If reaction conversion is specified, at least one '
            'reaction input component must be labeled is_limiting',
            ValidationError)
    if options.validate_ids:
        # The reaction_id suffix is a 32-character uuid4 hex string.
        if not re.fullmatch('^ord-[0-9a-f]{32}$', message.reaction_id):
            warnings.warn('Reaction ID is malformed', ValidationError)
    if options.require_provenance:
        if not message.HasField('provenance'):
            warnings.warn('Reaction requires provenance', ValidationError)


def validate_reaction_identifier(message):
    check_type_and_details(message)
    if not message.value:
        warnings.warn('value must be set', ValidationError)


def validate_reaction_input(message):
    if len(message.components) + len(message.crude_components) == 0:
        warnings.warn('Reaction inputs must have at least one component',
                      ValidationError)
    for component in message.components:
        if not component.WhichOneof('amount'):
            warnings.warn('Reaction input components require an amount',
                          ValidationError)


def validate_addition_device(message):
    check_type_and_details(message)


def validate_addition_speed(message):
    del message  # Unused.


def validate_crude_component(message):
    if not message.reaction_id:
        warnings.warn('CrudeComponents must specify a reaction_id',
                      ValidationError)
    if message.has_derived_amount and message.HasField('amount'):
        warnings.warn(
            'CrudeComponents with derived amounts cannot have their'
            ' mass or volume specified explicitly', ValidationError)
    if ((not message.HasField('has_derived_amount') or
         not message.has_derived_amount) and not message.HasField('amount')):
        warnings.warn(
            'Crude components should either have a derived amount or'
            ' a specified mass or volume', ValidationError)


def validate_compound(message):
    if len(message.identifiers) == 0:
        warnings.warn('Compounds must have at least one identifier',
                      ValidationError)
    if all(identifier.type == identifier.NAME
           for identifier in message.identifiers):
        warnings.warn(
            'Compounds should have more specific identifiers than '
            'NAME whenever possible', ValidationWarning)
    try:
        message_helpers.check_compound_identifiers(message)
    except ValueError as error:
        warnings.warn(str(error), ValidationWarning)


def validate_compound_feature(message):
    if not message.name:
        warnings.warn('Compound features must have names', ValidationError)


def validate_compound_preparation(message):
    check_type_and_details(message)
    if message.reaction_id and message.type != message.SYNTHESIZED:
        warnings.warn(
            'Reaction IDs should only be specified in compound'
            ' preparations when SYNTHESIZED', ValidationError)


def validate_compound_identifier(message):
    check_type_and_details(message)
    if not message.value:
        warnings.warn('value must be set', ValidationError)
    if Chem and message.type == message.SMILES:
        mol = Chem.MolFromSmiles(message.value)
        if mol is None:
            warnings.warn(
                f'RDKit {RDKIT_VERSION} could not validate'
                f' SMILES identifier {message.value}', ValidationError)
    elif message.type == message.INCHI:
        mol = Chem.MolFromInchi(message.value)
        if mol is None:
            warnings.warn(
                f'RDKit {RDKIT_VERSION} could not validate'
                f' InChI identifier {message.value}', ValidationError)
    elif message.type == message.MOLBLOCK:
        mol = Chem.MolFromMolBlock(message.value)
        if mol is None:
            warnings.warn(
                f'RDKit {RDKIT_VERSION} could not validate'
                ' MolBlock identifier', ValidationError)


def validate_vessel(message):
    del message  # Unused.


def validate_vessel_type(message):
    check_type_and_details(message)


def validate_vessel_material(message):
    check_type_and_details(message)


def validate_vessel_attachment(message):
    check_type_and_details(message)


def validate_vessel_preparation(message):
    check_type_and_details(message)


def validate_reaction_setup(message):
    del message  # Unused.


def validate_reaction_environment(message):
    check_type_and_details(message)


def validate_reaction_conditions(message):
    if message.conditions_are_dynamic and not message.details:
        warnings.warn(
            'Reaction conditions are dynamic, but no details'
            ' provided to explain how procedure deviates from'
            ' normal single-step reaction conditions.', ValidationError)
    if message.details and not message.conditions_are_dynamic:
        warnings.warn(
            'Reaction condition details provided but field '
            'conditions_are_dynamic is False. If the conditions '
            'cannot be fully captured by the schema, set to True.',
            ValidationWarning)


def validate_temperature_conditions(message):
    if not message.setpoint.value:
        warnings.warn(
            'Temperature setpoints should be specified; even if '
            'using ambient conditions, estimate room temperature and '
            'the precision of your estimate.', ValidationWarning)


def validate_temperature_control(message):
    check_type_and_details(message)


def validate_temperature_measurement(message):
    check_type_and_details(message)


def validate_pressure_conditions(message):
    del message


def validate_pressure_control(message):
    check_type_and_details(message)


def validate_atmosphere(message):
    check_type_and_details(message)


def validate_pressure_measurement(message):
    check_type_and_details(message)


def validate_stirring_conditions(message):
    del message  # Unused.


def validate_stirring_method(message):
    check_type_and_details(message)


def validate_stirring_rate(message):
    ensure_float_nonnegative(message, 'rpm')


def validate_illumination_conditions(message):
    del message  # Unused.


def validate_illumination_type(message):
    check_type_and_details(message)


def validate_electrochemistry_conditions(message):
    del message  # Unused.


def validate_electrochemistry_type(message):
    check_type_and_details(message)


def validate_electrochemistry_cell(message):
    check_type_and_details(message)


def validate_electrochemistry_measurement(message):
    del message  # Unused.


def validate_flow_conditions(message):
    del message  # Unused.


def validate_flow_type(message):
    check_type_and_details(message)


def validate_tubing(message):
    check_type_and_details(message)


def validate_reaction_notes(message):
    del message  # Unused.


def validate_reaction_observation(message):
    del message  # Unused.


def validate_reaction_workup(message):
    check_type_and_details(message)
    if (message.type == reaction_pb2.ReactionWorkup.WAIT and
            not message.duration.value):
        warnings.warn('"WAIT" workup steps require a defined duration',
                      ValidationError)
    if (message.type == reaction_pb2.ReactionWorkup.TEMPERATURE and
            not message.HasField('temperature')):
        warnings.warn(
            '"TEMPERATURE" workup steps require defined '
            'temperature conditions', ValidationError)
    if (message.type in (reaction_pb2.ReactionWorkup.EXTRACTION,
                         reaction_pb2.ReactionWorkup.FILTRATION) and
            not message.keep_phase):
        warnings.warn(
            'Workup step EXTRACTION or FILTRATION missing '
            'required field keep_phase', ValidationError)
    if (message.type in (reaction_pb2.ReactionWorkup.ADDITION,
                         reaction_pb2.ReactionWorkup.WASH,
                         reaction_pb2.ReactionWorkup.DRY_WITH_MATERIAL,
                         reaction_pb2.ReactionWorkup.SCAVENGING,
                         reaction_pb2.ReactionWorkup.DISSOLUTION,
                         reaction_pb2.ReactionWorkup.PH_ADJUST) and
            not message.input.components):
        warnings.warn('Workup step missing required inputs definition',
                      ValidationError)
    if (message.type == reaction_pb2.ReactionWorkup.STIRRING and
            not message.stirring):
        warnings.warn('Stirring workup step missing stirring definition',
                      ValidationError)
    if (message.type == reaction_pb2.ReactionWorkup.PH_ADJUST and
            not message.HasField('target_ph')):
        warnings.warn('pH adjustment workup missing target pH', ValidationError)
    if (message.type == reaction_pb2.ReactionWorkup.ALIQUOT and
            not message.WhichOneof('amount')):
        warnings.warn('Aliquot workup step missing volume/mass amount',
                      ValidationError)


def validate_reaction_outcome(message):
    # pylint: disable=singleton-comparison
    # Can only have one desired product
    if sum(product.is_desired_product for product in message.products) > 1:
        warnings.warn('Cannot have more than one desired product!',
                      ValidationError)
    # Check key values for product analyses
    # NOTE(ccoley): Could use any(), but using expanded loops for clarity
    analysis_keys = list(message.analyses.keys())
    for product in message.products:
        for field in [
                'analysis_identity', 'analysis_yield', 'analysis_purity',
                'analysis_selectivity'
        ]:
            for key in getattr(product, field):
                if key not in analysis_keys:
                    warnings.warn(
                        f'Undefined analysis key {key} '
                        'in ReactionProduct', ValidationError)
    # TODO(ccoley): While we do not currently check whether the parent Reaction
    # is *actually* used in a multistep reaction within a Dataset (i.e., in a
    # CrudeComponent); this is an additional check that could be added to the
    # submission pipeline.
    if not message.products and not message.HasField('conversion'):
        warnings.warn(
            'No products or conversion are specified for reaction;'
            ' this is permissible only for multistep reactions',
            ValidationWarning)


def validate_reaction_product(message):
    # pylint: disable=too-many-boolean-expressions
    if (message.compound.HasField('volume_includes_solutes') or
            message.compound.HasField('is_limiting') or
            message.compound.preparations or message.compound.vendor_source or
            message.compound.vendor_id or message.compound.vendor_lot):
        warnings.warn(
            'Compounds defined as reaction products should not have'
            ' any inapplicable fields defined.', ValidationError)


def validate_texture(message):
    check_type_and_details(message)


def validate_selectivity(message):
    ensure_float_nonnegative(message, 'precision')
    if message.type == message.EE:
        ensure_float_range(message, 'value', 0, 100)
        if 0 < message.value < 1:
            warnings.warn(
                'EE selectivity values are 0-100, not fractions '
                f'({message.value} used)', ValidationWarning)
    elif message.type in [message.ER, message.DR, message.EZ, message.ZE]:
        ensure_float_nonnegative(message, 'value')
    check_type_and_details(message)


def validate_date_time(message):
    if message.value:
        try:
            parser.parse(message.value).ctime()
        except parser.ParserError:
            warnings.warn(f'Could not parse DateTime string {message.value}',
                          ValidationError)


def validate_reaction_analysis(message):
    # TODO(ccoley): Will be lots to expand here if we add structured data.
    check_type_and_details(message)


def validate_reaction_provenance(message):
    # Prepare datetimes
    if not message.HasField('record_created'):
        warnings.warn('Reactions must have record_created defined',
                      ValidationError)
    experiment_start = None
    record_created = None
    record_modified = None
    try:
        if message.experiment_start.value:
            experiment_start = parser.parse(message.experiment_start.value)
        if message.record_created.time.value:
            record_created = parser.parse(message.record_created.time.value)
        for record in message.record_modified:
            # Use the last record as the most recent modification time.
            record_modified = parser.parse(record.time.value)
    except parser.ParserError:
        warnings.warn('Failed to parse DateTime string(s)')
    # Check signs of time differences
    if experiment_start and record_created:
        if (record_created - experiment_start).total_seconds() < 0:
            warnings.warn('Record creation time should be after experiment',
                          ValidationError)
    if record_modified and record_created:
        if (record_modified - record_created).total_seconds() < 0:
            warnings.warn('Record modified time should be after creation',
                          ValidationError)
    # TODO(ccoley) could check if publication_url is valid, etc.


def validate_record_event(message):
    if not message.time.value:
        warnings.warn('RecordEvent must have `time` specified', ValidationError)
    person = message.person
    if not (person.username or person.name or person.orcid):
        warnings.warn(
            'Person must have at least one of '
            '`username`, `name`, or `orcid` specified', ValidationError)
    if not person.email:
        warnings.warn('Person must have `email` specified', ValidationError)


def validate_person(message):
    # NOTE(ccoley): final character is checksum, but ignoring that for now
    if message.orcid:
        if not re.match('[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]',
                        message.orcid):
            warnings.warn('Invalid ORCID: Enter as 0000-0000-0000-0000',
                          ValidationError)
    if message.email:
        # Based on https://www.regular-expressions.info/email.html.
        # Added optional "[bot]" suffix to the username for GitHub actions.
        if not re.fullmatch(
                r'[a-zA-Z0-9._+-]+(?:\[bot\])?@[a-zA-Z0-9.-]+\.[a-z]{2,}',
                message.email):
            warnings.warn(f'Invalid email address: {message.email}',
                          ValidationError)


def validate_time(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_mass(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_moles(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_volume(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_concentration(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_pressure(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_temperature(message):
    check_value_and_units(message)
    if message.units == message.CELSIUS:
        ensure_float_range(message, 'value', min_value=-273.15)
    elif message.units == message.FAHRENHEIT:
        ensure_float_range(message, 'value', min_value=-459)
    elif message.units == message.KELVIN:
        ensure_float_range(message, 'value', min_value=0)
    ensure_float_nonnegative(message, 'precision')


def validate_current(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_voltage(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_length(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_wavelength(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_flow_rate(message):
    check_value_and_units(message)
    ensure_float_nonnegative(message, 'value')
    ensure_float_nonnegative(message, 'precision')


def validate_percentage(message):
    if 0 < message.value < 1:
        warnings.warn(
            'Percentage values are 0-100, not fractions '
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
        warnings.warn('Data format is required for bytes_data', ValidationError)


# pylint: enable=missing-function-docstring

_VALIDATOR_SWITCH = {
    dataset_pb2.Dataset:
        validate_dataset,
    dataset_pb2.DatasetExample:
        validate_dataset_example,
    reaction_pb2.Reaction:
        validate_reaction,
    # Basics
    reaction_pb2.ReactionIdentifier:
        validate_reaction_identifier,
    reaction_pb2.ReactionInput:
        validate_reaction_input,
    reaction_pb2.ReactionInput.AdditionDevice:
        validate_addition_device,
    reaction_pb2.ReactionInput.AdditionSpeed:
        validate_addition_speed,
    # Compounds
    reaction_pb2.CrudeComponent:
        validate_crude_component,
    reaction_pb2.Compound:
        validate_compound,
    reaction_pb2.Compound.Feature:
        validate_compound_feature,
    reaction_pb2.CompoundPreparation:
        validate_compound_preparation,
    reaction_pb2.CompoundIdentifier:
        validate_compound_identifier,
    # Setup
    reaction_pb2.Vessel:
        validate_vessel,
    reaction_pb2.VesselType:
        validate_vessel_type,
    reaction_pb2.VesselMaterial:
        validate_vessel_material,
    reaction_pb2.VesselAttachment:
        validate_vessel_attachment,
    reaction_pb2.VesselPreparation:
        validate_vessel_preparation,
    reaction_pb2.ReactionSetup:
        validate_reaction_setup,
    reaction_pb2.ReactionSetup.ReactionEnvironment:
        validate_reaction_environment,
    # Conditions
    reaction_pb2.ReactionConditions:
        validate_reaction_conditions,
    reaction_pb2.TemperatureConditions:
        validate_temperature_conditions,
    reaction_pb2.TemperatureConditions.TemperatureControl:
        validate_temperature_control,
    reaction_pb2.TemperatureConditions.Measurement:
        validate_temperature_measurement,
    reaction_pb2.PressureConditions:
        validate_pressure_conditions,
    reaction_pb2.PressureConditions.PressureControl:
        validate_pressure_control,
    reaction_pb2.PressureConditions.Atmosphere:
        validate_pressure_control,
    reaction_pb2.PressureConditions.Measurement:
        validate_pressure_measurement,
    reaction_pb2.StirringConditions:
        validate_stirring_conditions,
    reaction_pb2.StirringConditions.StirringMethod:
        validate_stirring_method,
    reaction_pb2.StirringConditions.StirringRate:
        validate_stirring_rate,
    reaction_pb2.IlluminationConditions:
        validate_illumination_conditions,
    reaction_pb2.IlluminationConditions.IlluminationType:
        validate_illumination_type,
    reaction_pb2.ElectrochemistryConditions:
        validate_electrochemistry_conditions,
    reaction_pb2.ElectrochemistryConditions.ElectrochemistryType:
        validate_electrochemistry_type,
    reaction_pb2.ElectrochemistryConditions.ElectrochemistryCell:
        validate_electrochemistry_cell,
    reaction_pb2.ElectrochemistryConditions.Measurement:
        validate_electrochemistry_measurement,
    reaction_pb2.FlowConditions:
        validate_flow_conditions,
    reaction_pb2.FlowConditions.FlowType:
        validate_flow_type,
    reaction_pb2.FlowConditions.Tubing:
        validate_tubing,
    # Annotations
    reaction_pb2.ReactionNotes:
        validate_reaction_notes,
    reaction_pb2.ReactionObservation:
        validate_reaction_observation,
    # Outcome
    reaction_pb2.ReactionWorkup:
        validate_reaction_workup,
    reaction_pb2.ReactionOutcome:
        validate_reaction_outcome,
    reaction_pb2.ReactionProduct:
        validate_reaction_product,
    reaction_pb2.ReactionProduct.Texture:
        validate_texture,
    reaction_pb2.Selectivity:
        validate_selectivity,
    reaction_pb2.DateTime:
        validate_date_time,
    reaction_pb2.ReactionAnalysis:
        validate_reaction_analysis,
    # Metadata
    reaction_pb2.ReactionProvenance:
        validate_reaction_provenance,
    reaction_pb2.RecordEvent:
        validate_record_event,
    reaction_pb2.Person:
        validate_person,
    # Units
    reaction_pb2.Time:
        validate_time,
    reaction_pb2.Mass:
        validate_mass,
    reaction_pb2.Moles:
        validate_moles,
    reaction_pb2.Volume:
        validate_volume,
    reaction_pb2.Concentration:
        validate_concentration,
    reaction_pb2.Pressure:
        validate_pressure,
    reaction_pb2.Temperature:
        validate_temperature,
    reaction_pb2.Current:
        validate_current,
    reaction_pb2.Voltage:
        validate_voltage,
    reaction_pb2.Length:
        validate_length,
    reaction_pb2.Wavelength:
        validate_wavelength,
    reaction_pb2.FlowRate:
        validate_flow_rate,
    reaction_pb2.Percentage:
        validate_percentage,
    reaction_pb2.Data:
        validate_data,
}
