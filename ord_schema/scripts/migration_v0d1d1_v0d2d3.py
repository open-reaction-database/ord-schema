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
"""Converts dataset protos generated using the ORD schema v0.1.1 to
files using the ORD schema v0.2.2.

Example usage:
$ python migration_v0d1d1_v0d2d2.py --input="old/*.pbtxt" --output="new/"
"""

import glob
import os
import datetime
import warnings
from google.protobuf import text_format

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.proto import dataset_old_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern for Dataset protos.')
flags.DEFINE_string('output', None, 'Output folder for new Dataset protos.')


# pylint: disable=invalid-name
def migrate_MergeFrom(repeated_new, repeated_old):
    """Mimics the "extend" method of a repeated message field when the field
    types are defined equivalently but are different classes. Roundtrips each
    message through a serialized pbtxt format. Modified repeated_new in-place.

    Args:
        repeated_new: A RepeatedCompositeContainer.
        repeated_old: A RepeatedCompositeContainer containing messages that
            can be serialized to text and parsed as a message within the
            repeated_new RepeatedCompositeContainer.
    """
    for message_old in repeated_old:
        migrate_CopyFrom(repeated_new.add(), message_old)


def migrate_CopyFrom(message_new, message_old):
    """Mimics the "CopyFrom" method of a message field when the field
    types are defined equivalently but are different classes. Roundtrips each
    message through a serialized pbtxt format. Modifies message_new in-place.

    Args:
        message_new: A message to write the value of message_old into.
        message_old: A message that can be serialized to text and parsed as a
            message with the same type as message_new.
    """
    serialized = text_format.MessageToString(message_old)
    text_format.Parse(serialized, message_new)


def migrate_MergeFrom_multiple(message_new, message_old, fields):
    """Iterates over submessage names in "fields" to repeated call
    migrate_MergeFrom.

    Args:
        message_new: A message containing at least one repeated field.
        message_old: A message containing at least one repeated field.
        fields: A list of field names corresponding to repeated fields that
            should be merged from the old message to the new one.
    """
    for field in fields:
        migrate_MergeFrom(getattr(message_new, field),
                          getattr(message_old, field))


def migrate_CopyFrom_multiple(message_new, message_old, fields):
    """Iterates over submessage names in "fields" to repeated call
    migrate_CopyFrom.

    Args:
        message_new: A message containing at least one submessage field.
        message_old: A message containing at least one submessage field.
        fields: A list of field names corresponding to submessage fields that
            should be copied from the old message to the new one.
    """
    for field in fields:
        if message_old.HasField(field):
            migrate_CopyFrom(getattr(message_new, field),
                             getattr(message_old, field))


def migrate_assign_multiple(message_new, message_old, fields):
    """Iterates over field names in "fields" to assign the values from the
     old message to the new one.

    Args:
        message_new: A message.
        message_old: A message.
        fields: A list of field names corresponding to fields whose values
            should be copied from the old message to the new one.
    """
    for field in fields:
        setattr(message_new, field, getattr(message_old, field))


def migrate_input(input_new, input_old):
    """Copies the contents from an old reaction input to the new reaction
    input format.

    Args:
        input_new: A reaction_pb2.ReactionInput message.
        input_old: A reaction_old_pb2.ReactionInput message.
    """
    if input_old.HasField('crude_components'):
        raise NotImplementedError('Crude component migration not defined')
    input_new.addition_order = input_old.addition_order
    migrate_CopyFrom_multiple(
        input_new, input_old,
        ('addition_time', 'addition_speed', 'addition_duration', 'flow_rate',
         'addition_device', 'addition_temperature'))
    # Changed input fields - components
    for component_old in input_old.components:
        component = input_new.components.add()
        migrate_MergeFrom_multiple(component, component_old,
                                   ('identifiers', 'preparations'))
        if len(component_old.features) > 0:
            raise NotImplementedError
        if component_old.HasField('is_limiting'):
            component.is_limiting = component_old.is_limiting
        # Vendor information is now in a submessage
        if (component_old.vendor_source or component_old.vendor_id or
                component_old.vendor_lot):
            component.source.vendor = component_old.vendor_source
            component.source.id = component_old.vendor_id
            component.source.lot = component_old.vendor_lot
        # Reaction role enum is relocated and renumbered
        role_name = component_old.ReactionRole.ReactionRoleType.DESCRIPTOR.\
            values_by_number[component_old.reaction_role].name
        component.reaction_role = reaction_pb2.ReactionRole.ReactionRoleType.\
            DESCRIPTOR.values_by_name[role_name].number
        # Amounts are nwo in a submessage
        kind = component_old.WhichOneof('amount')
        if kind:
            migrate_CopyFrom(getattr(component.amount, kind),
                             getattr(component_old, kind))
        if component_old.HasField('volume_includes_solutes'):
            component.amount.volume_includes_solutes = (
                component_old.volume_includes_solutes)


def migrate_workup(workup_new, workup_old):
    """Copies the contents from an old workup to the new workup format.

    Args:
        workup_new: A reaction_pb2.ReactionWorkup message.
        workup_old: A reaction_old_pb2.ReactionWorkup message.
    """
    type_name = workup_old.WorkupType.DESCRIPTOR.values_by_number[
        workup_old.type].name
    workup_new.type = workup_new.WorkupType.DESCRIPTOR.values_by_name[
        type_name].number
    if workup_old.HasField('input'):
        migrate_input(workup_new.input, workup_old.input)
    migrate_assign_multiple(workup_new, workup_old, ('details', 'keep_phase'))
    migrate_CopyFrom_multiple(workup_new, workup_old,
                              ('duration', 'temperature', 'stirring'))
    if workup_old.HasField('target_ph'):
        workup_new.target_ph = workup_old.target_ph
    if workup_old.HasField('is_automated'):
        workup_new.is_automated = workup_old.is_automated
    if workup_old.HasField('mass'):
        migrate_CopyFrom(workup_new.amount.mass, workup_old.mass)
    elif workup_old.HasField('volume'):
        migrate_CopyFrom(workup_new.amount.volume, workup_old.volume)


# pylint: disable=too-many-locals,too-many-nested-blocks
# pylint: disable=too-many-branches,too-many-statements
def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input, recursive=True)
    logging.info('Found %d datasets', len(filenames))
    for filename in filenames:
        logging.info('Processing %s', filename)
        dataset_old = message_helpers.load_message(filename,
                                                   dataset_old_pb2.Dataset)
        dataset = dataset_pb2.Dataset()

        # Copy unchanged dataset-level information
        migrate_assign_multiple(dataset, dataset_old,
                                ('name', 'description', 'dataset_id'))
        migrate_MergeFrom(dataset.reaction_ids, dataset_old.reaction_ids)

        for i, reaction_old in enumerate(dataset_old.reactions):
            logging.info(f'  reaction {i} of {len(dataset_old.reactions)}')
            reaction = dataset.reactions.add()
            # Copy over unchanged fields - all but inputs, workups, outputs
            migrate_MergeFrom_multiple(reaction, reaction_old,
                                       ('identifiers', 'observations'))
            migrate_CopyFrom_multiple(
                reaction, reaction_old,
                ('setup', 'notes', 'conditions', 'provenance'))
            reaction.reaction_id = reaction_old.reaction_id

            # Copy workups field, renamed from workup. Need to handle
            # the change in the "input" field.
            for workup_old in reaction_old.workup:
                workup = reaction.workups.add()
                migrate_workup(workup, workup_old)

            # Copy changed fields in input
            for input_key in reaction_old.inputs:
                migrate_input(reaction.inputs[input_key],
                              reaction_old.inputs[input_key])

            # Copy changed fields in output
            for outcome_old in reaction_old.outcomes:
                outcome = reaction.outcomes.add()
                migrate_CopyFrom_multiple(outcome, outcome_old,
                                          ('reaction_time', 'conversion'))
                for product_old in outcome_old.products:
                    product = outcome.products.add()
                    if product_old.HasField('is_desired_product'):
                        product.is_desired_product = (
                            product_old.is_desired_product)
                    migrate_MergeFrom(product.identifiers,
                                      product_old.compound.identifiers)
                    product.isolated_color = product_old.isolated_color
                    migrate_CopyFrom(product.texture, product_old.texture)
                    role_name = product_old.compound.ReactionRole.\
                        ReactionRoleType.DESCRIPTOR.values_by_number[
                            product_old.compound.reaction_role].name
                    product.reaction_role = reaction_pb2.ReactionRole.\
                        ReactionRoleType.DESCRIPTOR.values_by_name[
                            role_name].number

                    for analysis_key in product_old.analysis_identity:
                        analysis_old = outcome_old.analyses[analysis_key]
                        measurement = product.measurements.add(
                            analysis_key=analysis_key, type='IDENTITY')
                        if analysis_old.HasField('uses_internal_standard'):
                            measurement.uses_internal_standard = (
                                analysis_old.uses_internal_standard)
                        if analysis_old.HasField('uses_authentic_standard'):
                            measurement.uses_authentic_standard = (
                                analysis_old.uses_authentic_standard)

                    for analysis_key in product_old.analysis_yield:
                        analysis_old = outcome_old.analyses[analysis_key]
                        measurement = product.measurements.add(
                            analysis_key=analysis_key, type='YIELD')
                        migrate_CopyFrom(measurement.percentage,
                                         product_old.compound_yield)
                        if analysis_old.HasField('uses_internal_standard'):
                            measurement.uses_internal_standard = (
                                analysis_old.uses_internal_standard)
                        if analysis_old.HasField('uses_authentic_standard'):
                            measurement.uses_authentic_standard = (
                                analysis_old.uses_authentic_standard)

                    for analysis_key in product_old.analysis_purity:
                        analysis_old = outcome_old.analyses[analysis_key]
                        measurement = product.measurements.add(
                            analysis_key=analysis_key, type='PURITY')
                        migrate_CopyFrom(measurement.percentage,
                                         product_old.purity)
                        if analysis_old.HasField('uses_internal_standard'):
                            measurement.uses_internal_standard = (
                                analysis_old.uses_internal_standard)
                        if analysis_old.HasField('uses_authentic_standard'):
                            measurement.uses_authentic_standard = (
                                analysis_old.uses_authentic_standard)

                    for analysis_key in product_old.analysis_selectivity:
                        analysis_old = outcome_old.analyses[analysis_key]
                        measurement = product.measurements.add(
                            analysis_key=analysis_key, type='SELECTIVITY')
                        migrate_assign_multiple(measurement.selectivity,
                                                product_old.selectivity,
                                                ('type', 'details'))
                        if product_old.selectivity.HasField('value'):
                            measurement.float_value.value = (
                                product_old.selectivity.value)
                            measurement.float_value.precision = (
                                product_old.selectivity.precision)
                        if analysis_old.HasField('uses_internal_standard'):
                            measurement.uses_internal_standard = (
                                analysis_old.uses_internal_standard)
                        if analysis_old.HasField('uses_authentic_standard'):
                            measurement.uses_authentic_standard = (
                                analysis_old.uses_authentic_standard)

                    kind = product_old.compound.WhichOneof('amount')
                    if kind:
                        # Look for analysis_key of type weight to assume
                        # is the source of the amount
                        for analysis_key in outcome_old.analyses:
                            analysis_old = (outcome_old.analyses[analysis_key])
                            if analysis_old.type == analysis_old.WEIGHT:
                                measurement = product.measurements.add(
                                    analysis_key=analysis_key, type='AMOUNT')
                                migrate_CopyFrom(
                                    getattr(measurement.amount, kind),
                                    getattr(product_old.compound, kind))
                                warnings.warn(
                                    'Inferring that analysis_key '
                                    f'{analysis_key} was used for the amount')
                                break
                        else:
                            raise ValueError(
                                'Could not find a WEIGHT'
                                ' analysis that could have been used to'
                                ' calculate the product amount')

                    if len(product_old.compound.features) > 0:
                        raise NotImplementedError

                for analysis_key in outcome_old.analyses:
                    analysis_old = outcome_old.analyses[analysis_key]
                    analysis = outcome.analyses[analysis_key]
                    # Enums for analysis types didn't change
                    migrate_assign_multiple(analysis, analysis_old,
                                            ('type', 'chmo_id', 'details',
                                             'instrument_manufacturer'))
                    migrate_CopyFrom_multiple(analysis, analysis_old,
                                              ('instrument_last_calibrated',))

                    for data_key in analysis_old.processed_data:
                        migrate_CopyFrom(analysis.data[data_key],
                                         analysis_old.processed_data[data_key])
                    for data_key in analysis_old.raw_data:
                        if data_key in analysis.data:
                            data_key += '_raw'
                        migrate_CopyFrom(analysis.data[data_key],
                                         analysis_old.raw_data[data_key])

            # Add record_modified event
            reaction.provenance.record_modified.add(
                time=dict(value=str(datetime.datetime.now())),
                person=dict(username='connorcoley', email='ccoley@mit.edu'),
                details='Automatic migration from schema v0.1.1 to v0.2.2')

        # Save
        filename_new = os.path.join(FLAGS.output, os.path.basename(filename))
        logging.info('Saving new %s', filename_new)
        message_helpers.write_message(dataset, filename_new)


if __name__ == '__main__':
    flags.mark_flag_as_required('input')
    app.run(main)
