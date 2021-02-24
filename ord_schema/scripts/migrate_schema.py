# Copyright 2021 Open Reaction Database Project Authors
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
"""Migrates from dataset_old.proto (v0.2.15) to dataset.proto (v0.3.0).

See https://github.com/open-reaction-database/ord-schema/pull/557.
"""

import glob
import os
import datetime
from typing import Iterable

from absl import app
from absl import flags
from absl import logging

from google.protobuf import text_format  # pytype: disable=import-error

from ord_schema import message_helpers
from ord_schema.proto import dataset_old_pb2
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_old_pb2
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input_pattern', None, 'Input pattern for Dataset protos.')
flags.DEFINE_string('output_dir', None, 'Output folder for new Dataset protos.')


def migrate_message(old, new):
    """Migrates an unchanged message."""
    text_format.Parse(text_format.MessageToString(old), new)


def migrate_messages(old, new, fields: Iterable[str]):
    """Migrates a set of unchanged messages."""
    for field in fields:
        if old.HasField(field):
            migrate_message(getattr(old, field), getattr(new, field))


def migrate_repeated_message(old, new):
    """Migrates an unchanged repeated message."""
    for message in old:
        migrate_message(message, new.add())


def migrate_repeated_messages(old, new, fields: Iterable[str]):
    """Migrates a set of unchanged repeated messages."""
    for field in fields:
        migrate_repeated_message(getattr(old, field), getattr(new, field))


def migrate_scalars(old, new, fields: Iterable[str]):
    """Migrates scalar values with unchanged names."""
    for field in fields:
        setattr(new, field, getattr(old, field))


def migrate_optional_scalars(old, new, fields: Iterable[str]):
    """Migrates optional scalar values with unchanged names."""
    for field in fields:
        if old.HasField(field):
            setattr(new, field, getattr(old, field))


def migrate_vessel(old: reaction_old_pb2.Vessel, new: reaction_pb2.Vessel):
    """Migrates a Vessel message."""
    new.type = old.type.type
    new.details = old.type.details
    migrate_messages(old, new, ('material', 'volume'))
    migrate_repeated_messages(old, new, ('preparations', 'attachments'))


def migrate_stirring_conditions(old: reaction_old_pb2.StirringConditions,
                                new: reaction_pb2.StirringConditions):
    """Migrates a StirringConditions message."""
    new.type = old.method.type
    new.details = old.method.details
    migrate_messages(old, new, ('rate',))


def migrate_illumination_conditions(
        old: reaction_old_pb2.IlluminationConditions,
        new: reaction_pb2.IlluminationConditions):
    """Migrates an IlluminationConditions messsage."""
    new.type = old.type.type
    new.details = old.type.details
    migrate_messages(old, new, ('peak_wavelength', 'distance_to_vessel'))
    migrate_scalars(old, new, ('color',))


def migrate_electrochemistry_conditions(
        old: reaction_old_pb2.ElectrochemistryConditions,
        new: reaction_pb2.ElectrochemistryConditions):
    """Migrates an ElectrochemistryConditions messsage."""
    new.type = old.electrochemistry_type.type
    new.details = old.electrochemistry_type.details
    migrate_scalars(old, new, ('anode_material', 'cathode_material'))
    migrate_messages(old, new,
                     ('current', 'voltage', 'electrode_separation', 'cell'))
    migrate_repeated_messages(old, new, ('measurements',))


def migrate_flow_conditions(old: reaction_old_pb2.FlowConditions,
                            new: reaction_pb2.FlowConditions):
    """Migrates a FlowConditions messsage."""
    new.type = old.flow_type.type
    new.details = old.flow_type.details
    migrate_scalars(old, new, ('pump_type',))
    migrate_messages(old, new, ('tubing',))


def migrate_setup(old: reaction_old_pb2.ReactionSetup,
                  new: reaction_pb2.ReactionSetup):
    """Migrates a ReactionSetup messsage."""
    if old.HasField('vessel'):
        migrate_vessel(old.vessel, new.vessel)
    migrate_optional_scalars(old, new, ('is_automated',))
    migrate_scalars(old, new, ('automation_platform',))
    for key, value in old.automation_code.items():
        migrate_message(value, new.automation_code[key])
    migrate_messages(old, new, ('environment',))


def migrate_conditions(old: reaction_old_pb2.ReactionConditions,
                       new: reaction_pb2.ReactionConditions):
    """Migrates a ReactionConditions messsage."""
    migrate_messages(old, new, ('temperature', 'pressure'))
    if old.HasField('stirring'):
        migrate_stirring_conditions(old.stirring, new.stirring)
    if old.HasField('illumination'):
        migrate_illumination_conditions(old.illumination, new.illumination)
    if old.HasField('electrochemistry'):
        migrate_electrochemistry_conditions(old.electrochemistry,
                                            new.electrochemistry)
    if old.HasField('flow'):
        migrate_flow_conditions(old.flow, new.flow)
    migrate_optional_scalars(old, new, ('reflux', 'conditions_are_dynamic'))
    if old.HasField('pH'):
        new.ph = old.pH  # Lowercase now.
    migrate_scalars(old, new, ('details',))


def migrate_workup(old: reaction_old_pb2.ReactionWorkup,
                   new: reaction_pb2.ReactionWorkup):
    """Migrates a ReactionWorkup message."""
    migrate_messages(old, new, ('duration', 'input', 'amount', 'temperature'))
    migrate_scalars(old, new, ('type', 'details', 'keep_phase'))
    migrate_optional_scalars(old, new, ('target_ph', 'is_automated'))
    if old.HasField('stirring'):
        migrate_stirring_conditions(old.stirring, new.stirring)


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input_pattern)
    logging.info('Found %d datasets', len(filenames))
    for filename in filenames:
        logging.info(filename)
        old_dataset = message_helpers.load_message(filename,
                                                   dataset_old_pb2.Dataset)
        new_dataset = dataset_pb2.Dataset()

        # Copy unchanged dataset-level information.
        migrate_scalars(old_dataset, new_dataset,
                        ('name', 'description', 'dataset_id'))
        new_dataset.reaction_ids.extend(old_dataset.reaction_ids)

        num_changed = 0
        for old in old_dataset.reactions:
            new = new_dataset.reactions.add()
            # Copy over unchanged fields.
            migrate_repeated_messages(
                old, new, ('identifiers', 'observations', 'outcomes'))
            for key, value in old.inputs.items():
                migrate_message(value, new.inputs[key])
            migrate_messages(old, new, ('notes', 'provenance'))
            migrate_scalars(old, new, ('reaction_id',))

            # Migrate changed messages.
            if old.HasField('setup'):
                migrate_setup(old.setup, new.setup)
            if old.HasField('conditions'):
                migrate_conditions(old.conditions, new.conditions)
            for workup in old.workups:
                migrate_workup(workup, new.workups.add())

            if (new.SerializeToString(deterministic=True) !=
                    old.SerializeToString(deterministic=True)):
                num_changed += 1
                # Add a record_modified event.
                new.provenance.record_modified.add(
                    time=dict(value=str(datetime.datetime.now())),
                    person=dict(username='skearnes',
                                email='kearnes@google.com'),
                    details='Automatic migration from schema v0.2.15 to v0.3.0')

        # Save
        logging.info('Changes: %d/%d', num_changed, len(old_dataset.reactions))
        if num_changed:
            os.makedirs(FLAGS.output_dir, exist_ok=True)
            new_filename = os.path.join(FLAGS.output_dir,
                                        os.path.basename(filename))
            logging.info('Saving to %s', new_filename)
            message_helpers.write_message(new_dataset, new_filename)


if __name__ == '__main__':
    flags.mark_flags_as_required(['input_pattern', 'output_dir'])
    app.run(main)
