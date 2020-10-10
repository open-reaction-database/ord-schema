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
"""Creates a Dataset by enumerating a template with a spreadsheet.

Example usage:

$ python dataset_templating.py \
    --template=my_reaction.pbtxt \
    --spreadsheet=my_experiments.csv \
    --output=my_dataset.pbtxt
"""

import os

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema import templating

FLAGS = flags.FLAGS
flags.DEFINE_string('template', None,
                    'Path to a Reaction pbtxt file defining a template.')
flags.DEFINE_string(
    'spreadsheet', None,
    'Path to a spreadsheet file (with a header row) defining '
    'values to replace placeholders in the template.')
flags.DEFINE_string('output', None, 'Filename for output Dataset.')
flags.DEFINE_boolean('validate', True, 'If True, validate Reaction protos.')


def main(argv):
    del argv  # Only used by app.run().
    with open(FLAGS.template) as f:
        template_string = f.read()
    df = templating.read_spreadsheet(FLAGS.spreadsheet)
    logging.info('generating new Dataset from %s and %s', FLAGS.template,
                 FLAGS.spreadsheet)
    dataset = templating.generate_dataset(template_string,
                                          df,
                                          validate=FLAGS.validate)
    if FLAGS.output:
        output_filename = FLAGS.output
    else:
        basename, _ = os.path.splitext(FLAGS.spreadsheet)
        output_filename = os.path.join(f'{basename}_dataset.pbtxt')
    logging.info('writing new Dataset to %s', output_filename)
    message_helpers.write_message(dataset, output_filename)


if __name__ == '__main__':
    flags.mark_flags_as_required(['template', 'spreadsheet'])
    app.run(main)
