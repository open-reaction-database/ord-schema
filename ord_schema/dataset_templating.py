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
"""Script and functions for enumerating a dataset using a template reaction and
a spreadsheet file.

The templating code has specific expectations for how the reaction pbtxt and
spreadsheet are defined, namely that placeholder values in the pbtxt begin and
end with a "$" (dollar sign) and that these match a unique column header in the
spreadsheet file.

Example usage:

* For normal dataset generation from a template reaction:
  $ python dataset_templating.py --template=my_reaction.pbtxt
      --spreadsheet=my_experiments.csv --output=my_dataset.pbtxt
"""

import os
import re
import pandas as pd

from absl import app
from absl import flags
from absl import logging

from google.protobuf import text_format

from ord_schema import validations
from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('template', None,
                    'Path to a Reaction pbtxt file defining a template.')
flags.DEFINE_string(
    'spreadsheet', None,
    'Path to a spreadsheet file (with a header row) defining '
    'values to replace placeholders in the template.')
flags.DEFINE_string('output', None, 'Filename for output Dataset.')
flags.DEFINE_boolean('validate', True, 'If True, validate Reaction protos.')


def read_spreadsheet(file_name):
    """Reads a {csv, xls, xlsx} spreadsheet file."""
    _, suffix = os.path.splitext(file_name)
    if suffix in ['xls', 'xlsx']:
        return pd.read_excel(file_name, dtype=str)
    return pd.read_csv(file_name, dtype=str)


def _replace(string, substitutions):
    """Performs substring substitutions according to a dictionary.
    
    Inputs:
        string: A string whose contents should be modified.
        substitutions: A dictionary where each (key, value) pair defines
            a substring to replace and what its replacement should be.

    Returns:
        The modified string.
    """
    pattern = re.compile('|'.join(map(re.escape, substitutions.keys())))
    return pattern.sub(lambda match: substitutions[match.group(0)], string)


def generate_dataset(template_string, df, validate=True):
    """Generates a Dataset from a template reaction formatted as text according
    to entries in a pd.Dataframe.

    Inputs:
        template_string: The contents of a Reaction pbtxt where placeholder
            values to be replaced are defined between dollar sign. For example,
            a SMILES identifier value could be "$product_smiles$". Spaces
            are not allowed.
        df: Pandas Dataframe where each row corresponds to one reaction and
            column names match placeholders in the template_string.
        validate: Optional Boolean controlling whether Reaction messages should
            be validated as they are defined. Defaults to True.

    Returns:
        A Dataset message.

    Raises:
        ValueError: If there is no match for a placeholder string in df.
        ValueError: If validate is True and there are validation errors when
            validating an enumerated Reaction message.

    """
    placeholders = set(re.findall(r'\$\w+\$', template_string))
    for placeholder in placeholders:
        if placeholder not in df.columns:
            # Allow "$my_placeholder$" to match "my_placeholder" in df.
            if placeholder[1:-1] not in df.columns:
                raise ValueError(f'Placeholder {placeholder} not found as a'
                                 ' column in dataset spreadsheet')
            df.rename(columns={placeholder[1:-1]: placeholder}, inplace=True)

    reactions = []
    for _, substitutions in df[placeholders].iterrows():
        reaction_text = _replace(template_string, substitutions)
        reaction = text_format.Parse(reaction_text, reaction_pb2.Reaction())
        if validate:
            errors = validations.validate_message(reaction,
                                                  raise_on_error=False)
            if errors:
                raise ValueError(f'Enumerated Reaction is not valid: {errors}')
        reactions.append(reaction)

    return dataset_pb2.Dataset(reactions=reactions)


def main(argv):
    del argv  # Only used by app.run().
    flags.mark_flags_as_required(['template', 'spreadsheet'])
    
    with open(FLAGS.template, 'r') as fid:
        template_string = fid.read()
    df = read_spreadsheet(FLAGS.spreadsheet)

    logging.info('generating new Dataset from %s and %s', FLAGS.template,
                 FLAGS.spreadsheet)
    dataset = generate_dataset(template_string, df, validate=FLAGS.validate)

    if FLAGS.output:
        output_filename = FLAGS.output
    else:
        basename, _ = os.path.splitext(FLAGS.spreadsheet)
        output_filename = os.path.join(f'{basename}_dataset.pbtxt')
    logging.info('writing new Dataset to %s', output_filename)
    message_helpers.write_message(dataset, output_filename)


if __name__ == '__main__':
    app.run(main)
