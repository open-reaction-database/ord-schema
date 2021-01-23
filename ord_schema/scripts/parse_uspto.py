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
"""Parses CML from the NRD.

Data is at https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873.

See also https://depth-first.com/articles/2019/01/28/the-nextmove-patent-reaction-dataset/.
"""

import collections
import glob
import re
from xml.etree import ElementTree

from absl import app
from absl import flags
from absl import logging

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

FLAGS = flags.FLAGS
flags.DEFINE_string('input_pattern', None, 'Input pattern for CML files.')
flags.DEFINE_string('output', None, 'Output Dataset filename.')

# XML namespaces.
NAMESPACES = {
    'cml': 'http://www.xml-cml.org/schema',
    'dan': 'http://bitbucket.org/dan2097'
}


def parse_reaction(root):
    """Parses reaction CML into a Reaction message."""


def main(argv):
    del argv  # Only used by app.run().
    reactions = []
    for filename in glob.glob(FLAGS.input_pattern):
        logging.info(filename)
        tree = ElementTree.parse(filename)
        root = tree.getroot()
        for reaction_cml in root.iterfind('cml:reaction', namespaces=NAMESPACES):
            reaction = parse_reaction(reaction_cml)
            reactions.append(reaction)
    dataset = dataset_pb2.Dataset(reactions=reactions)
    if FLAGS.output:
        message_helpers.write_message(dataset, FLAGS.output)


if __name__ == '__main__':
    flags.mark_flag_as_required('input_pattern')
    app.run(main)
