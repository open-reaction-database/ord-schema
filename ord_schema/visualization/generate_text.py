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
"""Text generation script for Reaction messages.

This script is meant to convert a Reaction message into a full-fledged text
paragraph that would be appropriate for inclusion in a Supplemental Information
document for a publication.

Example usage:
* For normal operation from a pbtxt
  $ python generate_text.py --input_reaction=reaction.pb --input_format=pbtxt
      --type text
"""

import os
import re
import jinja2

from absl import app
from absl import flags

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import filters

FLAGS = flags.FLAGS
flags.DEFINE_string('input_file', None, 'File containing a Reaction message.')
flags.DEFINE_enum('type', 'text', ['text', 'html'],
                  'Text or HTML output format.')
flags.DEFINE_string('output', None, 'Filename for output Dataset.')

# pylint: disable=redefined-outer-name
with open(os.path.join(os.path.dirname(__file__), 'template.html'),
          'r') as fid:
    _HTML_TEMPLATE = fid.read()

with open(os.path.join(os.path.dirname(__file__), 'template.txt'), 'r') as fid:
    _TEXT_TEMPLATE = fid.read()


def _generate(reaction, template_string, line_breaks):
    """Renders a Jinja2 template string with a reaction message.

    Args:
        reaction: a Reaction message.
        template_string: Jinja template string.
        line_breaks: Whether to keep line breaks.

    Returns:
        Rendered string.
    """
    env = jinja2.Environment(loader=jinja2.BaseLoader())
    env.filters.update(filters.TEMPLATE_FILTERS)
    template = env.from_string(template_string)
    text = template.render(reaction=reaction)

    # Fix line breaks, extra spaces, "a" versus "an"
    if not line_breaks:
        text = ''.join(text.strip().splitlines())
    text = re.sub(r'[ ]{2,}', ' ', text)
    text = re.sub(r' a ([aeiouAEIOU])', r' an \1', text)
    text = re.sub(r'[ ]([\.\;\)\,]){1,2}', r'\1', text)
    return text


def generate_text(reaction):
    """Generates a textual reaction description."""
    return _generate(reaction, _TEXT_TEMPLATE, False)


def generate_html(reaction):
    return _generate(reaction, _HTML_TEMPLATE, True)


def main(argv):
    del argv  # Only used by app.run()
    reaction = message_helpers.load_message(FLAGS.input_file,
                                            reaction_pb2.Reaction)

    if FLAGS.type == 'html':
        text = generate_html(reaction)
    elif FLAGS.type == 'text':
        text = generate_text(reaction)

    if FLAGS.output:
        with open(FLAGS.output, 'w') as fid:
            fid.write(text)
    else:
        print(text)


if __name__ == '__main__':
    flags.mark_flag_as_required('input_file')
    app.run(main)
