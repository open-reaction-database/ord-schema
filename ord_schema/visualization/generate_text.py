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
"""Text generation for Reaction messages."""

import os
import re
import jinja2

from ord_schema.visualization import filters


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
    with open(os.path.join(os.path.dirname(__file__), 'template.txt'),
              'r') as f:
        template = f.read()
    return _generate(reaction, template_string=template, line_breaks=False)


def generate_html(reaction):
    """Generates an HTML reaction description."""
    with open(os.path.join(os.path.dirname(__file__), 'template.html'),
              'r') as f:
        template = f.read()
    return _generate(reaction, template, True)
