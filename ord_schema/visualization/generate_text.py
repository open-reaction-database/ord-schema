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
from typing import Optional

import jinja2

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import filters


def _generate(reaction: reaction_pb2.Reaction, template_string: str,
              line_breaks: bool, **kwargs) -> str:
    """Renders a Jinja2 template string with a reaction message.

    Args:
        reaction: a Reaction message.
        template_string: Jinja template string.
        line_breaks: Whether to keep line breaks.
        **kwargs: Additional keyword arguments to render().

    Returns:
        Rendered string.
    """
    env = jinja2.Environment(loader=jinja2.BaseLoader())
    env.filters.update(filters.TEMPLATE_FILTERS)
    template = env.from_string(template_string)
    text = template.render(reaction=reaction, **kwargs)

    # Fix line breaks, extra spaces, "a" versus "an"
    if not line_breaks:
        text = ''.join(text.strip().splitlines())
    text = re.sub(r'[ ]{2,}', ' ', text)
    text = re.sub(r' a ([aeiouAEIOU])', r' an \1', text)
    text = re.sub(r'[ ]([.;),]){1,2}', r'\1', text)
    return text


def generate_text(reaction: reaction_pb2.Reaction) -> str:
    """Generates a textual reaction description."""
    with open(os.path.join(os.path.dirname(__file__), 'template.txt'),
              'r') as f:
        template = f.read()
    return _generate(reaction, template_string=template, line_breaks=False)


def generate_html(reaction: reaction_pb2.Reaction,
                  compact=False,
                  bond_length: Optional[int] = None) -> str:
    """Generates an HTML reaction description."""
    # Special handling for e.g. USPTO reactions.
    reaction_smiles = message_helpers.get_reaction_smiles(reaction)
    if reaction_smiles and not reaction.inputs and not reaction.outcomes:
        reaction = message_helpers.reaction_from_smiles(reaction_smiles)
    with open(os.path.join(os.path.dirname(__file__), 'template.html'),
              'r') as f:
        template = f.read()
    if not bond_length:
        if compact:
            bond_length = 18
        else:
            bond_length = 25
    return _generate(reaction,
                     template_string=template,
                     line_breaks=True,
                     compact=compact,
                     bond_length=bond_length)
