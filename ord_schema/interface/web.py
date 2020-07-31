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

"""Query interface to the Open Reaction Database in Postgres.

This is a web interface. The client is stateless. Full information is in the
URL, so results can be linked.

There are three kinds of query.

  inputs and output SMILES or SMARTS

    These match against individual reagents and results. Each SMILES expression
    can be freely specified for "exact", "substructure", or "similarity"
    matching. SMARTS expressions allow no such specification.

  reaction ID

    Matches either zero or one reactions. The ID has the form,
      "ord-<32-bit hex>".

  reaction SMILES

    Matches against RDKit fingerprints.

Query parameters are communicated in URL GET params.

  input=<(smiles|smarts);(exact|substructure|similarity)>

    The token after the semicolon specifies the matching criterion. The default
    is "exact". The "input" param may be repeated.

  output=<(smiles|smarts);(exact|substructure|similarity)>

    The default matching criterion is "exact".

  reaction=<id>

  reaction_smiles=<smiles>

If multiple SMILES are given, then the query is interpreted as conjunction.

It is an error to specify more than one of "(input|output)", "reaction", and
"reaction_smiles" in the same URL.

All query parameters are assumed to be URL-encoded.
"""

import contextlib
import fcntl
import os
import random
import re
import string
import urllib
import uuid

import flask

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

from google.protobuf import text_format

import query

app = flask.Flask(__name__, template_folder='.')

def connect():
  return query.OrdPostgres(host='localhost', port=5430)

def get_reaction_by_id(reaction_id):
  db = connect()
  return db.get_reaction_by_id(reaction_id)

def get_reactions_by_smiles(reaction_smiles):
  pass

def get_reactions_by_reagents(inputs, output):
  pass

@app.route('/')
def show_root():
  """Shows the web form.

  If there are query params, then the query is executed and the form is
  populated with the results. The form fields are populated with the params.
  """
  reaction_id = flask.request.args.get('reaction')
  reaction_smiles = flask.request.args.get('reaction_smiles')
  inputs = flask.request.args.getlist("input")
  output = flask.request.args.get('output')
  if reaction_id is not None:
    if reaction_smiles is not None or output is not None or len(inputs) > 0:
      raise Exception('"reaction" can not be combined with other constraints')
    dataset = get_reaction_by_id(reaction_id)
  elif reaction_smiles is not None:
    if output is not None or len(inputs) > 0:
      raise Exception(
          '"reaction_smiles" can not be combined with other constraints')
    dataset = get_reactions_by_smiles(reaction_smiles)
  elif inputs is not None or output is not None:
    inputs = [i.split(';', 1) for i in inputs] if inputs else []
    output = output.split(';', 1) if output else None
    dataset = get_reactions_by_reagents(inputs, output)
  else:
    reactions = []
  return flask.render_template('web.html', dataset=dataset)
