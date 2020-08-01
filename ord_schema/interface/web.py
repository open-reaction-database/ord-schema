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

"""Query web interface to the Open Reaction Database in Postgres.

The client is stateless. Full information is in the URL so results can be
linked.

There are three kinds of query.

  inputs and output SMILES or SMARTS

    These match against individual reagents and results. Each SMILES expression
    can be freely specified for "exact", "substructure", or "similarity"
    matching. SMARTS expressions do not allow such a specification.

  reaction ID

    Matches either zero or one reactions. The ID has the form,
      "ord-<32-bit hex>".

  reaction SMILES

    Matches against RDKit fingerprints.

Query parameters are communicated in URL GET params.

  input=<pattern;(exact|substructure|similarity|smarts)>

  output=<pattern;(exact|substructure|similarity|smarts)>

    The token after the semicolon specifies the matching criterion. The default
    is "exact". The pattern is a SMILES string, unless the token is "smarts" in
    which case the pattern is a SMARTS string.

    The "input" param may be repeated. The "output" param may not.

  reaction_id=<id>

  reaction_smiles=<smiles>

If multiple conditions are given, then the query is interpreted as conjunction.

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


@app.route('/')
def show_root():
    """Shows the web form.

    If there are query params, then the query is executed and the form is
    populated with the results. The form fields are populated with the params.
    """
    reaction_id = flask.request.args.get('reaction_id')
    reaction_smiles = flask.request.args.get('reaction_smiles')
    inputs = flask.request.args.getlist("input")
    output = flask.request.args.get('output')

    predicate = query.Predicate()

    if reaction_id is not None:
        predicate.set_reaction_id(reaction_id)

    if reaction_smiles is not None:
        predicate.set_reaction_smiles(reaction_smiles)

    in_splits = [i.split(';', 1) for i in inputs] if inputs else []
    for smiles, mode_name in in_splits:
        mode = query.Predicate.MatchMode.from_name(mode_name)
        predicate.add_input(smiles, mode)

    if output is not None:
        smiles, mode_name = output.split(';', 1)
        mode = query.Predicate.MatchMode.from_name(mode_name)
        predicate.set_output(smiles, mode)

    if reaction_id or reaction_smiles or inputs or output:
        db = query.OrdPostgres(host='localhost', port=5430)
        reaction_ids = db.predicate_search_ids(predicate)
    else:
        reaction_ids = []
    return flask.render_template(
        'web.html', reaction_ids=reaction_ids, predicate=predicate.json())


@app.route('/id/<reaction_id>')
def show_id(reaction_id):
  """Returns the pbtxt of a single reaction as plain text."""
  predicate = query.Predicate()
  predicate.set_reaction_id(reaction_id)
  db = query.OrdPostgres(host='localhost', port=5430)
  dataset = db.predicate_search(predicate)
  if len(dataset.reactions) == 0:
      return flask.abort(404)
  return flask.Response(str(dataset.reactions[0]), mimetype='text/plain')
