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

import flask

from ord_schema.interface import query

app = flask.Flask(__name__, template_folder='.')


@app.route('/')
def show_root():
    """Shows the web form.

    If there are query params, then the query is executed and the form is
    populated with the results. The form fields are populated with the params.
    """
    reaction_ids = flask.request.args.get('reaction_ids')
    reaction_smarts = flask.request.args.get('reaction_smarts')
    inputs = flask.request.args.getlist('input')
    outputs = flask.request.args.getlist('output')
    do_chiral_sss = flask.request.args.get('use_stereochemistry')
    tanimoto_threshold = flask.request.args.get('similarity')

    if reaction_ids is not None:
        command = query.ReactionIdQuery(reaction_ids.split(','))
    elif reaction_smarts is not None:
        command = query.ReactionSmartsQuery(reaction_smarts)
    else:
        predicates = []
        for table, queries in [('inputs', inputs), ('outputs', outputs)]:
            splits = [i.split(';', 1) for i in queries] if queries else []
        for smiles, mode_name in splits:
            mode = query.ReactionComponentPredicate.MatchMode.from_name(
                mode_name)
            predicates.append(
                query.ReactionComponentPredicate(smiles, table, mode))
        command = query.ReactionComponentQuery(
            predicates, do_chiral_sss=do_chiral_sss,
            tanimoto_threshold=float(tanimoto_threshold))
    dataset = connect().run_query(command, return_ids=True)
    return flask.render_template('search.html',
                                 reaction_ids=dataset.reaction_ids,
                                 predicate=command.json())


@app.route('/id/<reaction_id>')
def show_id(reaction_id):
    """Returns the pbtxt of a single reaction as plain text."""
    dataset = connect().run_query(query.ReactionIdQuery([reaction_id]))
    if len(dataset.reactions) == 0:
        return flask.abort(404)
    return flask.Response(str(dataset.reactions[0]), mimetype='text/plain')


def connect():
    return query.OrdPostgres(dbname='ord',
                             user='ord-postgres',
                             password='ord-postgres',
                             host='localhost',
                             port=5432)
