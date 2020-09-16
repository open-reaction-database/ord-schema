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
linked. Query parameters are communicated in URL GET params:

  component=<pattern;source;(exact|substructure|similarity|smarts)>

    The second token specifies whether the predicate should match an input or
    an output.

    The last token specifies the matching criterion. The default is "exact".
    The pattern is a SMILES string, unless the token is "smarts" in which case
    the pattern is a SMARTS string.

    Component may be repeated any number of times.

  reaction_ids=<ids>

  reaction_smarts=<smarts>

These query types are mutually exclusive. All query parameters are assumed to
be URL-encoded.
"""

import os

import flask

from ord_schema.interface import query

app = flask.Flask(__name__, template_folder='.')
app.config['ORD_POSTGRES_HOST'] = os.getenv('ORD_POSTGRES_HOST', 'localhost')


@app.route('/')
def show_root():
    """Shows the web form.

    If there are query params, then the query is executed and the form is
    populated with the results. The form fields are populated with the params.
    """
    reaction_ids = flask.request.args.get('reaction_ids')
    reaction_smarts = flask.request.args.get('reaction_smarts')
    components = flask.request.args.getlist('component')
    use_stereochemistry = flask.request.args.get('use_stereochemistry')
    similarity = flask.request.args.get('similarity')

    if reaction_ids is not None:
        command = query.ReactionIdQuery(reaction_ids.split(','))
    elif reaction_smarts is not None:
        command = query.ReactionSmartsQuery(reaction_smarts)
    else:
        predicates = []
        for component in components:
            pattern, source, mode_name = component.split(';')
            table = query.ReactionComponentPredicate.SOURCE_TO_TABLE[source]
            mode = query.ReactionComponentPredicate.MatchMode.from_name(
                mode_name)
            predicates.append(
                query.ReactionComponentPredicate(pattern, table, mode))
        kwargs = {}
        if use_stereochemistry is not None:
            kwargs['do_chiral_sss'] = use_stereochemistry
        if similarity is not None:
            kwargs['tanimoto_threshold'] = float(similarity)
        command = query.ReactionComponentQuery(predicates, **kwargs)
    dataset = connect().run_query(command, return_ids=True)
    return flask.render_template('search.html',
                                 reaction_ids=dataset.reaction_ids,
                                 query=command.json())


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
                             host=app.config['ORD_POSTGRES_HOST'],
                             port=5432)
