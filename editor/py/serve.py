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
"""A web editor for Open Reaction Database structures."""

import collections
import contextlib
import datetime
import fcntl
import io
import json
import os
import re
import time
import uuid

import flask
import github
from google.protobuf import text_format
import psycopg2
import psycopg2.sql
import requests

from ord_schema import templating
from ord_schema import message_helpers
from ord_schema import resolvers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import generate_text
from ord_schema.visualization import drawing

# pylint: disable=invalid-name,no-member,inconsistent-return-statements
app = flask.Flask(__name__, template_folder='../html')

# For dataset merges operations like byte-value uploads and enumeration.
TEMP = '/tmp/ord-editor'
if not os.path.isdir(TEMP):
    os.mkdir(TEMP)

# Defaults for development, overridden in docker-compose.yml.
POSTGRES_HOST = os.getenv('POSTGRES_HOST', 'localhost')
POSTGRES_PORT = os.getenv('POSTGRES_PORT', '5432')
POSTGRES_USER = os.getenv('POSTGRES_USER', 'postgres')
POSTGRES_PASS = os.getenv('POSTGRES_PASSWORD', '')

# System user for immutable reactions imported from GitHub pull requests.
REVIEWER = '8df09572f3c74dbcb6003e2eef8e48fc'
# System user for automated testing.
TESTER = '680b0d9fe649417cb092d790907bd5a5'


@app.route('/')
def show_root():
    """The root path redirects to the "datasets" view."""
    return flask.redirect('/datasets')


@app.route('/datasets')
def show_datasets():
    """Lists all the user's pbtxts in the datasets table."""
    names = []
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'SELECT dataset_name FROM datasets WHERE user_id=%s')
        cursor.execute(query, [flask.g.user_id])
        for row in cursor:
            names.append(row[0])
    return flask.render_template('datasets.html',
                                 names=sorted(names),
                                 user_id=flask.g.user_id)


@app.route('/dataset/<name>')
def show_dataset(name):
    """Lists all Reactions contained in the named dataset."""
    dataset = get_dataset(name)
    reactions = []
    for reaction in dataset.reactions:
        reactions.append(reaction.identifiers)
    # Datasets belonging to the "review" user are immutable.
    freeze = flask.g.user_id == REVIEWER
    return flask.render_template('dataset.html',
                                 name=name,
                                 freeze=freeze,
                                 user_id=flask.g.user_id)


@app.route('/dataset/<name>/download')
def download_dataset(name):
    """Returns a pbtxt from the datasets table as an attachment."""
    pbtxt = get_pbtxt(name)
    data = io.BytesIO(pbtxt.encode('utf8'))
    return flask.send_file(data,
                           mimetype='application/protobuf',
                           as_attachment=True,
                           attachment_filename=f'{name}.pbtxt')


@app.route('/dataset/<name>/upload', methods=['POST'])
def upload_dataset(name):
    """Writes the request body to the datasets table without validation."""
    if exists_dataset(name):
        response = flask.make_response(f'dataset already exists: {name}', 409)
        flask.abort(response)
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL('INSERT INTO datasets VALUES (%s, %s, %s)')
        pbtxt = flask.request.get_data().decode('utf-8')
        dataset = dataset_pb2.Dataset()
        text_format.Parse(pbtxt, dataset)  # Validate.
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, name, pbtxt])
        flask.g.db.commit()
    return 'ok'


@app.route('/dataset/<name>/new', methods=['POST'])
def new_dataset(name):
    """Creates a new dataset."""
    if exists_dataset(name):
        response = flask.make_response(f'dataset already exists: {name}', 409)
        flask.abort(response)
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL('INSERT INTO datasets VALUES (%s, %s, %s)')
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, name, ''])
        flask.g.db.commit()
    return 'ok'


@app.route('/dataset/<name>/delete')
def delete_dataset(name):
    """Removes a Dataset."""
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'DELETE FROM datasets WHERE user_id=%s AND dataset_name=%s')
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, name])
        flask.g.db.commit()
    return flask.redirect('/datasets')


@app.route('/dataset/enumerate', methods=['POST'])
def enumerate_dataset():
    """Creates a new dataset based on a template reaction and a spreadsheet.

    Three pieces of information are expected to be POSTed in a json object:
        spreadsheet_name: the original filename of the uploaded spreadsheet.
        spreadsheet_data: a string containing the contents of the spreadsheet.
        template_string: a string containing a text-formatted Reaction proto,
            i.e., the contents of a pbtxt file.
    A new dataset is created from the template and spreadsheet using
    ord_schema.templating.generate_dataset.
    """
    # pylint: disable=broad-except
    data = flask.request.get_json(force=True)
    basename, suffix = os.path.splitext(data['spreadsheet_name'])
    spreadsheet_data = io.StringIO(data['spreadsheet_data'].lstrip('ï»¿'))
    dataframe = templating.read_spreadsheet(spreadsheet_data, suffix=suffix)
    dataset = None
    try:
        dataset = templating.generate_dataset(data['template_string'],
                                              dataframe,
                                              validate=False)
    except ValueError as error:
        flask.abort(flask.make_response(str(error), 400))
    except Exception as error:
        flask.abort(
            flask.make_response(
                f'Unexpected {error.__class__.__name__}: {error}', 400))
    put_dataset(f'{basename}_dataset', dataset)
    return 'ok'


@app.route('/dataset/<name>/reaction/<index>')
def show_reaction(name, index):
    """Render the page representing a single Reaction."""
    dataset = get_dataset(name)
    try:
        index = int(index)
    except ValueError:
        flask.abort(404)
    if len(dataset.reactions) <= index:
        flask.abort(404)
    # Reactions belonging to the "review" user are immutable.
    freeze = flask.g.user_id == REVIEWER
    return flask.render_template('reaction.html',
                                 name=name,
                                 index=index,
                                 freeze=freeze,
                                 user_id=flask.g.user_id)


@app.route('/reaction/download', methods=['POST'])
def download_reaction():
    """Returns a pbtxt file parsed from POST data as an attachment."""
    reaction = reaction_pb2.Reaction()
    reaction.ParseFromString(flask.request.get_data())
    data = io.BytesIO(text_format.MessageToBytes(reaction))
    return flask.send_file(data,
                           mimetype='application/protobuf',
                           as_attachment=True,
                           attachment_filename='reaction.pbtxt')


@app.route('/dataset/<name>/new/reaction')
def new_reaction(name):
    """Adds a new Reaction to the named Dataset and redirects to it."""
    dataset = get_dataset(name)
    dataset.reactions.add()
    put_dataset(name, dataset)
    return flask.redirect('/dataset/%s' % name)


@app.route('/dataset/<name>/clone/<index>')
def clone_reaction(name, index):
    """Copies a specific Reaction to the Dataset and view the Reaction."""
    dataset = get_dataset(name)
    try:
        index = int(index)
    except ValueError:
        flask.abort(404)
    if len(dataset.reactions) <= index:
        flask.abort(404)
    dataset.reactions.add().CopyFrom(dataset.reactions[index])
    index = len(dataset.reactions) - 1
    put_dataset(name, dataset)
    return flask.redirect('/dataset/%s/reaction/%s' % (name, index))


@app.route('/dataset/<name>/delete/reaction/<index>')
def delete_reaction(name, index):
    """Removes a specific Reaction from the Dataset and view the Dataset."""
    dataset = get_dataset(name)
    try:
        index = int(index)
    except ValueError:
        flask.abort(404)
    if len(dataset.reactions) <= index:
        flask.abort(404)
    del dataset.reactions[index]
    put_dataset(name, dataset)
    return flask.redirect('/dataset/%s' % name)


@app.route('/dataset/<name>/delete/reaction_id/<reaction_id>')
def delete_reaction_id(name, reaction_id):
    """Removes a Reaction reference from the Dataset and view the Dataset."""
    dataset = get_dataset(name)
    if reaction_id in dataset.reaction_ids:
        dataset.reaction_ids.remove(reaction_id)
        put_dataset(name, dataset)
        return flask.redirect('/dataset/%s' % name)
    flask.abort(404)


@app.route('/dataset/<name>/delete/reaction_id/')
def delete_reaction_id_blank(name):
    """Removes the first empty Reaction reference from the Dataset."""
    return delete_reaction_id(name, '')


@app.route('/dataset/proto/read/<name>')
def read_dataset(name):
    """Returns a Dataset as a serialized protobuf."""
    dataset = get_dataset(name)
    bites = dataset.SerializeToString(deterministic=True)
    response = flask.make_response(bites)
    response.headers.set('Content-Type', 'application/protobuf')
    return response


@app.route('/dataset/proto/write/<name>', methods=['POST'])
def write_dataset(name):
    """Inserts a protobuf including upload tokens into the datasets table."""
    dataset = dataset_pb2.Dataset()
    dataset.ParseFromString(flask.request.get_data())
    resolve_tokens(dataset)
    put_dataset(name, dataset)
    return 'ok'


@app.route('/dataset/proto/upload/<name>/<token>', methods=['POST'])
def write_upload(name, token):
    """Writes the POST body, names it <token>, and maybe updates the dataset.

    This is part of the upload mechanism. Fields named "bytes_value" can be
    populated only through browser file uploads (e.g. images), and the binary
    data in the files can not be integrated into the Dataset on the client
    because JS can not access the file system. Instead, JS populates the
    bytes_value with a random token and then uploads the file to this endpoint
    with the same token.

    Datasets protos and their bytes_values are reunited in resolve_tokens().

    Args:
        name: The dataset that owns the uploaded asset.
        token: The bytes_value placeholder used in pbtxt to reference the
            upload.

    Returns:
        A 200 response.
    """
    path = get_path(token)
    with open(path, 'wb') as upload:
        upload.write(flask.request.get_data())
    with lock(name):
        dataset = get_dataset(name)
        if resolve_tokens(dataset):
            put_dataset(name, dataset)
    return 'ok'


@app.route('/dataset/proto/download/<token>', methods=['POST'])
def read_upload(token):
    """Echoes a POST body back to the client as a file attachment.

    This is part of the upload mechanism. Since uploaded fields have no type
    information, their values can not be rendered in the browser. Instead, JS
    uses this endpoint to send a previously uploaded bytes_value back to the
    user as a download so they can access it again.

    Args:
        token: A placeholder name for the upload.

    Returns:
        The POST body from the request, after passing through a file.
    """
    data = io.BytesIO(flask.request.get_data())
    return flask.send_file(data,
                           mimetype='application/protobuf',
                           as_attachment=True,
                           attachment_filename=token)


@app.route('/dataset/proto/validate/<message_name>', methods=['POST'])
def validate_reaction(message_name):
    """Receives a serialized Reaction protobuf and runs validations."""
    message = message_helpers.create_message(message_name)
    message.ParseFromString(flask.request.get_data())
    options = validations.ValidationOptions(require_provenance=True)
    output = validations.validate_message(message,
                                          raise_on_error=False,
                                          options=options)
    return json.dumps({'errors': output.errors, 'warnings': output.warnings})


@app.route('/resolve/input', methods=['POST'])
def resolve_input():
    """Resolve an input string to a ReactionInput message."""
    string = flask.request.get_data().decode()
    try:
        reaction_input = resolvers.resolve_input(string)
        bites = reaction_input.SerializeToString(deterministic=True)
        response = flask.make_response(bites)
        response.headers.set('Content-Type', 'application/protobuf')
    except (ValueError, KeyError) as error:
        flask.abort(flask.make_response(str(error), 409))
    return response


@app.route('/resolve/<identifier_type>', methods=['POST'])
def resolve_compound(identifier_type):
    """Resolve a compound name to a SMILES string."""
    compound_name = flask.request.get_data()
    if not compound_name:
        return ''
    try:
        smiles, resolver = resolvers.name_resolve(identifier_type,
                                                  compound_name)
        return flask.jsonify((_canonicalize_smiles(smiles), resolver))
    except ValueError:
        return ''


@app.route('/canonicalize', methods=['POST'])
def canonicalize_smiles():
    """Canonicalizes a SMILES string from a POST request."""
    return flask.jsonify(_canonicalize_smiles(flask.request.get_data()))


def _canonicalize_smiles(smiles):
    """Canonicalizes a SMILES string."""
    try:
        return resolvers.canonicalize_smiles(smiles)
    except ValueError:
        return smiles  # Return the original SMILES on failure.


@app.route('/render/reaction', methods=['POST'])
def render_reaction():
    """Receives a serialized Reaction message and returns a block of HTML
    that contains a visual summary of the reaction."""
    reaction = reaction_pb2.Reaction()
    reaction.ParseFromString(flask.request.get_data())
    if not (reaction.inputs or reaction.outcomes):
        return ''
    try:
        html = generate_text.generate_html(reaction)
        return flask.jsonify(html)
    except (ValueError, KeyError):
        return ''


@app.route('/render/compound', methods=['POST'])
def render_compound():
    """Returns an HTML-tagged SVG for the given Compound."""
    compound = reaction_pb2.Compound()
    compound.ParseFromString(flask.request.get_data())
    try:
        mol = message_helpers.mol_from_compound(compound)
        return flask.jsonify(drawing.mol_to_svg(mol))
    except ValueError:
        return ''


@app.route('/dataset/proto/compare/<name>', methods=['POST'])
def compare(name):
    """For testing, compares a POST body to an entry in the datasets table.

    Returns HTTP status 200 if their pbtxt representations are equal as strings
    and 409 if they are not. See "make test".

    Args:
        name: The dataset record to compare.

    Returns:
        HTTP status 200 for a match and 409 if there is a difference.
    """
    remote = dataset_pb2.Dataset()
    remote.ParseFromString(flask.request.get_data())
    local = get_dataset(name)
    remote_ascii = text_format.MessageToString(remote)
    local_ascii = text_format.MessageToString(local)
    if remote_ascii != local_ascii:
        return 'differs', 409  # "Conflict"
    return 'equals'


@app.route('/js/<script>')
def js(script):
    """Accesses any built JS file by name from the Closure output directory."""
    path = flask.safe_join(
        os.path.join(os.path.dirname(__file__), '../gen/js/ord'), script)
    return flask.send_file(get_file(path), attachment_filename=script)


@app.route('/css/<sheet>')
def css(sheet):
    """Accesses any CSS file by name."""
    path = flask.safe_join(os.path.join(os.path.dirname(__file__), '../css'),
                           sheet)
    return flask.send_file(get_file(path), attachment_filename=sheet)


@app.route('/ketcher/iframe')
def ketcher_iframe():
    """Accesses a website serving Ketcher."""
    return flask.render_template('ketcher_iframe.html')


@app.route('/ketcher/info')
def indigo():
    """Dummy indigo endpoint to prevent 404 errors."""
    return '', 204


@app.route('/ketcher/<path:file>')
def ketcher(file):
    """Accesses any built Ketcher file by name."""
    path = flask.safe_join(
        os.path.join(os.path.dirname(__file__), '../ketcher/dist'), file)
    return flask.send_file(get_file(path), attachment_filename=file)


@app.route('/dataset/deps.js')
@app.route('/dataset/<file_name>/deps.js')
@app.route('/dataset/<file_name>/reaction/deps.js')
def deps(file_name=None):
    """Returns empty for deps table requests since this app doesn't use them."""
    del file_name  # Unused.
    return ''


@app.route('/ketcher/molfile', methods=['POST'])
def get_molfile():
    """Retrieves a POSTed Compound message string and returns a MolFile."""
    compound = reaction_pb2.Compound()
    compound.ParseFromString(flask.request.get_data())
    try:
        molblock = message_helpers.molblock_from_compound(compound)
        return flask.jsonify(molblock)
    except ValueError:
        return 'no existing structural identifier', 404


@app.route('/review')
def show_submissions():
    """For the review user only, render datasets with GitHub metadata."""
    if flask.g.user_id != REVIEWER:
        return flask.redirect('/')
    pull_requests = collections.defaultdict(list)
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'SELECT dataset_name FROM datasets WHERE user_id=%s')
        cursor.execute(query, [REVIEWER])
        for row in cursor:
            name = row[0]
            match = re.match('^PR_([0-9]+) ___(.*)___ (.*)', name)
            if match is None:
                continue
            number, title, short_name = match.groups()
            pull_requests[(number, title)].append((short_name, name))
    return flask.render_template('submissions.html',
                                 pull_requests=pull_requests)


@app.route('/review/sync')
def sync_reviews():
    """Import all current pull requests into the datasets table.

    These datasets have two extra pieces of metadata: a GitHub PR number and
    the PR title text. These are encoded into the dataset name in Postgres
    using delimiters."""
    if flask.g.user_id != REVIEWER:
        return flask.redirect('/')
    client = github.Github()
    repo = client.get_repo('Open-Reaction-Database/ord-data')
    user_id = flask.g.user_id
    with flask.g.db.cursor() as cursor:
        # First reset all datasets under review.
        query = psycopg2.sql.SQL('DELETE FROM datasets WHERE user_id=%s')
        cursor.execute(query, [REVIEWER])
        # Then import all pbtxts from open PR's.
        for pr in repo.get_pulls():
            for remote in pr.get_files():
                if not remote.filename.endswith('.pbtxt'):
                    continue
                pbtxt = requests.get(remote.raw_url).text
                name = 'PR_%d ___%s___ %s' % (pr.number, pr.title,
                                              remote.filename[:-6])
                query = psycopg2.sql.SQL(
                    'INSERT INTO datasets VALUES (%s, %s, %s)')
                cursor.execute(query, [user_id, name, pbtxt])
    flask.g.db.commit()
    return flask.redirect('/review')


@app.route('/revisions/<name>')
def show_revisions(name):
    """List all historical revisions of a dataset for download."""
    # Offset relative to UTC, for localized timestamp representation.
    hours = 0
    if 'tz' in flask.request.args:
        hours = int(flask.request.args.get('tz'))
    revisions = []
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL('SELECT timestamp FROM revisions '
                                 'WHERE user_id=%s and dataset_name=%s '
                                 'ORDER BY timestamp DESC')
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, name])
        for row in cursor:
            timestamp = int(row[0])
            localtime = (
                datetime.datetime.utcfromtimestamp(timestamp) -
                datetime.timedelta(hours=hours)).strftime('%Y-%m-%d %H:%M:%S')
            revisions.append((localtime, timestamp))
    return flask.render_template('revisions.html',
                                 name=name,
                                 revisions=revisions,
                                 user_id=user_id)


@app.route('/checkpoint/<name>', methods=['POST'])
def checkpoint(name):
    """Append a dataset to the "revisions" list."""
    dataset = dataset_pb2.Dataset()
    dataset.ParseFromString(flask.request.get_data())
    resolve_tokens(dataset)
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'INSERT INTO revisions VALUES (%s, %s, %s, %s)')
        timestamp = int(time.time())
        user_id = flask.g.user_id
        pbtxt = text_format.MessageToString(dataset, as_utf8=True)
        cursor.execute(query, [user_id, name, pbtxt, timestamp])
        flask.g.db.commit()
    return 'ok'


@app.route('/revision/<name>/<timestamp>')
def download_revision(name, timestamp):
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'SELECT pbtxt FROM revisions '
            'WHERE user_id=%s AND dataset_name=%s AND timestamp=%s')
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, name, int(timestamp)])
        if cursor.rowcount == 0:
            flask.abort(404)
        pbtxt = cursor.fetchone()[0]
    data = io.BytesIO(pbtxt.encode('utf8'))
    return flask.send_file(data,
                           mimetype='application/protobuf',
                           as_attachment=True,
                           attachment_filename=f'{name}-{timestamp}.pbtxt')


@app.after_request
def prevent_caching(response):
    """Prevents caching any of this app's resources on the client."""
    response.headers['Cache-Control'] = 'no-cache,no-store,must-revalidate'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    response.headers['Cache-Control'] = 'public, max-age=0'
    # Make the user ID accessible for logging.
    response.headers['User-Id'] = flask.g.get('user_id', 'unknown')
    return response


def get_dataset(name):
    """Reads a pbtxt proto from the datasets table and parses it."""
    pbtxt = get_pbtxt(name)
    dataset = dataset_pb2.Dataset()
    text_format.Parse(pbtxt, dataset)
    return dataset


def get_pbtxt(name):
    """Reads a pbtxt proto from the datasets."""
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'SELECT pbtxt FROM datasets WHERE user_id=%s AND dataset_name=%s')
        cursor.execute(query, [flask.g.user_id, name])
        if cursor.rowcount == 0:
            flask.abort(404)
        return cursor.fetchone()[0]


def put_dataset(name, dataset):
    """Write a dataset proto to the dataset table, clobbering if needed."""
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'INSERT INTO datasets VALUES (%s, %s, %s) '
            'ON CONFLICT (user_id, dataset_name) DO UPDATE SET pbtxt=%s')
        user_id = flask.g.user_id
        pbtxt = text_format.MessageToString(dataset, as_utf8=True)
        cursor.execute(query, [user_id, name, pbtxt, pbtxt])
        flask.g.db.commit()


@contextlib.contextmanager
def lock(file_name):
    """Blocks until an exclusive lock on the named file is obtained.

    This is part of the upload mechanism. Byte_value fields are populated
    through browser file uploads and so must be merged with Dataset protos on
    the server. The file uploads and the Dataset proto arrive asynchronously in
    concurrent connections. Locks ensure that only one request at a time
    accesses a Datset's .pbtxt file.

    Args:
        file_name: Name of the file to pass to the fcntl system call.

    Yields:
        The locked file descriptor.
    """
    path = get_path(file_name, suffix='.lock')
    with open(path, 'w') as lock_file:
        fcntl.lockf(lock_file, fcntl.LOCK_EX)
        try:
            yield lock_file
        finally:
            fcntl.lockf(lock_file, fcntl.LOCK_UN)


def resolve_tokens(proto):
    """Fills in bytes_value fields using client-generated placeholder tokens.

    This is part of the upload mechanism. It acts by recursion on the tree
    structure of the proto, hunting for fields named "bytes_value" and
    comparing the fields' values against uploaded files in the root directory.
    See write_dataset() and write_upload().

    Args:
        proto: A protobuf message that may contain a bytes_value somewhere.

    Returns:
        True if a bytes_value was matched anywhere in the tree.
    """
    matched = False
    if 'ListFields' in dir(proto):
        for descriptor, message in proto.ListFields():
            if descriptor.name == 'bytes_value':
                token = message[:20].decode('utf8')
                if token.startswith('upload_'):
                    path = get_path(token)
                    if os.path.isfile(path):
                        with open(path, 'rb') as f:
                            proto.bytes_value = f.read()
                        matched = True
            matched |= resolve_tokens(message)
    elif 'append' in dir(proto):
        for message in proto:
            matched |= resolve_tokens(message)
    elif 'keys' in dir(proto):
        for key in proto.keys():
            message = proto[key]
            matched |= resolve_tokens(message)
    return matched


def get_file(path):
    """Get file contents as a BytesIO object."""
    if not os.path.exists(path):
        flask.abort(404)
    with open(path, 'rb') as f:
        # NOTE(kearnes): Workaround for unclosed file warnings. See
        # https://github.com/pallets/werkzeug/issues/1785.
        data = io.BytesIO(f.read())
    return data


def get_user_path():
    """Returns the path of the current user's temp directory.

    Uses flask.safe_join to check for directory traversal attacks.

    Returns:
        Path to the user directory.
    """
    return flask.safe_join(TEMP, flask.g.user_id)


def get_path(file_name, suffix=''):
    """Returns a safe path in the user's temp directory.

    Uses flask.safe_join to check for directory traversal attacks.

    Args:
        file_name: Text filename.
        suffix: Text filename suffix. Defaults to '.pbtxt'.

    Returns:
        Path to the requested file.
    """
    return flask.safe_join(get_user_path(), f'{file_name}{suffix}')


def exists_dataset(name):
    """True if a dataset with the given name is defined for the current user."""
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'SELECT 1 FROM datasets WHERE user_id=%s AND dataset_name=%s')
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, name])
        return cursor.rowcount > 0


@app.route('/login')
def show_login():
    """Presents a form to set a new access token from a given user ID."""
    return flask.render_template('login.html')


@app.route('/authenticate', methods=['GET', 'POST'])
def authenticate():
    """Issue a new access token for a given user ID."""
    # GET authentications always login as the test user.
    if flask.request.method == 'GET':
        return issue_access_token(TESTER)
    user_id = flask.request.form.get('user_id')
    if user_id is None or re.match('^[0-9a-fA-F]{32}$', user_id) is None:
        return flask.redirect('/login')
    return issue_access_token(user_id)


def issue_access_token(user_id):
    """Login as the given user and set the access token in a response."""
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL('INSERT INTO logins VALUES (%s, %s, %s)')
        access_token = uuid.uuid4().hex
        timestamp = int(time.time())
        cursor.execute(query, [access_token, user_id, timestamp])
        flask.g.db.commit()
        response = flask.redirect('/')
        # Expires in a year.
        response.set_cookie('Access-Token', access_token, max_age=31536000)
        return response


def make_user(name='auto'):
    """Writes a new user ID and returns it

    Args:
        name: Hopefully a readable label for the user, not currently used in UI.

    Returns:
        The 32-character generated UUID of the user, currently used in the UI.
    """
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL('INSERT INTO users VALUES (%s, %s, %s)')
        user_id = uuid.uuid4().hex
        timestamp = int(time.time())
        cursor.execute(query, [user_id, name, timestamp])
        flask.g.db.commit()
    return user_id


@app.before_request
def init_user():
    """Connects to the DB and authenticates the user."""
    flask.g.db = psycopg2.connect(dbname='editor',
                                  user=POSTGRES_USER,
                                  password=POSTGRES_PASS,
                                  host=POSTGRES_HOST,
                                  port=int(POSTGRES_PORT))
    if flask.request.path in ('/login', '/authenticate'):
        return
    access_token = flask.request.cookies.get('Access-Token')
    if access_token is None:
        # Automatically login as a new user.
        user_id = make_user()
        return issue_access_token(user_id)
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            "SELECT user_id FROM logins WHERE access_token=%s")
        cursor.execute(query, [access_token])
        if cursor.rowcount == 0:
            # Automatically login as a new user.
            user_id = make_user()
            return issue_access_token(user_id)
        user_id = cursor.fetchone()[0]
    flask.g.user_id = user_id
    temp = get_user_path()
    if not os.path.isdir(temp):
        os.mkdir(temp)


@app.route('/logout')
def logout():
    """Clear the access token and redirect to /login."""
    response = flask.redirect('/login')
    response.set_cookie('Access-Token', '', expires=0)
    return response
