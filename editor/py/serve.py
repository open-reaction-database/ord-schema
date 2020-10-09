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

import contextlib
import difflib
import fcntl
import io
import json
import os
import pprint
import psycopg2
import psycopg2.sql
import re
import shutil
import subprocess
import time
import urllib
import uuid

import flask
import github
from google.protobuf import text_format

from ord_schema import dataset_templating
from ord_schema import message_helpers
from ord_schema import updates
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.visualization import generate_text
from ord_schema.visualization import drawing

# pylint: disable=invalid-name,no-member,inconsistent-return-statements
app = flask.Flask(__name__, template_folder='../html')
app.config['ORD_EDITOR_DB'] = os.path.abspath(os.getenv('ORD_EDITOR_DB', 'db'))
app.config['REVIEW_ROOT'] = flask.safe_join(app.config['ORD_EDITOR_DB'],
                                            '.review')
app.config['REVIEW_DATA_ROOT'] = flask.safe_join(app.config['REVIEW_ROOT'],
                                                 'ord-data')


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
    return flask.render_template(
        'datasets.html', names=sorted(names), user_id=flask.g.user_id)


@app.route('/dataset/<name>')
def show_dataset(name):
    """Lists all Reactions contained in the named dataset."""
    dataset = get_dataset(name)
    reactions = []
    for reaction in dataset.reactions:
        reactions.append(reaction.identifiers)
    # Reactions belong to the "review" user are immutable.
    freeze = flask.g.user_id == '8df09572f3c74dbcb6003e2eef8e48fc'
    return flask.render_template(
        'dataset.html', name=name, freeze=freeze, user_id=flask.g.user_id)


@app.route('/dataset/<name>/download')
def download_dataset(name):
    """Returns a pbtxt from the datasets table as an attachment."""
    pbtxt = get_pbtxt(name)
    data = io.BytesIO(pbtxt)
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
        pbtxt = flask.request.get_data()
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, pbtxt, name]) 
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
        pbtxt = ''
        cursor.execute(query, [user_id, name, pbtxt]) 
        flask.g.db.commit()
    return 'ok'


@app.route('/dataset/<name>/delete')
def delete_dataset(name):
    """Removes a Dataset."""
    with flask.g.db.cursor() as cursor:
        query = psycopg2.sql.SQL(
            'DELETE FROM datasets WHERE user_id=%s and name=%s')
        user_id = flask.g.user_id
        cursor.execute(query, [user_id, name])
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
    ord_schema.dataset_templating.generate_dataset.
    """
    # pylint: disable=broad-except
    data = flask.request.get_json(force=True)
    basename, suffix = os.path.splitext(data['spreadsheet_name'])
    spreadsheet_data = io.StringIO(data['spreadsheet_data'].lstrip('ï»¿'))
    dataframe = dataset_templating.read_spreadsheet(spreadsheet_data,
                                                    suffix=suffix)
    dataset = None
    try:
        dataset = dataset_templating.generate_dataset(data['template_string'],
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
    # Reactions belong to the "review" user are immutable.
    freeze = flask.g.user_id == '8df09572f3c74dbcb6003e2eef8e48fc'
    return flask.render_template(
        'reaction.html', name=name, index=index, freeze=freeze,
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


@app.route('/dataset/proto/upload/<file_name>/<token>', methods=['POST'])
def write_upload(file_name, token):
    """Writes the POST body, names it <token>, and maybe updates <file_name>.

    This is part of the upload mechanism. Fields named "bytes_value" can be
    populated only through browser file uploads (e.g. images), and the binary
    data in the files can not be integrated into the Dataset on the client
    because JS can not access the file system. Instead, JS populates the
    bytes_value with a random token and then uploads the file to this endpoint
    with the same token.

    Datasets protos and their bytes_values are reunited in resolve_tokens().

    Args:
        file_name: The .pbtxt containing the Dataset that owns the upload.
        token: The bytes_value placeholder used in pbtxt to reference the
            upload.

    Returns:
        A 200 response when the file was written.
    """
    path = get_path(token)
    with open(path, 'wb') as upload:
        upload.write(flask.request.get_data())
    with lock(file_name):
        dataset = get_dataset(name)
        if resolve_tokens(dataset):
            put_dataset(file_name, dataset)
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


@app.route('/resolve/<identifier_type>', methods=['POST'])
def resolve_compound(identifier_type):
    """Resolve a compound name to a SMILES string."""
    compound_name = flask.request.get_data()
    if not compound_name:
        return ''
    try:
        return flask.jsonify(
            updates.name_resolve(identifier_type, compound_name))
    except ValueError:
        return ''


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
    pbtxt = get_dataset(name)
    local = dataset_pb2.Dataset()
    text_format.Parse(pbtxt, local)
    remote_ascii = text_format.MessageToString(remote)
    local_ascii = text_format.MessageToString(local)
    if remote_ascii != local_ascii:
        app.logger.error(
            pprint.pformat(
                list(
                    difflib.context_diff(local_ascii.splitlines(),
                                         remote_ascii.splitlines()))))
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


@app.before_first_request
def init_submissions():
    """Clones ord-data for use in review mode."""
    # NOTE(kearnes): Use a local git clone so we aren't downloading the
    # uncompressed submission files.
    if not os.path.exists(app.config['REVIEW_DATA_ROOT']):
        subprocess.run([
            'git', 'clone',
            'https://github.com/Open-Reaction-Database/ord-data.git',
            app.config['REVIEW_DATA_ROOT']
        ],
                       check=True)


@app.route('/review')
def show_submissions():
    """Lists pending submissions to ord-data."""
    client = github.Github()
    repo = client.get_repo('Open-Reaction-Database/ord-data')
    pull_requests = {}
    for pull_request in repo.get_pulls():
        datasets = []
        for file in pull_request.get_files():
            if file.filename.endswith('.pbtxt'):
                datasets.append(file.filename[:-6])
        pull_requests[pull_request.number] = (pull_request.title, datasets)
    return flask.render_template('submissions.html',
                                 pull_requests=pull_requests)


@app.route('/review/<pull_request>/<file_name>')
def show_submission(pull_request, file_name):
    """Requests a dataset from a current pull request."""
    client = github.Github()
    repo = client.get_repo('Open-Reaction-Database/ord-data')
    pr = repo.get_pull(int(pull_request))
    with lock('review', user='.review'):
        subprocess.run(['git', 'checkout', 'main'],
                       cwd=app.config['REVIEW_DATA_ROOT'],
                       check=True)
        subprocess.run(['git', 'pull', 'origin', 'main'],
                       cwd=app.config['REVIEW_DATA_ROOT'],
                       check=True)
        subprocess.run([
            'git', 'fetch', 'origin',
            f'pull/{pull_request}/head:editor#{pull_request}'
        ],
                       cwd=app.config['REVIEW_DATA_ROOT'],
                       check=True)
        subprocess.run(['git', 'checkout', f'editor#{pull_request}'],
                       cwd=app.config['REVIEW_DATA_ROOT'],
                       check=True)
        subprocess.run(['git', 'pull', 'origin', f'pull/{pull_request}/head'],
                       cwd=app.config['REVIEW_DATA_ROOT'],
                       check=True)
        destination = flask.safe_join(app.config['REVIEW_ROOT'], pull_request)
        os.makedirs(destination, exist_ok=True)
        for file in pr.get_files():
            if file.filename.endswith('.pbtxt'):
                shutil.copy2(
                    flask.safe_join(app.config['REVIEW_DATA_ROOT'],
                                    file.filename), destination)
    parts = list(urllib.parse.urlparse(f'/dataset/{file_name}'))
    parts[4] = urllib.parse.urlencode({'user': f'.review/{pull_request}'})
    url = urllib.parse.urlunparse(parts)
    return flask.redirect(url, code=307)


@app.after_request
def prevent_caching(response):
    """Prevents caching any of this app's resources on the client."""
    response.headers['Cache-Control'] = 'no-cache,no-store,must-revalidate'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    response.headers['Cache-Control'] = 'public, max-age=0'
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


@contextlib.contextmanager
def lock(file_name, user=None):
    """Blocks until an exclusive lock on the named file is obtained.

    This is part of the upload mechanism. Byte_value fields are populated
    through browser file uploads and so must be merged with Dataset protos on
    the server. The file uploads and the Dataset proto arrive asynchronously in
    concurrent connections. Locks ensure that only one request at a time
    accesses a Datset's .pbtxt file.

    Args:
        file_name: Name of the file to pass to the fcntl system call.
        user: The current user; used to determine the lock file path.

    Yields:
        The locked file descriptor.
    """
    path = get_path(file_name, user=user, suffix='.lock')
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


def get_path(file_name, user=None, suffix='.pbtxt'):
    """Returns a safe path in the editor filesystem.

    Uses flask.safe_join to check for directory traversal attacks.

    Args:
        file_name: Text filename.
        user: The current user. If None, defaults to flask.g.user.
        suffix: Text filename suffix. Defaults to '.pbtxt'.

    Returns:
        Path to the requested file.
    """
    if not user:
        user = flask.g.user
    return flask.safe_join(app.config['ORD_EDITOR_DB'], user,
                           f'{file_name}{suffix}')


def get_user_path():
    """Checks that a username results in a valid path in the editor filesystem.

    Uses flask.safe_join to check for directory traversal attacks.

    Returns:
        Path to the user directory.
    """
    return flask.safe_join(app.config['ORD_EDITOR_DB'], flask.g.user)


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
    """Presents a form to set the access token."""
    user_id = flask.request.args.get('user_id', '')
    return flask.render_template('login.html', user_id=user_id)


@app.route('/authenticate', methods=['POST'])
def authenticate():
    """Issue a new access token for a given user ID."""
    user_id = flask.request.form.get('user_id')
    if user_id is None:
        return flask.redirect('/login')
    with psycopg2.connect(dbname='editor', port=5430) as db:
        with db.cursor() as cursor:
            query = psycopg2.sql.SQL(
                f"SELECT user_id FROM users WHERE user_id=%s")
            cursor.execute(query, [user_id])
            if cursor.rowcount == 0:
                return flask.redirect('/login')
            query = psycopg2.sql.SQL(f"INSERT INTO logins VALUES (%s, %s, %s)")
            access_token = uuid.uuid4().hex
            timestamp = int(time.time())
            cursor.execute(query, [access_token, user_id, timestamp])
            db.commit()
            response = flask.redirect('/')
            # Expires in a year.
            response.set_cookie('Access-Token', access_token, max_age=31536000)
            return response
    return flask.redirect('/login')


@app.before_request
def init_user():
    """Connects to the DB and authenticates the user."""
    if flask.request.path in ('/login', '/authenticate'):
        return
    access_token = flask.request.cookies.get('Access-Token')
    if access_token is None:
        return flask.redirect(f'/login')
    db = psycopg2.connect(dbname='editor', port=5430)
    cursor = db.cursor()
    query = psycopg2.sql.SQL(
        f"SELECT user_id FROM logins WHERE access_token=%s")
    cursor.execute(query, [access_token])
    if cursor.rowcount == 0:
        return flask.redirect(f'/login')
    user_id = cursor.fetchone()[0]
    flask.g.db = db
    flask.g.user_id = user_id


@app.route('/logout')
def logout():
    """Clear the access token and redirect to /login."""
    response = flask.redirect('/login')
    response.set_cookie('Access-Token', '', expires=0)
    return response
