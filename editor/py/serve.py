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
import fcntl
import os
import re
import json
import urllib
import uuid
import io
import tempfile
import urllib

import flask

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema import message_helpers
from ord_schema import validations
from ord_schema import dataset_templating
from ord_schema import updates
from ord_schema.visualization import generate_text
from ord_schema.visualization import drawing

from google.protobuf import text_format

# pylint: disable=invalid-name,no-member,inconsistent-return-statements
app = flask.Flask(__name__, template_folder='../html')


@app.route('/')
def show_root():
    """The root path redirects to the "datasets" view."""
    return flask.redirect('/datasets')


@app.route('/datasets')
def show_datasets():
    """Lists all .pbtxt files in the db/ directory."""
    redirect = ensure_user()
    if redirect:
        return redirect
    base_names = []
    for name in os.listdir(f'db/{flask.g.user}'):
        match = re.match(r'(.*).pbtxt', name)
        if match is not None:
            base_names.append(match.group(1))
    return flask.render_template('datasets.html', file_names=sorted(base_names))


@app.route('/dataset/<file_name>')
def show_dataset(file_name):
    """Lists all Reactions contained in the specified .pbtxt file."""
    redirect = ensure_user()
    if redirect:
        return redirect
    dataset = get_dataset(file_name)
    if dataset is None:
        flask.abort(404)
    reactions = []
    for reaction in dataset.reactions:
        reactions.append(reaction.identifiers)
    return flask.render_template('dataset.html', file_name=file_name)


@app.route('/dataset/<file_name>/download')
def download_dataset(file_name):
    """Returns a raw .pbtxt file from the db/ directory as an attachment."""
    path = f'db/{flask.g.user}/{file_name}.pbtxt'
    if not os.path.isfile(path):
        flask.abort(404)
    pbtxt = open(path, 'rb')
    return flask.send_file(pbtxt,
                           mimetype='application/protobuf',
                           as_attachment=True,
                           attachment_filename='%s.pbtxt' % file_name)


@app.route('/dataset/<file_name>/upload', methods=['POST'])
def upload_dataset(file_name):
    """Writes the request body to the db/ directory without validation."""
    path = f'db/{flask.g.user}/{file_name}.pbtxt'
    with open(path, 'wb') as upload:
        upload.write(flask.request.get_data())
    return 'ok'


@app.route('/dataset/<file_name>/new', methods=['POST'])
def new_dataset(file_name):
    """Creates a new dataset in the db/ directory."""
    path = f'db/{flask.g.user}/{file_name}.pbtxt'
    if os.path.isfile(path):
        flask.abort(404)
    with open(path, 'wb') as upload:
        upload.write(b'\n')
    return 'ok'


@app.route('/dataset/enumerate', methods=['POST'])
def enumerate_dataset():
    """Creates a new dataset in the db/ directory based on a template reaction
    pbtxt and a spreadsheet.

    Three pieces of information are expected to be POSTed in a json object:
        spreadsheet_name: the original filename of the uploaded spreadsheet.
        spreadsheet_data: a string containing the contents of the spreadsheet.
        template_string: a string containing a text-formatted Reaction proto,
            i.e., the contents of a pbtxt file.
    A new dataset is created from the template and spreadsheet using
    ord_schema.dataset_templating.generate_dataset."""
    data = flask.request.get_json(force=True)
    basename, suffix = os.path.splitext(data['spreadsheet_name'])
    with tempfile.NamedTemporaryFile(mode='w', suffix=suffix) as f:
        f.write(data['spreadsheet_data'].lstrip('ï»¿'))
        f.seek(0)
        dataframe = dataset_templating.read_spreadsheet(f.name)
    try:
        dataset = dataset_templating.generate_dataset(data['template_string'],
                                                      dataframe,
                                                      validate=False)
    except ValueError as e:
        flask.abort(flask.make_response(str(e), 400))
    message_helpers.write_message(
        dataset, f'db/{flask.g.user}/{basename}_dataset.pbtxt')
    return 'ok'


@app.route('/dataset/<file_name>/reaction/<index>')
def show_reaction(file_name, index):
    """Render the page representing a single Reaction."""
    redirect = ensure_user()
    if redirect:
        return redirect
    dataset = get_dataset(file_name)
    if dataset is None:
        flask.abort(404)
    try:
        index = int(index)
    except ValueError:
        flask.abort(404)
    if len(dataset.reactions) <= index:
        flask.abort(404)
    return flask.render_template('reaction.html',
                                 file_name=file_name,
                                 index=index)


@app.route('/reaction/download', methods=['POST'])
def download_reaction():
    """Returns a raw .pbtxt file taken from POST as an attachment."""
    reaction = reaction_pb2.Reaction()
    reaction.ParseFromString(flask.request.get_data())
    mem = io.BytesIO()
    mem.write(text_format.MessageToString(reaction).encode())
    mem.seek(0)
    return flask.send_file(mem,
                           mimetype='application/protobuf',
                           as_attachment=True,
                           attachment_filename='reaction.pbtxt')


@app.route('/dataset/<file_name>/new/reaction')
def new_reaction(file_name):
    """Adds a new Reaction to the given Dataset and view the Dataset."""
    dataset = get_dataset(file_name)
    if dataset is None:
        flask.abort(404)
    reaction = reaction_pb2.Reaction()
    dataset.reactions.append(reaction)
    put_dataset(file_name, dataset)
    return flask.redirect('/dataset/%s' % file_name)


@app.route('/dataset/<file_name>/clone/<index>')
def clone_reaction(file_name, index):
    """Copies a specific Reaction to the Dataset and view the Reaction."""
    dataset = get_dataset(file_name)
    if dataset is None:
        flask.abort(404)
    try:
        index = int(index)
    except ValueError:
        flask.abort(404)
    if len(dataset.reactions) <= index:
        flask.abort(404)
    original = dataset.reactions[index]
    clone = reaction_pb2.Reaction()
    clone.MergeFrom(original)
    dataset.reactions.append(clone)
    index = len(dataset.reactions) - 1
    put_dataset(file_name, dataset)
    return flask.redirect('/dataset/%s/reaction/%s' % (file_name, index))


@app.route('/dataset/<file_name>/delete/reaction/<index>')
def delete_reaction(file_name, index):
    """Removes a specific Reaction from the Dataset and view the Datset."""
    dataset = get_dataset(file_name)
    if dataset is None:
        flask.abort(404)
    try:
        index = int(index)
    except ValueError:
        flask.abort(404)
    if len(dataset.reactions) <= index:
        flask.abort(404)
    del dataset.reactions[index]
    put_dataset(file_name, dataset)
    return flask.redirect('/dataset/%s' % file_name)


@app.route('/dataset/<file_name>/delete/reaction_id/<reaction_id>')
def delete_reaction_id(file_name, reaction_id):
    """Removes a Reaction reference from the Dataset and view the Dataset."""
    dataset = get_dataset(file_name)
    if dataset is None:
        flask.abort(404)
    if reaction_id in dataset.reaction_ids:
        dataset.reaction_ids.remove(reaction_id)
        put_dataset(file_name, dataset)
        return flask.redirect('/dataset/%s' % file_name)
    flask.abort(404)


@app.route('/dataset/<file_name>/delete/reaction_id/')
def delete_reaction_id_blank(file_name):
    """Removes the first empty Reaction reference from the Dataset."""
    return delete_reaction_id(file_name, '')


@app.route('/dataset/proto/read/<file_name>')
def read_dataset(file_name):
    """Returns a Dataset as a serialized protobuf."""
    dataset = dataset_pb2.Dataset()
    with open(f'db/{flask.g.user}/{file_name}.pbtxt', 'rb') as pbtxt:
        text_format.Parse(pbtxt.read(), dataset)
    bites = dataset.SerializeToString(deterministic=True)
    response = flask.make_response(bites)
    response.headers.set('Content-Type', 'application/protobuf')
    return response


@app.route('/dataset/proto/write/<file_name>', methods=['POST'])
def write_dataset(file_name):
    """Receives a serialized Dataset protobuf and write it to a file."""
    dataset = dataset_pb2.Dataset()
    dataset.ParseFromString(flask.request.get_data())
    resolve_tokens(dataset)
    put_dataset(file_name, dataset)
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
    with open(f'db/{flask.g.user}/{token}', 'wb') as upload:
        upload.write(flask.request.get_data())
    with lock(file_name):
        dataset = get_dataset(file_name)
        if dataset is not None and resolve_tokens(dataset):
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
        token: A placeholder name for the upload in the db/ directory.

    Returns:
        The POST body from the request, after passing through a file.
    """
    path = f'db/{flask.g.user}/{token}'
    with open(path, 'wb') as upload:
        upload.write(flask.request.get_data())
    pbtxt = open(path, 'rb')
    return flask.send_file(pbtxt,
                           mimetype='application/protobuf',
                           as_attachment=True,
                           attachment_filename=token)


@app.route('/dataset/proto/validate/<message_name>', methods=['POST'])
def validate_reaction(message_name):
    """Receives a serialized Reaction protobuf and runs validations."""
    message = message_helpers.create_message(message_name)
    message.ParseFromString(flask.request.get_data())
    errors = validations.validate_message(message, raise_on_error=False)
    return json.dumps(errors)


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
    """Receives a serialized Compound message and returns base64-encoded png
    data corresponding to a line drawing of the molecule."""
    compound = reaction_pb2.Compound()
    compound.ParseFromString(flask.request.get_data())
    try:
        mol = message_helpers.mol_from_compound(compound)
        png_data = drawing.mol_to_png(mol)
        return flask.jsonify(png_data)
    except ValueError:
        return ''


@app.route('/dataset/proto/compare/<file_name>', methods=['POST'])
def compare(file_name):
    """For testing, compares a POST body to a Dataset in the db/ directory.

    Returns HTTP status 200 if their pbtxt representations are equal as strings
    and 409 if they are not. See "make test".

    Args:
        file_name: The .pbtxt file in the db/ directory to compare.

    Returns:
        HTTP status 200 for a match and 409 if there is a difference.
    """
    remote = dataset_pb2.Dataset()
    remote.ParseFromString(flask.request.get_data())
    with open(f'db/{flask.g.user}/{file_name}.pbtxt', 'rb') as pbtxt:
        local = dataset_pb2.Dataset()
        text_format.Parse(pbtxt.read(), local)
    remote_ascii = text_format.MessageToString(remote)
    local_ascii = text_format.MessageToString(local)
    if remote_ascii != local_ascii:
        return 'differs', 409  # "Conflict"
    return 'equals'


@app.route('/js/<script>')
def js(script):
    """Accesses any built JS file by name from the Closure output directory."""
    return flask.send_file('../gen/js/ord/%s' % script)


@app.route('/css/<sheet>')
def css(sheet):
    """Accesses any CSS file by name."""
    return flask.send_file('../css/%s' % sheet)


@app.route('/ketcher/iframe')
def ketcher_iframe():
    """Accesses a website serving Ketcher."""
    return flask.render_template('ketcher_iframe.html')


@app.route('/ketcher/info')
def indigo():
    """Dummy indigo endpoint to prevent 404 errors."""
    return ('', 204)


@app.route('/ketcher/<path:file>')
def ketcher(file):
    """Accesses any built Ketcher file by name."""
    return flask.send_file('../ketcher/dist/%s' % file)


@app.route('/dataset/deps.js')
@app.route('/dataset/<file_name>/deps.js')
@app.route('/dataset/<file_name>/reaction/deps.js')
def deps(file_name=None):
    """Returns empty for deps table requests since this app doesn't use them."""
    del file_name
    return ''


@app.route('/ketcher/molfile', methods=['POST'])
def get_molfile():
    """Retrieves a POSTed Compound message string and returns a MolFile."""
    compound = reaction_pb2.Compound()
    compound.ParseFromString(flask.request.get_data())
    try:
        molblock = message_helpers.molblock_from_compound(compound)
        return flask.jsonify(molblock)
    except ValueError as e:
        return 'no existing structural identifier', 404


@app.after_request
def prevent_caching(response):
    """Prevents caching any of this app's resources on the client."""
    response.headers['Cache-Control'] = 'no-cache,no-store,must-revalidate'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    response.headers['Cache-Control'] = 'public, max-age=0'
    return response


def get_dataset(file_name):
    """Reads a .pbtxt file from the db/ directory and parse it."""
    with lock(file_name):
        if (f'{file_name}.pbtxt') not in os.listdir(f'db/{flask.g.user}'):
            return None
        dataset = dataset_pb2.Dataset()
        with open(f'db/{flask.g.user}/{file_name}.pbtxt', 'rb') as pbtxt:
            text_format.Parse(pbtxt.read(), dataset)
        return dataset


def put_dataset(file_name, dataset):
    """Write a proto to a .pbtxt file in the db/ directory."""
    with lock(file_name):
        with open(f'db/{flask.g.user}/{file_name}.pbtxt', 'wb') as pbtxt:
            string = text_format.MessageToString(dataset, as_utf8=True)
            pbtxt.write(string.encode('utf8'))


@contextlib.contextmanager
def lock(file_name):
    """Blocks until an exclusive lock on the named file is obtained.

    This is part of the upload mechanism. Byte_value fields are populated
    through browser file uploads and so must be merged with Dataset protos on
    the server. The file uploads and the Dataset proto arrive asynchronously in
    concurrent connections. Locks ensure that only one request at a time
    accesses a Datset's .pbtxt file.

    Args:
        file_name: Name of the file in db/ to pass to the fcntl system call.

    Yields:
        The locked file descriptor.
    """
    path = f'db/{flask.g.user}/{file_name}.lock'
    lock_file = open(path, 'w')
    fcntl.lockf(lock_file, fcntl.LOCK_EX)
    try:
        yield lock_file
    finally:
        fcntl.lockf(lock_file, fcntl.LOCK_UN)


def resolve_tokens(proto):
    """Fills in bytes_value fields using client-generated placeholder tokens.

    This is part of the upload mechanism. It acts by recursion on the tree
    structure of the proto, hunting for fields named "bytes_value" and
    comparing the fields' values against uploaded files in the db/ directory.
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
                    path = f'db/{flask.g.user}/{token}'
                    if os.path.isfile(path):
                        proto.bytes_value = open(path, 'rb').read()
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


@app.before_request
def init_user():
    """Set the user cookie and initialize the db/{user} directory."""
    user = flask.request.cookies.get('ord-editor-user')
    if user is None:
        user = flask.request.args.get('user') or next_user()
        return redirect_to_user(user)
    flask.g.user = user
    if not os.path.isdir(f'db/{user}'):
        os.mkdir(f'db/{user}')
    if not os.listdir(f'db/{user}'):
        dataset = dataset_pb2.Dataset()
        put_dataset('dataset', dataset)


def ensure_user():
    """On page requests only, a user param in the URL overrides the cookie."""
    url_user = flask.request.args.get('user')
    cookie_user = flask.request.cookies.get('ord-editor-user')
    # If there is no user in the URL, set it from the cookie.
    if url_user is None:
        return redirect_to_user(cookie_user)
    # If the cookie and the URL disagree, the URL wins.
    if url_user != cookie_user:
        return redirect_to_user(url_user)
    # If the URL and the cookie are consistent, then do nothing.


def redirect_to_user(user):
    """Set the user cookie and return a redirect to the updated URL."""
    url = url_for_user(user)
    response = flask.redirect(url)
    # Expires in a year.
    response.set_cookie('ord-editor-user', user, max_age=31536000)
    return response


def url_for_user(user):
    """Replace the user in the request URL and return the updated URL."""
    parts = list(urllib.parse.urlparse(flask.request.url))
    parts[4] = urllib.parse.urlencode({'user': user})
    return urllib.parse.urlunparse(parts)


def next_user():
    """Return a user identifier that is not present in the db/ directory."""
    user = uuid.uuid4().hex
    while os.path.isdir(f'db/{user}'):
        user = uuid.uuid4().hex
    return user
