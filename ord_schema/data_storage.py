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
"""Utilities for offloading large `Data` values.

The original data is replaced with a URL pointing to the written data. This
keeps the protocol buffers small and avoids data capacity issues on GitHub.
"""

import hashlib
import os
import sys
import tempfile
import urllib

from absl import logging

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2

# Prefix for filenames that store reaction_pb2.Data values.
DATA_PREFIX = 'ord_data-'
DATA_REPO = 'ord-data'
DATA_URL_PREFIX = (
    f'https://github.com/Open-Reaction-Database/{DATA_REPO}/tree/main')


def write_data(message, dirname, min_size=0.0, max_size=1.0):
    """Writes a Data value to a file.

    If a value is a URL or is smaller than `min_size`, it is left unchanged.

    Args:
        message: Data message.
        dirname: Text output directory.
        min_size: Float minimum size of data before it will be written (in MB).
        max_size: Float maximum size of data to write (in MB).

    Returns:
        filename: Text filename containing the written data, or None if the
            value is a URL or is smaller than `min_size`.
        value_size: Float value size in MB, or None if no value was written.

    Raises:
        ValueError: if there is no value defined in `message` or if the value is
            larger than `max_size`.
    """
    if min_size > max_size:
        raise ValueError('min_size must be less than or equal to max_size')
    kind = message.WhichOneof('kind')
    if kind == 'string_value':
        value = message.string_value.encode()  # Convert to bytes.
    elif kind == 'float_value':
        return None, None
    elif kind == 'bytes_value':
        value = message.bytes_value
    elif kind == 'url':
        return None, None
    else:
        raise ValueError('no value to write')
    value_size = sys.getsizeof(value) / 1e6
    if value_size < min_size:
        return None, None
    if value_size > max_size:
        raise ValueError(
            f'value is larger than max_size ({value_size} vs {max_size}')
    value_hash = hashlib.sha256(value).hexdigest()
    suffix = message.format or 'txt'
    basename = f'{DATA_PREFIX}{value_hash}.{suffix}'
    filename = os.path.join(dirname, basename)
    with open(filename, 'wb') as f:
        f.write(value)
    return filename, value_size


def extract_data(message, root, min_size=0.0, max_size=1.0):
    """Replaces large Data values with pointers to offloaded data.

    Git LFS (https://git-lfs.github.com/) is convenient because it lives in the
    same repo as the associated Reaction records. However, it is not possible to
    get a permanent URL for the uploaded data because it is only committed to
    the PR branch. We have (at least) these options:

        1. Modify the URL just before or after the PR is merged to point to the
           correct branch.
        2. Modify the URL to point to its eventual destination (in the `main`
           branch) and deal with broken links during submission review, or
        3. Use relative paths (relative to the repository root). This means that
           users will have to traverse the repo manually to access referenced
           data instead of simply following a URL.
        4. Merge the data immediately in another repo so the URL is permanent.

    I think (2) is the best option because it yields URLs that will eventually
    work and it is simpler than (1). I don't like option (4) because it requires
    data to be committed and merged before review.

    Args:
        message: Protocol buffer message.
        root: Text root of the repository.
        min_size: Float minimum size of data before it will be written (in MB).
        max_size: Float maximum size of data to write (in MB).

    Returns:
        List of text filenames; the generated Data files.
    """
    dirname = tempfile.mkdtemp()
    data_messages = message_helpers.find_submessages(message,
                                                     reaction_pb2.Data)
    filenames = []
    for data_message in data_messages:
        data_filename, data_size = write_data(data_message,
                                              dirname,
                                              min_size=min_size,
                                              max_size=max_size)
        if data_filename:
            basename = os.path.basename(data_filename)
            output_filename = message_helpers.id_filename(basename)
            with_root = os.path.join(root, output_filename)
            os.makedirs(os.path.dirname(with_root), exist_ok=True)
            os.rename(data_filename, with_root)
            data_message.url = urllib.parse.urljoin(DATA_URL_PREFIX,
                                                    output_filename)
            logging.info('Created Data file (%g MB): %s', data_size, with_root)
            filenames.append(with_root)
    return filenames
