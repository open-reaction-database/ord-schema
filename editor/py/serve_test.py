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
"""Tests for editor.py.serve."""

import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import parameterized

import serve


class ServeTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        serve.app.config['ORD_EDITOR_DB'] = self.test_subdirectory
        serve.app.testing = True

    @parameterized.parameters([
        ('dataset', 200),
        ('../dataset', 404),  # flask.safe_join returns 404 on failure.
        ('other', 404),
    ])
    def test_check_path(self, file_name, expected):
        with serve.app.test_client() as client:
            response = client.get(f'/dataset/{file_name}/download',
                                  follow_redirects=True)
            self.assertEqual(response.status_code, expected)


if __name__ == '__main__':
    absltest.main()
