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

from absl.testing import absltest
from absl.testing import parameterized

import serve


class ServeTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        self.root = '/test'
        self.user = 'user'

    @parameterized.parameters(
        ('/test/user/foo.txt', '/test/user/foo.txt'),
        ('/test/user/bar/../foo.txt', '/test/user/foo.txt'))
    def test_check_path(self, path, expected):
        self.assertEqual(serve.check_path(path, root=self.root, user=self.user),
                         expected)

    @parameterized.parameters([
        'foo',
        '/test/../foo.txt',
        '/test/other-user/foo.txt',
        '/test/user/../other-user/foo.txt',
    ])
    def test_check_path_raises(self, path):
        with self.assertRaisesRegex(ValueError, 'path is not safe'):
            serve.check_path(path, root=self.root, user=self.user)


if __name__ == '__main__':
    absltest.main()
