#! /usr/bin/env python
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
"""Installer script."""

from distutils.command import build_py
from distutils import spawn
import glob
import subprocess

import setuptools


class BuildPyCommand(build_py.build_py):
    """Command that generates Python interface for protocol buffers."""

    def run(self):
        """Runs the protocol buffer compiler."""
        protoc = spawn.find_executable('protoc')
        if not protoc:
            raise RuntimeError('cannot find protoc')
        for source in glob.glob('ord_schema/proto/*.proto'):
            protoc_command = [
                protoc,
                # https://github.com/protocolbuffers/protobuf/blob/master/docs/field_presence.md
                '--experimental_allow_proto3_optional',
                '--python_out=.',
                source
            ]
            self.announce(f'running {protoc_command}')
            subprocess.check_call(protoc_command)
        # build_py.build_py is an old-style class, so super() doesn't work.
        build_py.build_py.run(self)

with open ("requirements.txt", "r") as f:
    requirements = f.readlines()

setuptools.setup(name='ord-schema',
                 description='Schema for the Open Reaction Database',
                 url='https://github.com/Open-Reaction-Database/ord-schema',
                 license='Apache License, Version 2.0',
                 packages=setuptools.find_packages(),
                 package_data={
                     'ord_schema.visualization': [
                         'reaction.html',
                         'template.html',
                         'template.txt',
                     ],
                 },
                 install_requires=requirements,
                 cmdclass={'build_py': BuildPyCommand}
)
