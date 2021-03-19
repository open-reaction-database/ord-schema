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
            subprocess.check_call(['pip', 'install', 'protoc-wheel-0>=3.14.0'])
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


if __name__ == '__main__':
    setuptools.setup(
        name='ord-schema',
        description='Schema for the Open Reaction Database',
        url='https://github.com/Open-Reaction-Database/ord-schema',
        license='Apache License, Version 2.0',
        packages=setuptools.find_packages(),
        package_data={
            'ord_schema.visualization': ['template.html', 'template.txt'],
        },
        install_requires=[
            'absl-py>=0.9.0',
            'flask>=1.1.2',
            'numpy>=1.18.1',
            'openpyxl>=3.0.5',
            'pandas>=1.0.4',
            'protobuf>=3.13.0',
            'protoc-wheel-0>=3.14.0',
            'pygithub>=1.51',
            'python-dateutil>=1.10.0',
            'jinja2>=2.0.0',
            'xlrd<2.0.0',
            'xlwt>=1.3.0',
            'joblib>=1.0.0',
        ],
        cmdclass={'build_py': BuildPyCommand})
