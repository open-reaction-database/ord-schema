#! /usr/bin/env python

from distutils.command import build_py
from distutils import spawn
import glob
import setuptools
import subprocess


class BuildPyCommand(build_py.build_py):
    """Command that generates Python interface for protocol buffers."""

    def run(self):
        """Runs the protocol buffer compiler."""
        protoc = spawn.find_executable('protoc')
        for source in glob.glob('proto/*.proto'):
            protoc_command = [protoc, '--python_out=ord_schema', source]
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
        install_requires=['absl-py>=0.9.0', 'protobuf>=3.11.0'],
        version="0.1",
        packages=setuptools.find_packages(),
        cmdclass={'build_py': BuildPyCommand})
