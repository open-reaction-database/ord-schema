#! /usr/bin/env python

import os

from distutils.command import build_py
from distutils import core
from distutils import spawn
import glob
import subprocess


class BuildPyCommand(build_py.build_py):
    """Command that generates Python interface for protocol buffers."""

    def run(self):
        """Runs the protocol buffer compiler."""
        protoc = spawn.find_executable('protoc')
        if not protoc:
            raise RuntimeError('cannot find protoc')
        for source in glob.glob('proto/*.proto'):
            protoc_command = [
                protoc,
                '--proto_path=..',
                '--python_out=.',
                os.path.join('ord-schema', source)
            ]
            self.announce(f'running {protoc_command}')
            subprocess.check_call(protoc_command)
        # build_py.build_py is an old-style class, so super() doesn't work.
        build_py.build_py.run(self)


if __name__ == '__main__':
    core.setup(
        name='ord-schema',
        description='Schema for the Open Reaction Database',
        url='https://github.com/Open-Reaction-Database/ord-schema',
        license='Apache License, Version 2.0',
        version="0.1",
        packages=setuptools.find_packages(),
        cmdclass={'build_py': BuildPyCommand})
