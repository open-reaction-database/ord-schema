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
import setuptools


with open("README.md") as f:
    long_description = f.read()

setuptools.setup(
    name="ord-schema",
    version="0.3.37",
    description="Schema for the Open Reaction Database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Open-Reaction-Database/ord-schema",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    package_data={
        "ord_schema.proto": [
            "dataset.proto",
            "reaction.proto",
            "test.proto",
        ],
    },
    python_requires=">=3.9",
    install_requires=[
        "docopt>=0.6.2",
        "flask>=1.1.2",
        "inflection>=0.5.1",
        "jinja2>=2.0.0",
        "joblib>=1.0.0",
        "numpy>=1.18.1",
        "openpyxl>=3.0.5",
        "pandas>=1.0.4",
        "protobuf<3.20,>=3.13.0",
        "psycopg2>=2.8.5",
        "pygithub>=1.51",
        "python-dateutil>=1.10.0",
        "rdkit>=2021.9.5",
        "sqlalchemy>=1.4.39",
        "xlrd>=2.0.1",
        "xlwt>=1.3.0",
    ],
    extras_require={
        "docs": [
            "ipython>=7.18.1",
            "Pygments>=2.7.2",
            "sphinx>=3.3.1",
            "sphinx-rtd-theme>=0.5.0",
            "sphinx-tabs>=1.3.0",
        ],
        "examples": [
            "glob2>=0.7",
            "matplotlib>=3.3.4",
            "scikit-learn>=0.24.1",
            "tensorflow>=2.4.1",
            "tqdm>=4.61.2",
            "wget>=3.2",
        ],
        "tests": [
            "black[jupyter]>=22.3.0",
            "coverage>=5.2.1",
            "pylint>=2.13.9",
            "pytest>=7.1.1",
            "pytest-cov>=3.0.0",
            "pytest-xdist>=3.0.2",
            "pytype>=2022.5.19",
            "testing-postgresql>=1.3.0",
            "treon>=0.1.3",
        ],
    },
)
