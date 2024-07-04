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
    version="0.3.81",
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
    package_data={"ord_schema.proto": ["*.pyi"]},
    python_requires=">=3.10",
    install_requires=[
        "docop",
        "flask",
        "inflection",
        "jinja2",
        "joblib",
        "numpy",
        "openpyxl",
        "pandas",
        "protobuf",
        "psycopg2-binary",
        "pygithub",
        "python-dateutil",
        "rdkit",
        "sqlalchemy",
    ],
    extras_require={
        "docs": [
            "ipython",
            "Pygments",
            "Sphinx",
            "sphinx-rtd-theme",
            "sphinx-tabs",
        ],
        "examples": [
            "glob2",
            "matplotlib",
            "scikit-learn",
            "tensorflow",
            "tqdm",
            "wget",
        ],
        "tests": [
            "black[jupyter]",
            "coverage",
            "pylint",
            "pytest",
            "pytest-cov",
            "pytest-xdist",
            "pytype",
            "testing-postgresql",
            "treon",
        ],
    },
)
