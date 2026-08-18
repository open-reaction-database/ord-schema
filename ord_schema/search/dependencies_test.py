# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Checks that the search extra declares everything the search path imports.

A development checkout installs every group, so an import reaching a distribution the
extra does not name still passes the rest of the suite; the failure lands on whoever
installs ``ord-schema[search]`` and cannot import the subpackage. This walks the import
graph the way that install does and asks the metadata which distribution provides each
module.
"""

import ast
import importlib.metadata
import pathlib
import re
import sys
import tomllib

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[2]
_ENTRY_POINTS = ("ord_schema.search.execute", "ord_schema.search.schema")


def _canonical(name: str) -> str:
    """Returns the PEP 503 normalized form of a distribution name."""
    return re.sub(r"[-_.]+", "-", name).lower()


def _declared(extra: str) -> set[str]:
    """Returns the distributions an install of ``ord-schema[extra]`` brings in."""
    with (_ROOT / "pyproject.toml").open("rb") as f:
        pyproject = tomllib.load(f)["project"]
    requirements = pyproject["dependencies"] + pyproject["optional-dependencies"][extra]
    # Everything up to the first version specifier, extra, or environment marker.
    return {_canonical(re.split(r"[<>=!~\[; ]", line)[0]) for line in requirements}


def _module_path(module: str) -> pathlib.Path | None:
    """Returns the file holding ``module``, or None if no file does.

    A name imported from a module -- a function, a class -- resolves to no file, which
    is the same answer as a module this package does not hold.
    """
    base = _ROOT / pathlib.Path(*module.split("."))
    for candidate in (base.with_suffix(".py"), base / "__init__.py"):
        if candidate.exists():
            return candidate
    return None


def _imports(path: pathlib.Path) -> set[str]:
    """Returns the modules ``path`` imports, as written."""
    found = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            # ``from a.b import c`` reaches a.b, and c when c is itself a module.
            found.add(node.module)
            found.update(f"{node.module}.{alias.name}" for alias in node.names)
    return found


def _reachable() -> set[str]:
    """Returns every module reachable from the search entry points."""
    seen, pending, third_party = set(), list(_ENTRY_POINTS), set()
    while pending:
        module = pending.pop()
        if module in seen:
            continue
        seen.add(module)
        path = _module_path(module)
        if path is None:
            # Anything this package does not hold is somebody else's to provide, except
            # the standard library, which no install declares.
            root = module.split(".")[0]
            if root not in sys.stdlib_module_names and root != "ord_schema":
                third_party.add(root)
            continue
        pending.extend(_imports(path))
    return third_party


def test_the_search_extra_declares_what_the_search_path_imports():
    declared = _declared("search")
    provided = importlib.metadata.packages_distributions()
    undeclared = {}
    for module in sorted(_reachable()):
        distributions = {_canonical(name) for name in provided.get(module, [])}
        if not distributions:
            pytest.fail(
                f"{module} is imported but no installed distribution provides it"
            )
        if not distributions & declared:
            undeclared[module] = sorted(distributions)
    assert not undeclared, (
        "reachable from ord_schema.search but not declared by the search extra: "
        + ", ".join(
            f"{module} ({'/'.join(names)})" for module, names in undeclared.items()
        )
    )


def test_the_walk_reaches_the_artifacts_the_search_path_reads():
    # A walk that stopped at the subpackage boundary would have missed the import that
    # motivated this test, which is in ord_schema.artifacts.pivot.
    assert "inflection" in _reachable()
