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

"""Checks that each install profile declares everything its modules import.

A development checkout installs every group, so a module reaching a distribution its own
install does not name still passes the rest of the suite; the failure lands on whoever
installs ``ord-schema`` or one extra of it and cannot import the subpackage. This walks
the import graph the way each of those installs does and asks the metadata which
distribution provides each module reached.
"""

import ast
import importlib.metadata
import pathlib
import re
import sys
import tomllib

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACKAGE = "ord_schema"

# What one install holds: the extras it names, and the subpackages those extras make
# importable. Every module belongs to one profile, the longest matching root winning,
# so a base install owns whatever no extra claims. A walk crosses freely into anything
# else the package holds, which is the point: a base module reaching into
# ord_schema.orm shows up here as sqlalchemy going undeclared.
_ROOTS: dict[str, tuple[str, ...]] = {
    "base": (),
    # The search subpackage reads through the artifacts, whose imports the search extra
    # therefore has to carry.
    "search": ("ord_schema.search", "ord_schema.artifacts"),
    "nl": ("ord_schema.search.nl",),
    "orm": ("ord_schema.orm",),
    "huggingface": ("ord_schema.huggingface",),
}
_EXTRAS: dict[str, tuple[str, ...]] = {
    "base": (),
    "search": ("search",),
    # The nl module reads through the search subpackage, so its install is both.
    "nl": ("search", "nl"),
    "orm": ("orm",),
    "huggingface": ("huggingface",),
}
_PROFILES = tuple(pytest.param(name, id=name) for name in _EXTRAS)


def _canonical(name: str) -> str:
    """Returns the PEP 503 normalized form of a distribution name."""
    return re.sub(r"[-_.]+", "-", name).lower()


def _declared(extras: tuple[str, ...]) -> set[str]:
    """Returns the distributions an install of ``ord-schema[extras]`` brings in.

    Args:
        extras: The optional-dependency groups named, if any.

    Returns:
        Canonical distribution names, the base dependencies included.
    """
    with (_ROOT / "pyproject.toml").open("rb") as f:
        project = tomllib.load(f)["project"]
    requirements = list(project["dependencies"])
    for extra in extras:
        requirements.extend(project["optional-dependencies"][extra])
    # Everything up to the first version specifier, extra, or environment marker.
    return {_canonical(re.split(r"[<>=!~\[; ]", line)[0]) for line in requirements}


def _is_installed(distribution: str) -> bool:
    """Returns whether ``distribution`` is installed in this environment."""
    try:
        importlib.metadata.distribution(distribution)
    except importlib.metadata.PackageNotFoundError:
        return False
    return True


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


def _modules_under(root: str) -> set[str]:
    """Returns the importable modules ``root`` holds, itself included.

    Args:
        root: A dotted module or package name inside this package.

    Returns:
        Dotted names, without the tests and fixtures an install does not run.
    """
    base = _ROOT / pathlib.Path(*root.split("."))
    if base.with_suffix(".py").exists():
        return {root}
    found = set()
    for path in base.rglob("*.py"):
        if path.name.endswith("_test.py") or path.name == "conftest.py":
            continue
        parts = path.relative_to(_ROOT).with_suffix("").parts
        found.add(".".join(parts[:-1] if parts[-1] == "__init__" else parts))
    return found


def _owner(module: str) -> str:
    """Returns the profile a module belongs to, the longest matching root winning."""
    matched = {
        root: profile
        for profile, roots in _ROOTS.items()
        for root in roots
        if module == root or module.startswith(f"{root}.")
    }
    longest = ""
    for root in matched:
        if len(root) > len(longest):
            longest = root
    return matched[longest] if longest else "base"


def _entry_points(profile: str) -> set[str]:
    """Returns the modules ``profile`` owns; its install has to import these."""
    return {module for module in _modules_under(_PACKAGE) if _owner(module) == profile}


def _imports(path: pathlib.Path) -> set[str]:
    """Returns the modules ``path`` requires, as written.

    An import inside a ``try`` is one the module states it can do without, so it is not
    something this module's install has to declare; the profile owning the guarded
    module is where that one is checked.
    """
    tree = ast.parse(path.read_text())
    guarded = {
        node
        for statement in ast.walk(tree)
        if isinstance(statement, ast.Try)
        for node in ast.walk(statement)
    }
    found = set()
    for node in ast.walk(tree):
        if node in guarded:
            continue
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            # ``from a.b import c`` reaches a.b, and c when c is itself a module.
            found.add(node.module)
            found.update(f"{node.module}.{alias.name}" for alias in node.names)
    return found


def _reachable(profile: str) -> set[str]:
    """Returns the third-party modules reachable from a profile's own modules.

    Args:
        profile: A key of ``_ROOTS``.

    Returns:
        The top-level name of every module imported from outside this package and
        outside the standard library, transitively.
    """
    seen: set[str] = set()
    pending = list(_entry_points(profile))
    third_party: set[str] = set()
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
            if root not in sys.stdlib_module_names and root != _PACKAGE:
                third_party.add(root)
            continue
        pending.extend(_imports(path))
    return third_party


@pytest.mark.parametrize("profile", _PROFILES)
def test_an_install_declares_what_its_modules_import(profile):
    extras = _EXTRAS[profile]
    install = f"ord-schema[{','.join(extras)}]" if extras else "ord-schema"
    declared = _declared(extras)
    missing = sorted(name for name in declared if not _is_installed(name))
    if missing:
        # Nothing here can say which distribution provides a module it cannot see, and
        # the heavy extras are deliberately absent from a development checkout.
        pytest.skip(f"{install} is not installed here, missing {sorted(missing)}")

    provided = importlib.metadata.packages_distributions()
    undeclared = {}
    for module in sorted(_reachable(profile)):
        distributions = {_canonical(name) for name in provided.get(module, [])}
        if not distributions:
            pytest.fail(
                f"{module} is imported but no installed distribution provides it"
            )
        if not distributions & declared:
            undeclared[module] = sorted(distributions)
    assert not undeclared, (
        f"imported by {install} but not declared by it: "
        + ", ".join(
            f"{module} ({'/'.join(names)})" for module, names in undeclared.items()
        )
    )


def test_the_walk_crosses_subpackages():
    # A walk that stopped at the subpackage boundary would have missed the import that
    # motivated this test, which the search path reaches through ord_schema.artifacts.
    assert "inflection" in _reachable("search")


@pytest.mark.parametrize("profile", _PROFILES)
def test_every_profile_owns_modules(profile):
    # A root that names nothing -- a typo, or a module that has moved -- would leave its
    # profile checking an empty set and passing for the wrong reason.
    assert _entry_points(profile)
