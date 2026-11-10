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
    "nl": ("ord_schema.search.nl", "ord_schema.search.nl_eval"),
    "orm": ("ord_schema.orm",),
    "huggingface": ("ord_schema.huggingface",),
    # Publishing checks a tree against the levels and paths the artifacts define.
    "upload": ("ord_schema.scripts.upload_artifacts",),
    "analysis": ("ord_schema.analysis",),
    # The review_dataset script is the review extra's command line.
    "review": ("ord_schema.review", "ord_schema.scripts.review_dataset"),
}
_EXTRAS: dict[str, tuple[str, ...]] = {
    "base": (),
    "search": ("search",),
    # The nl module reads through the search subpackage, so its install is both.
    "nl": ("search", "nl"),
    "orm": ("orm",),
    "huggingface": ("huggingface",),
    "upload": ("search", "huggingface"),
    "analysis": ("analysis",),
    # The review extra installs the analysis one, which its checks build on.
    "review": ("analysis", "review"),
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


def _all_imports(path: pathlib.Path, module: str) -> set[str]:
    """Returns every module ``path`` could import, relative and guarded ones included.

    Stricter than ``_imports``, which an install check can afford to be lenient in: a
    guarded import there is one the module can do without. A boundary between two
    subpackages is crossed by any import at all, so this resolves relative imports
    against ``module`` and keeps the ones inside a ``try``.

    Args:
        path: The file holding ``module``.
        module: The dotted name of the module in ``path``.

    Returns:
        Absolute dotted names. As in ``_imports``, ``from a import b`` yields both
        ``a`` and ``a.b``, since ``b`` may be a module.
    """
    package = module if path.name == "__init__.py" else module.rpartition(".")[0]
    found = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                parts = package.split(".")
                base = ".".join(parts[: len(parts) - (node.level - 1)])
                source = f"{base}.{node.module}" if node.module else base
            else:
                source = node.module or ""
            found.add(source)
            found.update(f"{source}.{alias.name}" for alias in node.names)
    return found


def _internal_reachable(root: str) -> set[str]:
    """Returns the modules of this package reachable from ``root``, transitively.

    Args:
        root: A dotted package name inside this package.

    Returns:
        Dotted names of the modules this package holds that ``root`` imports, directly
        or through another of its modules, by any of the routes ``_all_imports`` sees.
    """
    seen: set[str] = set()
    pending = list(_modules_under(root))
    while pending:
        module = pending.pop()
        if module in seen:
            continue
        path = _module_path(module)
        if path is None:
            continue
        seen.add(module)
        pending.extend(_all_imports(path, module))
    return seen


def _review_modules_reached_from_analysis() -> set[str]:
    """Returns the modules of ``ord_schema.review`` that ``analysis`` reaches."""
    return {
        module
        for module in _internal_reachable("ord_schema.analysis")
        if module == "ord_schema.review" or module.startswith("ord_schema.review.")
    }


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


def test_analysis_never_imports_review():
    # The review checks build on the analysis summaries. The reverse would make the
    # analysis extra need whatever review installs, and the two layers would no longer
    # be separable.
    assert not _review_modules_reached_from_analysis()


def test_the_internal_walk_follows_imports():
    # The check above passes on an empty walk too; this one shows the walk reaches
    # another subpackage when a module does import it.
    reached = _internal_reachable("ord_schema.search")
    assert any(module.startswith("ord_schema.artifacts") for module in reached)


def _write_package(root: pathlib.Path, analysis: dict[str, str]) -> None:
    """Writes a minimal ord_schema with a review subpackage and the given analysis.

    Args:
        root: The directory to write ``ord_schema/`` into.
        analysis: Source for each module of ``ord_schema.analysis``, by file name.
    """
    package = root / "ord_schema"
    (package / "review").mkdir(parents=True)
    (package / "analysis").mkdir()
    (package / "__init__.py").write_text("")
    (package / "review" / "__init__.py").write_text("")
    (package / "review" / "report.py").write_text("")
    files = {"__init__.py": ""} | analysis
    for name, source in files.items():
        (package / "analysis" / name).write_text(source)


_GUARDED = "try:\n    from ord_schema import review\nexcept ImportError:\n    pass\n"


@pytest.mark.parametrize(
    "analysis",
    [
        pytest.param({"__init__.py": "import ord_schema.review\n"}, id="absolute"),
        pytest.param({"summary.py": "from ord_schema import review\n"}, id="from"),
        pytest.param({"summary.py": "from .. import review\n"}, id="relative"),
        pytest.param(
            {"summary.py": "from ..review import report\n"}, id="relative-module"
        ),
        pytest.param({"__init__.py": "from .. import review\n"}, id="relative-init"),
        pytest.param({"summary.py": _GUARDED}, id="guarded"),
        pytest.param(
            {
                "summary.py": "from . import helpers\n",
                "helpers.py": "from ..review import report\n",
            },
            id="through-another-module",
        ),
    ],
)
def test_the_boundary_check_sees_every_route_into_review(
    tmp_path, monkeypatch, analysis
):
    # A route the check missed would let the rule be broken with the test still green.
    _write_package(tmp_path, analysis)
    monkeypatch.setattr(sys.modules[__name__], "_ROOT", tmp_path)
    assert _review_modules_reached_from_analysis()


def test_the_boundary_check_passes_a_clean_analysis(tmp_path, monkeypatch):
    # The routes above have to be caught without flagging an analysis that only reaches
    # its own modules and the rest of the package.
    _write_package(
        tmp_path,
        {
            "summary.py": "from . import helpers\nfrom ord_schema import validations\n",
            "helpers.py": "import collections\n",
        },
    )
    monkeypatch.setattr(sys.modules[__name__], "_ROOT", tmp_path)
    assert not _review_modules_reached_from_analysis()
