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

"""Rewrites notebooks into the one JSON form nbformat serializes, less cell timings.

A notebook is JSON, and the tools that write one do not agree on how: whether a
multi-line output is one string or a list of lines, how keys are ordered, how much a
line is indented. Editing a notebook through any library therefore reformats whatever
it touches, and the reformatting arrives mixed into the diff of the actual change.
Holding every notebook in nbformat's serialization is what makes that stop -- the
reformatting happens once, here, rather than in whichever pull request next opens one.

Executing a notebook additionally stamps every cell with ``metadata.execution``: four
wall-clock instants from the kernel, which change on every run whether or not the
output did. nbclient offers no way to decline them. Nothing reads them either, since
nbval compares outputs and ignores cell metadata, so they come out.
"""

import argparse
import io
import sys
from pathlib import Path

import nbformat


def format_notebook(path: Path) -> bool:
    """Rewrites one notebook if it is not already formatted.

    Args:
        path: The notebook to format.

    Returns:
        Whether the file was rewritten.
    """
    original = path.read_text()
    notebook = nbformat.read(io.StringIO(original), as_version=nbformat.NO_CONVERT)
    for cell in notebook.cells:
        cell.get("metadata", {}).pop("execution", None)
    formatted = nbformat.writes(notebook) + "\n"
    if formatted == original:
        return False
    path.write_text(formatted)
    return True


def main(argv: list[str] | None = None) -> int:
    """Formats every notebook named on the command line.

    Args:
        argv: Command-line arguments; ``sys.argv[1:]`` by default.

    Returns:
        1 if any notebook was rewritten, so a pre-commit run fails and the author
        re-stages the formatted file; 0 if they were all formatted already.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebooks", nargs="+", type=Path, help="Notebooks to format")
    args = parser.parse_args(argv)
    rewritten = [path for path in args.notebooks if format_notebook(path)]
    for path in rewritten:
        print(f"formatted {path}")
    return 1 if rewritten else 0


if __name__ == "__main__":
    sys.exit(main())
