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

"""Removes per-cell execution timings from notebooks.

Executing a notebook stamps every cell with ``metadata.execution`` -- four wall-clock
instants recording when the kernel started and finished it. nbclient offers no way to
turn that off, so the timings arrive with any regenerated output and change on every
run, whether or not the output did.

Nothing reads them. nbval compares outputs and ignores cell metadata, and a reader
learns nothing from the microsecond a cell ran on someone's laptop. Left in, they are
the bulk of a regeneration's diff and hide the output change that motivated it.
"""

import argparse
import sys
from pathlib import Path

import nbformat


def strip(path: Path) -> bool:
    """Removes the execution timings from one notebook, rewriting it if any were found.

    Args:
        path: The notebook to strip.

    Returns:
        Whether the file was rewritten.
    """
    notebook = nbformat.read(path, as_version=nbformat.NO_CONVERT)
    stripped = False
    for cell in notebook.cells:
        if cell.get("metadata", {}).pop("execution", None) is not None:
            stripped = True
    if stripped:
        # nbformat round-trips byte for byte, so the diff is the timings and nothing
        # else -- and a notebook that has none is never rewritten at all.
        nbformat.write(notebook, path)
    return stripped


def main(argv: list[str] | None = None) -> int:
    """Strips every notebook named on the command line.

    Args:
        argv: Command-line arguments; ``sys.argv[1:]`` by default.

    Returns:
        1 if any notebook was rewritten, so a pre-commit run fails and the author
        re-stages the stripped file; 0 if they were all clean.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebooks", nargs="+", type=Path, help="Notebooks to strip")
    args = parser.parse_args(argv)
    rewritten = [path for path in args.notebooks if strip(path)]
    for path in rewritten:
        print(f"stripped execution timings from {path}")
    return 1 if rewritten else 0


if __name__ == "__main__":
    sys.exit(main())
