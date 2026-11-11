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

"""Runs the dataset review over many datasets and tabulates what it finds.

Usage:
    python -m ord_schema.scripts.review_corpus --output_dir=corpus/
        (--input_pattern='data/*/*.parquet' | --from_huggingface)
        [--dataset_ids=ID,...] [--sample=N] [--profiles=uniform,varied]
        [--n_jobs=4] [--bundles=N] [--resume]

Each dataset is reviewed under each profile, and the results are written, as each
dataset finishes, to long-format tables in the output directory:

- ``summary.csv``: one row per dataset and profile: reactions, time taken, counts
  at each severity, completeness by level, shapes, compounds, varying columns, and
  how many distinct sources (DOIs and patents) the reactions cite.
- ``findings.csv``: every finding, structural, chemistry, validation and
  completeness, with its dataset and profile.
- ``completeness_items.csv``: for each rubric item, how many reactions it applied
  to and how many lacked it, per dataset and profile.
- ``failures.csv``: datasets that could not be read or reviewed, and why.

``--sample=N`` picks N datasets spread evenly across file sizes, so a sample spans
small screens and large bulk datasets alike. ``--bundles=N`` also writes the full
review bundle for the first N datasets reviewed. With ``--from_huggingface``,
datasets are fetched from the ord-data mirror one at a time and deleted once
reviewed, unless ``--keep_downloads``. ``--resume`` skips datasets already in
``summary.csv``. Needs the ``review`` and ``huggingface`` extras.
"""

import argparse
import csv
import dataclasses
import glob
import pathlib
import shutil
import sys
import tempfile
import time
import traceback
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

from huggingface_hub import HfApi

from ord_schema import datasets, huggingface
from ord_schema.logging import get_logger, silence_rdkit_logs
from ord_schema.review import report, structure
from ord_schema.validations import Severity

logger = get_logger(__name__)

_SUMMARY_FIELDS = [
    "dataset_id",
    "name",
    "source",
    "profile",
    "reactions",
    "seconds",
    "errors",
    "warnings",
    "notes",
    "structural_findings",
    "chemistry_findings",
    "validation_messages",
    "completeness_overall",
    "completeness_inputs",
    "completeness_conditions",
    "completeness_workups",
    "completeness_outcomes",
    "completeness_provenance",
    "shapes",
    "compounds",
    "varying_columns",
    "column_groups",
    "constant_columns",
    "distinct_sources",
]
_FINDING_FIELDS = [
    "dataset_id",
    "profile",
    "source",
    "check",
    "severity",
    "reactions",
    "message",
    "examples",
]
_ITEM_FIELDS = [
    "dataset_id",
    "profile",
    "item",
    "level",
    "weight",
    "applicable",
    "lacking",
]
_FAILURE_FIELDS = ["dataset_id", "source", "error"]


@dataclasses.dataclass(frozen=True)
class DatasetRef:
    """A dataset to review: a local file, or an ID to fetch from the mirror.

    Attributes:
        dataset_id: Its ID, as the file is named.
        path: A local file, or empty to fetch from the mirror.
        size: The file's size in bytes, for sampling.
    """

    dataset_id: str
    path: str
    size: int


def _local_refs(pattern: str) -> list[DatasetRef]:
    """Local files, named by their stems, or by folder and stem where stems repeat."""
    paths = [pathlib.Path(path) for path in sorted(glob.glob(pattern))]
    stems = [path.name.split(".")[0] for path in paths]
    refs = []
    for path, stem in zip(paths, stems, strict=True):
        dataset_id = stem if stems.count(stem) == 1 else f"{path.parent.name}/{stem}"
        refs.append(DatasetRef(dataset_id, str(path), path.stat().st_size))
    return refs


def _mirror_refs(revision: str) -> list[DatasetRef]:
    info = HfApi().dataset_info(
        huggingface.ORD_DATA_HF_REPO, revision=revision, files_metadata=True
    )
    refs = []
    for sibling in info.siblings or []:
        name = pathlib.Path(sibling.rfilename).name
        if sibling.rfilename.startswith("data/") and name.endswith(
            (".parquet", ".pb.gz")
        ):
            refs.append(DatasetRef(name.split(".")[0], "", sibling.size or 0))
    return sorted(refs, key=lambda ref: ref.dataset_id)


def spread_sample(refs: list[DatasetRef], n: int) -> list[DatasetRef]:
    """N datasets spread evenly across file sizes, smallest and largest included.

    Args:
        refs: The datasets to choose from.
        n: How many to choose.

    Returns:
        The chosen datasets, smallest first.
    """
    ordered = sorted(refs, key=lambda ref: (ref.size, ref.dataset_id))
    if n >= len(ordered):
        return ordered
    if n <= 1:
        return ordered[:n]
    step = (len(ordered) - 1) / (n - 1)
    return [ordered[round(i * step)] for i in range(n)]


def _distinct_sources(review: report.ReviewReport) -> str:
    """How many distinct DOIs and patents the reactions cite; "N+" past the cap."""
    values = review.evidence.summary.values
    untracked = review.evidence.summary.untracked
    count = 0
    capped = False
    for path in ("provenance.doi", "provenance.patent"):
        count += len(values.get(path, {}))
        capped = capped or bool(untracked.get(path))
    return f"{count}+" if capped else str(count)


def _round(value: float | None) -> str:
    return "" if value is None else f"{value:.3f}"


def _rows(
    review: report.ReviewReport, dataset_id: str, source: str, seconds: float
) -> dict[str, list[dict[str, Any]]]:
    """The table rows for one dataset reviewed under one profile."""
    profile = review.profile.value
    counts = review.severity_counts()
    levels = review.completeness.mean_levels()
    table = review.evidence.table
    summary = {
        "dataset_id": dataset_id,
        "name": review.name,
        "source": source,
        "profile": profile,
        "reactions": review.reaction_count,
        "seconds": f"{seconds:.1f}",
        "errors": counts[Severity.ERROR],
        "warnings": counts[Severity.WARNING],
        "notes": counts[Severity.INFO],
        "structural_findings": len(review.structure.findings),
        "chemistry_findings": len(review.chemistry),
        "validation_messages": len(review.validation),
        "completeness_overall": _round(review.completeness.mean_overall()),
        **{f"completeness_{level}": _round(value) for level, value in levels.items()},
        "shapes": len(review.evidence.summary.shapes),
        "compounds": len(review.evidence.summary.compounds),
        "varying_columns": len(table.varying.columns),
        "column_groups": len(review.column_stats),
        "constant_columns": len(table.constants),
        "distinct_sources": _distinct_sources(review),
    }
    findings = [
        {
            "dataset_id": dataset_id,
            "profile": profile,
            "source": kind,
            "check": f.check,
            "severity": f.severity.name,
            "reactions": "" if f.reactions is None else f.reactions,
            "message": f.message,
            "examples": " | ".join(f.examples),
        }
        for kind, group in (
            ("structure", review.structure.findings),
            ("chemistry", review.chemistry),
        )
        for f in group
    ]
    findings += [
        {
            "dataset_id": dataset_id,
            "profile": profile,
            "source": "validation",
            "check": "",
            "severity": g.severity.name,
            "reactions": g.reactions,
            "message": g.message,
            "examples": " | ".join(g.examples),
        }
        for g in review.validation
    ]
    items = [
        {
            "dataset_id": dataset_id,
            "profile": profile,
            "item": r.item.id,
            "level": r.item.level,
            "weight": r.weight,
            "applicable": r.applicable,
            "lacking": r.lacking,
        }
        for r in review.completeness.items.values()
    ]
    return {"summary": [summary], "findings": findings, "items": items}


def review_one(
    ref: DatasetRef,
    profiles: list[str],
    *,
    revision: str,
    bundle_dir: str | None,
    keep_downloads: bool,
    packet_budget_tokens: int,
) -> dict[str, list[dict[str, Any]]]:
    """Reviews one dataset under each profile; never raises.

    Returns:
        Rows for each table; a failure becomes a ``failures`` row.
    """
    silence_rdkit_logs()
    download_dir = None
    source = ref.path or f"{huggingface.ORD_DATA_HF_REPO}@{revision}"
    rows: dict[str, list[dict[str, Any]]] = {
        "summary": [],
        "findings": [],
        "items": [],
        "failures": [],
    }
    try:
        path = ref.path
        if not path:
            download_dir = tempfile.mkdtemp(prefix="review_corpus_")
            path = huggingface.fetch_dataset(
                ref.dataset_id, revision=revision, local_dir=download_dir
            )
        for profile in profiles:
            start = time.monotonic()
            dataset = datasets.load_dataset(path)
            review = report.build_report(dataset, structure.ValidationProfile(profile))
            for table, new in _rows(
                review, ref.dataset_id, source, time.monotonic() - start
            ).items():
                rows[table].extend(new)
            if bundle_dir:
                report.write_bundle(
                    review,
                    pathlib.Path(bundle_dir) / ref.dataset_id / profile,
                    packet_budget_tokens=packet_budget_tokens,
                )
    except Exception as error:  # noqa: BLE001 -- one bad dataset must not stop a corpus run
        rows["failures"].append(
            {
                "dataset_id": ref.dataset_id,
                "source": source,
                "error": f"{type(error).__name__}: {error}",
            }
        )
        logger.debug(traceback.format_exc())
    finally:
        if download_dir and not keep_downloads:
            shutil.rmtree(download_dir, ignore_errors=True)
    return rows


class _Tables:
    """The output tables, appended to as each dataset finishes."""

    def __init__(self, output_dir: pathlib.Path, *, resume: bool) -> None:
        self.paths = {
            "summary": (output_dir / "summary.csv", _SUMMARY_FIELDS),
            "findings": (output_dir / "findings.csv", _FINDING_FIELDS),
            "items": (output_dir / "completeness_items.csv", _ITEM_FIELDS),
            "failures": (output_dir / "failures.csv", _FAILURE_FIELDS),
        }
        for path, fields in self.paths.values():
            if not (resume and path.exists()):
                with path.open("w", newline="", encoding="utf-8") as f:
                    csv.DictWriter(f, fieldnames=fields).writeheader()

    def done(self) -> set[str]:
        """The datasets ``summary.csv`` already holds."""
        path, _ = self.paths["summary"]
        with path.open(newline="", encoding="utf-8") as f:
            return {row["dataset_id"] for row in csv.DictReader(f)}

    def append(self, rows: dict[str, list[dict[str, Any]]]) -> None:
        for table, new in rows.items():
            path, fields = self.paths[table]
            with path.open("a", newline="", encoding="utf-8") as f:
                csv.DictWriter(f, fieldnames=fields).writerows(new)


def _parse(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Runs the dataset review over many datasets."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input_pattern", help="Glob of local dataset files.")
    source.add_argument(
        "--from_huggingface",
        action="store_true",
        help="Fetch datasets from the ord-data Hugging Face mirror.",
    )
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--dataset_ids", help="Comma-separated IDs to review.")
    parser.add_argument(
        "--sample", type=int, help="Review N datasets spread across file sizes."
    )
    parser.add_argument(
        "--profiles",
        default="uniform,varied",
        help="Comma-separated profiles to review each dataset under.",
    )
    parser.add_argument("--n_jobs", type=int, default=1)
    parser.add_argument(
        "--bundles", type=int, default=0, help="Write full bundles for N datasets."
    )
    parser.add_argument("--revision", default="main")
    parser.add_argument("--keep_downloads", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--packet_budget_tokens", type=int, default=100_000, help="For the bundles."
    )
    args = parser.parse_args(argv)
    known = {p.value for p in structure.ValidationProfile}
    args.profile_list = [p.strip() for p in args.profiles.split(",") if p.strip()]
    unknown = sorted(set(args.profile_list) - known)
    if unknown:
        parser.error(f"unknown profiles {unknown}; choose from {sorted(known)}")
    return args


def _select(args: argparse.Namespace) -> list[DatasetRef]:
    refs = (
        _local_refs(args.input_pattern)
        if args.input_pattern
        else _mirror_refs(args.revision)
    )
    if args.dataset_ids:
        wanted = {i.strip() for i in args.dataset_ids.split(",") if i.strip()}
        refs = [ref for ref in refs if ref.dataset_id in wanted]
    if args.sample:
        refs = spread_sample(refs, args.sample)
    return refs


def _results(
    refs: list[DatasetRef], args: argparse.Namespace, output_dir: pathlib.Path
) -> Iterator[dict[str, list[dict[str, Any]]]]:
    def options(index: int) -> dict[str, Any]:
        return {
            "revision": args.revision,
            "bundle_dir": str(output_dir / "bundles") if index < args.bundles else None,
            "keep_downloads": args.keep_downloads,
            "packet_budget_tokens": args.packet_budget_tokens,
        }

    if args.n_jobs <= 1:
        for index, ref in enumerate(refs):
            yield review_one(ref, args.profile_list, **options(index))
        return
    with ProcessPoolExecutor(max_workers=args.n_jobs) as pool:
        futures = [
            pool.submit(review_one, ref, args.profile_list, **options(index))
            for index, ref in enumerate(refs)
        ]
        for future in as_completed(futures):
            yield future.result()


def main(argv: list[str] | None = None) -> int:
    """Reviews the chosen datasets; returns the process exit status."""
    args = _parse(argv)
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = _Tables(output_dir, resume=args.resume)
    refs = _select(args)
    if args.resume:
        done = tables.done()
        refs = [ref for ref in refs if ref.dataset_id not in done]
    logger.info(f"reviewing {len(refs)} datasets under {args.profile_list}")
    failed = 0
    for count, rows in enumerate(_results(refs, args, output_dir), start=1):
        tables.append(rows)
        failed += bool(rows["failures"])
        names = {row["dataset_id"] for row in rows["summary"] + rows["failures"]}
        logger.info(f"[{count}/{len(refs)}] {', '.join(sorted(names))}")
    logger.info(
        f"done: {len(refs) - failed} reviewed, {failed} failed; see {output_dir}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
