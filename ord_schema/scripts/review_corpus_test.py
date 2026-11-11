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

"""Tests for ord_schema.scripts.review_corpus."""

import csv
import pathlib
import shutil

import pytest

from ord_schema import datasets, message_helpers
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.scripts import review_corpus

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)


def _read(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _small_dataset(path: pathlib.Path) -> None:
    reaction = reaction_pb2.Reaction(reaction_id="r1")
    reaction.inputs["a"].components.append(
        message_helpers.build_compound(smiles="Nc1ccccc1", role="reactant")
    )
    product = reaction.outcomes.add().products.add()
    product.identifiers.add(type="SMILES", value="Brc1ccccc1N")
    dataset = dataset_pb2.Dataset(name="small", reactions=[reaction])
    datasets.save_dataset(dataset, path)


@pytest.fixture(name="corpus")
def corpus_fixture(tmp_path) -> pathlib.Path:
    folder = tmp_path / "corpus"
    folder.mkdir()
    shutil.copy(_NIELSEN, folder / "nielsen.pbtxt")
    _small_dataset(folder / "small.pbtxt")
    (folder / "broken.pbtxt").write_text("this is not a dataset")
    return folder


def test_each_dataset_is_tabulated_under_each_profile(corpus, tmp_path):
    output = tmp_path / "out"
    status = review_corpus.main(
        [f"--input_pattern={corpus}/*.pbtxt", f"--output_dir={output}", "--bundles=2"]
    )
    assert status == 0
    summary = _read(output / "summary.csv")
    assert sorted((r["dataset_id"], r["profile"]) for r in summary) == [
        ("nielsen", "uniform"),
        ("nielsen", "varied"),
        ("small", "uniform"),
        ("small", "varied"),
    ]
    small = next(r for r in summary if r["dataset_id"] == "small")
    assert small["reactions"] == "1"
    assert int(small["warnings"]) >= 1
    findings = _read(output / "findings.csv")
    assert ("small", "chemistry", "element-conservation") in {
        (r["dataset_id"], r["source"], r["check"]) for r in findings
    }
    items = _read(output / "completeness_items.csv")
    assert {r["dataset_id"] for r in items} == {"nielsen", "small"}
    (failure,) = _read(output / "failures.csv")
    assert failure["dataset_id"] == "broken"
    # The first two datasets in order get full bundles, under each profile; the
    # broken one has nothing to bundle.
    (bundled,) = (output / "bundles").iterdir()
    assert bundled.name == "nielsen"
    assert sorted(p.name for p in bundled.iterdir()) == ["uniform", "varied"]


def test_resume_skips_datasets_already_reviewed(corpus, tmp_path):
    output = tmp_path / "out"
    args = [f"--input_pattern={corpus}/*.pbtxt", f"--output_dir={output}"]
    review_corpus.main([*args, "--dataset_ids=small", "--profiles=uniform"])
    review_corpus.main([*args, "--resume", "--profiles=uniform", "--n_jobs=2"])
    summary = _read(output / "summary.csv")
    assert sorted(r["dataset_id"] for r in summary) == ["nielsen", "small"]
    assert [r["dataset_id"] for r in _read(output / "failures.csv")] == ["broken"]


def test_repeated_file_names_are_told_apart_by_folder(tmp_path):
    for folder in ("a", "b"):
        (tmp_path / folder).mkdir()
        shutil.copy(_NIELSEN, tmp_path / folder / "dataset.pbtxt")
    refs = review_corpus._local_refs(f"{tmp_path}/*/dataset.pbtxt")
    assert [ref.dataset_id for ref in refs] == ["a/dataset", "b/dataset"]


def test_the_sample_spreads_across_sizes():
    refs = [
        review_corpus.DatasetRef(f"d{i}", "", size) for i, size in enumerate(range(10))
    ]
    sample = review_corpus.spread_sample(refs, 4)
    assert [ref.size for ref in sample] == [0, 3, 6, 9]
    assert review_corpus.spread_sample(refs, 20) == refs


def test_an_unknown_profile_is_refused(tmp_path):
    with pytest.raises(SystemExit):
        review_corpus.main(
            ["--input_pattern=x", f"--output_dir={tmp_path}", "--profiles=strict"]
        )
