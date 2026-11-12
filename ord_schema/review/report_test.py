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

"""Tests for ord_schema.review.report and the review_dataset script."""

import csv
import pathlib
import re

import pytest
from rdkit import Chem

from ord_schema import datasets, message_helpers, parquet
from ord_schema.analysis.compounds import CompoundTable, StructureDrawings
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.review import report, structure
from ord_schema.scripts import review_dataset
from ord_schema.validations import Severity

_NIELSEN = (
    pathlib.Path(__file__).resolve().parents[2]
    / "examples/submissions/2_Nielsen_Deoxyfluorination_Screen"
    / "nielsen_fig1_dataset.pbtxt"
)


def _reaction(reaction_id: str, *, base: str, yield_: float) -> reaction_pb2.Reaction:
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.inputs["solvent"].components.append(
        message_helpers.build_compound(smiles="CS(C)=O", amount="1 mL", role="solvent")
    )
    reaction.inputs["base"].components.append(
        message_helpers.build_compound(smiles=base, amount="1 mmol", role="reagent")
    )
    outcome = reaction.outcomes.add()
    product = outcome.products.add()
    product.identifiers.add(type="SMILES", value="CCO")
    product.measurements.add(type="YIELD").percentage.value = yield_
    reaction.provenance.record_created.time.value = "2026-01-05"
    reaction.provenance.record_created.person.name = "A. Chemist"
    return reaction


def _dataset() -> dataset_pb2.Dataset:
    """A small screen with something for each part of the report to find."""
    reactions = [
        _reaction(f"r{i}", base=base, yield_=float(10 + i))
        for i, base in enumerate(["CCN(CC)CC", "CN(C)C", "C1CCNCC1"] * 4)
    ]
    # Yields as fractions: one validation warning, raised with different numbers.
    reactions[0].outcomes[0].products[0].measurements[0].percentage.value = 0.25
    reactions[1].outcomes[0].products[0].measurements[0].percentage.value = 0.75
    reactions[2].provenance.doi = "10.1000/example"  # Mixed sourcing.
    return dataset_pb2.Dataset(
        name="A <test> screen",
        description="Bases & yields",
        dataset_id="ord_dataset-0123456789abcdef0123456789abcdef",
        reactions=reactions,
    )


def test_messages_that_differ_only_in_numbers_are_one_group():
    review = report.build_report(_dataset())
    fractions = [g for g in review.validation if "not fractions" in g.message]
    assert len(fractions) == 1
    (group,) = fractions
    assert group.severity is Severity.WARNING
    assert group.reactions == 2
    assert group.examples == ["r0", "r1"]


def test_the_report_carries_the_structural_findings_and_profile():
    review = report.build_report(_dataset())
    assert review.profile is structure.ValidationProfile.UNIFORM
    assert review.structure.by_check()["mixed-sourcing"] == 1
    assert review.reaction_count == 12
    counts = review.severity_counts()
    assert counts[Severity.WARNING] >= 2


def test_reviewing_does_not_change_the_dataset():
    dataset = _dataset()
    before = dataset.SerializeToString(deterministic=True)
    report.build_report(dataset)
    assert dataset.SerializeToString(deterministic=True) == before


def test_the_bundle_holds_the_page_and_every_table(tmp_path):
    page = report.write_bundle(report.build_report(_dataset()), tmp_path)
    assert page == tmp_path / "ord_dataset-0123456789abcdef0123456789abcdef_review.html"
    tables = {path.name for path in (tmp_path / "tables").iterdir()}
    assert tables == {
        "findings.csv",
        "compounds.csv",
        "flat.csv",
        "constants.csv",
        "column_stats.csv",
        "exemplars.txtpb",
    }
    with (tmp_path / "tables" / "findings.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert {row["source"] for row in rows} == {"structure", "validation"}
    assert any(row["check"] == "mixed-sourcing" for row in rows)


def test_the_page_is_self_contained_and_escaped():
    page = report.render_html(report.build_report(_dataset()))
    assert page.startswith("<!DOCTYPE html>")
    assert "Review: A &lt;test&gt; screen" in page
    assert "<test>" not in page
    assert "Bases &amp; yields" in page
    assert "src=" not in page
    assert "<svg" in page
    for heading in (
        "Dataset-level findings",
        "Reaction validation",
        "How reactions are recorded",
        "Exemplar reactions",
        "What the dataset varies",
        "Compounds",
    ):
        assert heading in page


def test_a_varied_dataset_says_which_checks_it_skipped_and_leaves_out_the_charts():
    review = report.build_report(_dataset(), structure.ValidationProfile.VARIED)
    page = report.render_html(review)
    assert "Not applied under the <b>varied</b> profile" in page
    assert "tables/flat.csv" in page
    assert "<figure>" not in page


def test_long_lists_are_capped_on_the_page_but_not_in_the_tables(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "MAX_COMPOUNDS_SHOWN", 1)
    review = report.build_report(_dataset())
    page = report.write_bundle(review, tmp_path).read_text()
    assert "more REAGENT compounds, the least used, are in compounds.csv" in page
    with (tmp_path / "tables" / "compounds.csv").open(newline="") as f:
        assert len(list(csv.DictReader(f))) == len(review.evidence.summary.compounds)


def test_a_parquet_dataset_is_streamed(tmp_path):
    path = tmp_path / "dataset.parquet"
    parquet.save_dataset(_dataset(), path)
    view = datasets.load_dataset(path)
    assert isinstance(view, parquet.DatasetView)
    review = report.build_report(view)
    assert review.reaction_count == 12
    assert review.name == "A <test> screen"


def test_a_consistent_real_screen():
    review = report.build_report(datasets.load_dataset(_NIELSEN))
    assert review.structure.findings == []
    assert [g.reactions for g in review.validation] == [80]


def test_the_script_writes_the_bundle(tmp_path):
    status = review_dataset.main(
        [f"--input={_NIELSEN}", f"--output_dir={tmp_path}", "--profile=uniform"]
    )
    assert status == 0
    assert list(tmp_path.glob("*_review.html"))
    assert (tmp_path / "tables" / "findings.csv").exists()


def test_the_script_explains_a_git_lfs_pointer(tmp_path, caplog):
    pointer = tmp_path / "dataset.pb.gz"
    pointer.write_bytes(b"version https://git-lfs.github.com/spec/v1\noid sha256:0\n")
    status = review_dataset.main([f"--input={pointer}", f"--output_dir={tmp_path}"])
    assert status == 1
    assert "git lfs pull" in caplog.text


def test_the_script_refuses_an_unknown_profile(tmp_path):
    with pytest.raises(SystemExit):
        review_dataset.main(
            [f"--input={_NIELSEN}", f"--output_dir={tmp_path}", "--profile=bulk"]
        )


def test_a_compounds_structure_is_drawn_once_and_shown_on_hover():
    reactions = [
        _reaction(f"r{i}", base=base, yield_=float(10 + i))
        for i, base in enumerate(["CCN(CC)CC", "CN(C)C", "C1CCNCC1"] * 4)
    ]
    # Triethylamine as a solvent too, so it is in the cross-role index.
    reactions[0].inputs["solvent"].components[0].identifiers[0].value = "CCN(CC)CC"
    page = report.render_html(
        report.build_report(dataset_pb2.Dataset(name="x", reactions=reactions))
    )
    # One symbol per structure, defined before the content that uses it.
    symbols = re.findall(r'<symbol id="(mol-\d+)"', page)
    assert len(symbols) == len(set(symbols)) >= 3
    assert page.index("<symbol") < page.index("<main>")
    for symbol in symbols:
        assert f'href="#{symbol}"' in page
    # The base chart's bars, and the cross-role index.
    assert '<div class="chart"><svg' in page
    assert re.search(r'<g class="b\d+" tabindex="0">', page)
    cross = page[page.index("Compounds used in more than one role") :]
    assert '<span class="mol" tabindex="0"><code>CCN(CC)CC</code>' in cross
    assert ".moltip" in page


def test_a_name_that_parses_as_smiles_gets_no_drawing():
    drawings = StructureDrawings(CompoundTable())
    assert drawings.tip("CO") is None
    assert drawings.code("CO") == "<code>CO</code>"
    assert drawings.defs() == ""


def test_smiles_in_finding_examples_show_their_structure_on_hover():
    reactions = [
        _reaction(f"r{i}", base=base, yield_=float(10 + i))
        for i, base in enumerate(["CCN(CC)CC", "CN(C)C"] * 3)
    ]
    # One name for two structures: the finding's examples are their SMILES.
    for reaction in reactions:
        reaction.inputs["base"].components[0].identifiers.add(type="NAME", value="base")
    page = report.render_html(
        report.build_report(dataset_pb2.Dataset(name="x", reactions=reactions))
    )
    findings = page[
        page.index("Dataset-level findings") : page.index("Reaction validation")
    ]
    assert "name-structure-conflict" in findings
    assert '<span class="mol" tabindex="0"><code>CN(C)C</code>' in findings


def test_a_molblock_or_inchi_value_shows_its_structure_too():
    reactions = []
    for i, smiles in enumerate(["c1ccc2ccccc2c1", "N#Cc1ccc(C#N)cc1"] * 3):
        reaction = _reaction(f"r{i}", base="CCN(CC)CC", yield_=float(10 + i))
        catalyst = reaction.inputs["catalyst"].components.add()
        catalyst.identifiers.add(
            type="MOLBLOCK", value=Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))
        )
        catalyst.reaction_role = reaction_pb2.ReactionRole.CATALYST
        reactions.append(reaction)
    review = report.build_report(dataset_pb2.Dataset(name="x", reactions=reactions))
    drawings = StructureDrawings(review.evidence.summary.compounds)
    molblock = reactions[0].inputs["catalyst"].components[0].identifiers[0].value
    assert drawings.tip(molblock) is not None
    assert drawings.tip(Chem.MolToInchi(Chem.MolFromSmiles("c1ccc2ccccc2c1")))
    # A structure the dataset doesn't record gets none.
    assert drawings.tip(Chem.MolToMolBlock(Chem.MolFromSmiles("CCCC"))) is None
    page = report.render_html(review)
    assert re.search(r'<g class="b\d+" tabindex="0">', page)
