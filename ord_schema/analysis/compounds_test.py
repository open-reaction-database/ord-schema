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

"""Tests for ord_schema.analysis.compounds."""

import csv

from rdkit import Chem

from ord_schema import message_helpers
from ord_schema.analysis import compounds
from ord_schema.proto import reaction_pb2

_CAS = reaction_pb2.CompoundIdentifier.CAS_NUMBER


def _reaction(
    reaction_id: str, *inputs: reaction_pb2.Compound
) -> reaction_pb2.Reaction:
    """Returns a reaction with each of ``inputs`` as its own input."""
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    for index, compound in enumerate(inputs):
        reaction.inputs[f"input {index}"].components.append(compound)
    return reaction


def _with_cas(compound: reaction_pb2.Compound, cas: str) -> reaction_pb2.Compound:
    compound.identifiers.add(type=_CAS, value=cas)
    return compound


def _table(*reactions: reaction_pb2.Reaction) -> compounds.CompoundTable:
    table = compounds.CompoundTable()
    table.add_reactions(reactions)
    return table


def _record(table: compounds.CompoundTable, key: str) -> compounds.CompoundRecord:
    (record,) = [r for r in table.records() if r.key == key]
    return record


def test_compounds_are_identified_by_canonical_smiles():
    # The same structure written two ways is one compound, and the table says so.
    table = _table(
        _reaction("r1", message_helpers.build_compound(smiles="OCC", role="solvent")),
        _reaction("r2", message_helpers.build_compound(smiles="CCO", role="solvent")),
    )
    assert len(table) == 1
    record = _record(table, "smiles:CCO")
    assert record.reaction_count == 2
    assert record.smiles_variants == {"OCC", "CCO"}
    assert table.summary()["written_more_than_one_way"] == 1


def test_a_compound_without_a_structure_is_identified_by_name():
    table = _table(
        _reaction("r1", message_helpers.build_compound(name="Pd(PPh3)4")),
        _reaction("r2", message_helpers.build_compound(name="  pd(pph3)4 ")),
    )
    record = _record(table, "name:pd(pph3)4")
    assert record.reaction_count == 2
    assert record.smiles == ""


def test_a_compound_is_counted_once_per_reaction():
    solvent = message_helpers.build_compound(smiles="CCO", role="solvent")
    table = _table(_reaction("r1", solvent, solvent))
    record = _record(table, "smiles:CCO")
    assert record.reaction_count == 1
    assert record.roles["SOLVENT"] == 2
    assert record.example_reaction_ids == ["r1"]


def test_roles_come_from_the_proto_descriptor_and_group_in_reading_order():
    table = _table(
        _reaction(
            "r1",
            message_helpers.build_compound(smiles="CCO", role="solvent"),
            message_helpers.build_compound(smiles="c1ccccc1Br", role="reactant"),
            message_helpers.build_compound(name="catalyst X", role="catalyst"),
        )
    )
    assert list(table.by_role()) == ["REACTANT", "CATALYST", "SOLVENT"]


def test_a_compound_in_several_roles_is_indexed_across_roles():
    table = _table(
        _reaction("r1", message_helpers.build_compound(smiles="CCO", role="reactant")),
        _reaction("r2", message_helpers.build_compound(smiles="CCO", role="solvent")),
    )
    (record,) = table.cross_role()
    assert set(record.roles) == {"REACTANT", "SOLVENT"}
    groups = table.by_role()
    assert record in groups["REACTANT"]
    assert record in groups["SOLVENT"]


def test_workup_inputs_crude_components_and_products_are_included():
    reaction = _reaction("r1", message_helpers.build_compound(smiles="CCO"))
    workup = reaction.workups.add(type=reaction_pb2.ReactionWorkup.ADDITION)
    workup.input.components.append(
        message_helpers.build_compound(name="brine", role="workup")
    )
    reaction.inputs["crude"].crude_components.add(reaction_id="r0")
    product = reaction.outcomes.add().products.add()
    product.identifiers.add(type=reaction_pb2.CompoundIdentifier.SMILES, value="CC=O")
    table = _table(reaction)

    assert _record(table, "name:brine").locations == {compounds.WORKUP: 1}
    crude = _record(table, "crude:r0")
    assert set(crude.roles) == {compounds.CRUDE}
    # Labelled by the upstream reaction that made it; r1 is only where it was used.
    assert crude.label == "crude product of r0"
    assert crude.example_reaction_ids == ["r1"]
    # A product recorded with no role is still listed as a product.
    assert set(_record(table, "smiles:CC=O").roles) == {"PRODUCT"}
    assert _record(table, "smiles:CC=O").locations == {compounds.OUTCOME: 1}


def test_two_cas_numbers_for_one_structure_are_listed_not_resolved():
    table = _table(
        _reaction(
            "r1", _with_cas(message_helpers.build_compound(smiles="CCO"), "64-17-5")
        ),
        _reaction(
            "r2", _with_cas(message_helpers.build_compound(smiles="CCO"), "64-17-6")
        ),
    )
    record = _record(table, "smiles:CCO")
    assert table.conflicts(record) == ["more than one CAS number: 64-17-5, 64-17-6"]
    assert set(record.identifiers["CAS_NUMBER"]) == {"64-17-5", "64-17-6"}


def test_one_name_for_two_structures_is_a_conflict_on_both():
    table = _table(
        _reaction("r1", message_helpers.build_compound(smiles="CCO", name="spirit")),
        _reaction("r2", message_helpers.build_compound(smiles="CO", name="Spirit")),
    )
    assert table.name_conflicts() == {"spirit": ["CCO", "CO"]}
    assert table.conflicts(_record(table, "smiles:CCO")) == [
        "name 'spirit' is also used for: CO"
    ]
    assert table.conflicts(_record(table, "smiles:CO")) == [
        "name 'Spirit' is also used for: CCO"
    ]


def test_identifiers_that_disagree_within_one_compound_are_flagged():
    compound = message_helpers.build_compound(smiles="CCO")
    compound.identifiers.add(
        type=reaction_pb2.CompoundIdentifier.INCHI, value="InChI=1S/CH4O/c1-2/h2H,1H3"
    )
    table = _table(_reaction("r1", compound))
    (record,) = table.records()
    assert table.conflicts(record) == [
        "identifiers in one compound describe different structures: CCO, CO"
    ]


def test_an_unreadable_structure_is_flagged_and_falls_back_to_the_name():
    compound = message_helpers.build_compound(name="mystery")
    compound.identifiers.add(type=reaction_pb2.CompoundIdentifier.SMILES, value="C1CC")
    table = _table(_reaction("r1", compound))
    record = _record(table, "name:mystery")
    assert table.conflicts(record) == ["unreadable SMILES: C1CC"]


def test_tables_built_over_parts_of_a_dataset_merge():
    whole = _table(
        _reaction("r1", message_helpers.build_compound(smiles="CCO", role="solvent")),
        _reaction("r2", message_helpers.build_compound(smiles="CCO", role="reactant")),
    )
    first = _table(
        _reaction("r1", message_helpers.build_compound(smiles="CCO", role="solvent"))
    )
    second = _table(
        _reaction("r2", message_helpers.build_compound(smiles="CCO", role="reactant"))
    )
    first.merge(second)
    assert list(first.rows()) == list(whole.rows())
    assert first.reaction_count == 2
    # Merging copies, so adding to one table afterwards leaves the other alone.
    second.add_reaction(
        _reaction("r3", message_helpers.build_compound(smiles="CCO", role="reactant"))
    )
    assert _record(first, "smiles:CCO").reaction_count == 2


def test_the_csv_has_a_column_for_each_identifier_type_in_use(tmp_path):
    table = _table(
        _reaction(
            "r1",
            _with_cas(
                message_helpers.build_compound(smiles="CCO", name="ethanol"), "64-17-5"
            ),
        ),
        _reaction("r2", message_helpers.build_compound(name="brine", role="workup")),
    )
    path = tmp_path / "compounds.csv"
    table.write_csv(path)
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert {"SMILES", "NAME", "CAS_NUMBER"} <= set(rows[0])
    assert "INCHI" not in rows[0]
    ethanol = next(row for row in rows if row["key"] == "smiles:CCO")
    assert ethanol["CAS_NUMBER"] == "64-17-5"
    assert ethanol["roles"] == "UNSPECIFIED (1)"


def test_an_empty_table_still_writes_a_csv_header(tmp_path):
    path = tmp_path / "compounds.csv"
    compounds.CompoundTable().write_csv(path)
    assert path.read_text().splitlines() == [
        "key,label,smiles,roles,locations,reaction_count,"
        "smiles_variants,details,conflicts,example_reaction_ids"
    ]
    # A non-empty table has the same columns, plus one per identifier type in use.
    table = _table(
        _reaction("r1", message_helpers.build_compound(smiles="CCO", role="solvent"))
    )
    table.write_csv(path)
    header = path.read_text().splitlines()[0].split(",")
    assert [c for c in header if c != "SMILES"] == [
        "key",
        "label",
        "smiles",
        "roles",
        "locations",
        "reaction_count",
        "smiles_variants",
        "details",
        "conflicts",
        "example_reaction_ids",
    ]


def test_an_inchi_or_molblock_is_not_a_smiles_variant():
    ethanol = message_helpers.build_compound(smiles="CCO", role="solvent")
    ethanol.identifiers.add(type="INCHI", value="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    benzene = message_helpers.build_compound(smiles="c1ccccc1", role="solvent")
    benzene.identifiers.add(
        type="MOLBLOCK",
        value=Chem.MolToMolBlock(Chem.MolFromSmiles("c1ccccc1")),
    )
    table = _table(_reaction("r1", ethanol, benzene))
    for record in table.records():
        assert len(record.smiles_variants) == 1
    assert all(row["smiles_variants"] == "" for row in table.rows())


def test_the_html_page_is_self_contained_and_escapes_what_it_shows(tmp_path):
    table = _table(
        _reaction(
            "r1",
            message_helpers.build_compound(
                smiles="CCO", name="<b>ethanol</b>", role="solvent"
            ),
        ),
        _reaction("r2", message_helpers.build_compound(name="brine", role="workup")),
    )
    path = tmp_path / "compounds.html"
    table.write_html(path, title="Test dataset")
    page = path.read_text()
    assert page.startswith("<!DOCTYPE html>")
    assert "<title>Test dataset</title>" in page
    assert "&lt;b&gt;ethanol&lt;/b&gt;" in page
    assert "<b>ethanol</b>" not in page
    # A structure is drawn inline; nothing is fetched from elsewhere.
    assert "<svg" in page
    assert "src=" not in page
    assert "@media print" in page
    assert "<h2>SOLVENT (1)</h2>" in page
    assert "<h2>WORKUP (1)</h2>" in page


def test_the_cross_role_index_is_capped_like_the_role_tables():
    reactions = []
    for i, smiles in enumerate(["CCO", "CCCO", "CCCCO"]):
        reactions.append(
            _reaction(
                f"a{i}", message_helpers.build_compound(smiles=smiles, role="solvent")
            )
        )
        reactions.append(
            _reaction(
                f"b{i}", message_helpers.build_compound(smiles=smiles, role="reactant")
            )
        )
    sections = _table(*reactions).to_html_sections(max_rows=1)
    cross = sections.split("<section>")[1]
    assert cross.count("<tr><td>") == 1
    assert "2 more compounds used in more than one role" in cross
