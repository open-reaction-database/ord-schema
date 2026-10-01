# Copyright 2025 Open Reaction Database Project Authors
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
"""Tests for convert_udm_to_ord."""

import pathlib
import textwrap
from unittest.mock import patch

import pytest

from ord_schema import datasets as ds_io
from ord_schema import validations
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import convert_udm_to_ord as conv

TESTDATA = pathlib.Path(__file__).parent / "testdata"
SAMPLE_UDM = TESTDATA / "sample_udm.xml"

# sample_udm.xml has no SCIENTIST / CREATION_DATE; CLI fills ORD-required gaps.
SAMPLE_DEPOSITOR = {
    "person_name": "Test Scientist",
    "email": "test.scientist@example.com",
    "created_date": "2024-01-15",
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def dataset():
    """Converts sample_udm.xml once and reuses the result across tests."""
    return conv.convert(SAMPLE_UDM, **SAMPLE_DEPOSITOR)


# ---------------------------------------------------------------------------
# Basic structure
# ---------------------------------------------------------------------------


def test_basic_conversion_two_reactions(dataset):
    """One reaction with two variations → two ORD Reactions."""
    assert len(dataset.reactions) == 2


def test_dataset_name(dataset):
    assert dataset.name == "Test Dataset"


def test_dataset_description_contains_doi(dataset):
    assert "10.1000/test.doi" in dataset.description


def test_reaction_ids_are_canonical(dataset):
    """update_dataset should assign 'ord-' prefixed IDs."""
    for rxn in dataset.reactions:
        assert rxn.reaction_id.startswith("ord-"), rxn.reaction_id


# ---------------------------------------------------------------------------
# Identifiers
# ---------------------------------------------------------------------------


def test_reaction_identifier_rsmiles(dataset):
    rxn = dataset.reactions[0]
    types = [i.type for i in rxn.identifiers]
    assert reaction_pb2.ReactionIdentifier.REACTION_SMILES in types


# ---------------------------------------------------------------------------
# Inputs and roles
# ---------------------------------------------------------------------------


def test_reactant_role(dataset):
    """The reactant compound must carry the REACTANT role."""
    rxn = dataset.reactions[0]
    roles = [
        c.reaction_role
        for inp in rxn.inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.REACTANT
    ]
    assert len(roles) >= 1


def test_reagent_role(dataset):
    rxn = dataset.reactions[0]
    roles = [
        c.reaction_role
        for inp in rxn.inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.REAGENT
    ]
    assert len(roles) >= 1


def test_catalyst_role(dataset):
    rxn = dataset.reactions[0]
    roles = [
        c.reaction_role
        for inp in rxn.inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.CATALYST
    ]
    assert len(roles) >= 1


def test_solvent_role(dataset):
    rxn = dataset.reactions[0]
    roles = [
        c.reaction_role
        for inp in rxn.inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.SOLVENT
    ]
    assert len(roles) >= 1


# ---------------------------------------------------------------------------
# Amount and units
# ---------------------------------------------------------------------------


def test_amount_units_mass_gram(dataset):
    """Reactant amount unit 'g' should map to GRAM."""
    rxn = dataset.reactions[0]
    reactant_components = [
        c
        for inp in rxn.inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.REACTANT
    ]
    assert reactant_components
    comp = reactant_components[0]
    assert comp.amount.HasField("mass")
    assert comp.amount.mass.value == pytest.approx(1.5)
    assert comp.amount.mass.units == reaction_pb2.Mass.GRAM


def test_amount_units_moles_millimole(dataset):
    """Reagent amount unit 'mmol' should map to MILLIMOLE."""
    rxn = dataset.reactions[0]
    reagent_components = [
        c
        for inp in rxn.inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.REAGENT
    ]
    assert reagent_components
    comp = reagent_components[0]
    assert comp.amount.HasField("moles")
    assert comp.amount.moles.units == reaction_pb2.Moles.MILLIMOLE


def test_amount_units_volume_milliliter(dataset):
    """Solvent amount unit 'ml' should map to MILLILITER."""
    rxn = dataset.reactions[0]
    solvent_components = [
        c
        for inp in rxn.inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.SOLVENT
    ]
    assert solvent_components
    comp = solvent_components[0]
    assert comp.amount.HasField("volume")
    assert comp.amount.volume.units == reaction_pb2.Volume.MILLILITER


def test_combined_amount_string_is_parsed(tmp_path):
    """SURF-style <AMOUNT>0.3000 mmol</AMOUNT> splits into value + unit."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <REACTANT>
                  <MOLECULE MOL_ID="M1"/>
                  <AMOUNT>0.3000 mmol</AMOUNT>
                </REACTANT>
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "combined_amount.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    reactant = next(
        c
        for inp in dataset.reactions[0].inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.REACTANT
    )
    assert reactant.amount.moles.value == pytest.approx(0.3)
    assert reactant.amount.moles.units == reaction_pb2.Moles.MILLIMOLE


def test_section_unwrap_maps_inputs_and_products(tmp_path):
    """SURF VARIATION/SECTION chemistry is promoted to variation level."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES>
            <MOLECULE ID="M1"><NAME>Reactant</NAME></MOLECULE>
            <MOLECULE ID="M2"><NAME>Product</NAME></MOLECULE>
          </MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION>
                <SECTION>
                  <CONDITIONS>
                    <TEMPERATURE><exact>60</exact></TEMPERATURE>
                    <TIME><exact>0.5</exact></TIME>
                  </CONDITIONS>
                  <REACTANT>
                    <MOLECULE MOL_ID="M1"/>
                    <AMOUNT>1.0 mmol</AMOUNT>
                  </REACTANT>
                  <PRODUCT>
                    <MOLECULE MOL_ID="M2"/>
                    <YIELD>13.0</YIELD>
                  </PRODUCT>
                </SECTION>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "section.xml"
    p.write_text(xml)
    dataset = conv.convert(
        p,
        email="a@b.com",
        person_name="A",
        created_date="2024-01-15",
    )
    rxn = dataset.reactions[0]
    assert len(rxn.inputs) == 1
    reactant = next(iter(rxn.inputs.values())).components[0]
    assert reactant.amount.moles.value == pytest.approx(1.0)
    assert rxn.conditions.temperature.setpoint.value == pytest.approx(60.0)
    assert (
        rxn.conditions.temperature.setpoint.units == reaction_pb2.Temperature.CELSIUS
    )
    assert len(rxn.outcomes) == 1
    assert rxn.outcomes[0].reaction_time.value == pytest.approx(0.5)
    assert rxn.outcomes[0].reaction_time.units == reaction_pb2.Time.HOUR
    product = rxn.outcomes[0].products[0]
    assert product.measurements[0].percentage.value == pytest.approx(13.0)


def test_surf_template_converts_and_validates():
    """Real SURF-shaped fixture converts with depositor/dataset CLI fills."""
    surf = TESTDATA / "surf_template.xml"
    dataset = conv.convert(
        surf,
        name="SURF template",
        description="SURF Minisci template dataset",
        email="test@example.com",
        person_name="Test Scientist",
        created_date="2024-01-15",
    )
    assert len(dataset.reactions) == 5
    assert all(rxn.inputs for rxn in dataset.reactions)
    assert all(rxn.outcomes for rxn in dataset.reactions)
    assert dataset.reactions[0].provenance.doi == "10.1021/jo9010624"
    validations.validate_datasets({"_COMBINED": dataset})


# ---------------------------------------------------------------------------
# Conditions
# ---------------------------------------------------------------------------


def test_temperature_value(dataset):
    assert dataset.reactions[0].conditions.temperature.setpoint.value == pytest.approx(
        80.0
    )


def test_temperature_units(dataset):
    assert (
        dataset.reactions[0].conditions.temperature.setpoint.units
        == reaction_pb2.Temperature.CELSIUS
    )


def test_pressure_value(dataset):
    assert dataset.reactions[0].conditions.pressure.setpoint.value == pytest.approx(1.0)


def test_reflux(dataset):
    assert dataset.reactions[0].conditions.reflux is True


def test_ph(dataset):
    assert dataset.reactions[0].conditions.ph == pytest.approx(7.0)


def test_stirring_details(dataset):
    assert "600" in dataset.reactions[0].conditions.stirring.details


def test_stirring_type_inferred_from_details(dataset):
    """Structured UDM stirringRange has no method keyword; type is CUSTOM."""
    assert (
        dataset.reactions[0].conditions.stirring.type
        == reaction_pb2.StirringConditions.CUSTOM
    )


def test_stirring_rpm(dataset):
    assert dataset.reactions[0].conditions.stirring.rate.rpm == 600


def test_environment_fume_hood(dataset):
    env = dataset.reactions[0].setup.environment
    assert env.type == reaction_pb2.ReactionSetup.ReactionEnvironment.FUME_HOOD


# ---------------------------------------------------------------------------
# Outcomes
# ---------------------------------------------------------------------------


def test_reaction_time(dataset):
    rxn = dataset.reactions[0]
    assert rxn.outcomes[0].reaction_time.value == pytest.approx(2.0)
    assert rxn.outcomes[0].reaction_time.units == reaction_pb2.Time.HOUR


def test_yield_value(dataset):
    rxn = dataset.reactions[0]
    products = rxn.outcomes[0].products
    assert products
    measurements = [
        m
        for m in products[0].measurements
        if m.type == reaction_pb2.ProductMeasurement.YIELD
    ]
    assert measurements
    assert measurements[0].percentage.value == pytest.approx(85.0)


def test_second_variation_yield(dataset):
    rxn = dataset.reactions[1]
    products = rxn.outcomes[0].products
    assert products
    measurements = [
        m
        for m in products[0].measurements
        if m.type == reaction_pb2.ProductMeasurement.YIELD
    ]
    assert measurements
    assert measurements[0].percentage.value == pytest.approx(62.0)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def test_provenance_doi(dataset):
    assert dataset.reactions[0].provenance.doi == "10.1000/test.doi"


def test_provenance_organization(dataset):
    assert dataset.reactions[0].provenance.experimenter.organization == "Test Lab"


def test_provenance_scientist(dataset):
    assert dataset.reactions[0].provenance.experimenter.name == "Test Scientist"


def test_provenance_creation_date(dataset):
    assert "2024-01-15" in dataset.reactions[0].provenance.record_created.time.value


def test_is_not_mined(dataset):
    assert dataset.reactions[0].provenance.is_mined is False


def test_email_is_recorded_on_provenance(dataset):
    """CLI --email fills provenance when sample UDM has no SCIENTIST."""
    prov = dataset.reactions[0].provenance
    assert prov.record_created.person.email == "test.scientist@example.com"
    assert prov.experimenter.email == "test.scientist@example.com"


def test_email_is_recorded_on_modification_records(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <SCIENTIST>
                  <NAME>Test Scientist</NAME>
                  <EMAIL>test@example.com</EMAIL>
                </SCIENTIST>
                <MODIFICATION_DATE>2024-03-01</MODIFICATION_DATE>
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "moddate.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    prov = dataset.reactions[0].provenance
    modified = [e for e in prov.record_modified if "2024-03-01" in e.time.value]
    assert modified, "no record_modified event for the UDM modification date"
    assert modified[0].person.email == "test@example.com"


def test_plain_scientist_string_still_maps_name(tmp_path):
    """Legacy bare-name SCIENTIST strings keep working for the name field."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <SCIENTIST>Legacy Name</SCIENTIST>
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "plain_scientist.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert dataset.reactions[0].provenance.experimenter.name == "Legacy Name"
    assert dataset.reactions[0].provenance.experimenter.email == ""


def test_cli_person_flags_fill_missing_scientist(tmp_path):
    """CLI depositor flags populate provenance when UDM has no SCIENTIST."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "no_scientist.xml"
    p.write_text(xml)
    dataset = conv.convert(
        p,
        username="ada",
        person_name="Ada Lovelace",
        orcid="0000-0002-1825-0097",
        email="ada@example.com",
        created_date="2024-01-15",
    )
    person = dataset.reactions[0].provenance.record_created.person
    assert person.username == "ada"
    assert person.name == "Ada Lovelace"
    assert person.orcid == "0000-0002-1825-0097"
    assert person.email == "ada@example.com"
    assert "2024-01-15" in dataset.reactions[0].provenance.record_created.time.value
    experimenter = dataset.reactions[0].provenance.experimenter
    assert experimenter.username == "ada"
    assert experimenter.name == "Ada Lovelace"
    assert experimenter.orcid == "0000-0002-1825-0097"
    assert experimenter.email == "ada@example.com"


def test_udm_creation_date_wins_over_cli(tmp_path):
    """UDM CREATION_DATE is preferred over --created-date."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <CREATION_DATE>2020-06-01</CREATION_DATE>
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "creation_date.xml"
    p.write_text(xml)
    dataset = conv.convert(p, created_date="2024-01-15", email="a@b.com", person_name="A")
    assert "2020-06-01" in dataset.reactions[0].provenance.record_created.time.value


def test_validation_flag_hints_cover_common_gaps():
    """Validation failure text maps to the depositor CLI flags UDM often lacks."""
    email_hint = conv._validation_flag_hints("User email is required for record_created")
    assert "--email" in email_hint
    time_hint = conv._validation_flag_hints("RecordEvent must have `time` specified")
    assert "--created-date" in time_hint
    person_hint = conv._validation_flag_hints(
        "Person must have at least one of `username`, `name`, or `orcid` specified"
    )
    assert "--person-name" in person_hint
    assert conv._validation_flag_hints("unrelated error") == ""


def test_udm_scientist_wins_over_cli_person_flags(tmp_path):
    """UDM SCIENTIST name/email are preferred over CLI fallbacks."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <SCIENTIST>
                  <NAME>From UDM</NAME>
                  <EMAIL>udm@example.com</EMAIL>
                </SCIENTIST>
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "udm_wins.xml"
    p.write_text(xml)
    dataset = conv.convert(
        p,
        person_name="From CLI",
        email="cli@example.com",
        username="cliuser",
    )
    person = dataset.reactions[0].provenance.record_created.person
    assert person.name == "From UDM"
    assert person.email == "udm@example.com"
    assert person.username == "cliuser"


def test_cli_person_name_flag_is_distinct_from_dataset_name(tmp_path):
    """--person-name sets Person.name; --name still sets Dataset.name."""
    out = tmp_path / "out.pbtxt"
    args = conv.parse_args(
        [
            "--input",
            str(SAMPLE_UDM),
            "--output",
            str(out),
            "--name",
            "Dataset Title",
            "--person-name",
            "Person Display",
            "--no-validate",
        ]
    )
    assert args.name == "Dataset Title"
    assert args.person_name == "Person Display"
    conv.main(args)
    loaded = ds_io.load_dataset(out, as_dataset=True)
    assert loaded.name == "Dataset Title"


def test_sample_dataset_validates():
    """The sample UDM file converts to a dataset ORD validation accepts."""
    dataset = conv.convert(SAMPLE_UDM, **SAMPLE_DEPOSITOR)
    validations.validate_datasets({"_COMBINED": dataset})


def test_structureless_molecule_gets_name_identifier(dataset):
    """A CUSTOM identifier would need a details string; NAME is the right type."""
    solvent = next(
        c
        for inp in dataset.reactions[0].inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.SOLVENT
    )
    assert [i.type for i in solvent.identifiers] == [
        reaction_pb2.CompoundIdentifier.NAME
    ]
    assert solvent.identifiers[0].value == "Ethanol"


def test_sample_reactant_has_molblock(dataset):
    """Realistic MOLSTRUCTURE in the sample becomes a MOLBLOCK identifier."""
    reactant = next(
        c
        for inp in dataset.reactions[0].inputs.values()
        for c in inp.components
        if c.reaction_role == reaction_pb2.ReactionRole.REACTANT
    )
    types = [i.type for i in reactant.identifiers]
    assert reaction_pb2.CompoundIdentifier.MOLBLOCK in types
    assert any("M  END" in i.value for i in reactant.identifiers)


# ---------------------------------------------------------------------------
# Both variations are independent reactions
# ---------------------------------------------------------------------------


def test_two_variations_independent(dataset):
    """Each variation produces a fully independent Reaction (not sharing state)."""
    r0, r1 = dataset.reactions[0], dataset.reactions[1]
    assert r0.reaction_id != r1.reaction_id
    assert r0.conditions.temperature.setpoint.value == pytest.approx(80.0)
    assert r1.conditions.temperature.setpoint.value == pytest.approx(60.0)


# ---------------------------------------------------------------------------
# Robustness: missing / empty elements
# ---------------------------------------------------------------------------


def test_missing_molecules_does_not_crash(tmp_path):
    """A UDM file with no MOLECULES element should parse without crashing."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Empty</TITLE></LEGAL>
          <MOLECULES/>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1"/>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "empty.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert len(dataset.reactions) == 1


def test_missing_reactions_exits(tmp_path):
    """A UDM file with no <REACTIONS> element should raise SystemExit."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>NoReactions</TITLE></LEGAL>
          <MOLECULES/>
        </UDM>
    """)
    p = tmp_path / "no_reactions.xml"
    p.write_text(xml)
    with pytest.raises(SystemExit):
        conv.convert(p)


def test_invalid_xml_exits(tmp_path):
    """A non-XML file should raise SystemExit."""
    p = tmp_path / "bad.xml"
    p.write_text("not xml at all <<<")
    with pytest.raises(SystemExit):
        conv.convert(p)


def test_wrong_root_exits(tmp_path):
    """A well-formed XML file without <UDM> root should raise SystemExit."""
    xml = '<?xml version="1.0"?><ROOT/>'
    p = tmp_path / "wrong_root.xml"
    p.write_text(xml)
    with pytest.raises(SystemExit):
        conv.convert(p)


def test_missing_file_exits(tmp_path):
    with pytest.raises(SystemExit):
        conv.convert(tmp_path / "nonexistent.xml")


# ---------------------------------------------------------------------------
# CLI round-trip
# ---------------------------------------------------------------------------


def test_cli_roundtrip(tmp_path):
    """CLI converts sample_udm.xml to .pbtxt and the file can be loaded back."""
    out = tmp_path / "out.pbtxt"
    args = conv.parse_args(
        ["--input", str(SAMPLE_UDM), "--output", str(out), "--no-validate"]
    )
    conv.main(args)
    assert out.exists()
    loaded = ds_io.load_dataset(out, as_dataset=True)
    assert len(loaded.reactions) == 2


def test_cli_name_override(tmp_path):
    """--name flag overrides the dataset name from the UDM file."""
    out = tmp_path / "out.pbtxt"
    args = conv.parse_args(
        [
            "--input",
            str(SAMPLE_UDM),
            "--output",
            str(out),
            "--name",
            "Custom Name",
            "--no-validate",
        ]
    )
    conv.main(args)
    loaded = ds_io.load_dataset(out, as_dataset=True)
    assert loaded.name == "Custom Name"


# ---------------------------------------------------------------------------
# Failure 1: AMOUNT with inline XML attribute
# ---------------------------------------------------------------------------


def test_amount_inline_unit_attribute(tmp_path):
    """<AMOUNT units="g">2.5</AMOUNT> must be parsed, not silently dropped."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES>
            <MOLECULE ID="MOL-1"><NAME>Benzaldehyde</NAME></MOLECULE>
          </MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <REACTANT>
                  <MOLECULE MOL_ID="MOL-1"/>
                  <AMOUNT units="g">2.5</AMOUNT>
                </REACTANT>
                <PRODUCT>
                  <MOLECULE MOL_ID="MOL-1"/>
                </PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "inline.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    rxn = dataset.reactions[0]
    components = [c for inp in rxn.inputs.values() for c in inp.components]
    assert components, "no components found"
    comp = components[0]
    assert comp.amount.HasField("mass")
    assert comp.amount.mass.value == pytest.approx(2.5)
    assert comp.amount.mass.units == reaction_pb2.Mass.GRAM


# ---------------------------------------------------------------------------
# Failure 2: Multiple CONDITION_GROUPs
# ---------------------------------------------------------------------------


def test_multiple_condition_groups_does_not_crash(tmp_path):
    """Multiple CONDITION_GROUPs must not crash; first group's temperature is used."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
                <CONDITIONS>
                  <CONDITION_GROUP>
                    <TEMPERATURE units="c"><exact>80</exact></TEMPERATURE>
                  </CONDITION_GROUP>
                  <CONDITION_GROUP>
                    <TEMPERATURE units="c"><exact>100</exact></TEMPERATURE>
                  </CONDITION_GROUP>
                </CONDITIONS>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "multi_cg.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert len(dataset.reactions) == 1
    assert dataset.reactions[0].conditions.temperature.setpoint.value == pytest.approx(
        80.0
    )


# ---------------------------------------------------------------------------
# Failure 3: Unsafe filename sanitisation
# ---------------------------------------------------------------------------


def test_safe_filename_strips_path_separators():
    assert conv._safe_filename("Pd/C Hydrogenation") == "Pd_C Hydrogenation"


def test_safe_filename_traversal():
    result = conv._safe_filename("../../etc/cron.d/exploit")
    assert "/" not in result
    assert "\\" not in result


def test_safe_filename_special_chars():
    result = conv._safe_filename("Reaction: A → B")
    assert ":" not in result


def test_cli_title_with_slash_does_not_escape(tmp_path, monkeypatch):
    """Dataset name with '/' must not create a subdirectory."""
    monkeypatch.chdir(tmp_path)
    args = conv.parse_args(
        ["--input", str(SAMPLE_UDM), "--name", "A/B", "--no-validate"]
    )
    conv.main(args)
    assert (tmp_path / "A_B.pbtxt").exists()
    assert not (tmp_path / "A").is_dir()


# ---------------------------------------------------------------------------
# Failure 4: Modification date as plain string
# ---------------------------------------------------------------------------


def test_modification_date_string_is_preserved(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <MODIFICATION_DATE>2024-03-01</MODIFICATION_DATE>
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "moddate.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    prov = dataset.reactions[0].provenance
    assert any("2024-03-01" in e.time.value for e in prov.record_modified)


# ---------------------------------------------------------------------------
# Failure 5: No empty outcome for input-only variations
# ---------------------------------------------------------------------------


def test_no_outcome_when_no_products_or_duration(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <REACTANT>
                  <MOLECULE MOL_ID="M1"/>
                  <AMOUNT>1.0</AMOUNT>
                  <AMOUNT_UNIT>g</AMOUNT_UNIT>
                </REACTANT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "no_outcome.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert len(dataset.reactions[0].outcomes) == 0


# ---------------------------------------------------------------------------
# Failure 6: Same mol_id in multiple roles → separate input slots
# ---------------------------------------------------------------------------


def test_same_mol_in_two_roles_gets_separate_inputs(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="MOL-W"><NAME>Water</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <REACTANT>
                  <MOLECULE MOL_ID="MOL-W"/>
                  <AMOUNT>1.0</AMOUNT><AMOUNT_UNIT>g</AMOUNT_UNIT>
                </REACTANT>
                <SOLVENT>
                  <MOLECULE MOL_ID="MOL-W"/>
                  <AMOUNT>10</AMOUNT><AMOUNT_UNIT>ml</AMOUNT_UNIT>
                </SOLVENT>
                <PRODUCT>
                  <MOLECULE MOL_ID="MOL-W"/>
                </PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "dual_role.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    rxn = dataset.reactions[0]
    assert len(rxn.inputs) == 2
    roles = {c.reaction_role for inp in rxn.inputs.values() for c in inp.components}
    assert reaction_pb2.ReactionRole.REACTANT in roles
    assert reaction_pb2.ReactionRole.SOLVENT in roles


# ---------------------------------------------------------------------------
# Failure 7: Non-finite float values dropped
# ---------------------------------------------------------------------------


def test_infinite_amount_is_dropped(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <REACTANT>
                  <MOLECULE MOL_ID="M1"/>
                  <AMOUNT>1e999</AMOUNT>
                  <AMOUNT_UNIT>g</AMOUNT_UNIT>
                </REACTANT>
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "inf_amount.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    all_comps = [
        c for inp in dataset.reactions[0].inputs.values() for c in inp.components
    ]
    comp = next(iter(all_comps))
    assert not comp.amount.HasField("mass")
    assert not comp.amount.HasField("moles")
    assert not comp.amount.HasField("volume")


def test_nan_temperature_is_skipped(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
                <CONDITIONS>
                  <CONDITION_GROUP>
                    <TEMPERATURE units="c"><exact>nan</exact></TEMPERATURE>
                  </CONDITION_GROUP>
                </CONDITIONS>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "nan_temp.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    # Proto default for unset float is 0.0 — confirms value was never assigned.
    assert dataset.reactions[0].conditions.temperature.setpoint.value == pytest.approx(
        0.0
    )


# ---------------------------------------------------------------------------
# Failure 8: MOLSTRUCTURE with format attribute
# ---------------------------------------------------------------------------


_BENZALDEHYDE_MOLBLOCK = """
     RDKit          2D

  8  8  0  0  0  0  0  0  0  0999 V2000
   -3.0851    0.4695    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0677   -0.6327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6045   -0.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1586    1.1294    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3047    1.4594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3221    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8762   -1.0750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4129   -1.4050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  3  1  0
M  END
"""


def _udm_with_molstructure(molstructure: str) -> str:
    """Builds a one-variation UDM file whose molecule carries ``molstructure``."""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<UDM version="6.0.0">
  <LEGAL><TITLE>Test</TITLE></LEGAL>
  <MOLECULES>
    <MOLECULE ID="M1">
      <NAME>Benzaldehyde</NAME>
      <MOLSTRUCTURE format="mol">{molstructure}</MOLSTRUCTURE>
    </MOLECULE>
  </MOLECULES>
  <REACTIONS>
    <REACTION ID="R1">
      <VARIATION ID="V1">
        <REACTANT>
          <MOLECULE MOL_ID="M1"/>
          <AMOUNT>1.0</AMOUNT><AMOUNT_UNIT>g</AMOUNT_UNIT>
        </REACTANT>
        <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
      </VARIATION>
    </REACTION>
  </REACTIONS>
</UDM>
"""


def _first_component(dataset):
    """Returns the first input component of the first reaction."""
    return next(
        c for inp in dataset.reactions[0].inputs.values() for c in inp.components
    )


def test_molstructure_with_format_attribute(tmp_path):
    """A readable MolBlock becomes a MOLBLOCK identifier, attributes notwithstanding."""
    p = tmp_path / "molblock_attr.xml"
    p.write_text(_udm_with_molstructure(_BENZALDEHYDE_MOLBLOCK))
    comp = _first_component(conv.convert(p))
    molblock_ids = [
        i
        for i in comp.identifiers
        if i.type == reaction_pb2.CompoundIdentifier.MOLBLOCK
    ]
    assert molblock_ids, "no MOLBLOCK identifier found"
    assert "M  END" in molblock_ids[0].value
    assert "{" not in molblock_ids[0].value, "dict repr leaked into identifier value"


def test_indented_molstructure_is_recovered(tmp_path):
    """A MolBlock indented by an XML formatter is dedented back into shape."""
    indented = textwrap.indent(_BENZALDEHYDE_MOLBLOCK, "    ")
    p = tmp_path / "molblock_indented.xml"
    p.write_text(_udm_with_molstructure(indented))
    comp = _first_component(conv.convert(p))
    molblock_ids = [
        i
        for i in comp.identifiers
        if i.type == reaction_pb2.CompoundIdentifier.MOLBLOCK
    ]
    assert molblock_ids, "indented MolBlock was not recovered"
    assert conv._molblock_parses(molblock_ids[0].value)


def test_unreadable_molstructure_falls_back_to_name(tmp_path):
    """An unreadable MOLSTRUCTURE must not be recorded as a MOLBLOCK; ORD rejects it."""
    p = tmp_path / "molblock_bad.xml"
    p.write_text(_udm_with_molstructure("molblock content here"))
    comp = _first_component(conv.convert(p))
    types = [i.type for i in comp.identifiers]
    assert reaction_pb2.CompoundIdentifier.MOLBLOCK not in types
    assert types == [reaction_pb2.CompoundIdentifier.NAME]
    assert comp.identifiers[0].value == "Benzaldehyde"


# ---------------------------------------------------------------------------
# Failure 9: ValidationError exits cleanly
# ---------------------------------------------------------------------------


def test_validation_errors_exit_cleanly(tmp_path):
    """validate_datasets raising ValidationError must cause sys.exit(1)."""
    out = tmp_path / "out.pbtxt"
    args = conv.parse_args(["--input", str(SAMPLE_UDM), "--output", str(out)])
    with (
        patch(
            "ord_schema.scripts.convert_udm_to_ord.validations.validate_datasets",
            side_effect=validations.ValidationError("forced error"),
        ),
        pytest.raises(SystemExit) as exc_info,
    ):
        conv.main(args)
    assert exc_info.value.code == 1
    assert not out.exists()


def test_no_validate_writes_despite_mock_error(tmp_path):
    """--no-validate skips validation entirely; file is written regardless."""
    out = tmp_path / "out.pbtxt"
    args = conv.parse_args(
        ["--input", str(SAMPLE_UDM), "--output", str(out), "--no-validate"]
    )
    with patch(
        "ord_schema.scripts.convert_udm_to_ord.validations.validate_datasets",
        side_effect=validations.ValidationError("forced error"),
    ):
        conv.main(args)  # must NOT raise — validation is skipped
    assert out.exists()


# ---------------------------------------------------------------------------
# Failure 10: Empty <REACTIONS/> element
# ---------------------------------------------------------------------------


def test_empty_reactions_element_does_not_crash(tmp_path):
    """<REACTIONS/> (self-closing) must not raise AttributeError."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Empty</TITLE></LEGAL>
          <MOLECULES/>
          <REACTIONS/>
        </UDM>
    """)
    p = tmp_path / "empty_reactions.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert len(dataset.reactions) == 0


# ---------------------------------------------------------------------------
# Failure 11: Text elements with XML attributes
# ---------------------------------------------------------------------------


def test_producer_with_attribute_does_not_crash(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL>
            <TITLE>Test</TITLE>
            <PRODUCER id="P1">Test Lab</PRODUCER>
          </LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "producer_attr.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert dataset.reactions[0].provenance.experimenter.organization == "Test Lab"


def test_name_with_attribute_does_not_crash(tmp_path):
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES>
            <MOLECULE ID="M1">
              <NAME lang="en">Benzaldehyde</NAME>
            </MOLECULE>
          </MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "name_attr.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    rxn = dataset.reactions[0]
    product_ids = [i.value for p in rxn.outcomes[0].products for i in p.identifiers]
    assert any("Benzaldehyde" in v for v in product_ids)


# ---------------------------------------------------------------------------
# Failure 12: <PH>7.0</PH> plain string
# ---------------------------------------------------------------------------


def test_ph_plain_string_is_parsed(tmp_path):
    """<PH>7.0</PH> (no <exact> child) must be parsed, not silently dropped."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <VARIATION ID="V1">
                <PRODUCT><MOLECULE MOL_ID="M1"/></PRODUCT>
                <CONDITIONS>
                  <CONDITION_GROUP>
                    <PH>7.0</PH>
                  </CONDITION_GROUP>
                </CONDITIONS>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "ph_plain.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert dataset.reactions[0].conditions.ph == pytest.approx(7.0)


# ---------------------------------------------------------------------------
# Reaxys-style incomplete UDM policies
# ---------------------------------------------------------------------------


def test_reaction_product_id_creates_outcome(tmp_path):
    """REACTION/PRODUCT_ID fills outcomes when VARIATION has no PRODUCT block."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES>
            <MOLECULE ID="969129"><NAME>product-a</NAME></MOLECULE>
            <MOLECULE ID="R1"><NAME>reagent</NAME></MOLECULE>
          </MOLECULES>
          <REACTIONS>
            <REACTION ID="RXN">
              <PRODUCT_ID>969129</PRODUCT_ID>
              <VARIATION>
                <REAGENT><MOLECULE MOL_ID="R1"/></REAGENT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "product_id.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert len(dataset.reactions[0].outcomes) == 1
    assert dataset.reactions[0].outcomes[0].products[0].identifiers[0].value == (
        "product-a"
    )


def test_missing_amount_becomes_unmeasured_custom(tmp_path):
    """Inputs without AMOUNT get UnmeasuredAmount CUSTOM (ORD requires an amount)."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>water</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <PRODUCT_ID>M1</PRODUCT_ID>
              <VARIATION>
                <REAGENT><MOLECULE MOL_ID="M1"/></REAGENT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "no_amount.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    reagent = next(iter(dataset.reactions[0].inputs.values())).components[0]
    assert reagent.amount.WhichOneof("kind") == "unmeasured"
    assert reagent.amount.unmeasured.type == reaction_pb2.UnmeasuredAmount.CUSTOM
    assert "not reported" in reagent.amount.unmeasured.details


def test_pressure_without_unit_is_not_setpoint(tmp_path):
    """Unitless PRESSURE is not written to setpoint; raw value goes to details."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <PRODUCT_ID>M1</PRODUCT_ID>
              <VARIATION>
                <REAGENT><MOLECULE MOL_ID="M1"/></REAGENT>
                <CONDITIONS>
                  <CONDITION_GROUP>
                    <PRESSURE><exact>760.051</exact></PRESSURE>
                  </CONDITION_GROUP>
                </CONDITIONS>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "unitless_pressure.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    pressure = dataset.reactions[0].conditions.pressure
    assert not pressure.setpoint.HasField("value")
    assert "760.051" in dataset.reactions[0].conditions.details
    assert "unit omitted" in dataset.reactions[0].conditions.details


def test_empty_molecule_name_falls_back_to_mol_id(tmp_path):
    """Structureless empty-NAME molecules use the mol ID as NAME identifier."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES>
            <MOLECULE ID="28870379"><NAME/></MOLECULE>
            <MOLECULE ID="P1"><NAME>product</NAME></MOLECULE>
          </MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <PRODUCT_ID>P1</PRODUCT_ID>
              <VARIATION>
                <REAGENT>
                  <MOLECULE MOL_ID="28870379"/>
                  <NAME/>
                </REAGENT>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "empty_name.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    reagent = dataset.reactions[0].inputs["28870379_REAGENT"].components[0]
    assert reagent.identifiers[0].type == reaction_pb2.CompoundIdentifier.NAME
    assert reagent.identifiers[0].value == "28870379"


def test_reaction_reactant_id_creates_inputs(tmp_path):
    """REACTION/REACTANT_ID fills inputs when VARIATION has no role compounds."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES>
            <MOLECULE ID="1209228"><NAME>reactant-a</NAME></MOLECULE>
            <MOLECULE ID="P1"><NAME>product</NAME></MOLECULE>
          </MOLECULES>
          <REACTIONS>
            <REACTION ID="RXN">
              <REACTANT_ID>1209228</REACTANT_ID>
              <PRODUCT_ID>P1</PRODUCT_ID>
              <VARIATION>
                <CITATION CIT_ID="1"/>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "reactant_id.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    assert "1209228_REACTANT" in dataset.reactions[0].inputs
    reactant = dataset.reactions[0].inputs["1209228_REACTANT"].components[0]
    assert reactant.reaction_role == reaction_pb2.ReactionRole.REACTANT
    assert reactant.identifiers[0].value == "reactant-a"
    assert reactant.amount.WhichOneof("kind") == "unmeasured"


def test_free_text_preparation_sets_environment_custom(tmp_path):
    """Unknown PREPARATION text sets environment type=CUSTOM with details."""
    xml = textwrap.dedent("""\
        <?xml version="1.0" encoding="UTF-8"?>
        <UDM version="6.0.0">
          <LEGAL><TITLE>Test</TITLE></LEGAL>
          <MOLECULES><MOLECULE ID="M1"><NAME>A</NAME></MOLECULE></MOLECULES>
          <REACTIONS>
            <REACTION ID="R1">
              <PRODUCT_ID>M1</PRODUCT_ID>
              <REACTANT_ID>M1</REACTANT_ID>
              <VARIATION>
                <CONDITIONS>
                  <PREPARATION>Long literature procedure text without env keyword.</PREPARATION>
                </CONDITIONS>
              </VARIATION>
            </REACTION>
          </REACTIONS>
        </UDM>
    """)
    p = tmp_path / "prep_custom.xml"
    p.write_text(xml)
    dataset = conv.convert(p)
    env = dataset.reactions[0].setup.environment
    assert env.type == reaction_pb2.ReactionSetup.ReactionEnvironment.CUSTOM
    assert "literature procedure" in env.details
    validations.validate_message(env)
