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
"""Tests for ord_schema.scripts.convert_udm_to_ord."""

import getpass
import pathlib

import pytest

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.scripts import convert_udm_to_ord

_MOLECULES = """
  <MOLECULES>
    <MOLECULE ID="m1"><NAME>Reactant A</NAME></MOLECULE>
    <MOLECULE ID="m2"><NAME>Product A</NAME></MOLECULE>
  </MOLECULES>
"""

# SCIENTIST is UDM's authorDetails type (NAME/EMAIL/PHONE/ORGANISATION), not a
# bare string -- see udm_6_0_0.xsd.
_SCIENTIST = """
        <SCIENTIST>
          <NAME>Test Scientist</NAME>
          <EMAIL>scientist@example.com</EMAIL>
          <ORGANISATION>
            <NAME>Scientist Org</NAME>
            <ADDRESS>123 Lab St</ADDRESS>
          </ORGANISATION>
        </SCIENTIST>
"""

_ONE_VARIATION_REACTION = f"""
    <REACTION>
      <VARIATION>
        {_SCIENTIST}
        <CREATION_DATE>2020-01-01</CREATION_DATE>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
          <AMOUNT>1.0</AMOUNT>
        </REACTANT>
        <PRODUCT>
          <MOLECULE MOL_ID="m2"><NAME>Product A</NAME></MOLECULE>
          <YIELD><exact>85.0</exact></YIELD>
        </PRODUCT>
        <CONDITIONS>
          <PREPARATION>Stirred under nitrogen.</PREPARATION>
          <CONDITION_GROUP>
            <TEMPERATURE unit="K"><exact>298.0</exact></TEMPERATURE>
            <PRESSURE unit="bar"><exact>1.5</exact></PRESSURE>
            <STIRRING><exact>500</exact></STIRRING>
          </CONDITION_GROUP>
        </CONDITIONS>
        <COMMENT>A test reaction</COMMENT>
      </VARIATION>
    </REACTION>
"""


def _udm_xml(reactions_xml: str) -> str:
    """Wraps a <REACTION>...</REACTION> fragment in a minimal UDM document."""
    return f"""<UDM>
  <LEGAL>
    <TITLE>Test Dataset</TITLE>
    <PRODUCER>Data Vendor</PRODUCER>
    <DOI>10.0000/test-doi</DOI>
  </LEGAL>
  {_MOLECULES}
  <REACTIONS>
    {reactions_xml}
  </REACTIONS>
</UDM>
"""


def _write(tmp_path: pathlib.Path, filename: str, content: str) -> str:
    path = tmp_path / filename
    path.write_text(content)
    return str(path)


def test_simple(tmp_path):
    input_filename = _write(tmp_path, "input.xml", _udm_xml(_ONE_VARIATION_REACTION))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = [
        "--input",
        input_filename,
        "--output",
        output_filename,
        "--email",
        "submitter@example.com",
    ]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    assert dataset.name == "Test Dataset"
    assert dataset.description == "UDM dataset DOI: 10.0000/test-doi"
    assert len(dataset.reactions) == 1

    reaction = dataset.reactions[0]
    # The experimenter is the original UDM scientist who ran the reaction,
    # with their own organization/email/address taking precedence over the
    # dataset-level PRODUCER (a data vendor, not a research institution).
    assert reaction.provenance.experimenter.name == "Test Scientist"
    assert reaction.provenance.experimenter.email == "scientist@example.com"
    assert reaction.provenance.experimenter.organization == "Scientist Org"
    assert reaction.provenance.city == "123 Lab St"

    # record_created is whoever ran this conversion, not the UDM scientist.
    assert reaction.provenance.record_created.person.username == getpass.getuser()
    assert reaction.provenance.record_created.person.email == "submitter@example.com"
    assert reaction.provenance.record_created.person.name != "Test Scientist"
    assert reaction.provenance.record_created.time.value == "2020-01-01"

    assert reaction.conditions.temperature.setpoint.value == 298.0
    assert (
        reaction.conditions.temperature.setpoint.units
        == reaction.conditions.temperature.setpoint.KELVIN
    )
    assert reaction.conditions.pressure.setpoint.value == 1.5
    assert (
        reaction.conditions.pressure.setpoint.units
        == reaction.conditions.pressure.setpoint.BAR
    )
    assert reaction.conditions.stirring.rate.rpm == 500.0
    assert reaction.setup.environment.details == "Stirred under nitrogen."
    assert reaction.observations[0].comment == "A test reaction"

    [component] = reaction.inputs["m1"].components
    assert component.amount.moles.value == 1.0
    assert component.amount.moles.units == component.amount.moles.MOLE

    [outcome] = reaction.outcomes
    [product] = outcome.products
    assert product.identifiers[0].value == "Product A"
    [measurement] = product.measurements
    assert measurement.float_value.value == 85.0

    # Off by default: no raw UDM XML embedded without --include-udm-xml.
    assert "udm_reaction_xml" not in reaction.provenance.reaction_metadata
    assert "udm_parent_xml" not in reaction.provenance.reaction_metadata


def test_include_udm_xml(tmp_path):
    reaction_xml = """
    <REACTION>
      <VARIATION>
        <SCIENTIST><NAME>Scientist One</NAME></SCIENTIST>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
        </REACTANT>
      </VARIATION>
      <VARIATION>
        <SCIENTIST><NAME>Scientist Two</NAME></SCIENTIST>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
        </REACTANT>
      </VARIATION>
    </REACTION>
"""
    input_filename = _write(tmp_path, "input.xml", _udm_xml(reaction_xml))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = [
        "--input",
        input_filename,
        "--output",
        output_filename,
        "--no-validate",
        "--include-udm-xml",
    ]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    assert len(dataset.reactions) == 2

    for reaction in dataset.reactions:
        metadata = reaction.provenance.reaction_metadata
        assert metadata["udm_reaction_xml"].format == "xml"
        assert "<REACTION>" in metadata["udm_reaction_xml"].string_value
        # Both VARIATIONs came from the same <REACTION>, so they share
        # identical reaction-level XML, including both scientists.
        assert "Scientist One" in metadata["udm_reaction_xml"].string_value
        assert "Scientist Two" in metadata["udm_reaction_xml"].string_value
        assert metadata["udm_parent_xml"].format == "xml"
        assert "Test Dataset" in metadata["udm_parent_xml"].string_value
        assert "Data Vendor" in metadata["udm_parent_xml"].string_value
        # MOLECULES is deliberately excluded: it's a lookup table, not
        # reaction-specific context, and can be large.
        assert "MOLECULES" not in metadata["udm_parent_xml"].string_value

    # The two reactions' udm_reaction_xml differ only by REACTION-vs-VARIATION
    # scope, not by content -- both variations share one <REACTION> element.
    reaction_xmls = {
        r.provenance.reaction_metadata["udm_reaction_xml"].string_value
        for r in dataset.reactions
    }
    assert len(reaction_xmls) == 1


def test_scientist_organization_falls_back_to_producer(tmp_path):
    """When the scientist has no inline ORGANISATION, PRODUCER fills in."""
    reaction_xml = """
    <REACTION>
      <VARIATION>
        <SCIENTIST><NAME>No Org Scientist</NAME></SCIENTIST>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
        </REACTANT>
      </VARIATION>
    </REACTION>
"""
    input_filename = _write(tmp_path, "input.xml", _udm_xml(reaction_xml))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    reaction = dataset.reactions[0]
    assert reaction.provenance.experimenter.name == "No Org Scientist"
    assert reaction.provenance.experimenter.organization == "Data Vendor"


def test_record_created_defaults_to_submitter_without_email(tmp_path):
    """Without --email, record_created still identifies the submitter (by OS
    username) rather than the UDM scientist; the dataset just won't pass
    full validation without an explicit --email, since record_created still
    needs one even though the scientist's own email (if UDM supplies it)
    goes to experimenter, not record_created.
    """
    input_filename = _write(tmp_path, "input.xml", _udm_xml(_ONE_VARIATION_REACTION))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    reaction = dataset.reactions[0]
    assert reaction.provenance.record_created.person.username == getpass.getuser()
    assert not reaction.provenance.record_created.person.email
    assert reaction.provenance.experimenter.name == "Test Scientist"
    assert reaction.provenance.experimenter.email == "scientist@example.com"


def test_multiple_variations_become_separate_reactions(tmp_path):
    """Regression test: each <VARIATION> must produce its own Reaction.

    A prior version of this script reset the in-progress Reaction object on
    every VARIATION and only kept the last one, silently discarding every
    other variation's data.
    """
    reaction_xml = """
    <REACTION>
      <VARIATION>
        <SCIENTIST><NAME>Scientist One</NAME></SCIENTIST>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
          <AMOUNT>1.0</AMOUNT>
        </REACTANT>
        <PRODUCT>
          <MOLECULE MOL_ID="m2"><NAME>Product A</NAME></MOLECULE>
        </PRODUCT>
      </VARIATION>
      <VARIATION>
        <SCIENTIST><NAME>Scientist Two</NAME></SCIENTIST>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
          <AMOUNT>2.0</AMOUNT>
        </REACTANT>
        <PRODUCT>
          <MOLECULE MOL_ID="m2"><NAME>Product A</NAME></MOLECULE>
        </PRODUCT>
      </VARIATION>
    </REACTION>
"""
    input_filename = _write(tmp_path, "input.xml", _udm_xml(reaction_xml))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    assert len(dataset.reactions) == 2
    scientists = {r.provenance.experimenter.name for r in dataset.reactions}
    assert scientists == {"Scientist One", "Scientist Two"}
    amounts = {
        r.inputs["m1"].components[0].amount.moles.value for r in dataset.reactions
    }
    assert amounts == {1.0, 2.0}


def test_multiple_condition_groups_use_first_group_only(tmp_path):
    """Regression test: multiple <CONDITION_GROUP>s must not be merged, and
    must not multiply the Reaction.

    An earlier version of this code merged every group's fields into one
    ReactionConditions (fabricating composites, e.g. one group's temperature
    plus an unrelated group's pressure, that never existed in the source). A
    later version instead split each group into its own Reaction -- also
    wrong, since REACTANT/PRODUCT/YIELD/COMMENT/provenance are variation-
    level, not per-group, so splitting duplicated a single recorded outcome
    across multiple Reactions as if it were independently reproduced under
    each condition set. Only the first group is represented structurally;
    the existence of the rest is noted, not silently dropped.
    """
    reaction_xml = """
    <REACTION>
      <VARIATION>
        <SCIENTIST><NAME>Test Scientist</NAME></SCIENTIST>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
        </REACTANT>
        <PRODUCT>
          <MOLECULE MOL_ID="m2"><NAME>Product A</NAME></MOLECULE>
          <YIELD><exact>85.0</exact></YIELD>
        </PRODUCT>
        <CONDITIONS>
          <CONDITION_GROUP>
            <TEMPERATURE><exact>25.0</exact></TEMPERATURE>
          </CONDITION_GROUP>
          <CONDITION_GROUP>
            <PRESSURE unit="bar"><exact>2.0</exact></PRESSURE>
          </CONDITION_GROUP>
        </CONDITIONS>
      </VARIATION>
    </REACTION>
"""
    input_filename = _write(tmp_path, "input.xml", _udm_xml(reaction_xml))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    assert len(dataset.reactions) == 1
    reaction = dataset.reactions[0]

    # Only the outcome recorded once in the source, not duplicated.
    assert len(reaction.outcomes) == 1
    assert reaction.outcomes[0].products[0].measurements[0].float_value.value == 85.0

    # Only the first group's field is structurally represented.
    assert reaction.conditions.temperature.setpoint.value == 25.0
    assert reaction.conditions.pressure.setpoint.value == 0.0
    assert "2 CONDITION_GROUP" in reaction.conditions.details


def test_reaction_without_variation_still_emits_one_reaction(tmp_path):
    reaction_xml = """
    <REACTION>
      <RXNSTRUCTURE>
        <format>rsmiles</format>
        <value>CC.O&gt;&gt;CCO</value>
      </RXNSTRUCTURE>
    </REACTION>
"""
    input_filename = _write(tmp_path, "input.xml", _udm_xml(reaction_xml))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    assert len(dataset.reactions) == 1
    assert dataset.reactions[0].identifiers[0].value == "CC.O>>CCO"


def test_single_reactant_id_is_not_split_into_characters(tmp_path):
    """Regression test: a lone <REACTANT_ID> must not be iterated as a string.

    UDM's REACTANT_ID/PRODUCT_ID are maxOccurs="unbounded"; etree_to_dict
    collapses a single occurrence to a bare string rather than a one-element
    list, so naive iteration over it would previously produce one bogus
    identifier per character.
    """
    reaction_xml = """
    <REACTION>
      <REACTANT_ID>abc123</REACTANT_ID>
    </REACTION>
"""
    input_filename = _write(tmp_path, "input.xml", _udm_xml(reaction_xml))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    [component] = dataset.reactions[0].inputs["REACTANT_IDS"].components
    assert len(component.identifiers) == 1
    assert component.identifiers[0].value == "abc123"


def test_unknown_molecule_reference_does_not_crash(tmp_path):
    reaction_xml = """
    <REACTION>
      <VARIATION>
        <REACTANT>
          <MOLECULE MOL_ID="does_not_exist"><NAME>Missing</NAME></MOLECULE>
        </REACTANT>
      </VARIATION>
    </REACTION>
"""
    input_filename = _write(tmp_path, "input.xml", _udm_xml(reaction_xml))
    output_filename = str(tmp_path / "output.pbtxt")
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    [component] = dataset.reactions[0].inputs["does_not_exist"].components
    assert component.identifiers[0].type == component.identifiers[0].CUSTOM
    assert component.identifiers[0].value == "does_not_exist"


def test_missing_input_file(tmp_path):
    argv = ["--input", str(tmp_path / "does_not_exist.xml")]
    with pytest.raises(SystemExit):
        convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))


def test_not_udm_format(tmp_path):
    input_filename = _write(tmp_path, "input.xml", "<NOT_UDM></NOT_UDM>")
    argv = ["--input", input_filename]
    with pytest.raises(SystemExit):
        convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))
