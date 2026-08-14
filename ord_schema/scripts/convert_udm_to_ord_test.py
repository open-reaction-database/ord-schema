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

_ONE_VARIATION_REACTION = """
    <REACTION>
      <VARIATION>
        <SCIENTIST>Test Scientist</SCIENTIST>
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
          <CONDITION_GROUP>
            <TEMPERATURE><exact>25.0</exact></TEMPERATURE>
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
    <PRODUCER>Test Org</PRODUCER>
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
    # UDM's SCIENTIST field carries only a name, no email, so a person's
    # email -- required by ORD validation -- can never be filled in from UDM
    # source data alone. Validation is disabled here for that known reason.
    argv = ["--input", input_filename, "--output", output_filename, "--no-validate"]
    convert_udm_to_ord.main(convert_udm_to_ord.parse_args(argv))

    dataset = message_helpers.load_message(output_filename, dataset_pb2.Dataset)
    assert dataset.name == "Test Dataset"
    assert dataset.description == "UDM dataset DOI: 10.0000/test-doi"
    assert len(dataset.reactions) == 1

    reaction = dataset.reactions[0]
    assert reaction.provenance.experimenter.name == "Test Scientist"
    assert reaction.provenance.record_created.time.value == "2020-01-01"
    assert reaction.conditions.temperature.setpoint.value == 25.0
    assert (
        reaction.conditions.temperature.setpoint.units
        == reaction.conditions.temperature.setpoint.CELSIUS
    )
    assert reaction.observations[0].comment == "A test reaction"

    [component] = reaction.inputs["m1"].components
    assert component.amount.mass.value == 1.0
    assert component.amount.mass.units == component.amount.mass.GRAM

    [outcome] = reaction.outcomes
    [product] = outcome.products
    assert product.identifiers[0].value == "Product A"
    [measurement] = product.measurements
    assert measurement.float_value.value == 85.0


def test_multiple_variations_become_separate_reactions(tmp_path):
    """Regression test: each <VARIATION> must produce its own Reaction.

    A prior version of this script reset the in-progress Reaction object on
    every VARIATION and only kept the last one, silently discarding every
    other variation's data.
    """
    reaction_xml = """
    <REACTION>
      <VARIATION>
        <SCIENTIST>Scientist One</SCIENTIST>
        <REACTANT>
          <MOLECULE MOL_ID="m1"><NAME>Reactant A</NAME></MOLECULE>
          <AMOUNT>1.0</AMOUNT>
        </REACTANT>
        <PRODUCT>
          <MOLECULE MOL_ID="m2"><NAME>Product A</NAME></MOLECULE>
        </PRODUCT>
      </VARIATION>
      <VARIATION>
        <SCIENTIST>Scientist Two</SCIENTIST>
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
        r.inputs["m1"].components[0].amount.mass.value for r in dataset.reactions
    }
    assert amounts == {1.0, 2.0}


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
