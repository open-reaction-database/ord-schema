# Copyright 2020 Open Reaction Database Project Authors
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
"""Tests for ord_schema.validations."""

import sys
import warnings

import pytest
from google.protobuf import text_format

from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


# pylint: disable=too-many-public-methods


@pytest.fixture(autouse=True)
def setup():
    # Redirect warning messages to stdout so they can be filtered from the other test output.
    original_showwarning = warnings.showwarning

    # pylint: disable=too-many-arguments
    def _showwarning(message, category, filename, lineno, file=None, line=None):
        del file  # Unused.
        original_showwarning(
            message=message, category=category, filename=filename, lineno=lineno, file=sys.stdout, line=line
        )

    # pylint: enable=too-many-arguments
    warnings.showwarning = _showwarning
    yield
    # Restore the original showwarning.
    warnings.showwarning = original_showwarning


def _run_validation(message, **kwargs):
    original = type(message)()
    original.CopyFrom(message)
    output = validations.validate_message(message, **kwargs)
    # Verify that `message` is unchanged by the validation process.
    assert original == message
    return output


@pytest.mark.parametrize(
    "message",
    (
        reaction_pb2.ReactionNotes(),
        reaction_pb2.StirringConditions(type="UNSPECIFIED"),
        reaction_pb2.ReactionNotes(safety_notes=""),
    ),
)
def test_is_empty(message):
    assert validations.is_empty(message)


@pytest.mark.parametrize(
    "message",
    (
        reaction_pb2.StirringConditions(type="STIR_BAR"),
        reaction_pb2.ReactionNotes(is_heterogeneous=False),
        reaction_pb2.ReactionNotes(is_heterogeneous=True),
    ),
)
def test_is_not_empty(message):
    assert not validations.is_empty(message)


@pytest.mark.parametrize(
    "message",
    (
        reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER),
        reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR),
        reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM),
    ),
)
def test_units(message):
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    "message,expected",
    (
        (reaction_pb2.Volume(value=-15.0, units=reaction_pb2.Volume.MILLILITER), "non-negative"),
        (reaction_pb2.Time(value=-24, units=reaction_pb2.Time.HOUR), "non-negative"),
        (reaction_pb2.Mass(value=-32.1, units=reaction_pb2.Mass.GRAM), "non-negative"),
        (reaction_pb2.FlowRate(value=5), "units"),
        (reaction_pb2.Temperature(value=-5, units="KELVIN"), "between"),
        (reaction_pb2.Temperature(value=-500, units="CELSIUS"), "between"),
    ),
)
def test_units_should_fail(message, expected):
    with pytest.raises(validations.ValidationError, match=expected):
        _run_validation(message)


def test_orcid():
    message = reaction_pb2.Person(orcid="0000-0001-2345-678X")
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_orcid_should_fail():
    message = reaction_pb2.Person(orcid="abcd-0001-2345-678X")
    with pytest.raises(validations.ValidationError, match="Invalid"):
        _run_validation(message)


@pytest.mark.parametrize("email", ["ord+test@gmail.com", "student@alumni.mit.edu", "hypen-ated@gmail.com"])
def test_email(email):
    message = reaction_pb2.Person(email=email)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize("email", ["bad&character@gmail.com", "not-an-email", "bad@domain"])
def test_email_should_fail(email):
    message = reaction_pb2.Person(email=email)
    with pytest.raises(validations.ValidationError, match="Invalid"):
        _run_validation(message)


def test_synthesized_compound():
    message = reaction_pb2.CompoundPreparation(type="SYNTHESIZED", reaction_id="dummy_reaction_id")
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0
    message = reaction_pb2.CompoundPreparation(type="NONE", reaction_id="dummy_reaction_id")
    with pytest.raises(validations.ValidationError, match="only .* when SYNTHESIZED"):
        _run_validation(message)


def test_unmeasured_amount():
    message = reaction_pb2.ReactionInput()
    message.components.add().CopyFrom(
        reaction_pb2.Compound(
            identifiers=[dict(type="SMILES", value="c1ccccc1")],
            amount=dict(unmeasured=dict(type="SATURATED")),
        )
    )
    output = _run_validation(message)
    assert len(output.warnings) == 1
    message.components.add().CopyFrom(reaction_pb2.Compound(identifiers=[dict(type="SMILES", value="O")]))
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_crude_component():
    message = reaction_pb2.CrudeComponent()
    with pytest.raises(validations.ValidationError, match="reaction_id"):
        _run_validation(message)
    message = reaction_pb2.CrudeComponent(reaction_id="my_reaction_id")
    with pytest.raises(validations.ValidationError, match="amount"):
        _run_validation(message)
    message = reaction_pb2.CrudeComponent(reaction_id="my_reaction_id", has_derived_amount=False)
    with pytest.raises(validations.ValidationError, match="amount"):
        _run_validation(message)
    message = reaction_pb2.CrudeComponent(reaction_id="my_reaction_id", has_derived_amount=True)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_product_measurement():
    message = reaction_pb2.ProductMeasurement(analysis_key="my_analysis", type="YIELD", float_value=dict(value=60))
    output = _run_validation(message)
    assert "percentage" in output.warnings[0]
    message = reaction_pb2.ProductMeasurement(analysis_key="my_analysis", type="AREA", string_value="35.221")
    with pytest.raises(validations.ValidationError, match="numeric"):
        _run_validation(message)
    message = reaction_pb2.ProductMeasurement(
        analysis_key="my_analysis", type="YIELD", percentage=dict(value=60), selectivity=dict(type="EE")
    )
    with pytest.raises(validations.ValidationError, match="selectivity"):
        _run_validation(message)


def test_reaction():
    message = reaction_pb2.Reaction()
    with pytest.raises(validations.ValidationError, match="reaction input"):
        _run_validation(message)


def test_reaction_smiles():
    message = reaction_pb2.Reaction()
    message.identifiers.add(value="C>>C", type="REACTION_SMILES")
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0
    # Now disable the exception for reaction SMILES only.
    options = validations.ValidationOptions()
    options.allow_reaction_smiles_only = False
    with pytest.raises(validations.ValidationError, match="reaction input"):
        _run_validation(message, options=options)


def test_bad_reaction_smiles():
    message = reaction_pb2.Reaction()
    message.identifiers.add(value="test", type="REACTION_SMILES")
    with pytest.raises(validations.ValidationError, match="requires at least two > characters"):
        _run_validation(message)


@pytest.mark.parametrize(
    "workup_text",
    (
        "type: ALIQUOT amount {mass {value: 1.0 units: GRAM}}",
        "type: ADDITION input {components {"
        'identifiers {value: "CCO" type: SMILES} '
        "amount {mass {value: 10.0 units: GRAM}}}}",
    ),
)
def test_reaction_workup(workup_text):
    message = text_format.Parse(workup_text, reaction_pb2.ReactionWorkup())
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    "workup_text,error_msg",
    (
        ("type: ALIQUOT", "missing volume/mass"),
        (
            "type: ADDITION amount {mass {value: 1.0 units: GRAM}} "
            'input {components {identifiers {value: "CCO" type: SMILES} '
            "amount {mass {value: 10.0 units: GRAM}}}}",
            "should only be specified",
        ),
    ),
)
def test_bad_reaction_workup(workup_text, error_msg):
    message = text_format.Parse(workup_text, reaction_pb2.ReactionWorkup())
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert error_msg in output.warnings[0]


# pylint: disable=too-many-statements
def test_reaction_recursive():
    message = reaction_pb2.Reaction()
    # Reactions must have at least one input
    with pytest.raises(validations.ValidationError, match="reaction input"):
        _run_validation(message, recurse=False)
    dummy_input = message.inputs["dummy_input"]
    # Reactions must have at least one outcome
    with pytest.raises(validations.ValidationError, match="reaction outcome"):
        _run_validation(message, recurse=False)
    outcome = message.outcomes.add()
    output = _run_validation(message, recurse=False)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0
    # Inputs must have at least one component
    with pytest.raises(validations.ValidationError, match="component"):
        _run_validation(message)
    dummy_component = dummy_input.components.add()
    # Components must have at least one identifier
    with pytest.raises(validations.ValidationError, match="identifier"):
        _run_validation(message)
    dummy_component.identifiers.add(type="CUSTOM")
    # Custom identifiers must have details specified
    with pytest.raises(validations.ValidationError, match="details"):
        _run_validation(message)
    dummy_component.identifiers[0].details = "custom_identifier"
    dummy_component.identifiers[0].value = "custom_value"
    # Components of reaction inputs must have a defined amount
    with pytest.raises(validations.ValidationError, match="require an amount"):
        _run_validation(message)
    dummy_component.amount.mass.value = 1
    dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
    # Reactions don't require defined products or conversion because they
    # can be used to prepare crude for a second Reaction
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 1
    outcome.conversion.value = 75
    # If converseions are defined, must have limiting reagent flag
    with pytest.raises(validations.ValidationError, match="is_limiting"):
        _run_validation(message)
    dummy_component.is_limiting = True
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0

    # If a measurement uses an internal standard, a component must have
    # an INTERNAL_STANDARD reaction role
    outcome.analyses["dummy_analysis"].CopyFrom(reaction_pb2.Analysis(type="CUSTOM", details="test"))
    product = outcome.products.add(identifiers=[dict(type="SMILES", value="c1ccccc1")])
    product.measurements.add(
        analysis_key="dummy_analysis", type="YIELD", percentage=dict(value=75), uses_internal_standard=True
    )
    with pytest.raises(validations.ValidationError, match="INTERNAL_STANDARD"):
        _run_validation(message)
    # Assigning internal standard role to input should resolve the error
    message_input_istd = reaction_pb2.Reaction()
    message_input_istd.CopyFrom(message)
    message_input_istd.inputs["dummy_input"].components[0].reaction_role = reaction_pb2.ReactionRole.INTERNAL_STANDARD
    output = _run_validation(message_input_istd)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0
    # Assigning internal standard role to workup should resolve the error
    message_workup_istd = reaction_pb2.Reaction()
    message_workup_istd.CopyFrom(message)
    workup = message_workup_istd.workups.add(type="CUSTOM", details="test")
    istd = workup.input.components.add()
    istd.identifiers.add(type="SMILES", value="CCO")
    istd.amount.mass.value = 1
    istd.amount.mass.units = reaction_pb2.Mass.GRAM
    istd.reaction_role = reaction_pb2.ReactionRole.INTERNAL_STANDARD
    output = _run_validation(message_workup_istd)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_reaction_recursive_noraise_on_error():
    message = reaction_pb2.Reaction()
    message.inputs["dummy_input"].components.add()
    output = _run_validation(message, raise_on_error=False)
    expected = [
        'Reaction.inputs["dummy_input"].components[0]: ' "Compounds must have at least one identifier",
        "Reaction: Reactions should have at least 1 reaction outcome",
        "Reaction: All reaction input components require an amount",
    ]
    assert output.errors == expected
    assert len(output.warnings) == 1


def test_datetimes():
    message = reaction_pb2.ReactionProvenance()
    message.experiment_start.value = "2020-01-02"
    message.record_created.time.value = "2020-01-01"
    message.record_created.person.name = "test"
    message.record_created.person.email = "test@example.com"
    with pytest.raises(validations.ValidationError, match="after"):
        _run_validation(message)
    message.record_created.time.value = "2020-01-03"
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_reaction_id():
    message = reaction_pb2.Reaction()
    _ = message.inputs["test"]
    message.outcomes.add()
    message.reaction_id = "ord-c0bbd41f095a44a78b6221135961d809"
    options = validations.ValidationOptions(validate_ids=True)
    output = _run_validation(message, recurse=False, options=options)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    "reaction_id",
    (
        "ord-c0bbd41f095a4",
        "ord-c0bbd41f095a4c0bbd41f095a4c0bbd41f095a4",
        "foo-c0bbd41f095a44a78b6221135961d809",
        "ord-C0BBD41F095A44A78B6221135961D809",
        "ord-h0bbd41f095a44a78b6221135961d809",
        "ord-notARealId",
        "",
    ),
)
def test_bad_reaction_id(reaction_id):
    message = reaction_pb2.Reaction(reaction_id=reaction_id)
    _ = message.inputs["test"]
    message.outcomes.add()
    options = validations.ValidationOptions(validate_ids=True)
    with pytest.raises(validations.ValidationError, match="malformed"):
        _run_validation(message, recurse=False, options=options)


def test_missing_provenance():
    message = reaction_pb2.Reaction()
    _ = message.inputs["test"]
    message.outcomes.add()
    options = validations.ValidationOptions(require_provenance=True)
    with pytest.raises(validations.ValidationError, match="requires provenance"):
        _run_validation(message, recurse=False, options=options)


def test_data():
    message = reaction_pb2.Data()
    with pytest.raises(validations.ValidationError, match="requires one of"):
        _run_validation(message)
    message.bytes_value = b"test data"
    with pytest.raises(validations.ValidationError, match="format is required"):
        _run_validation(message)
    message.string_value = "test data"
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_dataset_bad_reaction_id():
    message = dataset_pb2.Dataset(reaction_ids=["foo"])
    options = validations.ValidationOptions(validate_ids=True)
    with pytest.raises(validations.ValidationError, match="malformed"):
        _run_validation(message, options=options)


def test_dataset_records_and_ids():
    message = dataset_pb2.Dataset(
        reactions=[reaction_pb2.Reaction()],
        reaction_ids=["ord-c0bbd41f095a44a78b6221135961d809"],
    )
    options = validations.ValidationOptions(validate_ids=True)
    with pytest.raises(validations.ValidationError, match="not both"):
        _run_validation(message, recurse=False, options=options)


def test_dataset_bad_id():
    message = dataset_pb2.Dataset(reactions=[reaction_pb2.Reaction()], dataset_id="foo")
    options = validations.ValidationOptions(validate_ids=True)
    with pytest.raises(validations.ValidationError, match="malformed"):
        _run_validation(message, recurse=False, options=options)


def test_dataset_example():
    message = dataset_pb2.DatasetExample()
    with pytest.raises(validations.ValidationError, match="description is required"):
        _run_validation(message)
    message.description = "test example"
    with pytest.raises(validations.ValidationError, match="url is required"):
        _run_validation(message)
    message.url = "example.com"
    with pytest.raises(validations.ValidationError, match="created is required"):
        _run_validation(message)
    message.created.time.value = "11 am"
    message.created.person.name = "test"
    message.created.person.email = "test@example.com"
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_dataset_crossreferences():
    message = dataset_pb2.Dataset()
    reaction1 = message.reactions.add()
    reaction2 = message.reactions.add()
    reaction3 = message.reactions.add()
    # Minimal reaction 1
    dummy_input = reaction1.inputs["dummy_input"]
    reaction1.outcomes.add()
    dummy_component = dummy_input.components.add()
    dummy_component.identifiers.add(type="CUSTOM")
    dummy_component.identifiers[0].details = "custom_identifier"
    dummy_component.identifiers[0].value = "custom_value"
    dummy_component.amount.mass.value = 1
    dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
    reaction2.CopyFrom(reaction1)
    reaction3.CopyFrom(reaction1)
    # It is not required to cross-reference SYNTHESIZED compounds.
    dummy_component.preparations.add(type="SYNTHESIZED")
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 3
    # References must refer to ID within this dataset.
    dummy_component.preparations[0].reaction_id = "placeholder_id"
    with pytest.raises(validations.ValidationError, match="undefined reaction_ids"):
        _run_validation(message)
    # Self-references aren't allowed
    reaction1.reaction_id = "placeholder_id"
    with pytest.raises(validations.ValidationError, match="own ID"):
        _run_validation(message)
    # Valid reference
    reaction1.reaction_id = ""
    reaction2.reaction_id = "placeholder_id"
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 3
    # Duplicate IDs not allowed
    reaction3.reaction_id = "placeholder_id"
    with pytest.raises(validations.ValidationError, match="same IDs"):
        _run_validation(message)
    reaction3.reaction_id = ""
    # Crude component also needs valid IDs
    dummy_input.crude_components.add(reaction_id="crude-making step", has_derived_amount=True)
    with pytest.raises(validations.ValidationError, match="undefined reaction_ids"):
        _run_validation(message)
    reaction3.reaction_id = "crude-making step"
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 3
