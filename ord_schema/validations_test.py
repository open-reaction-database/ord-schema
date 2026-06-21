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

from ord_schema import validations
from ord_schema.proto import dataset_pb2, reaction_pb2


@pytest.fixture(autouse=True)
def setup():
    # Redirect warning messages to stdout so they can be filtered from the other test output.
    original_showwarning = warnings.showwarning

    def _showwarning(message, category, filename, lineno, file=None, line=None):
        del file  # Unused.
        original_showwarning(
            message=message,
            category=category,
            filename=filename,
            lineno=lineno,
            file=sys.stdout,
            line=line,
        )

    warnings.showwarning = _showwarning  # ty: ignore[invalid-assignment]
    yield
    # Restore the original showwarning.
    warnings.showwarning = original_showwarning


def _run_validation(message, **kwargs):
    original = type(message)()
    original.CopyFrom(message)
    if "options" not in kwargs:
        kwargs["options"] = validations.ValidationOptions(require_provenance=False)
    output = validations.validate_message(message, **kwargs)
    # Verify that `message` is unchanged by the validation process.
    assert original == message
    return output


@pytest.mark.parametrize(
    "message",
    [
        reaction_pb2.ReactionNotes(),
        reaction_pb2.StirringConditions(type="UNSPECIFIED"),
        reaction_pb2.ReactionNotes(safety_notes=""),
    ],
)
def test_is_empty(message):
    assert validations.is_empty(message)


@pytest.mark.parametrize(
    "message",
    [
        reaction_pb2.StirringConditions(type="STIR_BAR"),
        reaction_pb2.ReactionNotes(is_heterogeneous=False),
        reaction_pb2.ReactionNotes(is_heterogeneous=True),
    ],
)
def test_is_not_empty(message):
    assert not validations.is_empty(message)


@pytest.mark.parametrize(
    "message",
    [
        reaction_pb2.Volume(value=15.0, units=reaction_pb2.Volume.MILLILITER),
        reaction_pb2.Time(value=24, units=reaction_pb2.Time.HOUR),
        reaction_pb2.Mass(value=32.1, units=reaction_pb2.Mass.GRAM),
        reaction_pb2.Temperature(value=25.0, units=reaction_pb2.Temperature.CELSIUS),
    ],
)
def test_units(message):
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    ("message", "expected"),
    [
        (
            reaction_pb2.Volume(value=-15.0, units=reaction_pb2.Volume.MILLILITER),
            "non-negative",
        ),
        (reaction_pb2.Time(value=-24, units=reaction_pb2.Time.HOUR), "non-negative"),
        (reaction_pb2.Mass(value=-32.1, units=reaction_pb2.Mass.GRAM), "non-negative"),
        (reaction_pb2.FlowRate(value=5), "units"),
        (reaction_pb2.Temperature(value=-5, units="KELVIN"), "between"),
        (reaction_pb2.Temperature(value=-500, units="CELSIUS"), "between"),
    ],
)
def test_units_should_fail(message, expected):
    with pytest.raises(validations.ValidationError, match=expected):
        _run_validation(message)


def test_orcid():
    message = reaction_pb2.Person(orcid="0000-0002-1825-0097")
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    "orcid",
    [
        "abcd-0001-2345-678X",  # Non-numeric.
        "0000-0001-2345-678X",  # Well-formed but incorrect checksum.
    ],
)
def test_orcid_should_fail(orcid):
    message = reaction_pb2.Person(orcid=orcid)
    with pytest.raises(validations.ValidationError, match="Invalid"):
        _run_validation(message)


@pytest.mark.parametrize(
    "email", ["ord+test@gmail.com", "student@alumni.mit.edu", "hypen-ated@gmail.com"]
)
def test_email(email):
    message = reaction_pb2.Person(email=email)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    "email", ["bad&character@gmail.com", "not-an-email", "bad@domain"]
)
def test_email_should_fail(email):
    message = reaction_pb2.Person(email=email)
    with pytest.raises(validations.ValidationError, match="Invalid"):
        _run_validation(message)


def test_synthesized_compound():
    message = reaction_pb2.CompoundPreparation(
        type="SYNTHESIZED", reaction_id="dummy_reaction_id"
    )
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0
    message = reaction_pb2.CompoundPreparation(
        type="NONE", reaction_id="dummy_reaction_id"
    )
    with pytest.raises(validations.ValidationError, match=r"only .* when SYNTHESIZED"):
        _run_validation(message)


def test_unmeasured_amount():
    message = reaction_pb2.ReactionInput()
    message.components.add().CopyFrom(
        reaction_pb2.Compound(
            identifiers=[{"type": "SMILES", "value": "c1ccccc1"}],
            amount={"unmeasured": {"type": "SATURATED"}},
        )
    )
    output = _run_validation(message)
    assert len(output.warnings) == 1
    message.components.add().CopyFrom(
        reaction_pb2.Compound(identifiers=[{"type": "SMILES", "value": "O"}])
    )
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_texture_in_reaction_input():
    def _make_dummy_reaction_input(component_texture_types, input_texture_type):
        message = reaction_pb2.ReactionInput(
            texture=reaction_pb2.Texture(type=input_texture_type)
        )
        for texture_type in component_texture_types:
            message.components.add().CopyFrom(
                reaction_pb2.Compound(
                    identifiers=[{"type": "SMILES", "value": "c1ccccc1"}],
                    texture=reaction_pb2.Texture(type=texture_type),
                )
            )
        return message

    # (1) a foam and a gas are unlikely to give a crystal
    c_texture_types = ["FOAM", "GAS"]
    i_texture_type = "CRYSTAL"
    output = _run_validation(
        _make_dummy_reaction_input(c_texture_types, i_texture_type)
    )
    assert len(output.warnings) == 1

    # (2) a wax and a liquid may give a liquid
    c_texture_types = ["WAX", "LIQUID"]
    i_texture_type = "LIQUID"
    output = _run_validation(
        _make_dummy_reaction_input(c_texture_types, i_texture_type)
    )
    assert len(output.warnings) == 0

    # (3) an oil and a liquid are unlikely to give a solid
    c_texture_types = ["OIL", "LIQUID"]
    i_texture_type = "SOLID"
    output = _run_validation(
        _make_dummy_reaction_input(c_texture_types, i_texture_type)
    )
    assert len(output.warnings) == 1

    # (4) a gas and a liquid may give a liquid
    c_texture_types = ["GAS", "LIQUID"]
    i_texture_type = "LIQUID"
    output = _run_validation(
        _make_dummy_reaction_input(c_texture_types, i_texture_type)
    )
    assert len(output.warnings) == 0


def test_crude_component():
    message = reaction_pb2.CrudeComponent()
    with pytest.raises(validations.ValidationError, match="reaction_id"):
        _run_validation(message)
    message = reaction_pb2.CrudeComponent(reaction_id="my_reaction_id")
    with pytest.raises(validations.ValidationError, match="amount"):
        _run_validation(message)
    message = reaction_pb2.CrudeComponent(
        reaction_id="my_reaction_id", has_derived_amount=False
    )
    with pytest.raises(validations.ValidationError, match="amount"):
        _run_validation(message)
    message = reaction_pb2.CrudeComponent(
        reaction_id="my_reaction_id", has_derived_amount=True
    )
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_product_measurement():
    message = reaction_pb2.ProductMeasurement(
        analysis_key="my_analysis", type="YIELD", float_value={"value": 60}
    )
    output = _run_validation(message)
    assert "percentage" in output.warnings[0]
    message = reaction_pb2.ProductMeasurement(
        analysis_key="my_analysis", type="AREA", string_value="35.221"
    )
    with pytest.raises(validations.ValidationError, match="numeric"):
        _run_validation(message)
    message = reaction_pb2.ProductMeasurement(
        analysis_key="my_analysis",
        type="YIELD",
        percentage={"value": 60},
        selectivity={"type": "EE"},
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
    with pytest.raises(
        validations.ValidationError, match="requires at least two > characters"
    ):
        _run_validation(message)


@pytest.mark.parametrize(
    ("identifier_type", "value"),
    [
        ("SMILES", "CO"),
        ("INCHI", "InChI=1S/CH4O/c1-2/h2H,1H3"),
        (
            "MOLBLOCK",
            "\n     RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n",
        ),
    ],
)
def test_compound_identifier(identifier_type, value):
    message = reaction_pb2.CompoundIdentifier(type=identifier_type, value=value)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    ("identifier_type", "value"),
    [
        (
            "SMILES",
            "[O-]1(C)[Ir+]234([O-](C)[Ir+]1567[CH]=8CC[CH]7=[CH]6CC[CH]85)[CH]=9CC[CH]4=[CH]3CC[CH]92",
        ),
        ("SMILES", "CO(C)(C)(C)C"),
        ("SMILES", "On1c(-c2ccccc2)c(-c2c(-c3ccccc3)nc3ccccc23)c2ccccc21"),
    ],
)
def test_bad_compound_identifier(identifier_type, value):
    message = reaction_pb2.CompoundIdentifier(type=identifier_type, value=value)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 1
    assert "could not sanitize" in output.warnings[0]


@pytest.mark.parametrize(
    ("identifier_type", "value"),
    [
        ("SMILES", "###"),
        ("INCHI", "###"),
        ("MOLBLOCK", "###"),
    ],
)
def test_invalid_compound_identifier(identifier_type, value):
    message = reaction_pb2.CompoundIdentifier(type=identifier_type, value=value)
    with pytest.raises(validations.ValidationError, match="could not validate"):
        _run_validation(message)


@pytest.mark.parametrize(
    ("identifier_type", "value"),
    [
        ("CAS_NUMBER", "64-17-5"),
        ("INCHI_KEY", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"),
        ("PUBCHEM_CID", "702"),
        ("CHEMSPIDER_ID", "682"),
    ],
)
def test_compound_identifier_format(identifier_type, value):
    message = reaction_pb2.CompoundIdentifier(type=identifier_type, value=value)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    ("identifier_type", "value", "expected"),
    [
        ("CAS_NUMBER", "64175", "CAS number"),
        ("INCHI_KEY", "not-a-key", "InChIKey"),
        ("PUBCHEM_CID", "abc", "PubChem CID"),
        ("CHEMSPIDER_ID", "12x", "ChemSpider ID"),
    ],
)
def test_bad_compound_identifier_format(identifier_type, value, expected):
    message = reaction_pb2.CompoundIdentifier(type=identifier_type, value=value)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 1
    assert expected in output.warnings[0]


@pytest.mark.parametrize("ph_value", [7.0, 0.0, 14.0])
def test_reaction_conditions_ph(ph_value):
    message = reaction_pb2.ReactionConditions(ph=ph_value)
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize("ph_value", [-1.0, 15.0])
def test_reaction_conditions_bad_ph(ph_value):
    message = reaction_pb2.ReactionConditions(ph=ph_value)
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "outside the expected range" in output.warnings[0]


def test_electrochemistry_missing_current():
    message = reaction_pb2.ElectrochemistryConditions(type="CONSTANT_CURRENT")
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "should specify the current" in output.warnings[0]


def test_electrochemistry_missing_voltage():
    message = reaction_pb2.ElectrochemistryConditions(type="CONSTANT_VOLTAGE")
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "should specify the voltage" in output.warnings[0]


def test_mass_spec_mz_range():
    message = reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails(
        type="TIC", tic_minimum_mz=100.0, tic_maximum_mz=50.0
    )
    with pytest.raises(validations.ValidationError, match="tic_minimum_mz"):
        _run_validation(message)


def test_mass_spec_eic_requires_masses():
    message = reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails(type="EIC")
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "should specify eic_masses" in output.warnings[0]


def test_mass_spec_negative_eic_mass():
    message = reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails(
        type="EIC", eic_masses=[-1.0]
    )
    with pytest.raises(
        validations.ValidationError, match="eic_masses must be non-negative"
    ):
        _run_validation(message)


def test_mass_spec_tic_with_eic_masses():
    message = reaction_pb2.ProductMeasurement.MassSpecMeasurementDetails(
        type="TIC", eic_masses=[100.0]
    )
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "should only be specified for EIC" in output.warnings[0]


def test_provenance_bad_url():
    message = reaction_pb2.ReactionProvenance(
        record_created={
            "time": {"value": "2021-01-01"},
            "person": {"username": "test", "email": "a@b.com"},
        },
        publication_url="not a url",
    )
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "valid URL" in output.warnings[0]


def test_illumination_dark_with_wavelength():
    message = reaction_pb2.IlluminationConditions(
        type="DARK", peak_wavelength={"value": 450.0, "units": "NANOMETER"}
    )
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "DARK or AMBIENT" in output.warnings[0]


@pytest.mark.parametrize(
    ("value", "is_mapped", "expected"),
    [
        # Atom maps present but the flag is not set.
        ("[CH4:1]>>[CH4:1]", False, "is_mapped is not set"),
        # Flag is set but the SMILES contains no atom maps.
        ("CO>>CO", True, "no atom maps"),
    ],
)
def test_reaction_identifier_atom_mapping(value, is_mapped, expected):
    message = reaction_pb2.ReactionIdentifier(
        type="REACTION_SMILES", value=value, is_mapped=is_mapped
    )
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert expected in output.warnings[0]


def test_data_bad_url():
    message = reaction_pb2.Data(url="not a url")
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert "valid URL" in output.warnings[0]


_ADDITION_INPUT = {
    "components": [
        {
            "identifiers": [{"value": "CCO", "type": "SMILES"}],
            "amount": {"mass": {"value": 10.0, "units": "GRAM"}},
        }
    ]
}


@pytest.mark.parametrize(
    "message",
    [
        reaction_pb2.ReactionWorkup(
            type="ALIQUOT", amount={"mass": {"value": 1.0, "units": "GRAM"}}
        ),
        reaction_pb2.ReactionWorkup(type="ADDITION", input=_ADDITION_INPUT),
    ],
)
def test_reaction_workup(message):
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    ("message", "error_msg"),
    [
        (reaction_pb2.ReactionWorkup(type="ALIQUOT"), "missing volume/mass"),
        (
            reaction_pb2.ReactionWorkup(
                type="ADDITION",
                amount={"mass": {"value": 1.0, "units": "GRAM"}},
                input=_ADDITION_INPUT,
            ),
            "should only be specified",
        ),
    ],
)
def test_bad_reaction_workup(message, error_msg):
    output = _run_validation(message)
    assert len(output.warnings) == 1
    assert error_msg in output.warnings[0]


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
    outcome.analyses["dummy_analysis"].CopyFrom(
        reaction_pb2.Analysis(type="CUSTOM", details="test")
    )
    product = outcome.products.add(
        identifiers=[{"type": "SMILES", "value": "c1ccccc1"}]
    )
    product.measurements.add(
        analysis_key="dummy_analysis",
        type="YIELD",
        percentage={"value": 75},
        uses_internal_standard=True,
    )
    with pytest.raises(validations.ValidationError, match="INTERNAL_STANDARD"):
        _run_validation(message)
    # Assigning internal standard role to input should resolve the error
    message_input_istd = reaction_pb2.Reaction()
    message_input_istd.CopyFrom(message)
    message_input_istd.inputs["dummy_input"].components[
        0
    ].reaction_role = reaction_pb2.ReactionRole.INTERNAL_STANDARD
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
        'Reaction.inputs["dummy_input"].components[0]: Compounds must have at least one identifier',
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
    options = validations.ValidationOptions(validate_ids=True, require_provenance=False)
    output = _run_validation(message, recurse=False, options=options)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


@pytest.mark.parametrize(
    "reaction_id",
    [
        "ord-c0bbd41f095a4",
        "ord-c0bbd41f095a4c0bbd41f095a4c0bbd41f095a4",
        "foo-c0bbd41f095a44a78b6221135961d809",
        "ord-C0BBD41F095A44A78B6221135961D809",
        "ord-h0bbd41f095a44a78b6221135961d809",
        "ord-notARealId",
        "",
    ],
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


def test_bad_doi():
    message = reaction_pb2.ReactionProvenance()
    message.record_created.time.value = "2023-07-01"
    message.record_created.person.email = "test@example.com"
    message.doi = "149"
    options = validations.ValidationOptions(require_provenance=True)
    with pytest.raises(validations.ValidationError, match="could not parse DOI"):
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


def test_dataset_missing_name():
    message = dataset_pb2.Dataset()
    options = validations.ValidationOptions(validate_ids=True)
    with pytest.raises(validations.ValidationError, match="name is required"):
        _run_validation(message, options=options)


def test_dataset_bad_reaction_id():
    message = dataset_pb2.Dataset(name="test", description="test", reaction_ids=["foo"])
    options = validations.ValidationOptions(validate_ids=True)
    with pytest.raises(validations.ValidationError, match="malformed"):
        _run_validation(message, options=options)


def test_dataset_records_and_ids():
    message = dataset_pb2.Dataset(
        name="test",
        description="test",
        reactions=[reaction_pb2.Reaction()],
        reaction_ids=["ord-c0bbd41f095a44a78b6221135961d809"],
    )
    options = validations.ValidationOptions(validate_ids=True)
    with pytest.raises(validations.ValidationError, match="not both"):
        _run_validation(message, recurse=False, options=options)


def test_dataset_bad_id():
    message = dataset_pb2.Dataset(
        name="test",
        description="test",
        reactions=[reaction_pb2.Reaction()],
        dataset_id="foo",
    )
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


def test_dataset_cross_references():
    message = dataset_pb2.Dataset(name="test", description="test")
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
    dummy_input.crude_components.add(
        reaction_id="crude-making step", has_derived_amount=True
    )
    with pytest.raises(validations.ValidationError, match="undefined reaction_ids"):
        _run_validation(message)
    reaction3.reaction_id = "crude-making step"
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 3


@pytest.mark.parametrize("message_cls", list(validations._VALIDATOR_SWITCH))
def test_validator_switch_dispatches(message_cls):
    """Default-empty instance of every dispatched type should validate without raising.

    Most validators emit warnings (e.g. missing required fields) on an empty proto;
    we only care that dispatch reaches the validator and returns without crashing.
    """
    message = message_cls()
    output = validations.validate_message(
        message,
        raise_on_error=False,
        options=validations.ValidationOptions(require_provenance=False),
    )
    assert isinstance(output, validations.ValidationOutput)


def _capture_warnings(func, *args, **kwargs):
    """Calls a validator directly and returns the recorded warnings.

    Used for branches that are shadowed by recursion in ``validate_message``
    (e.g. provenance-level DateTime parsing, which the DateTime validator would
    otherwise flag first).
    """
    with warnings.catch_warnings(record=True) as tape:
        warnings.simplefilter("always")
        func(*args, **kwargs)
    return tape


def test_type_and_details_missing_type():
    # Non-empty TypeDetailsMessage with an unset type.
    message = reaction_pb2.Texture(details="fine needles")
    with pytest.raises(validations.ValidationError, match="requires `type`"):
        _run_validation(message)


def test_type_and_details_custom_without_details():
    message = reaction_pb2.Texture(type="CUSTOM")
    with pytest.raises(validations.ValidationError, match="CUSTOM but details"):
        _run_validation(message)


def test_reaction_cxsmiles():
    message = reaction_pb2.ReactionIdentifier(
        type="REACTION_CXSMILES", value="CO>>CO |f:0|"
    )
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_amount_volume_includes_solutes_non_volume():
    message = reaction_pb2.Amount(
        mass=reaction_pb2.Mass(value=1.0, units="GRAM"),
        volume_includes_solutes=True,
    )
    with pytest.raises(validations.ValidationError, match="only be set for volume"):
        _run_validation(message)


@pytest.mark.parametrize(
    ("message", "expected"),
    [
        (
            reaction_pb2.CrudeComponent(
                reaction_id="r",
                has_derived_amount=True,
                amount={"mass": {"value": 1.0, "units": "GRAM"}},
            ),
            "cannot have their mass",
        ),
        (
            reaction_pb2.CrudeComponent(
                reaction_id="r", amount={"moles": {"value": 1.0, "units": "MOLE"}}
            ),
            "must be specified by mass or volume",
        ),
        (
            reaction_pb2.CrudeComponent(
                reaction_id="r",
                amount={
                    "volume": {"value": 1.0, "units": "LITER"},
                    "volume_includes_solutes": True,
                },
            ),
            "only be used for input Compounds",
        ),
    ],
)
def test_crude_component_amounts(message, expected):
    with pytest.raises(validations.ValidationError, match=expected):
        _run_validation(message)


def test_compound_inconsistent_identifiers():
    message = reaction_pb2.Compound(
        identifiers=[
            {"type": "SMILES", "value": "CO"},  # methanol
            {"type": "INCHI", "value": "InChI=1S/CH4/h1H4"},  # methane
        ]
    )
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert any(
        "could not validate that" in warning.lower() or "consistent" in warning.lower()
        for warning in output.warnings
    )


@pytest.mark.parametrize(
    ("message", "raises", "expected"),
    [
        (reaction_pb2.ReactionConditions(conditions_are_dynamic=True), True, "dynamic"),
        (reaction_pb2.ReactionConditions(details="ramped"), False, "details provided"),
    ],
)
def test_reaction_conditions_dynamic(message, raises, expected):
    if raises:
        with pytest.raises(validations.ValidationError, match=expected):
            _run_validation(message)
    else:
        output = _run_validation(message)
        assert len(output.warnings) == 1
        assert expected in output.warnings[0]


@pytest.mark.parametrize(
    ("message", "expected"),
    [
        (reaction_pb2.ReactionWorkup(type="WAIT"), "WAIT workup"),
        (reaction_pb2.ReactionWorkup(type="TEMPERATURE"), "temperature conditions"),
        (reaction_pb2.ReactionWorkup(type="EXTRACTION"), "keep_phase"),
        (reaction_pb2.ReactionWorkup(type="WASH"), "missing recommended inputs"),
        (reaction_pb2.ReactionWorkup(type="STIRRING"), "Stirring workup"),
        (
            reaction_pb2.ReactionWorkup(
                type="PH_ADJUST",
                input={
                    "components": [
                        {
                            "identifiers": [{"value": "O", "type": "SMILES"}],
                            "amount": {"volume": {"value": 1.0, "units": "LITER"}},
                        }
                    ]
                },
            ),
            "missing target pH",
        ),
        (
            reaction_pb2.ReactionWorkup(
                type="ALIQUOT", amount={"moles": {"value": 1.0, "units": "MOLE"}}
            ),
            "specified by mass or volume",
        ),
        (
            reaction_pb2.ReactionWorkup(
                type="ALIQUOT",
                amount={
                    "volume": {"value": 1.0, "units": "LITER"},
                    "volume_includes_solutes": True,
                },
            ),
            "only be used for input Compounds",
        ),
        (
            reaction_pb2.ReactionWorkup(
                type="CONCENTRATION", amount={"mass": {"value": 1.0, "units": "GRAM"}}
            ),
            "only be specified if workup type",
        ),
    ],
)
def test_reaction_workup_recommendations(message, expected):
    output = _run_validation(message)
    assert any(expected in warning for warning in output.warnings)


def test_reaction_outcome_bad_analysis_key():
    message = reaction_pb2.ReactionOutcome(
        products=[
            {
                "identifiers": [{"value": "CO", "type": "SMILES"}],
                "measurements": [{"type": "IDENTITY", "analysis_key": "missing"}],
            }
        ]
    )
    with pytest.raises(
        validations.ValidationError, match="does not match any known analysis"
    ):
        _run_validation(message)


def test_product_compound_side_product_desired():
    message = reaction_pb2.ProductCompound(
        identifiers=[{"value": "CO", "type": "SMILES"}],
        is_desired_product=True,
        reaction_role="SIDE_PRODUCT",
    )
    with pytest.raises(validations.ValidationError, match="SIDE_PRODUCT"):
        _run_validation(message)


def test_product_compound_inconsistent_identifiers():
    message = reaction_pb2.ProductCompound(
        identifiers=[
            {"type": "SMILES", "value": "CO"},  # methanol
            {"type": "INCHI", "value": "InChI=1S/CH4/h1H4"},  # methane
        ],
        is_desired_product=True,  # Also exercises the desired/not-SIDE_PRODUCT branch.
    )
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert any(
        "could not validate that" in warning.lower() or "consistent" in warning.lower()
        for warning in output.warnings
    )


@pytest.mark.parametrize(
    ("message", "raises", "expected"),
    [
        (
            reaction_pb2.ProductMeasurement(
                type="IDENTITY", analysis_key="x", percentage={"value": 50.0}
            ),
            True,
            "IDENTITY should not have",
        ),
        (
            reaction_pb2.ProductMeasurement(
                type="YIELD", analysis_key="x", float_value={"value": 5.0}
            ),
            False,
            "YIELD measurements",
        ),
        (
            reaction_pb2.ProductMeasurement(
                type="AREA", analysis_key="x", string_value="big"
            ),
            True,
            "must use numeric values",
        ),
    ],
)
def test_product_measurement_values(message, raises, expected):
    if raises:
        with pytest.raises(validations.ValidationError, match=expected):
            _run_validation(message)
    else:
        output = _run_validation(message)
        assert any(expected in warning for warning in output.warnings)


def test_bad_datetime():
    message = reaction_pb2.DateTime(value="definitely not a date")
    with pytest.raises(validations.ValidationError, match="Could not parse DateTime"):
        _run_validation(message)


def test_provenance_record_modified_before_created():
    message = reaction_pb2.ReactionProvenance(
        record_created={
            "time": {"value": "2021-06-01"},
            "person": {"username": "u", "email": "a@b.com"},
        },
        record_modified=[
            {
                "time": {"value": "2021-01-01"},
                "person": {"username": "u", "email": "a@b.com"},
            }
        ],
    )
    with pytest.raises(validations.ValidationError, match="after creation"):
        _run_validation(message)


def test_provenance_doi_should_be_trimmed():
    message = reaction_pb2.ReactionProvenance(
        record_created={
            "time": {"value": "2021-01-01"},
            "person": {"username": "u", "email": "a@b.com"},
        },
        doi="https://doi.org/10.1234/foo",
    )
    with pytest.raises(validations.ValidationError, match="trimmed"):
        _run_validation(message)


def test_provenance_unparseable_datetime_direct():
    # Shadowed by the DateTime validator under recursion, so call directly.
    message = reaction_pb2.ReactionProvenance()
    message.record_created.time.value = "garbage"
    message.record_created.person.username = "u"
    message.record_created.person.email = "a@b.com"
    tape = _capture_warnings(validations.validate_reaction_provenance, message)
    assert any("Failed to parse" in str(warning.message) for warning in tape)


def test_temperature_fahrenheit_below_absolute_zero():
    message = reaction_pb2.Temperature(value=-500.0, units="FAHRENHEIT")
    with pytest.raises(validations.ValidationError, match="between"):
        _run_validation(message)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (0.5, "not fractions"),
        (150.0, "outside the expected range"),
    ],
)
def test_percentage_ranges(value, expected):
    message = reaction_pb2.Percentage(value=value)
    output = _run_validation(message)
    assert any(expected in warning for warning in output.warnings)


def test_cross_ref_state_merge():
    state_a = validations.DatasetCrossRefState(defined_ids={"a", "shared"})
    state_b = validations.DatasetCrossRefState(
        defined_ids={"b", "shared"}, referenced_ids={"a"}, self_reference_count=1
    )
    state_a.merge(state_b)
    assert state_a.defined_ids == {"a", "b", "shared"}
    assert state_a.referenced_ids == {"a"}
    assert state_a.duplicate_count == 1  # "shared" defined in both.
    assert state_a.self_reference_count == 1


def test_validators_default_options():
    # Exercise the ``options is None`` default branches.
    validations.validate_reaction(reaction_pb2.Reaction(), options=None)
    validations.validate_dataset(dataset_pb2.Dataset(), options=None)
    validations.validate_dataset_streaming(
        name="n",
        description="d",
        dataset_id="",
        reaction_ids=[],
        has_reactions=True,
        state=validations.DatasetCrossRefState(),
        options=None,
    )


def test_workup_target_ph_out_of_range():
    message = reaction_pb2.ReactionWorkup(
        type="PH_ADJUST",
        target_ph=20.0,
        input={
            "components": [
                {
                    "identifiers": [{"value": "O", "type": "SMILES"}],
                    "amount": {"mass": {"value": 1.0, "units": "GRAM"}},
                }
            ]
        },
    )
    output = _run_validation(message)
    assert any("outside the expected range" in warning for warning in output.warnings)


def test_reaction_outcome_multiple_desired_products():
    product = {
        "identifiers": [{"value": "CO", "type": "SMILES"}],
        "is_desired_product": True,
        "reaction_role": "PRODUCT",
    }
    message = reaction_pb2.ReactionOutcome(products=[product, product])
    output = _run_validation(message)
    assert any("at most one" in warning for warning in output.warnings)


@pytest.mark.parametrize(
    ("measurement", "expects_warning"),
    [
        (
            reaction_pb2.ProductMeasurement(
                type="PURITY", analysis_key="x", float_value={"value": 5.0}
            ),
            True,
        ),
        (
            reaction_pb2.ProductMeasurement(
                type="AREA", analysis_key="x", float_value={"value": 5.0}
            ),
            False,
        ),
    ],
)
def test_product_measurement_purity_and_area(measurement, expects_warning):
    output = _run_validation(measurement)
    has_purity_warning = any(
        "PURITY measurements" in warning for warning in output.warnings
    )
    assert has_purity_warning == expects_warning


def test_provenance_clean():
    # A fully valid provenance: record_modified after record_created and a DOI
    # that needs no trimming (exercises the "no warning" branches).
    message = reaction_pb2.ReactionProvenance(
        record_created={
            "time": {"value": "2021-01-01"},
            "person": {"username": "u", "email": "a@b.com"},
        },
        record_modified=[
            {
                "time": {"value": "2021-06-01"},
                "person": {"username": "u", "email": "a@b.com"},
            }
        ],
        doi="10.1234/foo",
    )
    output = _run_validation(message)
    assert len(output.errors) == 0
    assert len(output.warnings) == 0


def test_provenance_record_modified_missing_email_direct():
    # Shadowed by validate_record_event under recursion, so call directly.
    message = reaction_pb2.ReactionProvenance()
    message.record_created.time.value = "2021-01-01"
    message.record_created.person.email = "a@b.com"
    record = message.record_modified.add()
    record.time.value = "2021-06-01"
    record.person.username = "u"
    tape = _capture_warnings(validations.validate_reaction_provenance, message)
    assert any(
        "email is required for record_modified" in str(warning.message)
        for warning in tape
    )


def test_validate_datasets(tmp_path):
    options = validations.ValidationOptions(require_provenance=False)
    good = dataset_pb2.Dataset(
        name="test",
        description="test",
        reactions=[
            reaction_pb2.Reaction(
                identifiers=[{"type": "REACTION_SMILES", "value": "CO>>CC"}]
            )
        ],
    )
    # A clean dataset does not raise.
    validations.validate_datasets({"good.pbtxt": good}, options=options)
    # This dataset has both a reaction-level error (an empty Reaction) and a
    # dataset-level error (reactions and reaction_ids are mutually exclusive);
    # with write_errors it also writes a sidecar file.
    bad = dataset_pb2.Dataset(
        name="test",
        description="test",
        reactions=[reaction_pb2.Reaction()],
        reaction_ids=["ord-c0bbd41f095a44a78b6221135961d809"],
    )
    bad_path = tmp_path / "bad.pbtxt"
    with pytest.raises(
        validations.ValidationError, match="validation encountered errors"
    ):
        validations.validate_datasets(
            {str(bad_path): bad}, write_errors=True, options=options
        )
    assert bad_path.with_suffix(".pbtxt.error").exists()
