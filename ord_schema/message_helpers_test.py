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
"""Tests for ord_schema.message_helpers."""

import tempfile
import time

import pandas as pd
import pytest
from google.protobuf import json_format
from google.protobuf import text_format
from rdkit import Chem

from ord_schema import message_helpers
from ord_schema.proto import reaction_pb2
from ord_schema.proto import test_pb2

_BENZENE_MOLBLOCK = """241
  -OEChem-07232015262D

 12 12  0     0  0  0  0  0  0999 V2000
    2.8660    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7320    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7320   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660    1.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4631    0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690    0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4631   -0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2690   -0.8100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -1.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  1  7  1  0  0  0  0
  2  4  1  0  0  0  0
  2  8  1  0  0  0  0
  3  5  2  0  0  0  0
  3  9  1  0  0  0  0
  4  6  2  0  0  0  0
  4 10  1  0  0  0  0
  5  6  1  0  0  0  0
  5 11  1  0  0  0  0
  6 12  1  0  0  0  0
M  END"""

# pylint: disable=no-self-use


class TestMessageHelpers:
    @pytest.mark.parametrize(
        "filename,expected",
        (
            ("ord-1234567890", "data/12/ord-1234567890"),
            ("test/ord-foo.pbtxt", "data/fo/ord-foo.pbtxt"),
            ("ord_dataset-f00.pbtxt", "data/f0/ord_dataset-f00.pbtxt"),
            ("ord_data-123456foo7.jpg", "data/12/ord_data-123456foo7.jpg"),
        ),
    )
    def test_id_filename(self, filename, expected):
        assert message_helpers.id_filename(filename) == expected

    @pytest.mark.parametrize(
        "value,identifier_type,expected",
        (
            ("c1ccccc1", "SMILES", "c1ccccc1"),
            ("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H", "INCHI", "c1ccccc1"),
            (_BENZENE_MOLBLOCK, "MOLBLOCK", "c1ccccc1"),
        ),
    )
    def test_smiles_from_compound(self, value, identifier_type, expected):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value=value, type=identifier_type)
        assert message_helpers.smiles_from_compound(compound) == expected

    @pytest.mark.parametrize(
        "value,identifier_type,expected",
        (
            ("c1ccccc1", "SMILES", "c1ccccc1"),
            ("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H", "INCHI", "c1ccccc1"),
            (_BENZENE_MOLBLOCK, "MOLBLOCK", "c1ccccc1"),
        ),
    )
    def test_mol_from_compound(self, value, identifier_type, expected):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value=value, type=identifier_type)
        mol = message_helpers.mol_from_compound(compound)
        assert Chem.MolToSmiles(mol) == expected
        mol, identifier = message_helpers.mol_from_compound(compound, return_identifier=True)
        assert Chem.MolToSmiles(mol) == expected
        assert identifier == compound.identifiers[0]

    @pytest.mark.parametrize("value,identifier_type", (("invalid_smiles", "SMILES"), ("invalid_inchi", "INCHI")))
    def test_mol_from_compound_failures(self, value, identifier_type):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value=value, type=identifier_type)
        with pytest.raises(ValueError, match="invalid structural identifier"):
            message_helpers.mol_from_compound(compound)

    def test_get_product_yield(self):
        product = reaction_pb2.ProductCompound()
        assert message_helpers.get_product_yield(product) is None
        product.measurements.add(type="YIELD", percentage=dict(value=23))
        assert 23 == message_helpers.get_product_yield(product)

    def test_check_compound_identifiers(self):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(value="c1ccccc1", type="SMILES")
        message_helpers.check_compound_identifiers(compound)
        compound.identifiers.add(value="InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H", type="INCHI")
        message_helpers.check_compound_identifiers(compound)
        compound.identifiers.add(value="c1ccc(O)cc1", type="SMILES")
        with pytest.raises(ValueError, match="inconsistent"):
            message_helpers.check_compound_identifiers(compound)

    def test_get_reaction_smiles(self):
        reaction = reaction_pb2.Reaction()
        reactant1 = reaction.inputs["reactant1"]
        reactant1.components.add(reaction_role="REACTANT").identifiers.add(value="c1ccccc1", type="SMILES")
        reactant1.components.add(reaction_role="SOLVENT").identifiers.add(value="N", type="SMILES")
        assert message_helpers.get_reaction_smiles(reaction, generate_if_missing=True) == "c1ccccc1>N>"
        reactant2 = reaction.inputs["reactant2"]
        reactant2.components.add(reaction_role="REACTANT").identifiers.add(value="Cc1ccccc1", type="SMILES")
        reactant2.components.add(reaction_role="SOLVENT").identifiers.add(value="N", type="SMILES")
        reaction.outcomes.add().products.add(reaction_role="PRODUCT").identifiers.add(value="O=C=O", type="SMILES")
        assert message_helpers.get_reaction_smiles(reaction, generate_if_missing=True) == "Cc1ccccc1.c1ccccc1>N>O=C=O"

    def test_get_reaction_smiles_failure(self):
        reaction = reaction_pb2.Reaction()
        reactant = reaction.inputs["reactant"]
        component = reactant.components.add(reaction_role="REACTANT")
        component.identifiers.add(value="benzene", type="NAME")
        with pytest.raises(ValueError, match="no valid reactants or products"):
            message_helpers.get_reaction_smiles(reaction, generate_if_missing=True)
        component.identifiers.add(value="c1ccccc1", type="SMILES")
        with pytest.raises(ValueError, match="must contain at least one"):
            message_helpers.get_reaction_smiles(reaction, generate_if_missing=True, allow_incomplete=False)
        reaction.outcomes.add().products.add(reaction_role="PRODUCT").identifiers.add(value="invalid", type="SMILES")
        with pytest.raises(ValueError, match="bad reaction SMILES"):
            message_helpers.get_reaction_smiles(reaction, generate_if_missing=True, allow_incomplete=False)

    def test_reaction_from_smiles(self):
        reaction_smiles = "[C:1].N>O>F.Cl"
        expected = reaction_pb2.Reaction()
        expected.identifiers.add(value=reaction_smiles, type="REACTION_SMILES")
        this_input = expected.inputs["from_reaction_smiles"]
        c_component = this_input.components.add()
        c_component.identifiers.add(value="C", type="SMILES", details="Extracted from reaction SMILES")
        c_component.reaction_role = reaction_pb2.ReactionRole.REACTANT
        c_component.amount.unmeasured.type = reaction_pb2.UnmeasuredAmount.CUSTOM
        c_component.amount.unmeasured.details = "Extracted from reaction SMILES"
        n_component = this_input.components.add()
        n_component.identifiers.add(value="N", type="SMILES", details="Extracted from reaction SMILES")
        n_component.reaction_role = reaction_pb2.ReactionRole.REACTANT
        n_component.amount.unmeasured.type = reaction_pb2.UnmeasuredAmount.CUSTOM
        n_component.amount.unmeasured.details = "Extracted from reaction SMILES"
        o_component = this_input.components.add()
        o_component.identifiers.add(value="O", type="SMILES", details="Extracted from reaction SMILES")
        o_component.reaction_role = reaction_pb2.ReactionRole.REAGENT
        o_component.amount.unmeasured.type = reaction_pb2.UnmeasuredAmount.CUSTOM
        o_component.amount.unmeasured.details = "Extracted from reaction SMILES"
        outcome = expected.outcomes.add()
        f_component = outcome.products.add()
        f_component.identifiers.add(value="F", type="SMILES", details="Extracted from reaction SMILES")
        f_component.reaction_role = reaction_pb2.ReactionRole.PRODUCT
        cl_component = outcome.products.add()
        cl_component.identifiers.add(value="Cl", type="SMILES", details="Extracted from reaction SMILES")
        cl_component.reaction_role = reaction_pb2.ReactionRole.PRODUCT
        assert message_helpers.reaction_from_smiles(reaction_smiles) == expected

    @pytest.mark.parametrize(
        "doi,expected",
        (
            (
                "https://dx.doi.org/10.1021/acscatal.0c02247",
                "10.1021/acscatal.0c02247",
            ),
        ),
    )
    def test_parse_doi(self, doi, expected):
        assert message_helpers.parse_doi(doi) == expected

    def test_fetch_dataset(self):
        dataset = message_helpers.fetch_dataset("ord_dataset-35a5a513f1dd44a3a97c88da99f81a00")
        assert len(dataset.reactions) == 7


class TestFindSubmessages:
    def test_scalar(self):
        message = test_pb2.Scalar(int32_value=5, float_value=6.7)
        assert len(message_helpers.find_submessages(message, test_pb2.Scalar)) == 0
        with pytest.raises(TypeError, match="must be a Protocol Buffer"):
            message_helpers.find_submessages(message, float)

    def test_nested(self):
        message = test_pb2.Nested()
        assert len(message_helpers.find_submessages(message, test_pb2.Nested.Child)) == 0
        message.child.value = 5.6
        submessages = message_helpers.find_submessages(message, test_pb2.Nested.Child)
        assert len(submessages) == 1
        # Show that the returned submessages work as references.
        submessages[0].value = 7.8
        assert round(abs(message.child.value - 7.8), 4) == 0

    def test_repeated_nested(self):
        message = test_pb2.RepeatedNested()
        message.children.add().value = 1.2
        message.children.add().value = 3.4
        assert len(message_helpers.find_submessages(message, test_pb2.RepeatedNested.Child)) == 2

    def test_map_nested(self):
        message = test_pb2.MapNested()
        message.children["one"].value = 1.2
        message.children["two"].value = 3.4
        assert len(message_helpers.find_submessages(message, test_pb2.MapNested.Child)) == 2

    def test_compounds(self):
        message = reaction_pb2.Reaction()
        message.inputs["test"].components.add().identifiers.add(type="NAME", value="aspirin")
        assert len(message_helpers.find_submessages(message, reaction_pb2.Compound)) == 1


class TestBuildData:
    def test_build_data(self, tmp_path):
        data = b"test data"
        filename = (tmp_path / "test.data").as_posix()
        with open(filename, "wb") as f:
            f.write(data)
        message = message_helpers.build_data(filename, description="binary data")
        assert message.bytes_value == data
        assert message.description == "binary data"
        assert message.format == "data"

    def test_bad_filename(self):
        with pytest.raises(ValueError, match="cannot deduce the file format"):
            message_helpers.build_data("testdata", "no description")


class TestBuildCompound:
    def test_smiles_and_name(self):
        compound = message_helpers.build_compound(smiles="c1ccccc1", name="benzene")
        expected = reaction_pb2.Compound(
            identifiers=[
                reaction_pb2.CompoundIdentifier(value="c1ccccc1", type="SMILES"),
                reaction_pb2.CompoundIdentifier(value="benzene", type="NAME"),
            ]
        )
        assert compound == expected

    def test_unmeasured_amount(self):
        compound = message_helpers.build_compound(smiles="CCO", amount="catalytic")
        expected = reaction_pb2.Compound(
            identifiers=[
                reaction_pb2.CompoundIdentifier(value="CCO", type="SMILES"),
            ],
            amount=dict(unmeasured=dict(type="CATALYTIC")),
        )
        assert compound == expected

    @pytest.mark.parametrize(
        "amount,expected",
        (
            ("1.2 g", reaction_pb2.Mass(value=1.2, units="GRAM")),
            ("3.4 mol", reaction_pb2.Moles(value=3.4, units="MOLE")),
            ("5.6 mL", reaction_pb2.Volume(value=5.6, units="MILLILITER")),
        ),
    )
    def test_amount(self, amount, expected):
        compound = message_helpers.build_compound(amount=amount)
        assert getattr(compound.amount, compound.amount.WhichOneof("kind")) == expected

    @pytest.mark.parametrize("amount", ("1.2", "-3.4 g"))
    def test_bad_amount(self, amount):
        with pytest.raises((KeyError, ValueError)):
            message_helpers.build_compound(amount=amount)

    def test_role(self):
        compound = message_helpers.build_compound(role="solvent")
        assert compound.reaction_role == reaction_pb2.ReactionRole.SOLVENT

    def test_bad_role(self):
        with pytest.raises(KeyError, match="not a supported type"):
            message_helpers.build_compound(role="flavorant")

    def test_is_limiting(self):
        assert message_helpers.build_compound(is_limiting=True).is_limiting
        assert not message_helpers.build_compound(is_limiting=False).is_limiting
        assert not message_helpers.build_compound().HasField("is_limiting")

    @pytest.mark.parametrize(
        "prep,details,expected",
        (
            ("dried", None, reaction_pb2.CompoundPreparation(type="DRIED")),
            (
                "dried",
                "in the fire of the sun",
                reaction_pb2.CompoundPreparation(type="DRIED", details="in the fire of the sun"),
            ),
            (
                "custom",
                "threw it on the ground",
                reaction_pb2.CompoundPreparation(type="CUSTOM", details="threw it on the ground"),
            ),
        ),
    )
    def test_prep(self, prep, details, expected):
        compound = message_helpers.build_compound(prep=prep, prep_details=details)
        assert compound.preparations[0] == expected

    def test_bad_prep(self):
        with pytest.raises(KeyError, match="not a supported type"):
            message_helpers.build_compound(prep="shaken")

    def test_prep_details_without_prep(self):
        with pytest.raises(ValueError, match="prep must be provided"):
            message_helpers.build_compound(prep_details="rinsed gently")

    def test_custom_prep_without_details(self):
        with pytest.raises(ValueError, match="prep_details must be provided"):
            message_helpers.build_compound(prep="custom")

    def test_vendor(self):
        assert message_helpers.build_compound(vendor="Sally").source.vendor == "Sally"


class TestSetSoluteMoles:
    def test_set_solute_moles_should_fail(self):
        solute = message_helpers.build_compound(name="Solute")
        solvent = message_helpers.build_compound(name="Solvent")
        with pytest.raises(ValueError, match="defined volume"):
            message_helpers.set_solute_moles(solute, [solvent], "10 mM")

        solute = message_helpers.build_compound(name="Solute", amount="1 mol")
        solvent = message_helpers.build_compound(name="Solvent", amount="1 L")
        with pytest.raises(ValueError, match="overwrite"):
            message_helpers.set_solute_moles(solute, [solvent], "10 mM")

    def test_set_solute_moles(self):
        solute = message_helpers.build_compound(name="Solute")
        solvent2 = message_helpers.build_compound(name="Solvent", amount="100 mL")
        message_helpers.set_solute_moles(solute, [solvent2], "1 molar")
        assert solute.amount.moles == reaction_pb2.Moles(units="MILLIMOLE", value=100)
        solvent3 = message_helpers.build_compound(name="Solvent", amount="75 uL")
        message_helpers.set_solute_moles(solute, [solvent3], "3 mM", overwrite=True)
        assert solute.amount.moles.units == reaction_pb2.Moles.NANOMOLE
        assert round(abs(solute.amount.moles.value - 225), 4) == 0
        solvent4 = message_helpers.build_compound(name="Solvent", amount="0.2 uL")
        message_helpers.set_solute_moles(solute, [solvent4], "30 mM", overwrite=True)
        assert solute.amount.moles == reaction_pb2.Moles(units="NANOMOLE", value=6)
        solvent5 = message_helpers.build_compound(name="Solvent", amount="0.8 uL")
        message_helpers.set_solute_moles(solute, [solvent4, solvent5], "30 mM", overwrite=True)
        assert solute.amount.moles == reaction_pb2.Moles(units="NANOMOLE", value=30)


class TestCompoundIdentifiers:
    def test_identifier_setters(self):
        compound = reaction_pb2.Compound()
        identifier = message_helpers.set_compound_name(compound, "water")
        assert identifier == reaction_pb2.CompoundIdentifier(type="NAME", value="water")
        assert compound.identifiers[0] == reaction_pb2.CompoundIdentifier(type="NAME", value="water")
        message_helpers.set_compound_smiles(compound, "O")
        assert compound.identifiers[1] == reaction_pb2.CompoundIdentifier(type="SMILES", value="O")
        identifier = message_helpers.set_compound_name(compound, "ice")
        assert identifier == reaction_pb2.CompoundIdentifier(type="NAME", value="ice")
        assert compound.identifiers[0] == reaction_pb2.CompoundIdentifier(type="NAME", value="ice")
        compound = reaction_pb2.Compound()
        _ = message_helpers.set_compound_molblock(compound, _BENZENE_MOLBLOCK)
        assert _BENZENE_MOLBLOCK == compound.identifiers[0].value

    def test_identifier_getters(self):
        compound = reaction_pb2.Compound()
        compound.identifiers.add(type="NAME", value="water")
        assert message_helpers.get_compound_name(compound) == "water"
        assert message_helpers.get_compound_smiles(compound) is None
        compound.identifiers.add(type="SMILES", value="O")
        assert message_helpers.get_compound_smiles(compound) == "O"
        assert message_helpers.smiles_from_compound(compound) == "O"
        compound = reaction_pb2.Compound()
        compound.identifiers.add(type="MOLBLOCK", value=_BENZENE_MOLBLOCK)
        assert message_helpers.get_compound_molblock(compound) == _BENZENE_MOLBLOCK
        assert message_helpers.molblock_from_compound(compound) == _BENZENE_MOLBLOCK


class TestSetDativeBonds:
    def test_has_transition_metal(self):
        assert not message_helpers.has_transition_metal(Chem.MolFromSmiles("P"))
        assert message_helpers.has_transition_metal(Chem.MolFromSmiles("Cl[Pd]Cl"))

    @pytest.mark.parametrize(
        "smiles,from_atoms,expected", (("[PH3][Pd](Cl)(Cl)[NH3]", ("N", "P"), "N->[Pd](<-P)(Cl)Cl"),)
    )
    def test_set_dative_bonds(self, smiles, from_atoms, expected):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        dative_mol = message_helpers.set_dative_bonds(mol, from_atoms=from_atoms)
        assert Chem.MolToSmiles(dative_mol) == expected


class TestLoadAndWriteMessage:
    @pytest.fixture
    def messages(self) -> list:
        yield [
            test_pb2.Scalar(int32_value=3, float_value=4.5),
            test_pb2.RepeatedScalar(values=[1.2, 3.4]),
            test_pb2.Enum(value="FIRST"),
            test_pb2.RepeatedEnum(values=["FIRST", "SECOND"]),
            test_pb2.Nested(child=test_pb2.Nested.Child(value=1.2)),
        ]

    @pytest.mark.parametrize("suffix", (".pbtxt", ".pb", ".json", ".pbtxt.gz", ".pb.gz", ".json.gz"))
    def test_round_trip(self, suffix, messages):
        for message in messages:
            with tempfile.NamedTemporaryFile(suffix=suffix) as f:
                message_helpers.write_message(message, f.name)
                f.flush()
                assert message == message_helpers.load_message(f.name, type(message))

    def test_gzip_reproducibility(self, messages, tmp_path):
        filename = (tmp_path / "test.pb.gz").as_posix()
        for message in messages:
            message_helpers.write_message(message, filename)
            with open(filename, "rb") as f:
                value = f.read()
            time.sleep(1)
            message_helpers.write_message(message, filename)
            with open(filename, "rb") as f:
                assert f.read() == value

    def test_bad_binary(self):
        with tempfile.NamedTemporaryFile(suffix=".pb") as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(message.SerializeToString())
            f.flush()
            # NOTE(kearnes): The decoder is not perfect; for example, it will
            # not be able to distinguish from a message with the same tags and
            # types (e.g. test_pb2.Scalar and test_pb2.RepeatedScalar).
            with pytest.raises(ValueError, match="parsing"):
                message_helpers.load_message(f.name, test_pb2.Nested)

    def test_bad_json(self):
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".json") as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(json_format.MessageToJson(message))
            f.flush()
            with pytest.raises(ValueError, match='no field named "values"'):
                message_helpers.load_message(f.name, test_pb2.Nested)

    def test_bad_pbtxt(self):
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".pbtxt") as f:
            message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
            f.write(text_format.MessageToString(message))
            f.flush()
            with pytest.raises(ValueError, match='no field named "values"'):
                message_helpers.load_message(f.name, test_pb2.Nested)

    def test_bad_suffix(self):
        message = test_pb2.RepeatedScalar(values=[1.2, 3.4])
        with pytest.raises(ValueError, match="not a valid MessageFormat"):
            message_helpers.write_message(message, "test.proto")


class TestCreateMessage:
    @pytest.mark.parametrize(
        "message_name,expected_class",
        (
            ("Reaction", reaction_pb2.Reaction),
            ("Temperature", reaction_pb2.Temperature),
            ("TemperatureConditions.Measurement", reaction_pb2.TemperatureConditions.Measurement),
        ),
    )
    def test_valid_messages(self, message_name, expected_class):
        message = message_helpers.create_message(message_name)
        assert isinstance(message, expected_class)

    @pytest.mark.parametrize("message_name", ("reaction", "aosdifjasdf"))
    def test_invalid_messages(self, message_name):
        with pytest.raises(ValueError, match="Cannot resolve"):
            message_helpers.create_message(message_name)


class TestMessagesToDataFrame:
    @pytest.mark.parametrize(
        "message,expected",
        (
            (test_pb2.Scalar(int32_value=3, float_value=4.5), {"int32_value": 3, "float_value": 4.5}),
            (
                test_pb2.Scalar(int32_value=0, float_value=0.0, bool_value=False),
                {"float_value": 0.0, "bool_value": False},
            ),
            (test_pb2.RepeatedScalar(values=[1.2, 3.4]), {"values[0]": 1.2, "values[1]": 3.4}),
            (test_pb2.Enum(value="FIRST"), {"value": "FIRST"}),
            (test_pb2.RepeatedEnum(values=["FIRST", "SECOND"]), {"values[0]": "FIRST", "values[1]": "SECOND"}),
            (test_pb2.Nested(child=test_pb2.Nested.Child(value=1.2)), {"child.value": 1.2}),
            (
                test_pb2.RepeatedNested(
                    children=[test_pb2.RepeatedNested.Child(value=1.2), test_pb2.RepeatedNested.Child(value=3.4)]
                ),
                {"children[0].value": 1.2, "children[1].value": 3.4},
            ),
            (test_pb2.Map(values={"a": 1.2, "b": 3.4}), {'values["a"]': 1.2, 'values["b"]': 3.4}),
            (
                test_pb2.MapNested(
                    children={"a": test_pb2.MapNested.Child(value=1.2), "b": test_pb2.MapNested.Child(value=3.4)}
                ),
                {'children["a"].value': 1.2, 'children["b"].value': 3.4},
            ),
        ),
    )
    def test_message_to_row(self, message, expected):
        row = message_helpers.message_to_row(message)
        pd.testing.assert_frame_equal(pd.DataFrame([row]), pd.DataFrame([expected]), check_like=True)

    def test_messages_to_dataframe(self):
        reaction1 = reaction_pb2.Reaction()
        input_test = reaction1.inputs["test"]
        outcome = reaction1.outcomes.add()
        component = input_test.components.add()
        component.identifiers.add(type="SMILES", value="CCO")
        component.is_limiting = True
        component.amount.mass.value = 1.2
        component.amount.mass.units = reaction_pb2.Mass.GRAM
        outcome.conversion.value = 3.4
        outcome.conversion.precision = 5.6
        reaction2 = reaction_pb2.Reaction()
        input_test = reaction2.inputs["test"]
        outcome = reaction2.outcomes.add()
        component = input_test.components.add()
        component.identifiers.add(type="SMILES", value="CCC")
        component.is_limiting = False
        component.amount.mass.value = 7.8
        component.amount.mass.units = reaction_pb2.Mass.GRAM
        outcome.conversion.value = 9.1
        outcome.conversion.precision = 2.3
        expected = pd.DataFrame(
            {
                'inputs["test"].components[0].identifiers[0].type': "SMILES",
                'inputs["test"].components[0].identifiers[0].value': ["CCO", "CCC"],
                'inputs["test"].components[0].is_limiting': [True, False],
                'inputs["test"].components[0].amount.mass.value': [1.2, 7.8],
                'inputs["test"].components[0].amount.mass.units': "GRAM",
                "outcomes[0].conversion.value": [3.4, 9.1],
                "outcomes[0].conversion.precision": [5.6, 2.3],
            }
        )
        pd.testing.assert_frame_equal(
            message_helpers.messages_to_dataframe([reaction1, reaction2]),
            expected,
            check_like=True,
        )
        # Drop constant columns and test again with drop_constant_columns=True.
        del expected['inputs["test"].components[0].identifiers[0].type']
        del expected['inputs["test"].components[0].amount.mass.units']
        pd.testing.assert_frame_equal(
            message_helpers.messages_to_dataframe([reaction1, reaction2], drop_constant_columns=True),
            expected,
            check_like=True,
        )
