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
"""Tests for ord_schema.templating."""

import pandas as pd
import pytest
from google.protobuf import text_format
from rdkit import Chem

from ord_schema import templating
from ord_schema.proto import reaction_pb2
from ord_schema.proto import dataset_pb2


@pytest.fixture
def valid_reaction() -> reaction_pb2.Reaction:
    message = reaction_pb2.Reaction()
    dummy_input = message.inputs["in"]
    outcome = message.outcomes.add()
    component = dummy_input.components.add()
    component.identifiers.add(type="SMILES", value="CCO")
    component.is_limiting = True
    component.amount.mass.value = 1
    component.amount.mass.units = reaction_pb2.Mass.GRAM
    outcome.conversion.value = 75
    outcome.conversion.precision = 99
    yield message


def test_valid_templating(valid_reaction):
    template_string = text_format.MessageToString(valid_reaction)
    template_string = template_string.replace('value: "CCO"', 'value: "$smiles$"')
    template_string = template_string.replace("value: 75", "value: $conversion$")
    df = pd.DataFrame.from_dict({"$smiles$": ["CCO", "CCCO", "CCCCO"], "$conversion$": [75, 50, 30]})
    dataset = templating.generate_dataset(template_string, df)
    expected_reactions = []
    for smiles, conversion in zip(["CCO", "CCCO", "CCCCO"], [75, 50, 30]):
        reaction = reaction_pb2.Reaction()
        reaction.CopyFrom(valid_reaction)
        reaction.inputs["in"].components[0].identifiers[0].value = smiles
        reaction.outcomes[0].conversion.value = conversion
        expected_reactions.append(reaction)
    expected_dataset = dataset_pb2.Dataset(reactions=expected_reactions)
    assert dataset == expected_dataset

    # Test without "$" in column names
    df = pd.DataFrame.from_dict({"smiles": ["CCO", "CCCO", "CCCCO"], "conversion": [75, 50, 30]})
    dataset = templating.generate_dataset(template_string, df)
    assert dataset == expected_dataset


def test_valid_templating_escapes(valid_reaction):
    smiles = ["CCO", "CCCO", "CCCCO"]
    mols = [Chem.MolFromSmiles(this_smiles) for this_smiles in smiles]
    molblocks = [Chem.MolToMolBlock(mol) for mol in mols]
    valid_reaction.inputs["in"].components[0].identifiers.add(type="MOLBLOCK", value="$molblock$")
    template_string = text_format.MessageToString(valid_reaction)
    df = pd.DataFrame.from_dict({"molblock": molblocks})
    dataset = templating.generate_dataset(template_string, df)
    expected_reactions = []
    for molblock in molblocks:
        reaction = reaction_pb2.Reaction()
        reaction.CopyFrom(valid_reaction)
        reaction.inputs["in"].components[0].identifiers[1].value = molblock
        expected_reactions.append(reaction)
    expected_dataset = dataset_pb2.Dataset(reactions=expected_reactions)
    assert dataset == expected_dataset


@pytest.mark.parametrize("suffix", [".csv", ".xls", ".xlsx"])
def test_read_spreadsheet(suffix, tmp_path):
    df = pd.DataFrame.from_dict({"smiles": ["CCO", "CCCO", "CCCCO"], "conversion": [75, 50, 30]})
    filename = (tmp_path / f"test{suffix}").as_posix()
    if suffix == ".csv":
        df.to_csv(filename, index=False)
    else:
        df.to_excel(filename, index=False)
    pd.testing.assert_frame_equal(templating.read_spreadsheet(filename), df)


def test_invalid_templating(valid_reaction):
    template_string = text_format.MessageToString(valid_reaction)
    template_string = template_string.replace('value: "CCO"', 'value: "$my_smiles$"')
    template_string = template_string.replace("precision: 99", "precision: $precision$")
    df = pd.DataFrame.from_dict({"$my_smiles$": ["CCO", "CCCO", "CCCCO"], "$precision$": [75, 50, -5]})
    expected_reactions = []
    for smiles, precision in zip(["CCO", "CCCO", "CCCCO"], [75, 50, -5]):
        reaction = reaction_pb2.Reaction()
        reaction.CopyFrom(valid_reaction)
        reaction.inputs["in"].components[0].identifiers[0].value = smiles
        reaction.outcomes[0].conversion.precision = precision
        expected_reactions.append(reaction)
    expected_dataset = dataset_pb2.Dataset(reactions=expected_reactions)
    with pytest.raises(ValueError, match="Enumerated Reaction is not valid"):
        templating.generate_dataset(template_string, df)
    dataset = templating.generate_dataset(template_string, df, validate=False)
    assert dataset == expected_dataset


def test_bad_placeholders(valid_reaction):
    template_string = text_format.MessageToString(valid_reaction)
    template_string = template_string.replace('value: "CCO"', 'value: "$my_smiles$"')
    template_string = template_string.replace("value: 75", "value: $conversion$")
    df = pd.DataFrame.from_dict({"$my_smiles$": ["CCO", "CCCO", "CCCCO"]})
    with pytest.raises(ValueError, match=r"\$conversion\$ not found"):
        templating.generate_dataset(template_string, df)


def test_missing_values(tmp_path):
    # pylint: disable=too-many-locals
    # Build a template reaction.
    reaction = reaction_pb2.Reaction()
    input1 = reaction.inputs["one"]
    input1_component1 = input1.components.add()
    input1_component1.identifiers.add(value="CCO", type="SMILES")
    input1_component1.amount.mass.value = 1.2
    input1_component1.amount.mass.units = reaction_pb2.Mass.GRAM
    input1_component2 = input1.components.add()
    input1_component2.is_limiting = True
    input1_component2.identifiers.add(value="c1ccccc1", type="SMILES")
    input1_component2.amount.volume.value = 3.4
    input1_component2.amount.volume.units = reaction_pb2.Volume.LITER
    input2 = reaction.inputs["two"]
    input2_component1 = input2.components.add()
    input2_component1.identifiers.add(value="COO", type="SMILES")
    input2_component1.amount.mass.value = 5.6
    input2_component1.amount.mass.units = reaction_pb2.Mass.GRAM
    outcome = reaction.outcomes.add()
    outcome.conversion.value = 75
    template_string = text_format.MessageToString(reaction)
    template_string = template_string.replace('value: "CCO"', 'value: "$smiles$"')
    template_string = template_string.replace("value: 5.6", "value: $mass$")
    # Build a spreadsheet and test for proper edits.
    filename = (tmp_path / "missing.csv").as_posix()
    with open(filename, "w") as f:
        f.write("smiles,mass\n")
        f.write("CN,\n")  # Missing mass.
        f.write(",1.5\n")  # Missing SMILES.
    df = pd.read_csv(filename)
    dataset = templating.generate_dataset(template_string, df)
    expected_dataset = dataset_pb2.Dataset()
    reaction1 = expected_dataset.reactions.add()
    reaction1.CopyFrom(reaction)
    reaction1.inputs["one"].components[0].identifiers[0].value = "CN"
    del reaction1.inputs["two"]  # No components after amount removal.
    reaction2 = expected_dataset.reactions.add()
    reaction2.CopyFrom(reaction)
    del reaction2.inputs["one"].components[0]  # No indentifiers.
    reaction2.inputs["two"].components[0].amount.mass.value = 1.5
    assert dataset == expected_dataset
