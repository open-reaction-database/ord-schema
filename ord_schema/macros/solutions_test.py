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
"""Tests for ord_schema.macros.solutions."""
import itertools

import pytest
from google.protobuf import text_format

from ord_schema.proto import reaction_pb2
from ord_schema.macros import solutions
from ord_schema import validations


def test_simple_solution():
    """Simple input/output pair test."""
    solvent_pbtxt = """
        identifiers {
            type: SMILES
            value: "O"
        }
        amount {
            volume {
                value: 1.0
                units: LITER
            }
            volume_includes_solutes: true
        }
    """
    solute_pbtxt = """
        identifiers {
            type: SMILES
            value: "[Na+].[Cl-]"
        }
        amount {
            moles {
                value: 1.0
                units: MOLE
            }
        }
    """
    expected_solvent = reaction_pb2.Compound()
    expected_solute = reaction_pb2.Compound()
    text_format.Parse(solvent_pbtxt, expected_solvent)
    text_format.Parse(solute_pbtxt, expected_solute)
    actual_compounds = solutions.simple_solution(
        solvent_smiles="O",
        solute_smiles="[Na+].[Cl-]",
        concentration="1M",
        volume="1L",
    )
    assert actual_compounds == [expected_solvent, expected_solute]


@pytest.mark.parametrize(
    "solute_smiles,concentration,volume", itertools.product(["[Na+].[Cl-]", None], ["1M", None], ["1L", None])
)
def test_simple_solution_legal(solute_smiles, concentration, volume):
    """Test that the macro never generates illegal solution configurations."""
    solution_components = solutions.simple_solution(
        solvent_smiles="O",
        solute_smiles=solute_smiles,
        concentration=concentration,
        volume=volume,
    )
    for component in solution_components:
        validations.validate_message(component)


@pytest.mark.parametrize("volume", ["1L", None])
def test_simple_saturated_solution(volume):
    """Test that the macro never generates illegal solution configurations."""
    solution_components = solutions.simple_solution(
        solvent_smiles="O",
        solute_smiles="[Na+].[Cl-]",
        volume=volume,
        saturated=True,
    )
    for component in solution_components:
        validations.validate_message(component)
    assert any(
        component.amount.unmeasured.type == reaction_pb2.UnmeasuredAmount.SATURATED for component in solution_components
    )


def test_simple_solution_concentration_and_saturated_illegal():
    with pytest.raises(ValueError, match="Cannot specify both"):
        solutions.simple_solution(
            solvent_smiles="O",
            solute_smiles="[Na+].[Cl-]",
            volume="1L",
            concentration="1M",
            saturated=True,
        )


def test_simple_solution_saturated_solute_missing():
    with pytest.raises(ValueError, match="Must specify"):
        solutions.simple_solution(solvent_smiles="O", solute_smiles=None, volume="1L", saturated=True)
