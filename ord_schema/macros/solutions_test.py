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

from absl.testing import absltest
from absl.testing import parameterized

from ord_schema.proto import reaction_pb2
from ord_schema.macros import solutions
from ord_schema import validations


class SolutionsTest(parameterized.TestCase, absltest.TestCase):

    @parameterized.parameters(
        itertools.product(['[Na+].[Cl-]', None], ['1M', None], ['1L', None]))
    def test_simple_solution(self, solute_smiles, concentration, volume):
        """Test that the macro never generates illegal solution configurations."""
        solution_components = solutions.simple_solution(
            solvent_smiles='O',
            solute_smiles=solute_smiles,
            concentration=concentration,
            volume=volume)
        for component in solution_components:
            validations.validate_message(component)


if __name__ == '__main__':
    absltest.main()
