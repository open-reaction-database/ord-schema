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

from absl.testing import absltest
from google.protobuf import text_format

from ord_schema import templating
from ord_schema.proto import reaction_pb2
from ord_schema.proto import dataset_pb2


class TemplatingTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        message = reaction_pb2.Reaction()
        dummy_input = message.inputs['in']
        outcome = message.outcomes.add()
        component = dummy_input.components.add()
        component.identifiers.add(type='SMILES', value='CCO')
        component.is_limiting = True
        component.mass.value = 1
        component.mass.units = reaction_pb2.Mass.GRAM
        outcome.conversion.value = 75
        outcome.conversion.precision = 99
        self.valid_reaction = message

    def test_valid_templating(self):
        template_string = text_format.MessageToString(self.valid_reaction)
        template_string = template_string.replace('value: "CCO"',
                                                  'value: "$my_smiles$"')
        template_string = template_string.replace('value: 75',
                                                  'value: $conversion$')
        df = pd.DataFrame.from_dict({
            '$my_smiles$': ['CCO', 'CCCO', 'CCCCO'],
            '$conversion$': ['75', '50', '30'],
        })
        dataset = templating.generate_dataset(template_string, df)
        expected_reactions = []
        for smiles, conversion in zip(['CCO', 'CCCO', 'CCCCO'], [75, 50, 30]):
            reaction = reaction_pb2.Reaction()
            reaction.CopyFrom(self.valid_reaction)
            reaction.inputs['in'].components[0].identifiers[0].value = smiles
            reaction.outcomes[0].conversion.value = conversion
            expected_reactions.append(reaction)
        expected_dataset = dataset_pb2.Dataset(reactions=expected_reactions)
        self.assertEqual(dataset, expected_dataset)

        # Test without "$" in column names
        df = pd.DataFrame.from_dict({
            'my_smiles': ['CCO', 'CCCO', 'CCCCO'],
            'conversion': ['75', '50', '30'],
        })
        dataset = templating.generate_dataset(template_string, df)
        self.assertEqual(dataset, expected_dataset)

    def test_invalid_templating(self):
        template_string = text_format.MessageToString(self.valid_reaction)
        template_string = template_string.replace('value: "CCO"',
                                                  'value: "$my_smiles$"')
        template_string = template_string.replace('precision: 99',
                                                  'precision: $precision$')
        df = pd.DataFrame.from_dict({
            '$my_smiles$': ['CCO', 'CCCO', 'CCCCO'],
            '$precision$': ['75', '50', '-5'],
        })
        expected_reactions = []
        for smiles, precision in zip(['CCO', 'CCCO', 'CCCCO'], [75, 50, -5]):
            reaction = reaction_pb2.Reaction()
            reaction.CopyFrom(self.valid_reaction)
            reaction.inputs['in'].components[0].identifiers[0].value = smiles
            reaction.outcomes[0].conversion.precision = precision
            expected_reactions.append(reaction)
        expected_dataset = dataset_pb2.Dataset(reactions=expected_reactions)
        with self.assertRaisesRegex(ValueError,
                                    'Enumerated Reaction is not valid'):
            templating.generate_dataset(template_string, df)
        dataset = templating.generate_dataset(template_string,
                                              df,
                                              validate=False)
        self.assertEqual(dataset, expected_dataset)

    def test_bad_placeholders(self):
        template_string = text_format.MessageToString(self.valid_reaction)
        template_string = template_string.replace('value: "CCO"',
                                                  'value: "$my_smiles$"')
        template_string = template_string.replace('value: 75',
                                                  'value: $conversion$')
        df = pd.DataFrame.from_dict({
            '$my_smiles$': ['CCO', 'CCCO', 'CCCCO'],
        })
        with self.assertRaisesRegex(ValueError, r'\$conversion\$ not found'):
            templating.generate_dataset(template_string, df)


if __name__ == '__main__':
    absltest.main()
