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
"""Tests for ord_schema.scripts.enumerate_dataset."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
import pandas as pd

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2
from ord_schema.scripts import enumerate_dataset


class EnumerateDatasetTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        template_string = """
        inputs {
            key: "test"
            value {
                components {
                    identifiers {
                        type: SMILES
                        value: "$input_smiles$"
                    }
                    mass {
                        value: $input_mass$
                        units: GRAM
                    }
                }
            }
        }
        outcomes {
            products {
                compound {
                    identifiers {
                        type: SMILES
                        value: "$product_smiles$"
                    }
                }
                compound_yield {
                    value: $product_yield$
                }
            }
        }
        """
        self.template = os.path.join(self.test_subdirectory, 'template.pbtxt')
        with open(self.template, 'w') as f:
            f.write(template_string)
        data = pd.DataFrame({
            'input_smiles': ['C', 'CC', 'CCC'],
            'input_mass': [1.2, 3.4, 5.6],
            'product_smiles': ['CO', 'CCO', 'CCCO'],
            'product_yield': [7.8, 9.0, 8.7],
        })
        self.spreadsheet = os.path.join(self.test_subdirectory,
                                        'spreadsheet.csv')
        data.to_csv(self.spreadsheet, index=False)
        self.expected = dataset_pb2.Dataset()
        r1 = self.expected.reactions.add()
        r1c1 = r1.inputs['test'].components.add()
        r1c1.identifiers.add(value='C', type='SMILES')
        r1c1.mass.CopyFrom(reaction_pb2.Mass(value=1.2, units='GRAM'))
        r1p1 = r1.outcomes.add().products.add()
        r1p1.compound.identifiers.add(value='CO', type='SMILES')
        r1p1.compound_yield.value = 7.8
        r2 = self.expected.reactions.add()
        r2c1 = r2.inputs['test'].components.add()
        r2c1.identifiers.add(value='CC', type='SMILES')
        r2c1.mass.CopyFrom(reaction_pb2.Mass(value=3.4, units='GRAM'))
        r2p1 = r2.outcomes.add().products.add()
        r2p1.compound.identifiers.add(value='CCO', type='SMILES')
        r2p1.compound_yield.value = 9.0
        r3 = self.expected.reactions.add()
        r3c1 = r3.inputs['test'].components.add()
        r3c1.identifiers.add(value='CCC', type='SMILES')
        r3c1.mass.CopyFrom(reaction_pb2.Mass(value=5.6, units='GRAM'))
        r3p1 = r3.outcomes.add().products.add()
        r3p1.compound.identifiers.add(value='CCCO', type='SMILES')
        r3p1.compound_yield.value = 8.7

    def test_main(self):
        output_filename = os.path.join(self.test_subdirectory, 'dataset.pbtxt')
        with flagsaver.flagsaver(template=self.template,
                                 spreadsheet=self.spreadsheet,
                                 output=output_filename):
            enumerate_dataset.main(())
        self.assertTrue(os.path.exists(output_filename))
        dataset = message_helpers.load_message(output_filename,
                                               dataset_pb2.Dataset)
        self.assertLen(dataset.reactions, 3)
        validations.validate_message(dataset, raise_on_error=True)
        self.assertEqual(dataset, self.expected)


if __name__ == '__main__':
    absltest.main()
