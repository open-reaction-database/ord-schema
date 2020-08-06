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
"""Tests for ord_schema.interface.build_database."""

import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
import pandas as pd

from ord_schema import message_helpers
from ord_schema.interface import build_database
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


class BuildDatabaseTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        reaction = reaction_pb2.Reaction()
        reaction.reaction_id = 'test'
        reaction.identifiers.add(value='reaction', type='REACTION_SMILES')
        input1 = reaction.inputs['input1']
        input1.components.add().identifiers.add(value='input1', type='SMILES')
        input2 = reaction.inputs['input2']
        input2.components.add().identifiers.add(value='input2a', type='SMILES')
        input2.components.add().identifiers.add(value='input2b', type='SMILES')
        outcome = reaction.outcomes.add()
        product = outcome.products.add()
        product.compound_yield.value = 2.5
        product.compound.identifiers.add(value='product', type='SMILES')
        self.dataset = dataset_pb2.Dataset(reactions=[reaction])
        message_helpers.write_message(
            self.dataset, os.path.join(self.test_subdirectory, 'test.pbtxt'))

    def test_main(self):
        input_pattern = os.path.join(self.test_subdirectory, '*.pbtxt')
        output_dir = os.path.join(self.test_subdirectory, 'tables')
        with flagsaver.flagsaver(input=input_pattern, output=output_dir,
                                 database=True):
            build_database.main(())
        with open(os.path.join(output_dir, 'reactions.csv')) as f:
            df = pd.read_csv(f)
        # NOTE(kearnes): Map keys are not always serialized in the same order.
        df['deserialized'] = df.serialized.apply(
            lambda x: reaction_pb2.Reaction.FromString(bytes.fromhex(x)))
        del df['serialized']
        pd.testing.assert_frame_equal(
            df,
            pd.DataFrame({
                'reaction_id': ['test'],
                'reaction_smiles': ['reaction'],
                'deserialized': [self.dataset.reactions[0]],
            }))
        with open(os.path.join(output_dir, 'inputs.csv')) as f:
            df = pd.read_csv(f)
        pd.testing.assert_frame_equal(
            df,
            pd.DataFrame({
                'reaction_id': ['test', 'test', 'test'],
                'smiles': ['input1', 'input2a', 'input2b'],
            }))
        with open(os.path.join(output_dir, 'outputs.csv')) as f:
            df = pd.read_csv(f)
        pd.testing.assert_frame_equal(
            df,
            pd.DataFrame({
                'reaction_id': ['test'],
                'smiles': ['product'],
                'yield': [2.5]
            }))


if __name__ == '__main__':
    absltest.main()
