"""Tests for ord_schema.update_reactions."""

import glob
import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from google.protobuf import text_format

from ord_schema import message_helpers
from ord_schema import update_reactions
from ord_schema.proto import reaction_pb2


class UpdateReactionsTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        reaction1 = reaction_pb2.Reaction()
        self.reaction1_filename = os.path.join(
            self.test_subdirectory, 'reaction1.pbtxt')
        with open(self.reaction1_filename, 'w') as f:
            f.write(text_format.MessageToString(reaction1))
        reaction2 = reaction_pb2.Reaction()
        reaction2.provenance.record_id = 'ord-test'
        self.reaction2_filename = os.path.join(
            self.test_subdirectory, 'reaction2.pbtxt')
        with open(self.reaction2_filename, 'w') as f:
            f.write(text_format.MessageToString(reaction2))

    def test_main(self):
        input_file = os.path.join(self.test_subdirectory, 'input_file.txt')
        with open(input_file, 'w') as f:
            f.write(f'A\t{self.reaction1_filename}\n')
            f.write(f'M\t{self.reaction2_filename}\n')
        output_dir = os.path.join(self.test_subdirectory, 'data')
        with flagsaver.flagsaver(input_file=input_file, root_dir=output_dir):
            update_reactions.main(())
        filenames = glob.glob(os.path.join(output_dir, '*', '*.pbtxt'))
        self.assertLen(filenames, 2)
        record_ids = []
        for filename in filenames:
            message = message_helpers.load_message(
                filename, reaction_pb2.Reaction, input_format='pbtxt')
            self.assertTrue(message.provenance.record_id)
            record_ids.append(message.provenance.record_id)
        self.assertIn('ord-test', record_ids)


if __name__ == '__main__':
    absltest.main()
