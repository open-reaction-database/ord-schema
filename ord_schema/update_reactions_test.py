"""Tests for ord_schema.update_reactions."""

import glob
import os
import tempfile

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver

from ord_schema import message_helpers
from ord_schema import update_reactions
from ord_schema.proto import reaction_pb2


class UpdateReactionsTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.test_subdirectory = tempfile.mkdtemp(dir=flags.FLAGS.test_tmpdir)
        reaction1 = reaction_pb2.Reaction()
        self.reaction1_filename = os.path.join(
            self.test_subdirectory, 'reaction1.pb')
        with open(self.reaction1_filename, 'wb') as f:
            f.write(reaction1.SerializeToString())
        # Reaction 2 already exists in the database.
        reaction2 = reaction_pb2.Reaction()
        reaction2.provenance.record_id = 'ord-test'
        self.reaction2_filename = os.path.join(
            self.test_subdirectory, 'data', 'te', 'test.pb')
        os.makedirs(os.path.dirname(self.reaction2_filename))
        with open(self.reaction2_filename, 'wb') as f:
            f.write(reaction2.SerializeToString())
        reaction3 = reaction_pb2.Reaction()
        reaction3.provenance.record_id = 'bad-id'
        self.reaction3_filename = os.path.join(
            self.test_subdirectory, 'reaction3.pb')
        with open(self.reaction3_filename, 'wb') as f:
            f.write(reaction3.SerializeToString())

    def test_main(self):
        input_file = os.path.join(self.test_subdirectory, 'input_file.txt')
        with open(input_file, 'w') as f:
            f.write(f'A\t{self.reaction1_filename}\n')
            f.write(f'M\t{self.reaction2_filename}\n')
        output_dir = os.path.join(self.test_subdirectory, 'data')
        with flagsaver.flagsaver(input_file=input_file, root_dir=output_dir):
            update_reactions.main(())
        filenames = glob.glob(os.path.join(output_dir, '*', '*.pb'))
        self.assertLen(filenames, 2)
        record_ids = []
        for filename in filenames:
            message = message_helpers.load_message(
                filename, reaction_pb2.Reaction, input_format='binary')
            self.assertTrue(message.provenance.record_id)
            record_ids.append(message.provenance.record_id)
        self.assertIn('ord-test', record_ids)
        # Only reaction3.pbtxt should be left since cleanup=True.
        self.assertLen(
            glob.glob(os.path.join(self.test_subdirectory, '*.pb')), 1)

    def test_main_bad_id(self):
        input_file = os.path.join(self.test_subdirectory, 'input_file.txt')
        with open(input_file, 'w') as f:
            f.write(f'M\t{self.reaction3_filename}\n')
        output_dir = os.path.join(self.test_subdirectory, 'data')
        with flagsaver.flagsaver(input_file=input_file, root_dir=output_dir):
            with self.assertRaisesRegex(ValueError, 'malformed record_id'):
                update_reactions.main(())

    def test_main_existing_id(self):
        input_file = os.path.join(self.test_subdirectory, 'input_file.txt')
        with open(input_file, 'w') as f:
            f.write(f'A\t{self.reaction2_filename}\n')
        output_dir = os.path.join(self.test_subdirectory, 'data')
        with flagsaver.flagsaver(input_file=input_file, root_dir=output_dir):
            with self.assertRaisesRegex(ValueError, 'record_id is already set'):
                update_reactions.main(())


if __name__ == '__main__':
    absltest.main()
