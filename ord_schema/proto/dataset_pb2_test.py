"""Tests for ord_schema.proto.dataset_pb2."""

from absl.testing import absltest

from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


class DatasetPb2Test(absltest.TestCase):

    def setUp(self):
        super().setUp()
        dataset = dataset_pb2.Dataset()
        dataset.name = 'test'
        dataset.description = 'test dataset'
        # Add a reaction directly to the dataset.
        reaction1 = dataset.reactions.add()
        reaction1.identifiers.add(value='C(C)Cl.Br>>C(C)Br.Cl',
                                  type='REACTION_SMILES')
        # Copy a reaction created elsewhere.
        reaction2 = reaction_pb2.Reaction()
        reaction2.identifiers.add(value='amide coupling', type='NAME')
        dataset.reactions.add().CopyFrom(reaction2)
        # Add an example.
        example = dataset.examples.add()
        example.description = 'test example'
        example.url = 'example.com'
        example.created.time.value = '11 am'
        self.dataset_pb = dataset.SerializeToString()

    def test_simple(self):
        dataset = dataset_pb2.Dataset.FromString(self.dataset_pb)
        self.assertEqual(dataset.name, 'test')
        self.assertEqual(dataset.description, 'test dataset')
        self.assertLen(dataset.reactions, 2)
        self.assertEqual(dataset.reactions[0].identifiers[0].type,
                         reaction_pb2.ReactionIdentifier.REACTION_SMILES)
        self.assertEqual(dataset.reactions[1].identifiers[0].type,
                         reaction_pb2.ReactionIdentifier.NAME)
        self.assertLen(dataset.examples, 1)
        self.assertEqual(dataset.examples[0].description, 'test example')
        self.assertEqual(dataset.examples[0].url, 'example.com')
        self.assertEqual(dataset.examples[0].created.time.value, '11 am')


if __name__ == '__main__':
    absltest.main()
