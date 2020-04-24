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
        dataset.version = 1
        # Add a reaction directly to the dataset.
        reaction1 = dataset.reactions['foo']
        reaction1.identifiers.add(value='C(C)Cl.Br>>C(C)Br.Cl',
                                  type='REACTION_SMILES')
        # Copy a reaction created elsewhere.
        reaction2 = reaction_pb2.Reaction()
        reaction2.identifiers.add(value='amide coupling', type='NAME')
        dataset.reactions['bar'].CopyFrom(reaction2)
        self.dataset_pb = dataset.SerializeToString()

    def test_simple(self):
        dataset = dataset_pb2.Dataset.FromString(self.dataset_pb)
        self.assertEqual(dataset.name, 'test')
        self.assertEqual(dataset.description, 'test dataset')
        self.assertLen(dataset.reactions, 2)
        self.assertEqual(dataset.reactions['foo'].identifiers[0].type,
                         reaction_pb2.ReactionIdentifier.REACTION_SMILES)
        self.assertEqual(dataset.reactions['bar'].identifiers[0].type,
                         reaction_pb2.ReactionIdentifier.NAME)


if __name__ == '__main__':
    absltest.main()
