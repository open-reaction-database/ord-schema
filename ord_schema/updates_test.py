"""Tests for ord_schema.updates."""

from absl.testing import absltest
from rdkit import Chem

from ord_schema import updates
from ord_schema.proto import reaction_pb2


class UpdatesTest(absltest.TestCase):

    def test_resolve_names(self):
        message = reaction_pb2.Reaction()
        message.inputs['test'].components.add().identifiers.add(
            type='NAME', value='aspirin')
        self.assertTrue(updates.resolve_names(message))
        self.assertEqual(
            message.inputs['test'].components[0].identifiers[1],
            reaction_pb2.CompoundIdentifier(type='SMILES',
                                            value='CC(=O)OC1=CC=CC=C1C(=O)O',
                                            details='NAME resolved by PubChem'))

    def test_add_binary_identifiers(self):
        smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        mol = Chem.MolFromSmiles(smiles)
        message = reaction_pb2.Reaction()
        message.inputs['test'].components.add().identifiers.add(
            type='SMILES', value=smiles)
        self.assertTrue(updates.add_binary_identifiers(message))
        self.assertEqual(
            message.inputs['test'].components[0].identifiers[1],
            reaction_pb2.CompoundIdentifier(
                type='RDKIT_BINARY',
                bytes_value=mol.ToBinary(),
                details='Generated from SMILES'))


class UpdateReactionTest(absltest.TestCase):

    def test_with_updates_simple(self):
        message = reaction_pb2.Reaction()
        updates.update_reaction(message)
        self.assertNotEqual(message, reaction_pb2.Reaction())

    def test_with_no_updates(self):
        message = reaction_pb2.Reaction()
        message.provenance.record_created.time.value = '2020-05-08'
        message.provenance.record_id = 'ord-test'
        copied = reaction_pb2.Reaction()
        copied.CopyFrom(message)
        updates.update_reaction(copied)
        self.assertEqual(copied, message)

    def test_with_resolve_names(self):
        reaction = reaction_pb2.Reaction()
        ethylamine = reaction.inputs['ethylamine']
        component = ethylamine.components.add()
        component.identifiers.add(type='NAME', value='ethylamine')
        updates.update_reaction(reaction)
        self.assertLen(component.identifiers, 2)
        self.assertEqual(
            component.identifiers[1],
            reaction_pb2.CompoundIdentifier(
                type='SMILES', value='CCN', details='NAME resolved by PubChem'))


if __name__ == '__main__':
    absltest.main()
