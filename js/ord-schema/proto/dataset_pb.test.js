/**
 * Copyright 2022 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

const ord_schema = require('..');

test('round-trip', () => {
    const dataset = new ord_schema.dataset_pb.Dataset();
    dataset.setName('test');
    dataset.setDescription('test dataset');
    // Add a reaction directly to the dataset.
    const reaction1 = dataset.addReactions();
    const identifier1 = reaction1.addIdentifiers();
    identifier1.setValue('C(C)Cl.Br>>C(C)Br.Cl');
    identifier1.setType(ord_schema.reaction_pb.ReactionIdentifier.ReactionIdentifierType.REACTION_SMILES);
    // Copy a reaction created elsewhere.
    const reaction2 = new ord_schema.reaction_pb.Reaction();
    const identifier2 = reaction2.addIdentifiers();
    identifier2.setValue('amide coupling');
    identifier2.setType(ord_schema.reaction_pb.ReactionIdentifier.ReactionIdentifierType.NAME);
    dataset.addReactions(reaction2);
    const serialized = dataset.serializeBinary();
    const other = ord_schema.dataset_pb.Dataset.deserializeBinary(serialized);
    expect(other.getName()).toBe('test');
    expect(other.getDescription()).toBe('test dataset');
    expect(other.getReactionsList()).toHaveLength(2);
    expect(other.getReactionsList()[0].getIdentifiersList()[0].getType()).toBe(ord_schema.reaction_pb.ReactionIdentifier.ReactionIdentifierType.REACTION_SMILES);
    expect(other.getReactionsList()[1].getIdentifiersList()[0].getType()).toBe(ord_schema.reaction_pb.ReactionIdentifier.ReactionIdentifierType.NAME);
});
