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
"""Tests for ord_schema.updates."""

import os
from collections.abc import Iterator

import pytest

from ord_schema import parquet_dataset, updates
from ord_schema.proto import dataset_pb2, reaction_pb2


class TestUpdateDataset:
    @pytest.fixture
    def dataset(self) -> Iterator[dataset_pb2.Dataset]:
        dataset = dataset_pb2.Dataset()
        reaction = dataset.reactions.add()
        # Minimal reaction.
        dummy_input = reaction.inputs["dummy_input"]
        reaction.outcomes.add()
        dummy_component = dummy_input.components.add()
        dummy_component.identifiers.add(type="CUSTOM")
        dummy_component.identifiers[0].details = "custom_identifier"
        dummy_component.identifiers[0].value = "custom_value"
        dummy_component.amount.mass.value = 1
        dummy_component.amount.mass.units = reaction_pb2.Mass.GRAM
        # Placeholders for referenced reactions.
        reaction2 = dataset.reactions.add()
        reaction3 = dataset.reactions.add()
        reaction2.CopyFrom(dataset.reactions[0])
        reaction3.CopyFrom(dataset.reactions[0])
        yield dataset

    def test_crossferences(self, dataset):
        dummy_input = dataset.reactions[0].inputs["dummy_input"]
        dummy_input.crude_components.add(reaction_id="crude", has_derived_amount=True)
        dataset.reactions[1].reaction_id = "crude"
        dummy_input.components[0].preparations.add(type="SYNTHESIZED", reaction_id="synthesized")
        dataset.reactions[2].reaction_id = "synthesized"
        # Check updated values.
        updates.update_dataset(dataset)
        assert dummy_input.crude_components[0].reaction_id == dataset.reactions[1].reaction_id
        assert dummy_input.crude_components[0].reaction_id != "crude"
        assert dummy_input.components[0].preparations[0].reaction_id == dataset.reactions[2].reaction_id
        assert dummy_input.components[0].preparations[0].reaction_id != "synthesized"

    def test_crossferences_are_proper_ids(self, dataset):
        dummy_input = dataset.reactions[0].inputs["dummy_input"]
        dummy_input.crude_components.add(reaction_id="ord-c0bbd41f095a44a78b6221135961d809", has_derived_amount=True)
        dataset.reactions[1].reaction_id = "ord-c0bbd41f095a44a78b6221135961d809"
        dummy_input.components[0].preparations.add(
            type="SYNTHESIZED", reaction_id="ord-d0bbd41f095a44a78b6221135961d809"
        )
        dataset.reactions[2].reaction_id = "ord-d0bbd41f095a44a78b6221135961d809"
        # Check updated values.
        updates.update_dataset(dataset)
        assert dummy_input.crude_components[0].reaction_id == dataset.reactions[1].reaction_id
        assert dummy_input.crude_components[0].reaction_id == "ord-c0bbd41f095a44a78b6221135961d809"
        assert dummy_input.components[0].preparations[0].reaction_id == dataset.reactions[2].reaction_id
        assert dummy_input.components[0].preparations[0].reaction_id == "ord-d0bbd41f095a44a78b6221135961d809"


class TestAssignDatasetId:
    def test_assigns_when_missing(self):
        dataset = dataset_pb2.Dataset()
        result = updates.assign_dataset_id(dataset)
        assert dataset.dataset_id == result
        assert dataset.dataset_id.startswith("ord_dataset-")

    def test_assigns_when_non_canonical(self):
        dataset = dataset_pb2.Dataset(dataset_id="placeholder")
        updates.assign_dataset_id(dataset)
        assert dataset.dataset_id != "placeholder"
        assert dataset.dataset_id.startswith("ord_dataset-")

    def test_keeps_canonical(self):
        dataset = dataset_pb2.Dataset(dataset_id="ord_dataset-c0bbd41f095a44a78b6221135961d809")
        updates.assign_dataset_id(dataset)
        assert dataset.dataset_id == "ord_dataset-c0bbd41f095a44a78b6221135961d809"


class TestAssignIdSubstitutions:
    def test_canonical_ids_get_none(self):
        new_ids, subs = updates.assign_id_substitutions(["ord-c0bbd41f095a44a78b6221135961d809"])
        assert new_ids == [None]
        assert subs == {}

    def test_placeholder_ids_get_substitutions(self):
        new_ids, subs = updates.assign_id_substitutions(["placeholder"])
        assert new_ids[0] is not None and new_ids[0].startswith("ord-")
        assert subs == {"placeholder": new_ids[0]}

    def test_empty_id_gets_new_id_but_no_substitution(self):
        # Empty old_id: nothing else could have referenced it, so no entry
        # in id_substitutions; the reaction still receives a fresh canonical ID.
        new_ids, subs = updates.assign_id_substitutions([""])
        assert new_ids[0] is not None and new_ids[0].startswith("ord-")
        assert subs == {}

    def test_parallel_to_input(self):
        old = ["ord-c0bbd41f095a44a78b6221135961d809", "p1", "", "p2"]
        new_ids, subs = updates.assign_id_substitutions(old)
        assert len(new_ids) == 4
        assert new_ids[0] is None
        assert subs == {"p1": new_ids[1], "p2": new_ids[3]}


class TestApplyReactionUpdates:
    def test_with_new_id(self):
        reaction = reaction_pb2.Reaction()
        modified = updates.apply_reaction_updates(reaction, new_id="ord-c0bbd41f095a44a78b6221135961d809")
        assert modified
        assert reaction.reaction_id == "ord-c0bbd41f095a44a78b6221135961d809"
        assert len(reaction.provenance.record_modified) == 1

    def test_no_new_id_no_other_changes_means_no_modification(self):
        reaction = reaction_pb2.Reaction(reaction_id="ord-c0bbd41f095a44a78b6221135961d809")
        modified = updates.apply_reaction_updates(reaction, new_id=None)
        assert not modified
        assert len(reaction.provenance.record_modified) == 0


class TestApplyCrossReferenceSubstitutions:
    def test_no_substitutions_is_noop(self):
        reaction = reaction_pb2.Reaction()
        reaction.inputs["x"].components.add().preparations.add(type="SYNTHESIZED", reaction_id="orig")
        updates.apply_cross_reference_substitutions(reaction, {})
        assert reaction.inputs["x"].components[0].preparations[0].reaction_id == "orig"

    def test_rewrites_preparations_and_crude_components(self):
        reaction = reaction_pb2.Reaction()
        reaction.inputs["x"].components.add().preparations.add(type="SYNTHESIZED", reaction_id="prep")
        reaction.inputs["x"].crude_components.add(reaction_id="crude", has_derived_amount=True)
        updates.apply_cross_reference_substitutions(reaction, {"prep": "ord-new1", "crude": "ord-new2"})
        assert reaction.inputs["x"].components[0].preparations[0].reaction_id == "ord-new1"
        assert reaction.inputs["x"].crude_components[0].reaction_id == "ord-new2"


class TestUpdateParquetDataset:
    def _make_dataset(self) -> dataset_pb2.Dataset:
        # Three reactions; r2 and r3 are referenced by r1 via placeholder IDs
        # that should be rewritten by the streaming update.
        r1 = reaction_pb2.Reaction(reaction_id="r1")
        comp = r1.inputs["x"].components.add()
        comp.identifiers.add(type="SMILES", value="CC")
        comp.preparations.add(type="SYNTHESIZED", reaction_id="r2")
        r1.inputs["x"].crude_components.add(reaction_id="r3", has_derived_amount=True)
        r2 = reaction_pb2.Reaction(reaction_id="r2")
        r3 = reaction_pb2.Reaction(reaction_id="r3")
        return dataset_pb2.Dataset(name="streaming", description="streaming test", reactions=[r1, r2, r3])

    def test_round_trip_assigns_ids_and_rewrites_cross_refs(self, tmp_path):
        original = self._make_dataset()
        input_path = os.path.join(tmp_path, "in.parquet")
        output_path = os.path.join(tmp_path, "out.parquet")
        parquet_dataset.write_dataset(original, input_path)
        updates.update_parquet_dataset(
            input_path, output_path, dataset_id="ord_dataset-c0bbd41f095a44a78b6221135961d809"
        )
        result = parquet_dataset.read_dataset(output_path)
        assert result.dataset_id == "ord_dataset-c0bbd41f095a44a78b6221135961d809"
        # Each reaction got a canonical ord- ID (the placeholders r1/r2/r3 were
        # all non-canonical, so all three should have been rewritten).
        new_ids = [r.reaction_id for r in result.reactions]
        assert all(rid.startswith("ord-") and len(rid) == 36 for rid in new_ids)
        # Cross-references on r1 should now point at the new IDs of r2 and r3.
        prep = result.reactions[0].inputs["x"].components[0].preparations[0]
        crude = result.reactions[0].inputs["x"].crude_components[0]
        assert prep.reaction_id == new_ids[1]
        assert crude.reaction_id == new_ids[2]
        # Each modified reaction got a record_modified event.
        assert all(len(r.provenance.record_modified) == 1 for r in result.reactions)

    def test_canonical_ids_are_preserved(self, tmp_path):
        canonical = "ord-c0bbd41f095a44a78b6221135961d809"
        reaction = reaction_pb2.Reaction(reaction_id=canonical)
        # An already-canonical reaction with a record_created date so it doesn't
        # trip any other update. Should round-trip unchanged.
        reaction.provenance.record_created.time.value = "2020-01-01"
        original = dataset_pb2.Dataset(name="x", description="x", reactions=[reaction])
        input_path = os.path.join(tmp_path, "in.parquet")
        output_path = os.path.join(tmp_path, "out.parquet")
        parquet_dataset.write_dataset(original, input_path)
        updates.update_parquet_dataset(
            input_path, output_path, dataset_id="ord_dataset-c0bbd41f095a44a78b6221135961d809"
        )
        result = parquet_dataset.read_dataset(output_path)
        assert result.reactions[0].reaction_id == canonical
        assert len(result.reactions[0].provenance.record_modified) == 0

    def test_streaming_matches_in_memory(self, tmp_path):
        # update_parquet_dataset (streaming) and update_dataset (in-memory)
        # should produce the same per-reaction shape modulo the random IDs.
        original = self._make_dataset()
        input_path = os.path.join(tmp_path, "in.parquet")
        output_path = os.path.join(tmp_path, "out.parquet")
        parquet_dataset.write_dataset(original, input_path)
        dataset_id = "ord_dataset-c0bbd41f095a44a78b6221135961d809"
        updates.update_parquet_dataset(input_path, output_path, dataset_id=dataset_id)
        streamed = parquet_dataset.read_dataset(output_path)
        in_memory = self._make_dataset()
        in_memory.dataset_id = dataset_id
        updates.update_dataset(in_memory)
        # Same number of reactions; cross-refs preserved by index.
        assert len(streamed.reactions) == len(in_memory.reactions)
        for s, m in zip(streamed.reactions, in_memory.reactions, strict=True):
            assert (len(s.inputs["x"].components) if "x" in s.inputs else 0) == (
                len(m.inputs["x"].components) if "x" in m.inputs else 0
            )
            assert len(s.provenance.record_modified) == len(m.provenance.record_modified)
