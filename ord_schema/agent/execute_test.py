# Copyright 2026 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for ord_schema.agent.execute."""

import pathlib
from collections.abc import Iterator

import pytest

from ord_schema import parquet, projection, structures
from ord_schema.agent import execute, query
from ord_schema.proto import dataset_pb2, reaction_pb2

_ROLE = reaction_pb2.ReactionRole


def _component(reaction_input, smiles, role):
    component = reaction_input.components.add(reaction_role=role)
    component.identifiers.add(type="SMILES", value=smiles)


def _reaction(reaction_id, *, components, product="Cc1ccccc1"):
    """A reaction with the given (smiles, role) input components."""
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    for smiles, role in components:
        _component(reaction.inputs["in"], smiles, role)
    outcome = reaction.outcomes.add()
    outcome.products.add().identifiers.add(type="SMILES", value=product)
    return reaction


# Two datasets, so that the second's structure ids only join correctly through a
# nonzero offset. Basenames sort "aa" before "bb".
_REACTIONS = {
    "aa": [
        # Pyridine as the solvent; benzene as the reactant.
        _reaction(
            "ord-aa01",
            components=[("c1ccncc1", _ROLE.SOLVENT), ("c1ccccc1", _ROLE.REACTANT)],
        ),
        # Pyridine as a reactant: contains pyridine, but not as a solvent.
        _reaction("ord-aa02", components=[("c1ccncc1", _ROLE.REACTANT)]),
    ],
    "bb": [
        # Ethanol: the only structure with a hydroxyl, and it lives in the second
        # dataset, so finding it proves the offset arithmetic.
        _reaction("ord-bb01", components=[("CCO", _ROLE.REACTANT)]),
    ],
}


@pytest.fixture(scope="module")
def corpus_dir(tmp_path_factory) -> pathlib.Path:
    """Builds sources, projections, and structures artifacts for both datasets."""
    root = tmp_path_factory.mktemp("corpus")
    for shard, reactions in _REACTIONS.items():
        source = root / "data" / f"ord_dataset-{shard}.parquet"
        source.parent.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            dataset_pb2.Dataset(
                dataset_id=f"ord_dataset-{shard}",
                name="test",
                description="test",
                reactions=reactions,
            ),
            str(source),
        )
        projected = root / "projections" / source.name
        projected.parent.mkdir(parents=True, exist_ok=True)
        projection.write_projection(source, projected)
        structured = root / "structures" / source.name
        structured.parent.mkdir(parents=True, exist_ok=True)
        structures.write_structures(projected, structured)
    return root


@pytest.fixture(scope="module")
def corpus(corpus_dir) -> Iterator[execute.Corpus]:
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={"pyridine": "c1ccncc1", "ethanol": "CCO"}.__getitem__,
    ) as value:
        yield value


def _exists(where):
    return {"op": "exists", "path": "inputs.components", "where": where}


def _search(corpus, where) -> set[str]:
    table = corpus.search(query.Query.model_validate({"where": where}))
    return set(table.column("reaction_id").to_pylist())


def test_substructure_binds_to_the_element(corpus):
    # Both aa reactions contain pyridine; only aa01 has it as the solvent. The
    # co-membership is the point: a reaction-granularity intersection would return
    # both.
    matched = _search(
        corpus,
        _exists(
            {
                "op": "and",
                "clauses": [
                    {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
                    {
                        "op": "eq",
                        "path": "reaction_role",
                        "value": {"literal": "SOLVENT"},
                    },
                ],
            }
        ),
    )
    assert matched == {"ord-aa01"}


def test_substructure_reaches_the_offset_dataset(corpus):
    # The hydroxyl is only in bb, whose ids are only correct through its offset.
    matched = _search(
        corpus,
        _exists({"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}),
    )
    assert matched == {"ord-bb01"}


def test_substructure_without_matches_returns_no_rows(corpus):
    matched = _search(
        corpus,
        _exists({"op": "substructure", "path": "smiles", "smarts": "[Pt]"}),
    )
    assert matched == set()


def test_a_compound_name_resolves_to_a_substructure_query(corpus):
    matched = _search(
        corpus,
        _exists({"op": "substructure", "path": "smiles", "compound": "pyridine"}),
    )
    assert matched == {"ord-aa01", "ord-aa02"}


def test_similarity_finds_the_identical_structure(corpus):
    matched = _search(
        corpus,
        _exists(
            {
                "op": "similarity",
                "path": "smiles",
                "smiles": "OCC",  # Ethanol, spelled differently than the corpus.
                "threshold": 0.99,
            }
        ),
    )
    assert matched == {"ord-bb01"}


def test_similarity_at_a_loose_threshold_stays_bounded_by_chemistry(corpus):
    # Pyridine and benzene are similar molecules but not Tanimoto-identical; a loose
    # threshold reaches them both, and ethanol never.
    matched = _search(
        corpus,
        _exists(
            {
                "op": "similarity",
                "path": "smiles",
                "smiles": "c1ccncc1",
                "threshold": 0.3,
            }
        ),
    )
    assert "ord-bb01" not in matched
    assert "ord-aa02" in matched


def test_a_plain_query_still_runs(corpus):
    matched = _search(
        corpus,
        _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}}),
    )
    assert matched == {"ord-bb01"}


def test_an_aggregate_with_a_structure_predicate(corpus):
    table = corpus.search(
        query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"}
                ),
                "aggregate": {"measures": [{"fn": "count", "name": "n"}]},
            }
        )
    )
    assert table.column("n").to_pylist() == [2]


def test_verification_prunes_screen_false_positives(corpus, monkeypatch):
    # Force the screen to pass everything: verification alone must still produce the
    # right answer, which is what makes the screen safe to loosen and never to skip.
    matched_rows = corpus._connection.execute(
        "SELECT count(*) FROM corpus_structures"
    ).fetchone()[0]
    assert matched_rows > 0
    original = execute._pattern_blob

    def _passes_everything(molecule):
        blob, _ = original(molecule)
        return b"\x00" * len(blob), 0

    monkeypatch.setattr(execute, "_pattern_blob", _passes_everything)
    matched = _search(
        corpus,
        _exists({"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}),
    )
    assert matched == {"ord-bb01"}


def test_unpaired_artifacts_are_refused(corpus_dir, tmp_path):
    with pytest.raises(execute.PairingError, match="only one side"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "ord_dataset-aa.parquet"),
        )


def test_a_mismatched_pair_is_refused(corpus_dir, tmp_path):
    # Rebuild bb's projection from a different source; its structures artifact now
    # describes molecules the projection does not hold, so the pair must not open.
    changed = tmp_path / "ord_dataset-bb.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-bb",
            name="test",
            description="changed",
            reactions=[
                _reaction("ord-bb01", components=[("CCN", _ROLE.REACTANT)]),
                _reaction("ord-bb02", components=[("CCO", _ROLE.REACTANT)]),
            ],
        ),
        str(changed),
    )
    projections = tmp_path / "projections"
    projections.mkdir()
    for shard in ("aa", "bb"):
        source = corpus_dir / "projections" / f"ord_dataset-{shard}.parquet"
        (projections / source.name).write_bytes(source.read_bytes())
    projection.write_projection(changed, projections / "ord_dataset-bb.parquet")
    with pytest.raises(execute.PairingError, match="different sources"):
        execute.Corpus(
            str(projections / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
        )


def test_nothing_matched_is_refused(tmp_path):
    with pytest.raises(execute.PairingError, match="no projections"):
        execute.Corpus(
            str(tmp_path / "absent" / "*.parquet"),
            str(tmp_path / "absent" / "*.parquet"),
        )
