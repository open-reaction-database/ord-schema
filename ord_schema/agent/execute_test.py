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

import contextlib
import logging
import pathlib
import threading
import time
from collections.abc import Iterator

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary

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


# Two datasets, so that the second's structure IDs only join correctly through a
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
def deep_corpus(tmp_path_factory) -> Iterator[execute.Corpus]:
    """A corpus holding a structure at every indexed path, and a reaction at none.

    A workup's components and a product measurement's authentic standard are reached by
    the longest traversals the index generates, and the main fixture has neither, so
    every assertion about them would otherwise compare nothing to nothing. The inputs
    and the products carry structures too, which is what lets a query at one deep path
    assert that a structure living at another is not found there.

    The second reaction has inputs and nothing else, so the deep levels are what the
    projection writes for a level the source never recorded: NULL rather than an empty
    list. A question about those levels has to answer for it too.
    """
    reaction = reaction_pb2.Reaction(reaction_id="ord-cc01")
    _component(reaction.inputs["in"], "CCO", _ROLE.REACTANT)
    _component(reaction.workups.add(type="EXTRACTION").input, "c1ccncc1", _ROLE.SOLVENT)
    outcome = reaction.outcomes.add()
    product = outcome.products.add()
    product.identifiers.add(type="SMILES", value="Cc1ccccc1")
    standard = product.measurements.add(analysis_key="nmr").authentic_standard
    standard.identifiers.add(type="SMILES", value="c1ccccc1")
    shallow = reaction_pb2.Reaction(reaction_id="ord-cc02")
    _component(shallow.inputs["in"], "CCO", _ROLE.REACTANT)
    root = tmp_path_factory.mktemp("deep")
    source = root / "data" / "ord_dataset-cc.parquet"
    source.parent.mkdir(parents=True, exist_ok=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-cc",
            name="test",
            description="test",
            reactions=[reaction, shallow],
        ),
        str(source),
    )
    projected = root / "projections" / source.name
    projected.parent.mkdir(parents=True, exist_ok=True)
    projection.write_projection(source, projected)
    structured = root / "structures" / source.name
    structured.parent.mkdir(parents=True, exist_ok=True)
    structures.write_structures(projected, structured)
    with execute.Corpus(
        str(projected), str(structured), resolver={}.__getitem__
    ) as value:
        yield value


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


def _reactions(table) -> set[str]:
    return set(table.column("reaction_id").to_pylist())


def _search(corpus, where) -> set[str]:
    return _reactions(corpus.search(query.Query.model_validate({"where": where})))


def _role_and_structure(smarts, role):
    return _exists(
        {
            "op": "and",
            "clauses": [
                {"op": "substructure", "path": "smiles", "smarts": smarts},
                {"op": "eq", "path": "reaction_role", "value": {"literal": role}},
            ],
        }
    )


def test_substructure_binds_to_the_element(corpus):
    # aa01 holds benzene as a REACTANT and pyridine as its SOLVENT. Bound, that is one
    # component required to be both, which nothing satisfies. The unbound reading --
    # "contains benzene" intersected with "contains a solvent" -- would return aa01, so
    # this is the assertion that fails if element binding is ever lost.
    assert _search(corpus, _role_and_structure("c1ccccc1", "SOLVENT")) == set()
    # The same query against the component that really is the solvent does match, so
    # the empty set above is binding at work rather than a predicate that never fires.
    assert _search(corpus, _role_and_structure("c1ccncc1", "SOLVENT")) == {"ord-aa01"}


def test_substructure_reaches_the_offset_dataset(corpus):
    # The hydroxyl is only in bb, whose IDs are only correct through its offset.
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


def test_the_screen_only_accelerates_the_answer(corpus):
    # The fingerprint screen is an accelerator, never part of the answer: a library
    # built with no screen at all has to match exactly the same structures. That is
    # the property which makes the screen safe to loosen and never safe to trust.
    screened = corpus._library()
    unscreened = rdSubstructLibrary.SubstructLibrary(
        rdSubstructLibrary.CachedMolHolder()
    )
    for index in range(len(screened)):
        unscreened.AddMol(screened.GetMol(index))
    for smarts in ("[OX2H]", "c1ccncc1", "c1ccccc1"):
        pattern = Chem.MolFromSmarts(smarts)
        with_screen = set(screened.GetMatches(pattern, maxResults=len(screened)))
        without_screen = set(unscreened.GetMatches(pattern, maxResults=len(unscreened)))
        assert with_screen == without_screen, smarts


def test_an_unparseable_structure_holds_its_slot(corpus_dir, tmp_path):
    # A structure RDKit could not read still occupies its ID. Skipping it instead would
    # shift every later ID down one, and the bitmap would then name molecules that are
    # not the ones matched -- wrong chemistry, corpus-wide, and silent. The nulled row
    # is aa's first structure, so every remaining structure in the corpus sits after it
    # and answers correctly only if the empty slot is really there.
    root = _copy_corpus(corpus_dir, tmp_path)
    target = root / "structures" / "ord_dataset-aa.parquet"
    with pq.ParquetFile(target) as artifact:
        table = artifact.read()
        schema = artifact.schema_arrow
    assert table.column("smiles")[0].as_py() == "c1ccncc1", "expected pyridine first"
    columns = {field.name: table.column(field.name) for field in schema}
    for name in ("mol_binary", "pattern_fp"):
        values = table.column(name).to_pylist()
        values[0] = None
        columns[name] = pa.array(values, type=schema.field(name).type)
    pq.write_table(pa.table(columns, schema=schema), target)
    with execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        # The library holds distinct molecules, so its length is not the corpus's. What
        # has to hold is that every structure ID resolves to exactly one entry.
        library = value._library()
        assert len(value._starts) - 1 == len(library)
        assert len(value._members) == value._total
        assert sorted(value._members) == list(range(value._total))
        # Benzene sits after the hole in aa, ethanol after it in bb. Both still name
        # their own reactions, which holds only while the empty slot is there.
        assert _search(
            value,
            _exists({"op": "substructure", "path": "smiles", "smarts": "c1ccccc1"}),
        ) == {"ord-aa01"}
        assert _search(
            value, _exists({"op": "substructure", "path": "smiles", "smarts": "[OX2H]"})
        ) == {"ord-bb01"}
        # The unreadable structure matches nothing at all, rather than matching whatever
        # molecule happens to land on its index.
        assert (
            _search(
                value,
                _exists({"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"}),
            )
            == set()
        )


def test_unpaired_artifacts_are_refused(corpus_dir, tmp_path):
    with pytest.raises(execute.PairingError, match="no counterpart"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "ord_dataset-aa.parquet"),
        )


def test_a_mismatched_pair_is_refused(corpus_dir, tmp_path):
    # Rebuild bb's projection from a different source, leaving its structures
    # artifact describing molecules the projection does not hold: the pair must not
    # open, or IDs would join to the wrong chemistry.
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
    with pytest.raises(execute.PairingError, match="no counterpart"):
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


def _copy_corpus(corpus_dir, tmp_path) -> pathlib.Path:
    """Copies the built corpus so a test can damage one artifact in isolation."""
    for name in ("projections", "structures"):
        (tmp_path / name).mkdir()
        for source in (corpus_dir / name).glob("*.parquet"):
            (tmp_path / name / source.name).write_bytes(source.read_bytes())
    return tmp_path


def test_offsets_place_each_dataset_in_its_own_slice(corpus):
    # Toluene is the product of every reaction and lives in BOTH datasets, at local ID
    # 2 in aa and 1 in bb. Reaching all three reactions therefore requires each ID to
    # be offset by its own file's base: swapping the two offsets, or shifting both,
    # produces a different answer.
    matched = corpus.search(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "outcomes.products",
                    "where": {
                        "op": "substructure",
                        "path": "smiles",
                        "smarts": "Cc1ccccc1",
                    },
                }
            }
        )
    )
    assert set(matched.column("reaction_id").to_pylist()) == {
        "ord-aa01",
        "ord-aa02",
        "ord-bb01",
    }


def test_a_structures_artifact_short_of_its_projection_is_refused(corpus_dir, tmp_path):
    # The dangerous desynchronization: an artifact rederived from a rewritten
    # projection keeps the same source hash, so the stamps still agree. A short one
    # would pair cleanly and then alias the next dataset's molecules -- in range for
    # get_bit, wrong about the chemistry, and silent.
    root = _copy_corpus(corpus_dir, tmp_path)
    target = root / "structures" / "ord_dataset-aa.parquet"
    with pq.ParquetFile(target) as artifact:
        table = artifact.read()
        schema = artifact.schema_arrow
    pq.write_table(table.slice(0, table.num_rows - 1).cast(schema), target)
    with pytest.raises(execute.PairingError, match="would join to another dataset"):
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
        )


def _rewrite_structures(root, shard, **columns):
    """Replaces whole columns of one structures artifact, keeping its schema."""
    target = root / "structures" / f"ord_dataset-{shard}.parquet"
    with pq.ParquetFile(target) as artifact:
        table = artifact.read()
        schema = artifact.schema_arrow
    replaced = {field.name: table.column(field.name) for field in schema}
    for name, values in columns.items():
        replaced[name] = pa.array(values, type=schema.field(name).type)
    pq.write_table(pa.table(replaced, schema=schema), target)


@pytest.mark.parametrize(
    ("label", "identifiers"),
    [("a duplicate", [0, 1, 1]), ("a gap", [0, 1, 3])],
)
def test_structure_ids_that_are_not_one_unbroken_run_are_refused(
    corpus_dir, tmp_path, label, identifiers
):
    # A library index is a structure ID only while the IDs are every integer from zero,
    # each exactly once. A duplicate or a gap slides later structures onto a neighbor's
    # index, and the row count -- unchanged by either -- cannot see it.
    root = _copy_corpus(corpus_dir, tmp_path)
    _rewrite_structures(root, "aa", structure_id=identifiers)
    with (
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        ) as value,
        pytest.raises(execute.PairingError, match="not one unbroken run"),
    ):
        value._library()


@pytest.mark.parametrize("missing", ["mol_binary", "pattern_fp"])
def test_half_a_derived_structure_is_refused(corpus_dir, tmp_path, missing):
    # The artifact writes a molecule and its fingerprint together or not at all, and
    # the library relies on that. Left to itself, one without the other goes two
    # different wrong ways: a missing molecule discards a live fingerprint, so the
    # structure quietly matches nothing while still counting as searchable, and a
    # missing fingerprint reaches RDKit as a type error naming no row at all.
    root = _copy_corpus(corpus_dir, tmp_path)
    target = root / "structures" / "ord_dataset-bb.parquet"
    with pq.ParquetFile(target) as artifact:
        values = artifact.read().column(missing).to_pylist()
    values[0] = None
    _rewrite_structures(root, "bb", **{missing: values})
    with (
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        ) as value,
        pytest.raises(execute.PairingError, match="writes both or neither"),
    ):
        value._library()


def test_two_artifacts_of_one_dataset_are_refused(corpus_dir, tmp_path):
    # Which of them answers a query would be arbitrary, so neither may.
    root = _copy_corpus(corpus_dir, tmp_path)
    source = root / "projections" / "ord_dataset-aa.parquet"
    (root / "projections" / "copy.parquet").write_bytes(source.read_bytes())
    with pytest.raises(execute.PairingError, match="same source dataset"):
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
        )


def test_artifacts_pair_across_differing_directory_layouts(corpus_dir, tmp_path):
    # Pairing is by source dataset, not by filename, so a corpus whose two trees are
    # laid out differently -- or whose files share a basename across directories --
    # still pairs. Keying on the basename would drop a dataset here without a word.
    root = _copy_corpus(corpus_dir, tmp_path)
    for shard, name in (
        ("one", "ord_dataset-aa.parquet"),
        ("two", "ord_dataset-bb.parquet"),
    ):
        (root / "sharded" / shard).mkdir(parents=True)
        (root / "sharded" / shard / "data.parquet").write_bytes(
            (root / "projections" / name).read_bytes()
        )
    with execute.Corpus(
        str(root / "sharded" / "*" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        matched = value.search(
            query.Query.model_validate(
                {
                    "where": _exists(
                        {"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}
                    )
                }
            )
        )
    assert matched.column("reaction_id").to_pylist() == ["ord-bb01"]


def test_a_stale_artifact_is_refused(corpus_dir, tmp_path):
    root = _copy_corpus(corpus_dir, tmp_path)
    target = root / "structures" / "ord_dataset-aa.parquet"
    with pq.ParquetFile(target) as artifact:
        table = artifact.read()
        schema = artifact.schema_arrow
    metadata = dict(schema.metadata)
    metadata[b"ord.rdkit_version"] = b"0000.00.0"
    pq.write_table(table.cast(schema.with_metadata(metadata)), target)
    with pytest.raises(execute.PairingError, match="stale"):
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
        )
    # The same corpus opens when the caller takes responsibility for the mismatch.
    with execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        require_current=False,
        resolver={}.__getitem__,
    ) as value:
        assert value.search(query.Query.model_validate({})).num_rows == 3


def test_similarity_bounds_the_match_set_at_the_threshold(corpus):
    # Tanimoto against pyridine in this fixture: benzene 0.3333, toluene 0.1765,
    # ethanol 0.0. Exact sets on either side of benzene pin the popcount band and the
    # comparison in both directions -- a band that dropped smaller candidates, or a
    # strict >, would change one of these.
    def at(threshold):
        return _search(
            corpus,
            _exists(
                {
                    "op": "similarity",
                    "path": "smiles",
                    "smiles": "c1ccncc1",
                    "threshold": threshold,
                }
            ),
        )

    assert at(0.3) == {"ord-aa01", "ord-aa02"}  # Reaches benzene, in aa01.
    assert at(0.34) == {"ord-aa01", "ord-aa02"}  # Benzene excluded; pyridine remains.
    assert at(1.0) == {"ord-aa01", "ord-aa02"}  # Identity: >= must not become >.


def test_similarity_at_one_is_an_exact_structure_query(corpus):
    # threshold 1.0 is legal and is the exact-structure query; a strict comparison
    # would return nothing for it.
    matched = _search(
        corpus,
        _exists(
            {
                "op": "similarity",
                "path": "smiles",
                "smiles": "OCC",
                "threshold": 1.0,
            }
        ),
    )
    assert matched == {"ord-bb01"}


def test_an_explicit_hydrogen_smarts_runs_after_the_rewrite(corpus):
    # The rewritten pattern has to survive the screen and the verification, not just
    # RDKit: [H]OC matches no stored molecule, [O&!H0]C matches ethanol's hydroxyl.
    with pytest.warns(UserWarning, match="rewritten"):
        request = query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "smarts": "[H]OC"}
                )
            }
        )
    matched = set(corpus.search(request).column("reaction_id").to_pylist())
    assert matched == {"ord-bb01"}


def test_forall_and_not_compose_with_a_structure_predicate(corpus):
    pyridine = {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"}
    # aa02's only input is pyridine, so it is the one reaction where every component
    # satisfies the predicate.
    every = corpus.search(
        query.Query.model_validate(
            {"where": {"op": "forall", "path": "inputs.components", "where": pyridine}}
        )
    )
    assert "ord-aa02" in set(every.column("reaction_id").to_pylist())
    negated = _search(corpus, _exists({"op": "not", "clause": pyridine}))
    assert negated == {"ord-aa01", "ord-bb01"}


def test_matching_runs_across_threads(corpus_dir):
    # RDKit matches across its own threads with the GIL released; the answer must not
    # depend on how many it uses.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        threads=2,
        resolver={}.__getitem__,
    ) as value:
        matched = value.search(
            query.Query.model_validate(
                {
                    "where": _exists(
                        {"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}
                    )
                }
            )
        )
    assert matched.column("reaction_id").to_pylist() == ["ord-bb01"]


def test_a_timeout_does_not_disturb_a_later_search(corpus):
    # Timer.cancel only sets a flag, so a timer past its own check fires anyway. An
    # interrupt outliving the search that armed it would fail a later query and report
    # it as having timed out.
    request = query.Query.model_validate(
        {"where": _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}})}
    )
    # A timeout this short may or may not fire; either outcome is fine for the query
    # that armed it. What must never happen is the interrupt landing on a later query.
    for _ in range(20):
        with contextlib.suppress(TimeoutError):
            corpus.search(request, timeout_seconds=0.001)
        assert _search(
            corpus, _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}})
        ) == {"ord-bb01"}


def test_a_timeout_interrupts_a_query_that_outlasts_it():
    # Nothing in the fixture corpus is slow enough to time out, so the interrupt itself
    # is pinned here against a query that is. Were it ever to stop reaching the running
    # statement, every timeout would quietly become no timeout at all.
    connection = duckdb.connect()
    try:
        cursor = connection.cursor()
        started = time.perf_counter()
        with pytest.raises(TimeoutError, match=r"exceeded 0\.2 seconds"):
            execute._run_with_timeout(
                cursor, "SELECT count(*) FROM range(100000000000)", {}, 0.2
            )
        # The bound is what ends it, rather than the query finishing on its own.
        assert time.perf_counter() - started < 30
    finally:
        connection.close()


def test_an_interrupt_reaches_only_the_cursor_it_was_called_on():
    # Timeouts are per search because interrupt() is per cursor. That is DuckDB's
    # behavior rather than this module's, and the whole timeout design rests on it: were
    # an interrupt ever to reach the connection's other cursors, one search timing out
    # would abort every search in flight and report each as having timed out. Pinned
    # here so the change is caught on a version bump instead of in production.
    connection = duckdb.connect()
    try:
        victim, other = connection.cursor(), connection.cursor()
        slow = "SELECT count(*) FROM range(30000000000)"
        # The query has to be long enough that an interrupt lands mid-flight, or the
        # assertions below would hold for a query that simply finished first.
        timer = threading.Timer(0.2, victim.interrupt)
        timer.start()
        try:
            with pytest.raises(duckdb.InterruptException):
                victim.execute(slow).fetchall()
        finally:
            timer.cancel()
        # Same query, same duration, interrupted through a sibling cursor instead.
        outcome: list[str] = []

        def run():
            try:
                other.execute(slow).fetchall()
                outcome.append("completed")
            except duckdb.InterruptException:
                outcome.append("interrupted")

        thread = threading.Thread(target=run)
        thread.start()
        time.sleep(0.2)  # Well inside the query, as the self-interrupt above showed.
        victim.interrupt()
        connection.interrupt()
        time.sleep(0.3)  # Longer than the self-interrupt took to land.
        # Neither reached it: the query is still running. Asserting that it is alive is
        # the whole test, since interrupting it below would end it either way.
        still_running = thread.is_alive()
        other.interrupt()  # Ends the query, so the test does not run for minutes.
        thread.join(timeout=_WAIT)
    finally:
        connection.close()
    assert still_running, "a sibling cursor or the connection interrupted the query"
    assert outcome == ["interrupted"]  # Its own cursor did stop it.


def test_a_timeout_does_not_disturb_a_concurrent_search(corpus_dir):
    # A search arming a timer must interrupt only itself. Interrupting anything wider
    # would abort whatever else is in flight and report it, wrongly, as a timeout.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value._library()
        timed = query.Query.model_validate(
            {"where": _exists({"op": "substructure", "path": "smiles", "smarts": "*"})}
        )
        steady = query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "eq", "path": "smiles", "value": {"literal": "CCO"}}
                )
            }
        )
        failures: list[BaseException] = []
        wrong: list[list[str]] = []
        stop = threading.Event()

        def timing_out():
            try:
                while not stop.is_set():
                    with contextlib.suppress(TimeoutError):
                        value.search(timed, timeout_seconds=0.001)
            except BaseException as error:  # noqa: BLE001
                failures.append(error)

        armer = threading.Thread(target=timing_out)
        armer.start()
        try:
            for _ in range(40):
                matched = value.search(steady).column("reaction_id").to_pylist()
                if matched != ["ord-bb01"]:
                    wrong.append(matched)
        finally:
            stop.set()
            armer.join(timeout=_WAIT)
    assert failures == []
    assert wrong == []


def test_a_generous_timeout_returns_the_answer(corpus):
    matched = corpus.search(
        query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}
                )
            }
        ),
        timeout_seconds=60,
    )
    assert matched.column("reaction_id").to_pylist() == ["ord-bb01"]


# Long enough that a loaded machine does not trip it, short enough that a regression
# reports as a failure instead of parking the suite forever. Nothing here waits on real
# work: every barrier is released by threads that have already been started.
_WAIT = 60


class _CursorCounter:
    """A connection that counts cursors and refuses to be queried directly.

    Forwards everything else, so a search that takes a cursor works and one that runs
    its own query on the shared connection fails.
    """

    def __init__(self, connection):
        self.connection = connection
        self.cursors = 0

    def cursor(self):
        self.cursors += 1
        return self.connection.cursor()

    def execute(self, *args, **kwargs):
        raise AssertionError("the search queried the shared connection")

    def __getattr__(self, name):
        return getattr(self.connection, name)


def test_a_search_runs_on_a_cursor_rather_than_the_shared_connection(corpus_dir):
    # A DuckDBPyConnection carries the pending result of its last execute, so searches
    # sharing one read each other's rows instead of raising. Counting cursors alone
    # would pass an implementation that took one and then queried the connection
    # anyway, so the connection is made unusable for queries: the search has to run on
    # what cursor() handed it.
    request = query.Query.model_validate(
        {"where": _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}})}
    )
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value._connection = counter = _CursorCounter(value._connection)
        # The first search also materializes its columns, which is a cursor of its own.
        assert value.search(request).num_rows == 1
        # Steady state, with that table already built: one search, one cursor.
        taken = counter.cursors
        assert value.search(request).num_rows == 1
        assert counter.cursors == taken + 1


def test_concurrent_searches_return_their_own_answers(corpus_dir):
    # Distinct queries with distinct answers, run at once against one corpus: each
    # thread has to come back with the reactions its own query selected. Sharing one
    # connection also makes searches raise "No open result set" as often as it makes
    # them lie, so the exceptions are collected too -- a thread that dies takes its
    # traceback to threading.excepthook, where an assertion on results alone sees
    # nothing wrong.
    # The similarity thread covers the one path that runs two statements on a single
    # cursor -- its screen, then the compiled query -- so reuse within one search shows
    # up here and nowhere else.
    expected = {
        "hydroxyl": ["ord-bb01"],
        "pyridine": ["ord-aa01", "ord-aa02"],
        "benzene": ["ord-aa01"],
        "like ethanol": ["ord-bb01"],
    }
    predicates = {
        "hydroxyl": {"op": "substructure", "path": "smiles", "smarts": "[OX2H]"},
        "pyridine": {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
        "benzene": {"op": "substructure", "path": "smiles", "smarts": "c1ccccc1"},
        "like ethanol": {
            "op": "similarity",
            "path": "smiles",
            "smiles": "OCC",
            "threshold": 0.99,
        },
    }
    rounds = 20
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value._library()  # Build once up front, so the threads race on searching.
        # Built before any thread starts: validation that raised inside a thread would
        # leave the others parked on a barrier nobody else reaches.
        requests = {
            label: query.Query.model_validate({"where": _exists(predicate)})
            for label, predicate in predicates.items()
        }
        wrong: list[tuple[str, list[str]]] = []
        failures: list[BaseException] = []
        completed: list[str] = []
        ready = threading.Barrier(len(expected))

        def search(label):
            try:
                ready.wait(timeout=_WAIT)
                for _ in range(rounds):
                    matched = sorted(
                        value.search(requests[label]).column("reaction_id").to_pylist()
                    )
                    if matched != expected[label]:
                        wrong.append((label, matched))
                completed.append(label)
            except BaseException as error:  # noqa: BLE001
                failures.append(error)

        threads = [threading.Thread(target=search, args=(label,)) for label in expected]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=_WAIT)
        alive = [thread for thread in threads if thread.is_alive()]
    assert failures == []
    assert alive == []
    assert wrong == []
    assert sorted(completed) == sorted(expected)


def test_concurrent_first_searches_build_one_library(corpus_dir, monkeypatch):
    # The library is about 1.5 GB at corpus scale. First searches arriving together must
    # not each build their own copy, so the build is serialized and whoever waits takes
    # the finished one.
    threads_count = 4
    builds = 0
    counting = threading.Lock()
    original = execute.rdSubstructLibrary.SubstructLibrary

    def counted(*args, **kwargs):
        nonlocal builds
        with counting:
            builds += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(execute.rdSubstructLibrary, "SubstructLibrary", counted)
    ready = threading.Barrier(threads_count)

    def resolver(name):
        # Every thread is inside search and about to want the library; releasing them
        # together is what makes the race reachable at all. This assumes each search
        # resolves the name itself -- were that cache ever hoisted to the corpus, the
        # threads that skipped it would leave this barrier one short, so it times out
        # rather than parking the suite.
        ready.wait(timeout=_WAIT)
        return {"pyridine": "c1ccncc1"}[name]

    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver=resolver,
    ) as value:
        request = query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "compound": "pyridine"}
                )
            }
        )
        results: list[list[str]] = []
        failures: list[BaseException] = []

        def search():
            try:
                matched = value.search(request).column("reaction_id").to_pylist()
                results.append(sorted(matched))
            except BaseException as error:  # noqa: BLE001
                failures.append(error)

        threads = [threading.Thread(target=search) for _ in range(threads_count)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=_WAIT)
        alive = [thread for thread in threads if thread.is_alive()]
    assert failures == []
    assert alive == []
    assert builds == 1
    assert results == [["ord-aa01", "ord-aa02"]] * threads_count


_SUBSTRUCTURE = {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"}
_SOLVENT = {"op": "eq", "path": "reaction_role", "value": {"literal": "SOLVENT"}}


def _index_spent(body) -> bool:
    """Returns whether the occurrence index takes any clause of this query."""
    spent = False

    def index(path, fields, allocate):
        nonlocal spent
        condition = execute._index_condition(path, fields, allocate)
        spent = spent or condition is not None
        return condition

    query.compile_query(query.Query.model_validate(body), index=index)
    return spent


def _no_index_condition(path, fields, allocate):
    """Stands in for _index_condition so every quantifier compiles over the elements."""
    del path, fields, allocate  # Unused.


# Every shape where the index takes a clause. The index answers one quantifier, not the
# query, so everything around that quantifier -- aggregates, orderings, limits, scalars,
# negation, disjunction, a second quantifier -- composes and belongs here.
_ROUTABLE = {
    "structure alone": {"where": _exists(_SUBSTRUCTURE)},
    "structure and role": {
        "where": _exists({"op": "and", "clauses": [_SUBSTRUCTURE, _SOLVENT]})
    },
    "role and structure, reversed": {
        "where": _exists({"op": "and", "clauses": [_SOLVENT, _SUBSTRUCTURE]})
    },
    "similarity": {
        "where": _exists(
            {"op": "similarity", "path": "smiles", "smiles": "OCC", "threshold": 0.4}
        )
    },
    "a compound name": {
        "where": _exists(
            {"op": "substructure", "path": "smiles", "compound": "pyridine"}
        )
    },
    "products rather than inputs": {
        "where": {
            "op": "exists",
            "path": "outcomes.products",
            "where": {"op": "substructure", "path": "smiles", "smarts": "Cc1ccccc1"},
        }
    },
    "matches nothing": {
        "where": _exists({"op": "substructure", "path": "smiles", "smarts": "[Pt]"})
    },
    # Toluene is every reaction's product and nobody's input, so asking for it among
    # the inputs must come back empty. An index that lost the path it was asked about
    # would answer with the products and return all three reactions.
    "a structure that lives only at another path": {
        "where": _exists(
            {"op": "substructure", "path": "smiles", "smarts": "Cc1ccccc1"}
        )
    },
    # One aromatic carbon matches pyridine and benzene both, which are aa01's two
    # inputs: one reaction reached through two occurrence rows -- and a semi-join names
    # the reaction once however many rows carry it.
    "several components of one reaction": {
        "where": _exists({"op": "substructure", "path": "smiles", "smarts": "c"})
    },
    # The shapes below wrap the indexed quantifier in the rest of the grammar, which is
    # the point of answering a clause rather than the query.
    "an aggregate over the quantifier": {
        "where": _exists(_SUBSTRUCTURE),
        "aggregate": {"measures": [{"fn": "count", "name": "n"}]},
    },
    "an ordering": {
        "where": _exists(_SUBSTRUCTURE),
        "order_by": [{"key": "reaction_id"}],
    },
    "an ordering and a limit": {
        "where": _exists(_SUBSTRUCTURE),
        "order_by": [{"key": "reaction_id"}],
        "limit": 1,
    },
    "a scalar beside the quantifier": {
        "where": {
            "op": "and",
            "clauses": [
                _exists(_SUBSTRUCTURE),
                {"op": "not_null", "path": "reaction_id"},
            ],
        }
    },
    # No element holds a match, which the semi-join states as absent membership.
    "a negation of the quantifier": {
        "where": {"op": "not", "clause": _exists(_SUBSTRUCTURE)}
    },
    "a disjunction of two quantifiers": {
        "where": {
            "op": "or",
            "clauses": [
                _exists(_SUBSTRUCTURE),
                _exists({"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}),
            ],
        }
    },
    # Two quantifiers are two semi-joins, each with its own bitmap.
    "two structure predicates": {
        "where": {
            "op": "and",
            "clauses": [
                _exists({"op": "and", "clauses": [_SUBSTRUCTURE, _SOLVENT]}),
                _exists({"op": "substructure", "path": "smiles", "smarts": "c1ccccc1"}),
            ],
        }
    },
}

# Element predicates an occurrence row cannot carry, so the quantifier compiles over
# the elements and the projection answers it -- whatever surrounds the quantifier.
_NOT_ROUTABLE = {
    # No structure means no bitmap, and the index is only worth its scan with one.
    "a role alone": {"where": _exists(_SOLVENT)},
    # amount lives on the element but not in the index.
    "an unindexed element field": {
        "where": _exists(
            {
                "op": "and",
                "clauses": [
                    _SUBSTRUCTURE,
                    {"op": "gt", "path": "amount.mass_grams", "value": {"literal": 1}},
                ],
            }
        )
    },
    # An equality on an element field the index does not carry. Routed, the condition
    # would be dropped and the answer would be the unconditioned one.
    "an equality on another element field": {
        "where": _exists(
            {
                "op": "and",
                "clauses": [
                    _SUBSTRUCTURE,
                    {
                        "op": "eq",
                        "path": "smiles",
                        "value": {"literal": "c1ccncc1"},
                    },
                ],
            }
        )
    },
    # The index carries the role a row has, never the roles it does not have, so an
    # inequality routed as an equality would invert the answer.
    "a role inequality": {
        "where": _exists(
            {
                "op": "and",
                "clauses": [
                    _SUBSTRUCTURE,
                    {
                        "op": "ne",
                        "path": "reaction_role",
                        "value": {"literal": "SOLVENT"},
                    },
                ],
            }
        )
    },
    # A role naming a compound binds a parameter, and the index SQL binds none: it
    # carries the structure parameter and nothing else.
    "a role given as a compound": {
        "where": _exists(
            {
                "op": "and",
                "clauses": [
                    _SUBSTRUCTURE,
                    {
                        "op": "eq",
                        "path": "reaction_role",
                        "value": {"compound": "pyridine"},
                    },
                ],
            }
        )
    },
    # Two roles on one element is a contradiction the projection answers with no rows.
    # An index keeping only the last would answer with the rows holding that one.
    "two roles at once": {
        "where": _exists(
            {
                "op": "and",
                "clauses": [
                    _SUBSTRUCTURE,
                    _SOLVENT,
                    {
                        "op": "eq",
                        "path": "reaction_role",
                        "value": {"literal": "REACTANT"},
                    },
                ],
            }
        )
    },
    # Two structure predicates on one element are two bitmaps, and one occurrence row
    # carries one structure: keeping either alone would drop the other's condition.
    "two structure predicates on one element": {
        "where": _exists(
            {
                "op": "and",
                "clauses": [
                    _SUBSTRUCTURE,
                    {"op": "substructure", "path": "smiles", "smarts": "c1ccccc1"},
                ],
            }
        )
    },
    # A negation inside the element is about what a row does not say, and a disjunction
    # there needs the element's own fields; both are element logic, not row membership.
    "a negation inside the quantifier": {
        "where": _exists({"op": "not", "clause": _SUBSTRUCTURE})
    },
    "a disjunction inside the quantifier": {
        "where": _exists({"op": "or", "clauses": [_SUBSTRUCTURE, _SOLVENT]})
    },
    # forall is not exists: the index shows the elements that match, never that every
    # element does.
    "a universal": {
        "where": {"op": "forall", "path": "inputs.components", "where": _SUBSTRUCTURE}
    },
    # isolated_color is an element field no occurrence row carries, like any of the
    # dozens beside the one indexed field.
    "an equality on an element field of the products": {
        "where": {
            "op": "exists",
            "path": "outcomes.products",
            "where": {
                "op": "and",
                "clauses": [
                    {"op": "substructure", "path": "smiles", "smarts": "Cc1ccccc1"},
                    {
                        "op": "eq",
                        "path": "isolated_color",
                        "value": {"literal": "white"},
                    },
                ],
            },
        }
    },
}


@pytest.mark.parametrize("label", sorted(_ROUTABLE))
def test_the_index_answers_exactly_what_the_projection_would(corpus, label):
    # The index is a second way to reach the same reactions, so the only thing that
    # makes it safe is that it never reaches different ones. Each query runs twice
    # through search itself -- once with the index, once with the hook answering None
    # everywhere, which is a corpus that never built one -- so both sides are the
    # production path. An ordered query is compared in order, since there the order is
    # the answer; the rest as multisets, since neither relation promises one.
    body = _ROUTABLE[label]
    assert _index_spent(body), "expected the index to take a clause"
    request = query.Query.model_validate(body)
    with_index = corpus.search(request).to_pylist()
    with pytest.MonkeyPatch.context() as patcher:
        patcher.setattr(execute, "_index_condition", _no_index_condition)
        direct = corpus.search(request).to_pylist()
    if request.order_by:
        assert list(map(repr, with_index)) == list(map(repr, direct)), label
    else:
        assert sorted(map(repr, with_index)) == sorted(map(repr, direct)), label


@pytest.mark.parametrize("label", sorted(_NOT_ROUTABLE))
def test_the_index_declines_what_an_occurrence_row_cannot_carry(label):
    # Taking one of these clauses would answer a different question than was asked, so
    # the hook has to decline rather than approximate; the quantifier then compiles over
    # the elements and the projection answers it.
    assert not _index_spent(_NOT_ROUTABLE[label])


@pytest.mark.parametrize(
    ("label", "expected"),
    [
        # A pyridine that is not the solvent, which is aa02's reactant. Routed as an
        # equality, this would come back as aa01 -- the one reaction it excludes.
        ("a role inequality", {"ord-aa02"}),
        # One component required to be two roles at once, which nothing is. Routed with
        # only the last role kept, this would come back as every pyridine reactant.
        ("two roles at once", set()),
        # Declined by the executor rather than by the shape, so the search exercises
        # the decline that happens after the terms were read.
        ("an equality on another element field", {"ord-aa01", "ord-aa02"}),
        # No molecule is a pyridine ring and a benzene ring at once. Kept as one of the
        # two clauses -- either one -- this would return that clause's matches instead.
        ("two structure predicates on one element", set()),
    ],
)
def test_a_declined_query_is_answered_rather_than_only_declined(
    corpus, label, expected
):
    # The planner declining is half the property; the other half is that the projection
    # answers, and answers what the comments beside these cases say routing would get
    # wrong. Asserted for the two whose routed answer would not merely differ but invert
    # or over-return.
    assert (
        _reactions(corpus.search(query.Query.model_validate(_NOT_ROUTABLE[label])))
        == expected
    )


def test_the_index_names_each_reaction_once(corpus):
    # A reaction holding two matching components is two rows in the index and one row in
    # the projection. A caller gets the table rather than a set, so a duplicate would
    # reach them; this asserts none does, which comparing the two paths as sets cannot.
    matched = corpus.search(
        query.Query.model_validate(
            {
                "where": _exists(
                    # Both of aa01's inputs hold aromatic carbon, so both are rows.
                    {"op": "substructure", "path": "smiles", "smarts": "c"}
                )
            }
        )
    ).column("reaction_id")
    assert len(matched) == len(set(matched.to_pylist()))


def test_a_limited_query_takes_some_of_the_match_set(corpus):
    # A limit with no ordering asks for some rows rather than particular ones, so
    # equality against the unindexed answer is not the property; the count and the
    # membership are. The semi-join leaves the limit on the outer query, where a dropped
    # one would return the whole set and the count catches it.
    body = {"where": _exists(_SUBSTRUCTURE), "limit": 1}
    assert _index_spent(body), "expected the index to take the clause"
    limited = corpus.search(query.Query.model_validate(body))
    matched = limited.column("reaction_id").to_pylist()
    assert len(matched) == 1
    assert set(matched) <= {"ord-aa01", "ord-aa02"}


def test_a_role_literal_cannot_close_its_quote(corpus):
    # The role is the one thing a query supplies that reaches the index SQL as text
    # rather than as a parameter. Unescaped, this literal closes the string and leaves
    # `OR '1'='1` as a condition, which matches every row at every path -- so the
    # assertion that bites is that nothing comes back, not that nothing raises.
    matched = _search(corpus, _role_and_structure("c1ccncc1", "SOLVENT' OR '1'='1"))
    assert matched == set()


def test_the_index_is_built_once_and_only_when_wanted(corpus_dir, caplog):
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        # An element predicate with no structure in it cannot use the index, so it must
        # not pay to build one.
        value.search(
            query.Query.model_validate(
                {
                    "where": _exists(
                        {"op": "eq", "path": "smiles", "value": {"literal": "CCO"}}
                    )
                }
            )
        )
        assert not value._occurrences_built
        caplog.set_level(logging.INFO, logger="ord_schema.agent.execute")
        request = query.Query.model_validate(_ROUTABLE["structure alone"])
        assert value.search(request).num_rows == 2
        assert value._occurrences_built
        # A second search reuses it rather than rebuilding, which the answer alone
        # cannot show -- a rebuild answers identically and costs a pass per path.
        assert value.search(request).num_rows == 2
    builds = [
        record for record in caplog.records if record.message.startswith("indexed ")
    ]
    assert len(builds) == 1


@pytest.mark.parametrize(
    ("path", "smarts", "expected"),
    [
        # A workup component: a structure at a path the inputs do not reach.
        ("workups.input.components", "c1ccncc1", {"ord-cc01"}),
        ("workups.input.components", "CCO", set()),
        # An authentic standard: the deepest path, reached by flattening twice.
        ("outcomes.products.measurements.authentic_standard", "c1ccccc1", {"ord-cc01"}),
        ("outcomes.products.measurements.authentic_standard", "CCO", set()),
    ],
)
def test_the_index_reaches_the_deep_paths_the_projection_does(
    deep_corpus, path, smarts, expected
):
    # Four paths are indexed and the main fixture holds structures at two of them, so
    # the other two are asserted here: the same query, routed and declined, over a
    # corpus that actually has something at each.
    where = {"op": "substructure", "path": "smiles", "smarts": smarts}
    body = {"where": {"op": "exists", "path": path, "where": where}}
    assert _index_spent(body), "expected the index to take the clause"
    request = query.Query.model_validate(body)
    assert _reactions(deep_corpus.search(request)) == expected
    # The same question compiled over the elements, so the projection answers what the
    # index just answered.
    with pytest.MonkeyPatch.context() as patcher:
        patcher.setattr(execute, "_index_condition", _no_index_condition)
        assert _reactions(deep_corpus.search(request)) == expected


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        # ord-cc02 records no workups and no outcomes, so the projection writes NULL at
        # both. Nothing pyridine-shaped sits at a level with no elements, which makes
        # the negation true of it -- and true is what the index says, since it holds no
        # occurrence for that reaction. A level left as NULL would answer neither way
        # and drop the reaction from the projection's answer alone.
        ("workups.input.components", {"ord-cc02"}),
        # The deepest path: ord-cc01 records one and it is benzene, so neither reaction
        # has pyridine there -- one by holding something else, the other by holding
        # nothing at all.
        (
            "outcomes.products.measurements.authentic_standard",
            {"ord-cc01", "ord-cc02"},
        ),
        # The level both reactions do record, so the same negation over a level that is
        # there says the same thing on either path.
        ("inputs.components", {"ord-cc01", "ord-cc02"}),
    ],
)
def test_a_negation_answers_the_same_over_a_level_the_source_never_recorded(
    deep_corpus, path, expected
):
    body = {
        "where": {
            "op": "not",
            "clause": {
                "op": "exists",
                "path": path,
                "where": {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
            },
        }
    }
    assert _index_spent(body), "expected the index to take the clause"
    request = query.Query.model_validate(body)
    assert _reactions(deep_corpus.search(request)) == expected
    with pytest.MonkeyPatch.context() as patcher:
        patcher.setattr(execute, "_index_condition", _no_index_condition)
        assert _reactions(deep_corpus.search(request)) == expected


def test_a_universal_holds_of_a_level_the_source_never_recorded(deep_corpus):
    # A forall is the absence of a counterexample, and a level with no elements holds
    # none. The projection writes NULL there, so this is the one place the answer would
    # otherwise be neither true nor false -- and a reaction would vanish from a question
    # it plainly satisfies.
    request = query.Query.model_validate(
        {
            "where": {
                "op": "forall",
                "path": "workups.input.components",
                "where": {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
            }
        }
    )
    assert _reactions(deep_corpus.search(request)) == {"ord-cc01", "ord-cc02"}


@pytest.mark.parametrize(
    ("role", "expected"), [("SOLVENT", {"ord-cc01"}), ("REACTANT", set())]
)
def test_a_role_binds_to_the_element_at_a_deep_path(deep_corpus, role, expected):
    # The index writes reaction_role for all four paths, and every other role condition
    # here asks about inputs.components. A role that came out right for the inputs and
    # NULL elsewhere would answer "pyridine as the solvent in the workup" with nothing,
    # silently, while the same query without the role still worked.
    request = query.Query.model_validate(
        {
            "where": {
                "op": "exists",
                "path": "workups.input.components",
                "where": {
                    "op": "and",
                    "clauses": [
                        {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
                        {
                            "op": "eq",
                            "path": "reaction_role",
                            "value": {"literal": role},
                        },
                    ],
                },
            }
        }
    )
    assert _reactions(deep_corpus.search(request)) == expected


def test_a_routed_query_is_bounded_by_the_timeout(corpus, monkeypatch):
    # The timeout is the only bound on what a query costs, and the routed path is a
    # different statement reached by a different branch. The fixture is too small to
    # outlast any timeout worth setting, so what is asserted is that the bound is
    # carried: a routed query dropping it would run unbounded on a real corpus.
    seen: list[tuple[str, float]] = []
    original = execute._run_with_timeout

    def recording(cursor, sql, parameters, timeout_seconds):
        seen.append((sql, timeout_seconds))
        return original(cursor, sql, parameters, timeout_seconds)

    monkeypatch.setattr(execute, "_run_with_timeout", recording)
    body = _ROUTABLE["structure and role"]
    assert _index_spent(body), "expected the index to take the clause"
    request = query.Query.model_validate(body)
    assert _reactions(corpus.search(request, timeout_seconds=60)) == {"ord-aa01"}
    assert [timeout for _, timeout in seen] == [60]
    assert "FROM occurrences" in seen[0][0]


def test_a_spent_clause_reads_the_index_and_a_declined_one_does_not(corpus_dir):
    # Everything else here compares the two paths, which agree -- so a search that
    # quietly stopped spending the index would keep passing while the index cost a
    # build and answered nothing. A row only the index holds separates them: the
    # semi-join filters reactions by membership, so redirecting an existing reaction
    # through a planted occurrence row changes the indexed answer and cannot change the
    # projection's.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        spent = query.Query.model_validate(_ROUTABLE["structure and role"])
        assert _reactions(value.search(spent)) == {"ord-aa01"}
        # aa01's solvent is its pyridine, so this copies the row the query just matched
        # on and attributes it to bb01, which holds no pyridine at all.
        value._connection.execute(
            "INSERT INTO occurrences SELECT 'ord-bb01', path, global_id, "
            "reaction_role FROM occurrences "
            "WHERE reaction_id = 'ord-aa01' AND reaction_role = 'SOLVENT'"
        )
        assert _reactions(value.search(spent)) == {"ord-aa01", "ord-bb01"}
        # The same question about an element field no occurrence row carries, so the
        # quantifier compiles over the elements and the planted row is invisible.
        declined_body = {
            "where": _exists(
                {
                    "op": "and",
                    "clauses": [
                        _SUBSTRUCTURE,
                        _SOLVENT,
                        {"op": "is_null", "path": "amount.mass_grams"},
                    ],
                }
            )
        }
        assert not _index_spent(declined_body)
        declined = query.Query.model_validate(declined_body)
        assert _reactions(value.search(declined)) == {"ord-aa01"}


def test_concurrent_first_searches_build_one_index(corpus_dir, caplog, monkeypatch):
    # The index is a materialized table, several passes over the projections. First
    # searches arriving together must not each build their own; whoever waits takes the
    # finished one. The barrier sits at the door of the build rather than in the
    # resolver, which a search reaches only afterwards, so the race is arranged rather
    # than hoped for -- and the build's own log line is what counts them, since
    # concurrent CREATEs would also collide in DuckDB's catalog and raise.
    caplog.set_level(logging.INFO, logger="ord_schema.agent.execute")
    threads_count = 4
    ready = threading.Barrier(threads_count)
    building = execute.Corpus._occurrences

    def synchronized(self):
        ready.wait(timeout=_WAIT)
        return building(self)

    monkeypatch.setattr(execute.Corpus, "_occurrences", synchronized)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={"pyridine": "c1ccncc1"}.__getitem__,
    ) as value:
        request = query.Query.model_validate(_ROUTABLE["a compound name"])
        results: list[int] = []
        failures: list[BaseException] = []

        def search():
            try:
                results.append(value.search(request).num_rows)
            except BaseException as error:  # noqa: BLE001
                failures.append(error)

        threads = [threading.Thread(target=search) for _ in range(threads_count)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=_WAIT)
        alive = [thread for thread in threads if thread.is_alive()]
    builds = [
        record for record in caplog.records if record.message.startswith("indexed ")
    ]
    assert failures == []
    assert alive == []
    assert len(builds) == 1
    assert results == [2] * threads_count


def test_an_index_that_reaches_nothing_is_refused(corpus_dir, monkeypatch):
    # An index that reaches no structure at all answers every structure query with no
    # matches, which reads like an answer rather than like a corpus nobody can search.
    # The fixture has no authentic standards, so indexing that path alone is a traversal
    # that reaches nothing -- what a path binding against the wrong schema looks like.
    path = "outcomes.products.measurements.authentic_standard"
    monkeypatch.setattr(execute, "INDEXED_PATHS", {path: execute.INDEXED_PATHS[path]})
    request = query.Query.model_validate(
        {
            "where": {
                "op": "exists",
                "path": path,
                "where": {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
            }
        }
    )
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        with pytest.raises(execute.PairingError, match="reached 0 of the corpus's 5"):
            value.search(request)
        # The table is in the catalog, holding nothing; what was not set is the flag, so
        # nobody reads it.
        assert not value._occurrences_built
        value._connection = counter = _CursorCounter(value._connection)
        with pytest.raises(execute.PairingError, match="reached 0 of the corpus's 5"):
            value.search(request)
        # A corpus does not change under an open Corpus, so the second refusal is the
        # first one's reason rather than several more passes to reach it again.
        assert counter.cursors == 0


def test_a_reaction_two_files_both_state_is_refused(tmp_path):
    # An occurrence names its reaction by ID, so a semi-join answers for every row
    # carrying that ID. Two files stating one reaction pair cleanly -- what pairing
    # checks is that each projection has its structures -- and each copy's structures
    # would then answer the other's queries: an original dataset globbed beside a
    # corrected re-release returns reactions whose own components hold nothing.
    for shard, smiles in (("aa", "c1ccncc1"), ("bb", "CCO")):
        source = tmp_path / "data" / f"ord_dataset-{shard}.parquet"
        source.parent.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            dataset_pb2.Dataset(
                dataset_id=f"ord_dataset-{shard}",
                name="test",
                description="test",
                reactions=[_reaction("ord-x01", components=[(smiles, _ROLE.SOLVENT)])],
            ),
            str(source),
        )
        projected = tmp_path / "projections" / source.name
        projected.parent.mkdir(parents=True, exist_ok=True)
        projection.write_projection(source, projected)
        structured = tmp_path / "structures" / source.name
        structured.parent.mkdir(parents=True, exist_ok=True)
        structures.write_structures(projected, structured)
    with (
        execute.Corpus(
            str(tmp_path / "projections" / "*.parquet"),
            str(tmp_path / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        ) as value,
        pytest.raises(execute.PairingError, match="ord-x01 is stated by 2 rows"),
    ):
        value.search(query.Query.model_validate(_ROUTABLE["structure alone"]))


def test_an_index_that_misses_one_path_is_refused(corpus_dir, monkeypatch):
    # The dangerous shape is not an index that reaches nothing -- it is one that reaches
    # almost everything. A single traversal that binds and finds no element leaves the
    # other paths carrying the total, so a count of occurrences looks healthy while
    # every query at the dead path answers with silence. Counting the structures reached
    # is what tells the two apart. Dropping a path is what a schema walk that fell out
    # of step with the projection would do.
    kept = {
        path: expression
        for path, expression in execute.INDEXED_PATHS.items()
        if path != "outcomes.products"
    }
    monkeypatch.setattr(execute, "INDEXED_PATHS", kept)
    # Toluene is every reaction's product and nobody's input, so dropping the products
    # path leaves exactly its two structures unreached.
    with (
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        ) as value,
        pytest.raises(execute.PairingError, match="reached 3 of the corpus's 5"),
    ):
        value.search(query.Query.model_validate(_ROUTABLE["structure alone"]))


def test_a_path_the_index_does_not_cover_is_left_to_the_projection(
    corpus_dir, monkeypatch
):
    # The planner's path check is unfalsifiable against the real schema, since every
    # path a structure can sit at is indexed. It becomes load-bearing the moment that
    # stops being true, and routing a query to a path the index does not carry answers
    # it with the empty set.
    kept = {
        path: expression
        for path, expression in execute.INDEXED_PATHS.items()
        if path != "outcomes.products"
    }
    monkeypatch.setattr(execute, "INDEXED_PATHS", kept)
    body = _ROUTABLE["products rather than inputs"]
    assert not _index_spent(body)
    request = query.Query.model_validate(body)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        # Answered from the projection, which needs no index and so never builds one.
        assert _reactions(value.search(request)) == {"ord-aa01", "ord-aa02", "ord-bb01"}
        assert not value._occurrences_built


def test_the_indexed_paths_are_the_ones_a_structure_can_sit_at():
    # Derived from the projection schema, so what the walk accepts and refuses is worth
    # pinning against schemas built for the purpose: today's projection exercises only
    # the accepting branch, and the refusals exist for the schema change that will.
    assert sorted(execute.INDEXED_PATHS) == [
        "inputs.components",
        "outcomes.products",
        "outcomes.products.measurements.authentic_standard",
        "workups.input.components",
    ]
    # A structure at a path with no role is one the index cannot carry -- and one it
    # would then report as a structure it failed to reach.
    with pytest.raises(ValueError, match="no reaction_role"):
        execute._indexed_paths(
            pa.schema(
                [
                    pa.field(
                        "inputs",
                        pa.list_(pa.struct([pa.field("structure_id", pa.uint32())])),
                    )
                ]
            )
        )
    # A structure on a level that is not repeated is one the build cannot unnest.
    with pytest.raises(ValueError, match="rather than a repeated expression"):
        execute._indexed_paths(
            pa.schema(
                [
                    pa.field(
                        "setup",
                        pa.struct(
                            [
                                pa.field("structure_id", pa.uint32()),
                                pa.field("reaction_role", pa.string()),
                            ]
                        ),
                    )
                ]
            )
        )
    # A schema with no structure anywhere would assemble a build with no SQL in it.
    with pytest.raises(ValueError, match="no structure at any path"):
        execute._indexed_paths(pa.schema([pa.field("reaction_id", pa.string())]))


def test_one_molecule_in_two_datasets_is_matched_once_and_found_twice(corpus):
    # Structures are deduplicated per dataset, so toluene holds a row in aa and another
    # in bb. The library matches it once; the answer still has to name both structures,
    # and through them every reaction that made it.
    library = corpus._library()
    assert len(library) < corpus._total, "expected the corpus to hold a duplicate"
    entries = {
        corpus._members[corpus._starts[entry]]: corpus._starts[entry + 1]
        - corpus._starts[entry]
        for entry in range(len(library))
    }
    assert max(entries.values()) == 2, "toluene should cover two structure IDs"
    assert sum(entries.values()) == corpus._total
    matched = corpus.search(
        query.Query.model_validate(
            {
                "where": {
                    "op": "exists",
                    "path": "outcomes.products",
                    "where": {
                        "op": "substructure",
                        "path": "smiles",
                        "smarts": "Cc1ccccc1",
                    },
                }
            }
        )
    )
    assert set(matched.column("reaction_id").to_pylist()) == {
        "ord-aa01",
        "ord-aa02",
        "ord-bb01",
    }


def _substructure(smarts):
    return query.Query.model_validate(
        {"where": _exists({"op": "substructure", "path": "smiles", "smarts": smarts})}
    )


def test_a_repeated_predicate_is_matched_once(corpus_dir, monkeypatch):
    # Matching is a pass over every molecule in the corpus and depends on nothing but
    # the query, so asking twice must cost once -- and must answer the same.
    matched = 0
    original = execute.Corpus._substructure_ids

    def counted(self, parameter, resolve):
        nonlocal matched
        matched += 1
        return original(self, parameter, resolve)

    monkeypatch.setattr(execute.Corpus, "_substructure_ids", counted)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={"pyridine": "c1ccncc1", "azabenzene": "c1ccncc1"}.__getitem__,
    ) as value:
        first = set(
            value.search(_substructure("c1ccncc1")).column("reaction_id").to_pylist()
        )
        second = set(
            value.search(_substructure("c1ccncc1")).column("reaction_id").to_pylist()
        )
        assert first == second == {"ord-aa01", "ord-aa02"}
        assert matched == 1
        # A different pattern is a different question, cache or no cache.
        assert set(
            value.search(_substructure("[OX2H]")).column("reaction_id").to_pylist()
        ) == {"ord-bb01"}
        assert matched == 2
        # Two names for one molecule are one match, since the key is what they resolve
        # to rather than what they were called.
        by_name = query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "compound": "pyridine"}
                )
            }
        )
        other_name = query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "compound": "azabenzene"}
                )
            }
        )
        assert value.search(by_name).num_rows == 2
        assert value.search(other_name).num_rows == 2
        # One match between them, and one more than before: a name resolves to SMILES,
        # and the identical text stated as a pattern was read as SMARTS.
        assert matched == 3


def test_the_cache_does_not_confuse_two_predicates(corpus):
    # One molecule, several questions. Toluene's Tanimoto against the inputs is benzene
    # 0.273, pyridine 0.176, ethanol 0.062, so each threshold below selects a different
    # set of reactions -- a key that ignored the threshold would answer one with
    # another. The same SMILES read as a substructure is a different operation again.
    def similarity(threshold):
        return set(
            corpus.search(
                query.Query.model_validate(
                    {
                        "where": _exists(
                            {
                                "op": "similarity",
                                "path": "smiles",
                                "smiles": "Cc1ccccc1",
                                "threshold": threshold,
                            }
                        )
                    }
                )
            )
            .column("reaction_id")
            .to_pylist()
        )

    tight = similarity(0.25)
    assert tight == {"ord-aa01"}  # Benzene only.
    assert similarity(0.15) == {"ord-aa01", "ord-aa02"}  # Pyridine joins it.
    assert similarity(0.05) == {"ord-aa01", "ord-aa02", "ord-bb01"}  # Ethanol too.
    # Substructure with that SMILES read as a SMARTS matches no input at all.
    assert (
        set(corpus.search(_substructure("Cc1ccccc1")).column("reaction_id").to_pylist())
        == set()
    )
    assert similarity(0.25) == tight  # Still itself after the others went through.


def test_one_predicate_asked_at_once_is_matched_once(corpus_dir, monkeypatch):
    # Several clients asking for the same scaffold at once is what a popular query looks
    # like, and the cache alone does not cover it: it holds finished answers, so callers
    # arriving while the first pass runs would each start their own. The pass below is
    # held open long enough that they all arrive during it.
    threads_count = 4
    matched = 0
    counting = threading.Lock()
    original = execute.Corpus._substructure_ids

    def held(self, parameter, resolve):
        nonlocal matched
        with counting:
            matched += 1
        time.sleep(0.5)
        return original(self, parameter, resolve)

    monkeypatch.setattr(execute.Corpus, "_substructure_ids", held)
    ready = threading.Barrier(threads_count)

    def resolver(name):
        # Releases every thread together, which is what puts them in one another's way.
        ready.wait(timeout=_WAIT)
        return {"pyridine": "c1ccncc1"}[name]

    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver=resolver,
    ) as value:
        value._library()  # Built up front, so the threads race on matching alone.
        request = query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "compound": "pyridine"}
                )
            }
        )
        results: list[list[str]] = []
        failures: list[BaseException] = []

        def search():
            try:
                matches = value.search(request).column("reaction_id").to_pylist()
                results.append(sorted(matches))
            except BaseException as error:  # noqa: BLE001
                failures.append(error)

        threads = [threading.Thread(target=search) for _ in range(threads_count)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=_WAIT)
        alive = [thread for thread in threads if thread.is_alive()]
        waiting = dict(value._matching)
    assert failures == []
    assert alive == []
    assert matched == 1
    assert results == [["ord-aa01", "ord-aa02"]] * threads_count
    assert waiting == {}  # The bookkeeping is dropped with the answer published.


def test_a_failed_match_is_not_left_in_the_way(corpus_dir, monkeypatch):
    # A pass that raises publishes nothing, so the predicate has to be left askable: the
    # error belongs to the caller that hit it, not to the corpus.
    failing = True
    original = execute.Corpus._substructure_ids

    def sometimes(self, parameter, resolve):
        if failing:
            raise RuntimeError("the library gave up")
        return original(self, parameter, resolve)

    monkeypatch.setattr(execute.Corpus, "_substructure_ids", sometimes)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        with pytest.raises(RuntimeError, match="gave up"):
            value.search(_substructure("c1ccncc1"))
        assert value._matching == {}
        failing = False
        assert set(
            value.search(_substructure("c1ccncc1")).column("reaction_id").to_pylist()
        ) == {"ord-aa01", "ord-aa02"}


def test_the_cache_stays_bounded(corpus_dir, monkeypatch):
    # Each entry is a bitmap over the whole corpus, so an unbounded cache would grow
    # without limit on a server asked many different questions.
    monkeypatch.setattr(execute, "_CACHED_MATCHES", 2)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        for smarts in ("c1ccncc1", "[OX2H]", "c1ccccc1", "[Pt]"):
            value.search(_substructure(smarts))
        assert len(value._matched) == 2
        # The two most recent survive; the earliest are gone.
        assert [key[2] for key in value._matched] == ["c1ccccc1", "[Pt]"]


@contextlib.contextmanager
def _no_narrow_table(self, columns: frozenset[str]) -> Iterator[str | None]:
    """Stands in for _narrowed_table so a search reads the projection directly."""
    del self, columns  # Unused.
    yield None


def test_a_narrow_table_answers_what_the_projection_would(corpus_dir):
    # Materializing the columns a query names is a second relation to answer from, so
    # what makes it safe is that it answers identically. Each of these is run against
    # the narrow table and against the projection, and compared as multisets: a query
    # with no ordering has none, and the two relations really do return the same rows
    # in different orders at corpus scale. What must never differ is which rows.
    requests = {
        "a scalar leaf": {
            "where": {
                "op": "not_null",
                "path": "conditions.temperature.setpoint_kelvin",
            }
        },
        "a nested quantifier": {
            "where": _exists(
                {"op": "eq", "path": "smiles", "value": {"literal": "CCO"}}
            )
        },
        "an aggregate": {
            "where": None,
            "aggregate": {
                "group_by": ["conditions.stirring.type"],
                "measures": [{"fn": "count", "name": "n"}],
            },
        },
        "an ordering and a limit": {
            "where": None,
            "order_by": [{"key": "reaction_id"}],
            "limit": 2,
        },
        "a structure predicate the index declines": {
            "where": {
                "op": "and",
                "clauses": [
                    _exists(
                        {"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}
                    ),
                    {"op": "not_null", "path": "reaction_id"},
                ],
            }
        },
    }
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        for label, body in requests.items():
            request = query.Query.model_validate(body)
            narrowed = value.search(request).to_pylist()
            # The same request with the materialization turned off.
            with pytest.MonkeyPatch.context() as patcher:
                patcher.setattr(execute.Corpus, "_narrowed_table", _no_narrow_table)
                direct = value.search(request).to_pylist()
            if request.order_by:
                # An ordered query is the one case where the order is the answer, so
                # this is the one comparison that may not sort first.
                assert list(map(repr, narrowed)) == list(map(repr, direct)), label
            else:
                assert sorted(map(repr, narrowed)) == sorted(map(repr, direct)), label


def test_the_columns_a_query_names_are_the_ones_materialized(corpus_dir):
    # Too few and the query fails to resolve; too many and the corpus holds columns
    # nobody asked about. A name inside a string literal costs a column, nothing more.
    assert execute._mentioned(
        "SELECT reaction_id FROM reactions WHERE conditions.temperature.x > 1"
    ) == frozenset({"reaction_id", "conditions"})
    assert "provenance" in execute._mentioned("SELECT provenance.doi FROM reactions")
    # reaction_id comes along whether or not the query said so.
    assert "reaction_id" in execute._mentioned("SELECT 1 FROM reactions")
    # A substring of a real column is not that column.
    assert "conditions" not in execute._mentioned(
        "SELECT reaction_id FROM reactions WHERE notes.preconditions_x = 'a'"
    )
    # An element's field is not the reaction's column of the same name: a query
    # filtering components by their SMILES reads no reaction-level smiles column, and
    # materializing one costs the largest column in the projection.
    assert "smiles" not in execute._mentioned(
        "SELECT reaction_id FROM reactions WHERE len(list_filter("
        "flatten(list_transform(map_values(inputs), x -> x.components)), "
        "e0 -> e0.smiles = 'CCO')) > 0"
    )
    # Where a path starts is where a top-level column is named.
    assert "smiles" in execute._mentioned("SELECT smiles FROM reactions")


def test_a_materialized_column_set_is_reused(corpus_dir):
    request = query.Query.model_validate(
        {"where": {"op": "not_null", "path": "conditions.temperature.setpoint_kelvin"}}
    )
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value.search(request)
        assert len(value._narrowed) == 1
        built = dict(value._narrowed)
        value.search(request)
        assert dict(value._narrowed) == built  # Reused, not rebuilt under a new name.
        # A query naming other columns is a different set, materialized separately.
        value.search(
            query.Query.model_validate(
                {"where": {"op": "not_null", "path": "provenance.doi"}}
            )
        )
        assert len(value._narrowed) == 2


def _stated_bytes(*readings):
    """Stands in for _memory_bytes, reading the given figures in order.

    Stated rather than measured so a budget test admits exactly the tables it means to,
    whatever DuckDB's allocator does with three rows.
    """
    figures = iter(readings)
    return lambda cursor: next(figures)


def _tables(value) -> set[str]:
    """Returns the names of the tables the corpus connection holds."""
    cursor = value._connection.cursor()
    try:
        rows = cursor.execute("SELECT table_name FROM duckdb_tables()").fetchall()
    finally:
        cursor.close()
    return {row[0] for row in rows}


def test_a_materialization_is_measured_in_bytes(corpus_dir):
    # duckdb_tables reports a row *count* under a name that reads like a size, and
    # spending it as one makes the budget a number that can never be crossed and the
    # eviction below dead code. Three rows, and the table is tens of kilobytes.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value.search(
            query.Query.model_validate(
                {"where": {"op": "not_null", "path": "provenance.doi"}}
            )
        )
        (entry,) = value._narrowed.values()
        assert entry.held > 1024


def test_the_cache_evicts_the_least_recently_used_and_drops_its_table(
    corpus_dir, monkeypatch
):
    # An unevicted cache is one table per column set a server is ever asked for, none of
    # them ever freed. Sizes are stated here rather than measured so the budget admits
    # exactly one table whatever DuckDB's allocator does with three rows.
    monkeypatch.setattr(execute, "_NARROW_BUDGET_BYTES", 150)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        monkeypatch.setattr(
            execute, "_memory_bytes", _stated_bytes(0, 100, 0, 100, 0, 100)
        )
        first = frozenset({"reaction_id", "provenance"})
        second = frozenset({"reaction_id", "conditions"})
        with value._narrowed_table(first) as name:
            dropped = name
        assert dropped in _tables(value)
        with value._narrowed_table(second) as name:
            kept = name
        assert list(value._narrowed) == [second]
        assert dropped not in _tables(value)  # Dropped, not merely forgotten.
        assert kept in _tables(value)


def test_a_table_a_search_is_reading_is_not_evicted(corpus_dir, monkeypatch):
    # A search takes the table's name, then resolves compounds and matches structures --
    # seconds during which it holds nothing but a string. Evicting there would fail a
    # query that had already been answered correctly, with a catalog error naming a
    # table the caller has never heard of.
    monkeypatch.setattr(execute, "_NARROW_BUDGET_BYTES", 150)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        monkeypatch.setattr(
            execute, "_memory_bytes", _stated_bytes(0, 100, 0, 100, 0, 100)
        )
        first = frozenset({"reaction_id", "provenance"})
        second = frozenset({"reaction_id", "conditions"})
        with value._narrowed_table(first) as reading:
            with value._narrowed_table(second):
                pass
            # Over budget, and the only candidate is being read: it stays, and the
            # cache stays over its budget until the reader is done.
            assert reading in _tables(value)
            assert list(value._narrowed) == [first, second]
        # Released, so the next materialization can take it.
        with value._narrowed_table(frozenset({"reaction_id", "notes"})):
            pass
        assert first not in value._narrowed
        assert reading not in _tables(value)


@contextlib.contextmanager
def _stalled_build(value, monkeypatch) -> Iterator[threading.Event]:
    """Runs a materialization that stops mid-build, and yields what releases it.

    The stall stands in for the pass over the corpus a real build spends its seconds
    in. Whatever the caller does while it holds is what a search does while another
    search is materializing.
    """
    building = threading.Event()
    finish = threading.Event()
    original = execute._memory_bytes

    def blocking(cursor):
        # The reading taken after the CREATE, so the table exists and the build is
        # holding whatever it holds while nothing else has recorded it.
        if building.is_set():
            return original(cursor)
        building.set()
        finish.wait()
        return original(cursor)

    monkeypatch.setattr(execute, "_memory_bytes", blocking)
    stalled: list[BaseException] = []

    def build() -> None:
        try:
            with value._narrowed_table(frozenset({"reaction_id", "conditions"})):
                pass
        except BaseException as error:  # noqa: BLE001 - reported below, not handled.
            stalled.append(error)

    builder = threading.Thread(target=build)
    builder.start()
    try:
        assert building.wait(timeout=30), "the build never reached its stall"
        yield finish
    finally:
        finish.set()
        builder.join(timeout=30)
    assert not builder.is_alive()
    assert stalled == []


def _in_a_thread(work) -> tuple[threading.Thread, list]:
    """Runs ``work`` in a thread, and returns it with the list its result lands in."""
    done: list = []
    thread = threading.Thread(target=lambda: done.append(work()))
    thread.start()
    return thread, done


def test_a_materialized_table_is_handed_out_while_another_is_being_built(
    corpus_dir, monkeypatch
):
    # A build is a pass over the corpus, and a Corpus is shared between searches. A
    # search whose columns are already materialized has nothing to wait for, and making
    # it wait turns one slow first query into a stall for everyone. Asked from a thread
    # rather than here, so a hit that does wait fails this test rather than hanging it.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        cached = frozenset({"reaction_id", "provenance"})
        with value._narrowed_table(cached):
            pass
        with _stalled_build(value, monkeypatch):

            def read() -> str | None:
                with value._narrowed_table(cached) as name:
                    return name

            reader, answered = _in_a_thread(read)
            reader.join(timeout=10)
            assert answered, "a cache hit waited on an unrelated build"
            assert answered[0] is not None
        assert not reader.is_alive()


def test_one_column_set_is_built_once_however_many_searches_want_it(
    corpus_dir, monkeypatch
):
    # Two searches naming the same columns arrive together on a cold cache. Building
    # twice costs two passes over the corpus and leaves two tables of one column set,
    # of which only the last is remembered -- the other is memory nothing will free.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        wanted = frozenset({"reaction_id", "conditions"})
        with _stalled_build(value, monkeypatch):

            def read() -> str | None:
                with value._narrowed_table(wanted) as name:
                    return name

            # The stalled build is materializing these very columns.
            second, answered = _in_a_thread(read)
            second.join(timeout=1)
            assert not answered, "the second ask did not wait for the build in flight"
        second.join(timeout=30)
        assert answered, "the second ask never got its table"
        assert answered[0] is not None
        assert value._narrow_serial == 1  # One name taken, so one table built.
        assert [name for name in _tables(value) if name.startswith("narrow_")] == [
            answered[0]
        ]


def test_a_table_a_second_search_took_from_the_cache_is_not_evicted(
    corpus_dir, monkeypatch
):
    # A hit hands back a name the search will read seconds later, exactly as a build
    # does, so it has to take the read the same way. Held at zero readers, the table is
    # the first thing eviction takes, and the search fails on a name that no longer
    # resolves after it had already been answered correctly.
    monkeypatch.setattr(execute, "_NARROW_BUDGET_BYTES", 150)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        monkeypatch.setattr(
            execute, "_memory_bytes", _stated_bytes(0, 100, 0, 100, 0, 100)
        )
        first = frozenset({"reaction_id", "provenance"})
        with value._narrowed_table(first) as built:
            pass
        # Taken from the cache this time, and read while the budget forces an eviction.
        with value._narrowed_table(first) as reading:
            assert reading == built
            with value._narrowed_table(frozenset({"reaction_id", "conditions"})):
                pass
            assert reading in _tables(value)


def test_a_failed_materialization_leaves_nothing_behind(corpus_dir, monkeypatch):
    # A table nobody tracks is memory nobody frees, and under a name a later attempt
    # would collide with. The failure is raised after the CREATE, which is the window.
    request = query.Query.model_validate(
        {"where": {"op": "not_null", "path": "provenance.doi"}}
    )
    failing = True
    original = execute._memory_bytes
    reads = 0

    def sometimes(cursor):
        # The second reading is the one taken after the CREATE, so failing there is the
        # window where a table exists and nothing has recorded it.
        nonlocal reads
        reads += 1
        if failing and reads == 2:
            raise duckdb.Error("no memory accounting today")
        return original(cursor)

    monkeypatch.setattr(execute, "_memory_bytes", sometimes)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        with pytest.raises(duckdb.Error, match="no memory accounting"):
            value.search(request)
        assert not value._narrowed
        assert not [name for name in _tables(value) if name.startswith("narrow_")]
        failing = False
        # The column set is still materializable, under a name of its own.
        value.search(request)
        assert len(value._narrowed) == 1


def test_a_literal_naming_the_relation_is_not_a_rewrite(tmp_path):
    # The columns are read off the compiled SQL, which carries string literals inline --
    # so a query searching for the text "FROM reactions" puts that text in the SQL. The
    # narrow table is reached by compiling again against it, not by editing that SQL: an
    # edit would rewrite the literal too, and this reaction would stop being findable by
    # what its own notes say.
    note = "distilled FROM reactions overnight"
    reaction = _reaction("ord-dd01", components=[("CCO", _ROLE.REACTANT)])
    reaction.notes.procedure_details = note
    source = tmp_path / "data" / "ord_dataset-dd.parquet"
    source.parent.mkdir(parents=True, exist_ok=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-dd",
            name="test",
            description="test",
            reactions=[reaction],
        ),
        str(source),
    )
    projected = tmp_path / "projections" / source.name
    projected.parent.mkdir(parents=True, exist_ok=True)
    projection.write_projection(source, projected)
    structured = tmp_path / "structures" / source.name
    structured.parent.mkdir(parents=True, exist_ok=True)
    structures.write_structures(projected, structured)
    request = query.Query.model_validate(
        {
            "where": {
                "op": "eq",
                "path": "notes.procedure_details",
                "value": {"literal": note},
            }
        }
    )
    with execute.Corpus(
        str(projected), str(structured), resolver={}.__getitem__
    ) as value:
        assert _reactions(value.search(request)) == {"ord-dd01"}
        assert len(value._narrowed) == 1


def test_a_column_set_too_large_to_keep_is_not_kept(corpus_dir, monkeypatch):
    # Materializing is only worth it if the table survives to answer the next query.
    monkeypatch.setattr(execute, "_NARROW_BUDGET_BYTES", 1)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        matched = value.search(
            query.Query.model_validate(
                {
                    "where": _exists(
                        {"op": "eq", "path": "smiles", "value": {"literal": "CCO"}}
                    )
                }
            )
        )
        # Answered anyway, straight from the projection.
        assert matched.column("reaction_id").to_pylist() == ["ord-bb01"]
        assert not value._narrowed


def test_a_column_set_too_large_to_keep_is_not_built_twice(corpus_dir, monkeypatch):
    # What the refusal costs is a pass over the corpus, and the projection it reads
    # does not change while the corpus is open. Asking again is the common case -- a
    # notebook loop, a server serving one shape of question -- and rebuilding a table
    # only to drop it again holds the narrow lock against every other search each time.
    monkeypatch.setattr(execute, "_NARROW_BUDGET_BYTES", 1)
    columns = frozenset({"reaction_id", "provenance"})
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        with value._narrowed_table(columns) as name:
            assert name is None
        built = value._narrow_serial
        with value._narrowed_table(columns) as name:
            assert name is None
        assert value._narrow_serial == built  # No second CREATE to name.


def test_a_library_that_lost_a_molecule_is_refused(corpus_dir, monkeypatch):
    # The holders and the entry table are filled in one branch, so a divergence means
    # one of them missed a row -- after which every entry above it maps to a neighbor's
    # structures: in range, wrong molecule, and invisible to a count of either alone.
    class _Short:
        """A library that holds fewer molecules than the corpus has distinct SMILES."""

        def __init__(self, molecules, patterns):
            del molecules, patterns  # Unused: nothing gets as far as searching it.

        def __len__(self):
            return 1

    monkeypatch.setattr(execute.rdSubstructLibrary, "SubstructLibrary", _Short)
    with (
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        ) as value,
        pytest.raises(execute.PairingError, match="distinct SMILES"),
    ):
        value._library()


def test_a_pattern_and_a_name_that_read_alike_are_not_one_answer(corpus_dir):
    # A stated substructure pattern is SMARTS; a resolved compound is SMILES. The same
    # text is two query molecules across those parsers -- C1=CC=CC=C1 is benzene as
    # SMILES and six aliphatic carbons as SMARTS -- and resolvers answer in exactly that
    # Kekule form, so a key holding only the text would answer one with the other.
    kekule = "C1=CC=CC=C1"
    by_name = query.Query.model_validate(
        {"where": _exists({"op": "substructure", "path": "smiles", "compound": "b"})}
    )

    def answers(first, second):
        with execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={"b": kekule}.__getitem__,
        ) as value:
            return _reactions(value.search(first)), _reactions(value.search(second))

    # The corpus stores benzene aromatic, so the name matches it and the pattern, read
    # as six aliphatic carbons, matches nothing. Whichever runs first, both hold.
    assert answers(_substructure(kekule), by_name) == (set(), {"ord-aa01"})
    assert answers(by_name, _substructure(kekule)) == ({"ord-aa01"}, set())


def test_a_structure_parameter_naming_nothing_is_an_error(corpus):
    # The grammar guarantees one of the two, so this cannot come from a query -- but
    # resolving "" would ask the external resolver for a compound with no name and cache
    # whatever came back, which is a long way from where the mistake was.
    parameter = query.StructureParameter(
        name="p0", op="substructure", pattern=None, compound=None, threshold=None
    )
    cursor = corpus._connection.cursor()
    try:
        with pytest.raises(ValueError, match="neither a pattern nor a compound"):
            corpus._matches(cursor, parameter, {}.__getitem__)
    finally:
        cursor.close()


def test_the_cache_keeps_what_was_asked_for_most_recently(corpus_dir, monkeypatch):
    # Least recently *used*, not least recently computed: a scaffold asked for over and
    # over is exactly what the cache is for, and it would be the first thing dropped if
    # a hit did not count as a use.
    monkeypatch.setattr(execute, "_CACHED_MATCHES", 2)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value.search(_substructure("c1ccncc1"))
        value.search(_substructure("[OX2H]"))
        value.search(_substructure("c1ccncc1"))  # Asked again: now the newest.
        value.search(_substructure("[Pt]"))
        assert [key[2] for key in value._matched] == ["c1ccncc1", "[Pt]"]


def test_a_waiter_behind_a_failed_match_answers_for_itself(corpus_dir, monkeypatch):
    # A pass that raises publishes nothing, so the callers queued behind it have to take
    # the matching over rather than inherit an error they had no part in -- and none of
    # them may be left waiting on an event that will never be set.
    threads_count = 4
    failing = True
    original = execute.Corpus._substructure_ids
    ready = threading.Barrier(threads_count)

    def sometimes(self, parameter, resolve):
        if failing:
            time.sleep(0.5)  # Long enough that the others are waiting on this one.
            raise RuntimeError("the library gave up")
        return original(self, parameter, resolve)

    monkeypatch.setattr(execute.Corpus, "_substructure_ids", sometimes)
    released = threading.Event()

    def resolver(name):
        # The barrier is one-shot: every thread reaches it before any gets through, so
        # the searches after them do not wait for a crowd that has already dispersed.
        if not released.is_set():
            ready.wait(timeout=_WAIT)
            released.set()
        return {"pyridine": "c1ccncc1"}[name]

    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver=resolver,
    ) as value:
        value._library()
        request = query.Query.model_validate(_ROUTABLE["a compound name"])
        failures: list[BaseException] = []
        results: list[int] = []

        def search():
            try:
                results.append(value.search(request).num_rows)
            except BaseException as error:  # noqa: BLE001
                failures.append(error)

        threads = [threading.Thread(target=search) for _ in range(threads_count)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=_WAIT)
        alive = [thread for thread in threads if thread.is_alive()]
        assert alive == []
        assert results == []
        assert [type(error) for error in failures] == [RuntimeError] * threads_count
        assert value._matching == {}
        # Nothing was cached, so the predicate is still askable once matching works.
        failing = False
        assert value.search(request).num_rows == 2


def test_the_match_limit_has_to_be_passed_or_matches_are_dropped():
    # GetMatches defaults to 1000 results and truncates in silence, which at corpus
    # scale would turn a broad pattern into a wrong answer rather than a slow one. The
    # fixture has five structures, so no test of it can reach the limit; the default is
    # pinned here instead, against a library large enough to cross it.
    benzene = Chem.MolFromSmiles("c1ccccc1")
    holder = rdSubstructLibrary.CachedMolHolder()
    for _ in range(1500):
        holder.AddMol(benzene)
    library = rdSubstructLibrary.SubstructLibrary(holder)
    pattern = Chem.MolFromSmarts("c1ccccc1")
    assert len(library.GetMatches(pattern)) == 1000
    assert len(library.GetMatches(pattern, maxResults=len(library))) == 1500


def test_the_artifact_fingerprint_is_the_one_the_library_screens_with():
    # The library is fed the pattern_fp the artifact already stores rather than
    # recomputing it. That is only sound while RDKit's own PatternHolder fingerprint
    # is the same function -- if it ever diverges, the screen would start rejecting
    # true matches, so the equality is pinned here rather than assumed.
    molecule = Chem.MolFromSmiles("Cc1ccncc1O")
    holder = rdSubstructLibrary.PatternHolder(structures.PATTERN_FP_SIZE)
    holder.AddMol(molecule)
    stored = Chem.PatternFingerprint(molecule, fpSize=structures.PATTERN_FP_SIZE)
    assert list(holder.GetFingerprint(0).GetOnBits()) == list(stored.GetOnBits())
