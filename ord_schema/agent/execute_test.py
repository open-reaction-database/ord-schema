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
        assert len(value._library()) == value._total
        # Benzene sits after the hole inside aa, ethanol after it over in bb. Both
        # still name their own reactions, which holds only while the slot is held.
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
            execute.Corpus._run_with_timeout(
                cursor, "SELECT count(*) FROM range(100000000000)", {}, 0.2
            )
        # The bound is what ends it, rather than the query finishing on its own.
        assert time.perf_counter() - started < 30
    finally:
        connection.close()


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
        assert value.search(request).num_rows == 1
        assert counter.cursors == 1


def test_concurrent_searches_return_their_own_answers(corpus_dir):
    # Distinct queries with distinct answers, run at once against one corpus: each
    # thread has to come back with the reactions its own query selected. Sharing one
    # connection also makes searches raise "No open result set" as often as it makes
    # them lie, so the exceptions are collected too -- a thread that dies takes its
    # traceback to threading.excepthook, where an assertion on results alone sees
    # nothing wrong.
    # The similarity thread earns its place: that path runs its screen and then the
    # compiled query as two statements on one cursor, so it is where reuse within a
    # single search would show up.
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
    # The library is about 2 GB at corpus scale. First searches arriving together must
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
