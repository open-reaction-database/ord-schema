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

"""Tests for ord_schema.search.execute."""

import contextlib
import dataclasses
import logging
import pathlib
import re
import shutil
import threading
import time
from collections.abc import Iterator, Sequence

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary

from ord_schema import parquet
from ord_schema.artifacts import (
    base,
    occurrences,
    pivot,
    projection,
    structures,
)
from ord_schema.proto import dataset_pb2, reaction_pb2
from ord_schema.search import execute, query

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
def deep_root(tmp_path_factory) -> pathlib.Path:
    """Sources, projections, and structures holding a structure at every indexed path.

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
    return root


@pytest.fixture(scope="module")
def deep_corpus(deep_root) -> Iterator[execute.Corpus]:
    """The corpus over ``deep_root``, holding no pivot artifacts."""
    with execute.Corpus(
        str(deep_root / "projections" / "*.parquet"),
        str(deep_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
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


def _exists_at(path, where):
    return {"op": "exists", "path": path, "where": where}


def _exists(where):
    return _exists_at("inputs.components", where)


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


def _capped(monkeypatch, tmp_path, *caps: str, relative: str = "") -> None:
    """Builds a cgroup v2 tree this process sits at the bottom of.

    Args:
        monkeypatch: The fixture, which points the readers at the tree.
        tmp_path: Where to build it.
        caps: What each level states, root first, one per level of ``relative`` plus the
            root itself. A level with no cap given states none at all.
        relative: The cgroup this process occupies, relative to the mount.
    """
    root = tmp_path / "cgroup"
    directories = [root]
    for part in relative.split("/") if relative else []:
        directories.append(directories[-1] / part)
    for directory, cap in zip(directories, caps, strict=False):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "memory.max").write_text(cap)
    directories[-1].mkdir(parents=True, exist_ok=True)

    proc = tmp_path / "cgroup.proc"
    proc.write_text(f"0::/{relative}\n")
    monkeypatch.setattr(execute, "_CGROUP_V2_ROOT", root)
    monkeypatch.setattr(execute, "_CGROUP_V1_MEMORY_ROOT", tmp_path / "absent")
    monkeypatch.setattr(execute, "_PROC_SELF_CGROUP", proc)


def _open(corpus_dir, **kwargs) -> execute.Corpus:
    return execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        **kwargs,
    )


def test_a_memory_limit_is_given_to_duckdb(corpus_dir):
    with _open(corpus_dir, memory_limit="512MiB") as corpus:
        setting = corpus._connection.execute(
            "SELECT current_setting('memory_limit')"
        ).fetchone()
    assert setting is not None
    assert execute._setting_bytes(setting[0]) == 512 * 1024**2


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("512.0 MiB", 512 * 1024**2),
        ("6.3 GiB", int(6.3 * 1024**3)),
        ("1024 B", 1024),
        # What the setting says when nothing bounds it.
        ("-1", None),
    ],
)
def test_a_size_setting_reads_as_bytes(value, expected):
    assert execute._setting_bytes(value) == expected


@pytest.mark.parametrize("cap", ["max", str(2**63 - 4096)])
def test_an_unbounded_cgroup_is_no_cap(monkeypatch, tmp_path, cap):
    # Both versions spell "unlimited" as something readable, so a corpus on a machine
    # with cgroups but no limit must not warn about the room it has.
    _capped(monkeypatch, tmp_path, cap)
    assert execute._container_cap_bytes() is None


def test_a_cap_above_this_process_is_found(monkeypatch, tmp_path):
    # A container inside a slice sits several levels below the mount, and the level
    # holding the cap is not the one it occupies: reading only the root finds nothing
    # here, and only the leaf finds the wrong thing.
    _capped(
        monkeypatch,
        tmp_path,
        "max",
        str(8 * 1024**3),
        "max",
        relative="system.slice/docker-abc.scope",
    )
    assert execute._container_cap_bytes() == 8 * 1024**3


def test_the_smallest_cap_that_applies_is_the_one_that_counts(monkeypatch, tmp_path):
    # Every level's cap binds, so the process reaches the smallest of them first --
    # which is not always the innermost.
    _capped(
        monkeypatch,
        tmp_path,
        str(4 * 1024**3),
        str(12 * 1024**3),
        relative="nested",
    )
    assert execute._container_cap_bytes() == 4 * 1024**3


def test_a_cgroup_line_for_another_controller_is_not_read_as_this_one(
    monkeypatch, tmp_path
):
    # A v1 line names its controllers where the v2 line names none, so matching by
    # position rather than by name would read one hierarchy's path against the other's
    # mount.
    root = tmp_path / "cgroup"
    (root / "leaf").mkdir(parents=True)
    (root / "memory.max").write_text("max")
    (root / "leaf" / "memory.max").write_text(str(2 * 1024**3))
    proc = tmp_path / "cgroup.proc"
    proc.write_text("5:memory:/other\n0::/leaf\n")
    monkeypatch.setattr(execute, "_CGROUP_V2_ROOT", root)
    monkeypatch.setattr(execute, "_CGROUP_V1_MEMORY_ROOT", tmp_path / "absent")
    monkeypatch.setattr(execute, "_PROC_SELF_CGROUP", proc)
    assert execute._container_cap_bytes() == 2 * 1024**3


def test_a_cap_leaving_the_process_nothing_is_warned_about(
    corpus_dir, tmp_path, monkeypatch, caplog
):
    _capped(monkeypatch, tmp_path, str(2 * 1024**3))
    caplog.set_level(logging.WARNING, logger="ord_schema.search.execute")
    caplog.clear()
    with _open(corpus_dir, memory_limit="1GiB"):
        pass
    assert [
        record for record in caplog.records if "1.0 GB for everything" in record.message
    ]


def test_a_cap_leaving_the_process_room_is_not_warned_about(
    corpus_dir, tmp_path, monkeypatch, caplog
):
    _capped(monkeypatch, tmp_path, str(12 * 1024**3))
    caplog.set_level(logging.WARNING, logger="ord_schema.search.execute")
    caplog.clear()
    with _open(corpus_dir, memory_limit="1GiB"):
        pass
    assert not [
        record for record in caplog.records if "for everything" in record.message
    ]


def test_no_cap_to_read_is_not_warned_about(corpus_dir, tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(execute, "_CGROUP_V2_ROOT", tmp_path / "absent")
    monkeypatch.setattr(execute, "_CGROUP_V1_MEMORY_ROOT", tmp_path / "absent")
    monkeypatch.setattr(execute, "_PROC_SELF_CGROUP", tmp_path / "absent")
    caplog.set_level(logging.WARNING, logger="ord_schema.search.execute")
    caplog.clear()
    with _open(corpus_dir, memory_limit="1GiB"):
        pass
    assert not [
        record for record in caplog.records if "for everything" in record.message
    ]


def test_a_resolved_name_is_compared_in_the_corpus_s_own_spelling(corpus_dir):
    # PubChem answers with a Kekule form and the projection stores what RDKit
    # canonicalizes, so comparing the two as strings matches nothing -- and says so by
    # returning no rows, which is the worst way for a query to be wrong.
    with _open(corpus_dir, resolver={"pyridine": "C1=CC=NC=C1"}.__getitem__) as corpus:
        table = corpus.search(
            query.Query.model_validate(
                {
                    "where": {
                        "op": "exists",
                        "path": "inputs.components",
                        "where": {
                            "op": "eq",
                            "path": "smiles",
                            "value": {"compound": "pyridine"},
                        },
                    }
                }
            )
        )
    assert table.num_rows > 0


def test_a_name_resolving_to_nonsense_is_left_alone(corpus_dir):
    # Canonicalizing cannot fix an unparseable answer, and swallowing it here would
    # hide a resolver fault behind an empty result.
    with _open(
        corpus_dir, resolver={"nonsense": "not a molecule"}.__getitem__
    ) as corpus:
        table = corpus.search(
            query.Query.model_validate(
                {
                    "where": {
                        "op": "exists",
                        "path": "inputs.components",
                        "where": {
                            "op": "eq",
                            "path": "smiles",
                            "value": {"compound": "nonsense"},
                        },
                    }
                }
            )
        )
    assert table.num_rows == 0


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
            warm=False,
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
            warm=False,
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


def test_the_footer_cache_reaches_the_cursor_a_search_runs_on(corpus):
    # A cursor is its own session, and a session-scoped setting would leave the parse
    # charged to every search while the corpus's own connection -- which runs no
    # queries -- reported the cache as on.
    cursor = corpus._connection.cursor()
    try:
        setting = cursor.execute(
            "SELECT current_setting('parquet_metadata_cache')"
        ).fetchone()
    finally:
        cursor.close()
    assert setting == (True,)


def test_a_rewritten_artifact_is_not_answered_from_its_cached_footer(
    corpus_dir, tmp_path
):
    # Holding the footers means holding something the file can go on to contradict. The
    # rewrite lands in the same second as the read that cached it, which is the case a
    # cache keyed on a coarse timestamp would answer stale. Nothing is materialized, so
    # both searches reach the Parquet files rather than the second one reading a table
    # the first one built.
    root = _copy_corpus(corpus_dir, tmp_path)
    request = query.Query.model_validate(
        {
            "where": _exists(
                {
                    "op": "eq",
                    "path": "reaction_role",
                    "value": {"literal": "SOLVENT"},
                }
            )
        }
    )
    with execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=0,
    ) as value:
        assert _reactions(value.search(request)) == {"ord-aa01"}
        projected = root / "projections" / "ord_dataset-aa.parquet"
        table = pq.read_table(projected)
        kept = table.filter(
            pa.compute.not_equal(table.column("reaction_id"), "ord-aa01")
        )
        pq.write_table(kept, projected)
        assert _reactions(value.search(request)) == set()


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
        warm=False,
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
    "same compound": {
        "where": _exists(
            {"op": "same_compound", "path": "smiles", "smiles": "c1ccncc1"}
        )
    },
    "same parent": {
        "where": _exists({"op": "same_parent", "path": "smiles", "smiles": "c1ccncc1"})
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
    # The shape the index answers half of: over the corpus, what the other half costs
    # is set by whether its column was materialized -- 0.43s against 3.34s for the same
    # question, which is why the budget is what it is.
    "a clause the projection has to answer beside it": {
        "where": {
            "op": "and",
            "clauses": [
                _exists(_SUBSTRUCTURE),
                {
                    "op": "exists",
                    "path": "outcomes.products",
                    "where": {
                        "op": "eq",
                        "path": "isolated_color",
                        "value": {"literal": "white"},
                    },
                },
            ],
        }
    },
    "a negated clause beside the indexed one": {
        "where": {
            "op": "and",
            "clauses": [
                _exists(_SUBSTRUCTURE),
                {"op": "not", "clause": {"op": "not_null", "path": "provenance.doi"}},
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
        warm=False,
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
        caplog.set_level(logging.INFO, logger="ord_schema.search.execute")
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


def test_the_indexed_paths_and_field_come_from_the_artifact():
    # One walk, in the artifact, and this pins that it stays one. A path routed here
    # that no artifact carries is a quantifier answered from nothing, and a field name
    # spelled here rather than read is a column the artifact does not hold.
    assert set(execute.INDEXED_PATHS) == set(occurrences.PATHS)
    assert execute._INDEXED_FIELD == occurrences.INDEXED_FIELD


def test_a_corpus_builds_the_occurrence_index_at_open(corpus_dir):
    # The default, so the first structure query answers from an index that is already
    # there rather than paying for one -- over ORD a pass per uncovered path.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        assert value._occurrences_built


def test_a_cold_corpus_leaves_the_index_to_the_first_query(corpus_dir):
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        warm=False,
    ) as value:
        assert not value._occurrences_built


def test_max_rows_bounds_what_a_search_returns(corpus_dir, caplog):
    # The corpus holds three reactions and the query asks for none of them in
    # particular, which is the shape that returns everything at corpus scale.
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            max_rows=2,
        ) as value,
    ):
        assert value.search(query.Query.model_validate({})).num_rows == 2
    assert any("bounds this query at 2 rows" in m for m in _messages(caplog))


def test_an_answer_that_reaches_the_bound_says_it_may_be_cut(corpus_dir, caplog):
    # Cut rather than refused: the answer is nearer what the caller wanted than an
    # error, so the log says it may be short, and the table carries TRUNCATED for a
    # reader who never sees the log.
    with (
        caplog.at_level(logging.WARNING, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            max_rows=1,
        ) as value,
    ):
        found = value.search(query.Query.model_validate({"limit": 3}))
    assert found.num_rows == 1
    assert any("there may be matches it did not return" in m for m in _messages(caplog))


def test_an_answer_the_bound_never_reached_says_nothing(corpus_dir, caplog):
    # The ordinary case, and the reason the warning waits for the result: warning
    # wherever the bound was merely applied is how a reader learns to skip the line.
    with (
        caplog.at_level(logging.WARNING, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            max_rows=100,
        ) as value,
    ):
        assert value.search(query.Query.model_validate({"limit": 5000})).num_rows == 3
    assert not _messages(caplog)


def test_a_truncated_aggregate_says_the_distribution_is_partial(corpus_dir, caplog):
    # A truncated list of reactions is a sample of the matching ones; a truncated list
    # of groups is part of a distribution read as the whole of it, and ordered by
    # nothing it is an arbitrary part.
    request = {
        "aggregate": {
            "group_by": ["reaction_id"],
            "measures": [{"fn": "count", "name": "n"}],
        }
    }
    with (
        caplog.at_level(logging.WARNING, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            max_rows=2,
        ) as value,
    ):
        assert value.search(query.Query.model_validate(request)).num_rows == 2
    assert any("an arbitrary part of the distribution" in m for m in _messages(caplog))


def test_a_search_under_a_smaller_limit_than_max_rows_is_untouched(corpus_dir):
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        max_rows=100,
    ) as value:
        assert value.search(query.Query.model_validate({"limit": 2})).num_rows == 2


def test_a_max_rows_no_query_can_satisfy_is_refused(corpus_dir):
    with pytest.raises(ValueError, match="which no query can return"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            max_rows=0,
        )


def test_a_corpus_builds_the_substructure_library_at_open(corpus_dir):
    # Warming covers both builds a substructure query needs, so what an open corpus
    # holds is what a container has to be sized for rather than a floor a later query
    # raises it to.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        assert value._substructure_library is not None


def test_a_cold_corpus_leaves_the_library_to_the_first_query(corpus_dir):
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        warm=False,
    ) as value:
        assert value._substructure_library is None


def test_projections_that_do_not_join_to_their_offsets_are_refused(
    corpus_dir, tmp_path, monkeypatch
):
    # read_parquet globs every path it is handed, so anything making DuckDB resolve one
    # differently than Python did drops that file's rows from the join rather than
    # failing, leaving reactions no query can find. Staged by pointing the view at
    # copies filed elsewhere, which is what such a resolution looks like from here. The
    # structures side has been counted against its footers all along; the reactions
    # side was covered only by the occurrence index reaching every structure, and a
    # path read from a pivot artifact reaches them whatever the view holds.
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    for projected in sorted((corpus_dir / "projections").glob("*.parquet")):
        shutil.copy(projected, elsewhere / projected.name)
    original = execute._sql_paths

    def redirected(paths):
        listed = [str(path) for path in paths]
        if all("projections" in path for path in listed):
            return original(str(elsewhere / pathlib.Path(path).name) for path in listed)
        return original(listed)

    monkeypatch.setattr(execute, "_sql_paths", redirected)
    with pytest.raises(execute.PairingError, match="the projections hold"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        )


def test_deriving_pivots_beside_warming_is_refused(corpus_dir, tmp_path):
    # Warming reads the levels the index covers, and a set covering some of the
    # projections is what that read refuses -- so the two together would make an
    # interrupted derivation permanent, the pass that completes it sitting behind the
    # refusal.
    with pytest.raises(ValueError, match="derive_pivots was set beside warm"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(tmp_path / "pivots"),
            derive_pivots=True,
        )


def test_warming_refuses_at_open_what_the_index_would_refuse_on_a_query(tmp_path):
    # Two datasets stating one reaction pair cleanly and only the index refuses them,
    # so this is the refusal warming is able to move. Opening such a corpus and failing
    # later leaves one that looks searchable and is not.
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
    with pytest.raises(execute.PairingError, match="ord-x01 is stated by 2 rows"):
        execute.Corpus(
            str(tmp_path / "projections" / "*.parquet"),
            str(tmp_path / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        )


def test_checking_the_index_builds_it_and_counts_every_path(corpus_dir):
    # The point of the check is that it happens at startup rather than inside whichever
    # query first wants the index, so it has to leave the index built. Every indexed
    # path is reported, including the ones this corpus records nothing at: a zero is
    # ordinary, and omitting it would read as a path the walk lost.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        counts = value.check_index()
        assert value._occurrences_built
        assert set(counts) == set(execute.INDEXED_PATHS)
        # Four components across the three reactions, and nothing anywhere else.
        assert counts["inputs.components"] == 4
        assert counts["workups.input.components"] == 0


def test_checking_a_corpus_the_index_refuses_fails_the_check(corpus_dir, tmp_path):
    # The refusal a deployment wants at startup: a corpus stating one reaction twice
    # would have each copy's structures answering for the other. Raised here rather
    # than at whatever hour the first structure query arrives.
    root = _copy_corpus(corpus_dir, tmp_path)
    doubled = root / "projections" / "ord_dataset-aa.parquet"
    table = pq.read_table(doubled)
    pq.write_table(pa.concat_tables([table, table]), doubled)
    with (
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            require_current=False,
            warm=False,
        ) as value,
        pytest.raises(execute.PairingError, match="is stated by 2 rows"),
    ):
        value.check_index()


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
    # What reaches the query is what the earlier phases left, not the whole bound: the
    # screen and the index build come out of the same budget.
    (spent,) = [timeout for _, timeout in seen]
    assert 0 < spent <= 60
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
    caplog.set_level(logging.INFO, logger="ord_schema.search.execute")
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
        warm=False,
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
        warm=False,
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
            warm=False,
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
            warm=False,
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
        warm=False,
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
    # A structure on a level that is not repeated is one no pivot holds the elements of,
    # so no artifact can carry it and the build would have to unnest a level that is not
    # one. Refused by the artifact's walk, which is the walk this reads.
    with pytest.raises(ValueError, match="no repeated level reaches it"):
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


def _over_products(where):
    return {"op": "exists", "path": "outcomes.products", "where": where}


_WHITE_PRODUCT = _over_products(
    {"op": "eq", "path": "isolated_color", "value": {"literal": "white"}}
)


def test_a_materialized_pivot_is_reused(corpus_dir):
    request = query.Query.model_validate(
        {"where": _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}})}
    )
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value.search(request)
        assert list(value._pivoted) == ["inputs.components"]
        built = dict(value._pivoted)
        value.search(request)
        assert dict(value._pivoted) == built  # Reused, not rebuilt under a new name.
        # A quantifier over another level is another pivot, materialized separately.
        value.search(query.Query.model_validate({"where": _WHITE_PRODUCT}))
        assert list(value._pivoted) == ["inputs.components", "outcomes.products"]


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


def test_a_budget_no_memory_could_answer_to_is_refused(corpus_dir):
    # Negative is not an amount of memory, and taken as one it makes every pivot too
    # large and every eviction immediate, which reads as a corpus that cannot cache.
    with pytest.raises(ValueError, match="not an amount of memory"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivot_budget_bytes=-1,
        )
    # Zero is an amount of memory, and a caller asking to materialize nothing means it
    # -- including the build that would otherwise measure a pivot before dropping it,
    # which is the peak a machine has to have room for.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=0,
    ) as value:
        value.search(query.Query.model_validate({"where": _WHITE_PRODUCT}))
        assert not value._pivoted
        assert value._pivot_serial == 0  # No table was built to find that out.
        assert not [name for name in _tables(value) if name.startswith("pivoted_")]


def test_a_materialization_is_measured_in_bytes(corpus_dir):
    # duckdb_tables reports a row *count* under a name that reads like a size, and
    # spending it as one makes the budget a number that can never be crossed and the
    # eviction below dead code. Three rows, and the table is tens of kilobytes.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        value.search(query.Query.model_validate({"where": _WHITE_PRODUCT}))
        (entry,) = value._pivoted.values()
        assert entry.held > 1024


def test_the_cache_evicts_the_least_recently_used_and_drops_its_table(
    corpus_dir, monkeypatch
):
    # An unevicted cache is one table per level a server is ever asked about, none of
    # them ever freed. Sizes are stated here rather than measured so the budget admits
    # exactly one table whatever DuckDB's allocator does with three rows.
    budget = 150
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=budget,
    ) as value:
        monkeypatch.setattr(
            execute, "_memory_bytes", _stated_bytes(0, 100, 0, 100, 0, 100)
        )
        with value._pivoted_table("inputs.components") as name:
            dropped = name
        assert dropped in _tables(value)
        with value._pivoted_table("outcomes.products") as name:
            kept = name
        assert list(value._pivoted) == ["outcomes.products"]
        assert dropped not in _tables(value)  # Dropped, not merely forgotten.
        assert kept in _tables(value)


def test_a_table_a_search_is_reading_is_not_evicted(corpus_dir, monkeypatch):
    # A search takes the table's name, then resolves compounds and matches structures --
    # seconds during which it holds nothing but a string. Evicting there would fail a
    # query that had already been answered correctly, with a catalog error naming a
    # table the caller has never heard of.
    budget = 150
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=budget,
    ) as value:
        monkeypatch.setattr(
            execute, "_memory_bytes", _stated_bytes(0, 100, 0, 100, 0, 100)
        )
        with value._pivoted_table("inputs.components") as reading:
            with value._pivoted_table("outcomes.products"):
                pass
            # Over budget, and the only candidate is being read: it stays, and the
            # cache stays over its budget until the reader is done.
            assert reading in _tables(value)
            assert list(value._pivoted) == ["inputs.components", "outcomes.products"]
        # Released, so the next materialization can take it.
        with value._pivoted_table("identifiers"):
            pass
        assert "inputs.components" not in value._pivoted
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
            with value._pivoted_table("outcomes.products"):
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
    # search over a level already materialized has nothing to wait for, and making it
    # wait turns one slow first query into a stall for everyone. Asked from a thread
    # rather than here, so a hit that does wait fails this test rather than hanging it.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        with value._pivoted_table("inputs.components"):
            pass
        with _stalled_build(value, monkeypatch):

            def read() -> str | None:
                with value._pivoted_table("inputs.components") as name:
                    return name

            reader, answered = _in_a_thread(read)
            reader.join(timeout=10)
            assert answered, "a cache hit waited on an unrelated build"
            assert answered[0] is not None
        assert not reader.is_alive()


def test_one_pivot_is_built_once_however_many_searches_want_it(corpus_dir, monkeypatch):
    # Two searches over the same level arrive together on a cold cache. Building twice
    # costs two passes over the corpus and leaves two tables of one level, of which only
    # the last is remembered -- the other is memory nothing will free.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        with _stalled_build(value, monkeypatch):

            def read() -> str | None:
                # The stalled build is materializing this very level.
                with value._pivoted_table("outcomes.products") as name:
                    return name

            second, answered = _in_a_thread(read)
            second.join(timeout=1)
            assert not answered, "the second ask did not wait for the build in flight"
        second.join(timeout=30)
        assert answered, "the second ask never got its table"
        assert answered[0] is not None
        assert value._pivot_serial == 1  # One name taken, so one table built.
        assert [name for name in _tables(value) if name.startswith("pivoted_")] == [
            answered[0]
        ]


def test_a_table_a_second_search_took_from_the_cache_is_not_evicted(
    corpus_dir, monkeypatch
):
    # A hit hands back a name the search will read seconds later, exactly as a build
    # does, so it has to take the read the same way. Held at zero readers, the table is
    # the first thing eviction takes, and the search fails on a name that no longer
    # resolves after it had already been answered correctly.
    budget = 150
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=budget,
    ) as value:
        monkeypatch.setattr(
            execute, "_memory_bytes", _stated_bytes(0, 100, 0, 100, 0, 100)
        )
        with value._pivoted_table("inputs.components") as built:
            pass
        # Taken from the cache this time, and read while the budget forces an eviction.
        with value._pivoted_table("inputs.components") as reading:
            assert reading == built
            with value._pivoted_table("outcomes.products"):
                pass
            assert reading in _tables(value)


def test_a_failed_materialization_leaves_nothing_behind(corpus_dir, monkeypatch):
    # A table nobody tracks is memory nobody frees, and under a name a later attempt
    # would collide with. The failure is raised after the CREATE, which is the window.
    request = query.Query.model_validate({"where": _WHITE_PRODUCT})
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
        assert not value._pivoted
        assert not [name for name in _tables(value) if name.startswith("pivoted_")]
        failing = False
        # The level is still materializable, under a name of its own.
        value.search(request)
        assert len(value._pivoted) == 1


def test_a_pivot_too_large_to_keep_is_not_kept(corpus_dir, monkeypatch):
    # Materializing is only worth it if the table survives to answer the next query.
    budget = 1
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=budget,
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
        assert not value._pivoted


def test_a_pivot_too_large_to_keep_says_so_and_says_what_to_change(corpus_dir, caplog):
    # The slowest thing a corpus does is unnest the projection because a pivot did not
    # fit, and nothing about the answer says that is what happened. Whoever is asking
    # why a query takes seconds is reading the log, so the log has to name the level,
    # the two figures, and what changes it.
    caplog.set_level(logging.INFO, logger="ord_schema.search.execute")
    request = query.Query.model_validate({"where": _WHITE_PRODUCT})

    def warnings() -> list[str]:
        return [
            record.getMessage()
            for record in caplog.records
            if record.levelno >= logging.WARNING
            and record.name == "ord_schema.search.execute"
        ]

    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=1,
    ) as value:
        value.search(request)
        assert len(warnings()) == 1
        said = warnings()[0]
        assert "outcomes.products" in said
        assert "pivot_budget_bytes" in said
        assert execute._format_bytes(value._pivot_budget) in said
        # Both figures have to be figures. Stated in gigabytes, a fixture costing
        # kilobytes against a budget of one byte reads "0.0 GB over 0.0 GB", which
        # names the problem in units that hide it.
        assert "1 B" in said
        assert re.search(r"takes \d+(\.\d+)? (kB|MB|GB)", said), said
        # Said again for the next query, since that query is slow for the same reason
        # and its asker was not necessarily watching when the corpus first found out.
        value.search(request)
        assert len(warnings()) == 2


def test_two_searches_wanting_one_refused_pivot_build_it_once(corpus_dir, monkeypatch):
    # The refusal is settled under the build lock as well as before it, so a search
    # that waited out another's build finds the answer rather than repeating a pass
    # over the corpus to reach the same refusal.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=1,
    ) as value:
        with _stalled_build(value, monkeypatch):

            def read() -> str | None:
                with value._pivoted_table("outcomes.products") as name:
                    return name

            second, answered = _in_a_thread(read)
            second.join(timeout=1)
            assert not answered, "the second ask did not wait for the build in flight"
        second.join(timeout=30)
        assert answered == [None]  # Refused, and told so rather than building again.
        assert value._pivot_serial == 1


def test_a_pivot_too_large_to_keep_is_not_built_twice(corpus_dir, monkeypatch):
    # What the refusal costs is a pass over the corpus, and the projection it reads
    # does not change while the corpus is open. Asking again is the common case -- a
    # notebook loop, a server serving one shape of question -- and rebuilding a table
    # only to drop it again holds the build lock against every other search each time.
    budget = 1
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivot_budget_bytes=budget,
    ) as value:
        with value._pivoted_table("inputs.components") as name:
            assert name is None
        built = value._pivot_serial
        with value._pivoted_table("inputs.components") as name:
            assert name is None
        assert value._pivot_serial == built  # No second CREATE to name.


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
            warm=False,
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


def _wide_reactions() -> list[reaction_pb2.Reaction]:
    """Reactions exercising the element shapes the ORD corpus barely holds.

    Measured over the whole corpus, no reaction records a NULL or empty ``outcomes``
    and only 119 record a NULL ``products`` under one, so corpus agreement is weak
    evidence about levels that are absent, empty, or carry a NULL leaf. These say it
    outright. Two outcomes on one reaction matter most: ORD is effectively
    single-outcome, which is exactly why a correlation that dropped the outcome
    ordinal answered identically and looked correct.
    """
    reactions = []

    # Two outcomes. The high yield belongs to the product that is not desired, so a
    # correlation blind to either ordinal pairs them and answers yes.
    split = reaction_pb2.Reaction(reaction_id="ord-wd01")
    first = split.outcomes.add()
    undesired = first.products.add(is_desired_product=False, isolated_color="white")
    undesired.measurements.add(type="YIELD").percentage.value = 90
    second = split.outcomes.add()
    desired = second.products.add(is_desired_product=True, isolated_color="yellow")
    desired.measurements.add(type="YIELD").percentage.value = 10
    reactions.append(split)

    # One outcome, one product, both conditions on the same element.
    together = reaction_pb2.Reaction(reaction_id="ord-wd02")
    product = together.outcomes.add().products.add(
        is_desired_product=True, isolated_color="white"
    )
    product.measurements.add(type="YIELD").percentage.value = 90
    reactions.append(together)

    # An outcome carrying no products at all: the level is present and empty.
    empty = reaction_pb2.Reaction(reaction_id="ord-wd03")
    empty.outcomes.add()
    reactions.append(empty)

    # No outcomes at all: the projection writes NULL rather than an empty list.
    reactions.append(reaction_pb2.Reaction(reaction_id="ord-wd04"))

    # A product carrying no isolated_color, so the leaf a body reads is NULL.
    blank = reaction_pb2.Reaction(reaction_id="ord-wd05")
    blank.outcomes.add().products.add(is_desired_product=True)
    reactions.append(blank)

    # Several products where only the last matches, so a filter that stopped early or
    # read only the first element would miss it.
    several = reaction_pb2.Reaction(reaction_id="ord-wd06")
    outcome = several.outcomes.add()
    outcome.products.add(is_desired_product=False, isolated_color="red")
    outcome.products.add(is_desired_product=False, isolated_color="green")
    outcome.products.add(is_desired_product=True, isolated_color="white")
    reactions.append(several)

    # One outcome holding both, so the *product* ordinal is what separates them. This
    # is the shape of the 942 reactions a product-blind correlation over-returned on
    # ORD, which an outcome-blind one could not see: the corpus is single-outcome.
    sibling = reaction_pb2.Reaction(reaction_id="ord-wd07")
    outcome = sibling.outcomes.add()
    wanted = outcome.products.add(is_desired_product=True, isolated_color="blue")
    wanted.measurements.add(type="YIELD").percentage.value = 10
    other = outcome.products.add(is_desired_product=False, isolated_color="blue")
    other.measurements.add(type="YIELD").percentage.value = 90
    reactions.append(sibling)

    return reactions


@pytest.fixture(scope="module")
def wide_corpus(tmp_path_factory) -> Iterator[execute.Corpus]:
    """A corpus of element shapes, for comparing the pivot route against the level."""
    root = tmp_path_factory.mktemp("wide")
    source = root / "data" / "ord_dataset-wd.parquet"
    source.parent.mkdir(parents=True, exist_ok=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-wd",
            name="test",
            description="test",
            reactions=_wide_reactions(),
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


def _white(op):
    return {
        "op": op,
        "path": "outcomes.products",
        "where": {
            "op": "eq",
            "path": "isolated_color",
            "value": {"literal": "white"},
        },
    }


_DIFFERENTIAL = {
    "exists a white product": _white("exists"),
    "every product is white": _white("forall"),
    "no white product": {"op": "not", "clause": _white("exists")},
    "not every product is white": {"op": "not", "clause": _white("forall")},
    # Two conditions on one element: the pivot answers this from one row, which is the
    # co-membership the nested form gets by binding an element.
    "a desired white product": {
        "op": "exists",
        "path": "outcomes.products",
        "where": {
            "op": "and",
            "clauses": [
                {
                    "op": "eq",
                    "path": "isolated_color",
                    "value": {"literal": "white"},
                },
                {
                    "op": "eq",
                    "path": "is_desired_product",
                    "value": {"literal": True},
                },
            ],
        },
    },
    # The leaf is NULL for ord-wd05, so the comparison is NULL rather than false and
    # both routes have to fold it the same way.
    "a product that is not white": {
        "op": "exists",
        "path": "outcomes.products",
        "where": {
            "op": "ne",
            "path": "isolated_color",
            "value": {"literal": "white"},
        },
    },
    "every measurement is a yield": {
        "op": "forall",
        "path": "outcomes.products.measurements",
        "where": {"op": "eq", "path": "type", "value": {"literal": "YIELD"}},
    },
    # The nested-correlation case, and the one the wide corpus was built for: ord-wd01
    # carries the 90% yield on the product that is *not* desired, in a different
    # outcome from the one that is. A correlation joining on anything short of the
    # whole ordinal prefix pairs them and answers yes.
    "a desired product with a yield above 50%": {
        "op": "exists",
        "path": "outcomes.products",
        "where": {
            "op": "and",
            "clauses": [
                {
                    "op": "eq",
                    "path": "is_desired_product",
                    "value": {"literal": True},
                },
                {
                    "op": "exists",
                    "path": "measurements",
                    "where": {
                        "op": "and",
                        "clauses": [
                            {"op": "eq", "path": "type", "value": {"literal": "YIELD"}},
                            {
                                "op": "gt",
                                "path": "percentage.value",
                                "value": {"literal": 50},
                            },
                        ],
                    },
                },
            ],
        },
    },
    "a yield above 50%": {
        "op": "exists",
        "path": "outcomes.products.measurements",
        "where": {
            "op": "and",
            "clauses": [
                {"op": "eq", "path": "type", "value": {"literal": "YIELD"}},
                {
                    "op": "gt",
                    "path": "percentage.value",
                    "value": {"literal": 50},
                },
            ],
        },
    },
}


@pytest.mark.parametrize("where", list(_DIFFERENTIAL.values()), ids=list(_DIFFERENTIAL))
def test_a_pivot_answers_what_the_level_answers(wide_corpus, monkeypatch, where):
    # The guard on two routes to one answer. A pivot that lost an ordinal, folded a
    # NULL level the other way, or correlated on the wrong prefix shows up here as a
    # set that differs -- which is how the 942-reaction over-match was found.
    pivoted = _search(wide_corpus, where)
    monkeypatch.setattr(
        execute.Corpus, "_pivoted_table", _no_pivoted_table, raising=True
    )
    assert _search(wide_corpus, where) == pivoted


@contextlib.contextmanager
def _no_pivoted_table(self, path: str) -> Iterator[str | None]:
    """Stands in for _pivoted_table so a quantifier compiles over the elements."""
    del self, path  # Unused.
    yield None


def test_the_differential_cases_are_not_all_the_same_answer(wide_corpus):
    # A differential test comparing two routes proves nothing if every case answers
    # with the whole corpus, which is what a fixture that lost its variety would do.
    answers = {
        name: frozenset(_search(wide_corpus, where))
        for name, where in _DIFFERENTIAL.items()
    }
    assert len(set(answers.values())) >= 5
    assert any(answer for answer in answers.values())
    assert any(len(answer) < 6 for answer in answers.values())


def test_a_pivot_row_says_which_element_it_was(wide_corpus):
    # Nothing in phase 1 joins on the ordinals, so the differential test cannot see
    # them; they are what a correlation across levels will need, and a pivot that
    # numbered its elements wrongly would look correct until then. ord-wd01 carries two
    # outcomes of one product each, and ord-wd06 one outcome of three.
    with wide_corpus._pivoted_table("outcomes.products") as name:
        assert name is not None
        rows = wide_corpus._connection.execute(
            f"SELECT reaction_id, outcome_index, product_index, "  # noqa: S608
            f"element.isolated_color FROM {name} "
            "WHERE reaction_id IN ('ord-wd01', 'ord-wd06') "
            "ORDER BY reaction_id, outcome_index, product_index"
        ).fetchall()
    assert rows == [
        ("ord-wd01", 1, 1, "white"),
        ("ord-wd01", 2, 1, "yellow"),
        ("ord-wd06", 1, 1, "red"),
        ("ord-wd06", 1, 2, "green"),
        ("ord-wd06", 1, 3, "white"),
    ]


def test_a_pivot_holds_every_element_including_an_empty_one(wide_corpus):
    # Completeness is what lets a pivot answer forall, so an element whose fields are
    # all NULL still gets a row: ord-wd05 records a product with no isolated_color.
    with wide_corpus._pivoted_table("outcomes.products") as name:
        assert name is not None
        blank = wide_corpus._connection.execute(
            f"SELECT count(*) FROM {name} WHERE reaction_id = 'ord-wd05'"  # noqa: S608
        ).fetchone()
        # A level that is present but empty contributes no rows, and one the source
        # never recorded contributes none either -- which is what makes forall
        # vacuously true of both.
        absent = wide_corpus._connection.execute(
            f"SELECT count(*) FROM {name} "  # noqa: S608
            "WHERE reaction_id IN ('ord-wd03', 'ord-wd04')"
        ).fetchone()
    assert blank == (1,)
    assert absent == (0,)


def _write_pivots(
    root: pathlib.Path,
    levels: Sequence[str],
    into: pathlib.Path | None = None,
) -> pathlib.Path:
    """Derives pivot artifacts for ``levels`` from every projection under ``root``.

    Written under ``<level>/<shard>/`` rather than directly under the level, because
    that is where ``derive_pivots`` puts them: it mirrors the projections' own layout,
    and a helper that flattened it would test a tree nothing produces.

    Args:
        root: Holds the projections to derive from.
        levels: Which repeated levels to write.
        into: Where to write, defaulting to ``<root>/pivots``. A test wanting a tree
            nothing else sees passes its own, since ``root`` is shared by the module.

    Returns:
        The directory written into.
    """
    pivots = into if into is not None else root / "pivots"
    for level in levels:
        for projected in sorted((root / "projections").glob("*.parquet")):
            shard = pivots / level / projected.stem[-2:]
            shard.mkdir(parents=True, exist_ok=True)
            pivot.write_pivot(projected, shard / projected.name, level_path=level)
    return pivots


@pytest.fixture(scope="module")
def wide_root(tmp_path_factory) -> pathlib.Path:
    """The wide corpus on disk, so pivot artifacts can be derived beside it."""
    root = tmp_path_factory.mktemp("wide_artifacts")
    source = root / "data" / "ord_dataset-wd.parquet"
    source.parent.mkdir(parents=True, exist_ok=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-wd",
            name="test",
            description="test",
            reactions=_wide_reactions(),
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


@pytest.mark.parametrize("where", list(_DIFFERENTIAL.values()), ids=list(_DIFFERENTIAL))
def test_a_pivot_artifact_answers_what_building_one_answers(wide_root, where):
    # The artifact is the same rows by another route, so the whole point is that it
    # cannot answer differently from the pivot the executor would have built.
    pivots = _write_pivots(wide_root, ("outcomes.products",))
    projected = str(wide_root / "projections" / "*.parquet")
    structured = str(wide_root / "structures" / "*.parquet")
    with execute.Corpus(projected, structured, resolver={}.__getitem__) as built:
        expected = _search(built, where)
    with execute.Corpus(
        projected, structured, resolver={}.__getitem__, pivots_dir=str(pivots)
    ) as read:
        assert _search(read, where) == expected


def test_a_pivot_artifact_is_read_rather_than_built(wide_root, caplog):
    pivots = _write_pivots(wide_root, ("outcomes.products",))
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
        ) as corpus,
    ):
        _search(corpus, _white("exists"))
    messages = [record.message for record in caplog.records]
    assert any("read 1 pivot artifacts" in message for message in messages)
    assert not any("building the pivot" in message for message in messages)


# The levels the occurrence index's paths range within: three are indexed paths of
# their own, and the authentic standard is one compound per measurement, so its rows
# come from the measurements' pivot.
_OCCURRENCE_LEVELS = (
    "inputs.components",
    "workups.input.components",
    "outcomes.products",
    "outcomes.products.measurements",
)


def _indexed_corpus(root, pivots, **kwargs) -> execute.Corpus:
    return execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        resolver={"pyridine": "c1ccncc1", "ethanol": "CCO"}.__getitem__,
        pivots_dir=None if pivots is None else str(pivots),
        **kwargs,
    )


@pytest.mark.parametrize(
    "smarts",
    ["c1ccncc1", "[OX2H]", "c1ccccc1", "[#6]"],
)
def test_an_index_read_from_pivots_answers_what_unnesting_answers(
    tmp_path, corpus_dir, smarts
):
    # The index is the same occurrences by another route, so the routes cannot disagree.
    # Run over both datasets, because an occurrence's ID is its element's plus its
    # dataset's offset and only a corpus with a second dataset has a nonzero one: a
    # build that dropped the offset answers identically on the first dataset alone.
    pivots = _write_pivots(corpus_dir, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    where = _exists({"op": "substructure", "path": "smiles", "smarts": smarts})
    with _indexed_corpus(corpus_dir, None) as unnested:
        expected = _search(unnested, where)
    # Every pattern here matches something, so neither route is compared against the
    # empty set it would also return had the build produced no rows at all.
    assert expected
    with _indexed_corpus(corpus_dir, pivots) as read:
        assert _search(read, where) == expected


def test_an_index_read_from_pivots_keeps_a_structure_in_its_own_dataset(
    tmp_path, corpus_dir
):
    # The hydroxyl is only in bb, whose IDs are only correct through its offset. An
    # offset dropped or taken from the wrong dataset still lands inside the corpus's
    # structures, so it returns a reaction rather than raising -- the wrong one.
    pivots = _write_pivots(corpus_dir, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    with _indexed_corpus(corpus_dir, pivots) as corpus:
        matched = _search(
            corpus,
            _exists({"op": "substructure", "path": "smiles", "smarts": "[OX2H]"}),
        )
    assert matched == {"ord-bb01"}


def test_an_index_read_from_pivots_keeps_the_element_s_role(tmp_path, corpus_dir):
    # The role travels beside the structure, and it is what binds a query to one
    # element: aa01 holds benzene as a reactant and pyridine as its solvent. A role read
    # wrong from the artifact answers "pyridine as the solvent" with silence and
    # "benzene as the solvent" with aa01, neither of which raises.
    pivots = _write_pivots(corpus_dir, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    with _indexed_corpus(corpus_dir, pivots) as corpus:
        assert _search(corpus, _role_and_structure("c1ccncc1", "SOLVENT")) == {
            "ord-aa01"
        }
        assert _search(corpus, _role_and_structure("c1ccccc1", "SOLVENT")) == set()


# One structure per indexed path, with the reactions holding it. The two deepest paths
# are why this runs over deep_root: the main corpus records no workups and no authentic
# standards, so a differential there compares nothing to nothing at exactly the paths
# whose traversal is longest -- and the authentic standard is the one whose level is not
# the indexed path, reached from the measurements' artifact through the remainder.
_INDEXED_OCCURRENCES = (
    ("inputs.components", "CCO", {"ord-cc01", "ord-cc02"}),
    ("workups.input.components", "c1ccncc1", {"ord-cc01"}),
    ("outcomes.products", "Cc1ccccc1", {"ord-cc01"}),
    ("outcomes.products.measurements.authentic_standard", "c1ccccc1", {"ord-cc01"}),
)


@pytest.mark.parametrize(("path", "smarts", "expected"), _INDEXED_OCCURRENCES)
def test_an_index_read_from_pivots_answers_at_every_indexed_path(
    tmp_path, deep_root, path, smarts, expected
):
    pivots = _write_pivots(deep_root, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    where = _exists_at(path, {"op": "substructure", "path": "smiles", "smarts": smarts})
    with _indexed_corpus(deep_root, None) as unnested:
        assert _search(unnested, where) == expected
    with _indexed_corpus(deep_root, pivots) as read:
        assert _search(read, where) == expected


def test_an_index_reads_the_pivots_rather_than_unnesting_the_projection(
    tmp_path, corpus_dir, caplog
):
    pivots = _write_pivots(corpus_dir, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        _indexed_corpus(corpus_dir, pivots) as corpus,
    ):
        _search(
            corpus, _exists({"op": "substructure", "path": "smiles", "smarts": "[Pt]"})
        )
    assert any(
        f"and {len(execute.INDEXED_PATHS)} from pivot artifacts" in record.message
        for record in caplog.records
    )


def test_a_path_whose_level_has_no_pivot_is_still_unnested(
    tmp_path, corpus_dir, caplog
):
    # A directory holding some levels and not others is a partial speedup rather than a
    # missing index, so the answer has to be the whole answer either way.
    pivots = _write_pivots(corpus_dir, ("inputs.components",), into=tmp_path / "pivots")
    where = _exists({"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"})
    with _indexed_corpus(corpus_dir, None) as unnested:
        expected = _search(unnested, where)
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        _indexed_corpus(corpus_dir, pivots) as corpus,
    ):
        assert _search(corpus, where) == expected
    assert any(
        "and 1 from pivot artifacts" in record.message for record in caplog.records
    )


def test_requiring_pivots_refuses_a_level_with_no_artifacts(wide_root):
    # The silent case: a missing level is otherwise unnested in process on the query
    # that first wants it, minutes charged to whoever asked with no error anywhere.
    pivots = _write_pivots(wide_root, ("outcomes.products",))
    with pytest.raises(execute.PairingError, match="holds no pivot artifacts"):
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
            require_pivots=True,
        )


def test_requiring_pivots_without_a_directory_says_so(wide_root):
    with pytest.raises(ValueError, match="nowhere to look"):
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            require_pivots=True,
        )


def test_requiring_pivots_and_budgeting_for_one_is_refused(wide_root, tmp_path):
    # Requiring every level leaves nothing to build, so a budget beside it is money set
    # aside for something that cannot happen -- more likely a misunderstanding of one
    # of the two than a deliberate pairing.
    pivots = _write_pivots(wide_root, ("outcomes.products",), into=tmp_path / "pivots")
    with pytest.raises(ValueError, match="pass one or the other"):
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
            require_pivots=True,
            pivot_budget_bytes=0,
        )


def test_requiring_pivots_opens_a_corpus_holding_every_level(wide_root, tmp_path):
    # Its own directory: this is the one tree holding every level, and the tests below
    # are about what happens when a level is missing from the shared one.
    pivots = _write_pivots(wide_root, sorted(pivot.LEVELS), into=tmp_path / "pivots")
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
        require_pivots=True,
    ) as corpus:
        assert _search(corpus, _white("exists"))


def test_a_level_without_artifacts_is_still_built(wide_root, caplog):
    # A partial set of artifacts is a partial speedup, not a missing answer.
    pivots = _write_pivots(wide_root, ("outcomes.products",))
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
        ) as corpus,
    ):
        found = _search(
            corpus,
            {
                "op": "exists",
                "path": "outcomes.products.measurements",
                "where": {
                    "op": "eq",
                    "path": "type",
                    "value": {"literal": "YIELD"},
                },
            },
        )
    assert found == {"ord-wd01", "ord-wd02", "ord-wd07"}
    assert any(
        "building the pivot over outcomes.products.measurements" in record.message
        for record in caplog.records
    )


def _messages(caplog) -> list[str]:
    """Returns what the executor logged."""
    return [
        record.message
        for record in caplog.records
        if record.name == "ord_schema.search.execute"
    ]


def test_a_level_with_no_artifacts_is_derived_rather_than_built(
    wide_root, tmp_path, caplog
):
    # The pass costs what building in memory costs, and the difference is what survives
    # it: a file the next process reads, against a table that dies with this one.
    pivots = tmp_path / "pivots"
    # The answer the level itself gives, so a derivation that lost rows shows up below.
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as built:
        expected = _search(built, _white("exists"))
    caplog.clear()
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
            derive_pivots=True,
            warm=False,
        ) as corpus,
    ):
        assert _search(corpus, _white("exists")) == expected
    assert sorted(path.name for path in pivots.glob("outcomes.products/*.parquet")) == [
        "ord_dataset-wd.parquet"
    ]
    said = _messages(caplog)
    assert any("deriving the pivot artifacts for outcomes.products" in m for m in said)
    assert not any("building the pivot" in m for m in said)


def test_a_derived_pivot_is_read_by_the_next_corpus_over_the_directory(
    wide_root, tmp_path, caplog
):
    # Deriving into a directory is only worth more than building in memory if what it
    # wrote is what a later process finds there -- including one that derives nothing.
    pivots = tmp_path / "pivots"
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
        derive_pivots=True,
        warm=False,
    ) as deriving:
        expected = _search(deriving, _white("exists"))
    caplog.clear()
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
        ) as reading,
    ):
        assert _search(reading, _white("exists")) == expected
    said = _messages(caplog)
    assert any("read 1 pivot artifacts" in m for m in said)
    assert not any("deriving" in m or "building the pivot" in m for m in said)


def test_a_level_already_derived_is_not_derived_again(wide_root, tmp_path, caplog):
    # The second query wants the same level, and re-globbing is what tells it there is
    # nothing left to write. A derivation per query would be minutes per query at scale.
    pivots = tmp_path / "pivots"
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
            derive_pivots=True,
            warm=False,
        ) as corpus,
    ):
        _search(corpus, _white("exists"))
        _search(corpus, _white("exists"))
    said = _messages(caplog)
    assert sum("deriving the pivot artifacts" in m for m in said) == 1


def test_a_derivation_interrupted_partway_is_finished_by_the_next(
    corpus_dir, tmp_path, caplog
):
    # A set covering some of the projections is exactly what the pairing check refuses,
    # so reading the level before deriving it would make an interrupted run permanent:
    # the pass that would complete the set is the one behind the refusal.
    pivots = tmp_path / "pivots" / "inputs.components"
    pivots.mkdir(parents=True)
    projections = sorted((corpus_dir / "projections").glob("*.parquet"))
    assert len(projections) > 1, "one projection cannot be a partial set"
    pivot.write_pivot(
        projections[0], pivots / projections[0].name, level_path="inputs.components"
    )
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(tmp_path / "pivots"),
            derive_pivots=True,
            warm=False,
        ) as corpus,
    ):
        found = _search(
            corpus, _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}})
        )
    assert found == {"ord-bb01"}
    assert len(list(pivots.glob("*.parquet"))) == len(projections)
    # The one already there was left alone rather than written again.
    assert any("1 already current" in m for m in _messages(caplog))


def test_a_stranger_s_pivots_are_still_refused_when_deriving(
    wide_root, corpus_dir, tmp_path
):
    # Deriving before reading must not turn the pairing check into a formality: a pass
    # over this corpus's projections adds its own artifacts beside the stranger's rather
    # than displacing them, and answering from the union would be wrong twice over.
    stranger = tmp_path / "pivots" / "outcomes.products"
    stranger.mkdir(parents=True)
    for projected in sorted((corpus_dir / "projections").glob("*.parquet")):
        pivot.write_pivot(
            projected, stranger / projected.name, level_path="outcomes.products"
        )
    with (
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(tmp_path / "pivots"),
            derive_pivots=True,
            warm=False,
        ) as corpus,
        pytest.raises(execute.PairingError, match="another corpus"),
    ):
        _search(corpus, _white("exists"))


def test_checking_the_pivots_derives_none_of_them(wide_root, tmp_path):
    # check_pivots reaches every level the schema has, and deriving there would be
    # 39 unnests of the projection at startup for a deployment that asked about two.
    pivots = tmp_path / "pivots"
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
        derive_pivots=True,
        warm=False,
    ) as corpus:
        assert corpus.check_pivots() == {}
    assert not pivots.exists()


def test_deriving_with_nowhere_to_derive_into_is_refused(wide_root):
    # Silently building in memory instead would answer every query the slow way for a
    # deployment that believed it had asked for artifacts.
    with pytest.raises(ValueError, match="without a pivots_dir"):
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            derive_pivots=True,
        )


def test_a_pivot_from_another_corpus_is_refused(wide_root, corpus_dir, tmp_path):
    # The failure this pins is silent: a pivot names reactions by ID, so artifacts
    # derived elsewhere answer a quantifier against rows this corpus does not hold.
    stranger = _write_pivots(corpus_dir, ("outcomes.products",))
    with (
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(stranger),
            warm=False,
        ) as corpus,
        pytest.raises(execute.PairingError, match="another corpus"),
    ):
        _search(corpus, _white("exists"))


def test_a_stranger_s_pivot_over_an_indexed_level_is_refused_at_open(
    wide_root, corpus_dir
):
    # Warming reads the levels the index covers, so a pivots_dir belonging to another
    # corpus fails the deployment rather than whichever query first quantified over one
    # of those levels. Filed at an indexed level for that reason: the levels the index
    # does not cover are still found by the query that wants them.
    stranger = _write_pivots(corpus_dir, ("inputs.components",))
    with pytest.raises(execute.PairingError, match="another corpus"):
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(stranger),
        )


def test_two_pivots_of_one_dataset_under_a_level_are_refused(corpus_dir, tmp_path):
    # A stray copy under the level satisfies the set comparison over source hashes, and
    # both artifacts are then read. What that costs is in _check_pivots; what makes it
    # worth a refusal here is that the corpus opens and answers, so nothing else says.
    pivots = _write_pivots(corpus_dir, ("inputs.components",), into=tmp_path / "pivots")
    stray = pivots / "inputs.components" / "again"
    stray.mkdir(parents=True)
    for artifact in sorted((pivots / "inputs.components").rglob("*.parquet")):
        shutil.copy(artifact, stray / artifact.name)
    with pytest.raises(execute.PairingError, match="the same source dataset"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
        )


def test_a_stale_pivot_is_refused(corpus_dir, tmp_path):
    # Warming reads the levels the index covers, so what a stale artifact fails is the
    # deployment. The stamps are the only thing saying an artifact was written by the
    # library reading it; a pivot derived by an older one can hold columns this schema
    # walk does not find.
    pivots = _write_pivots(corpus_dir, ("inputs.components",), into=tmp_path / "pivots")
    target = next((pivots / "inputs.components").rglob("*.parquet"))
    with pq.ParquetFile(target) as artifact:
        table = artifact.read()
        schema = artifact.schema_arrow
    metadata = dict(schema.metadata)
    metadata[b"ord.rdkit_version"] = b"0000.00.0"
    pq.write_table(table.cast(schema.with_metadata(metadata)), target)
    with pytest.raises(execute.PairingError, match="is a stale"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
        )


def test_a_pivot_filed_under_the_wrong_level_is_refused(wide_root, tmp_path):
    pivots = tmp_path / "pivots"
    (pivots / "outcomes.products").mkdir(parents=True)
    projected = next((wide_root / "projections").glob("*.parquet"))
    pivot.write_pivot(
        projected,
        pivots / "outcomes.products" / projected.name,
        level_path="outcomes.products.measurements",
    )
    with (
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(pivots),
            warm=False,
        ) as corpus,
        pytest.raises(execute.PairingError, match="wrong level"),
    ):
        _search(corpus, _white("exists"))


def test_a_structure_predicate_inside_a_pivoted_quantifier_agrees(corpus, monkeypatch):
    # The occurrence index declines this: its shape is one structure predicate and
    # string equalities, and a not_null clause is neither. The pivot takes it instead,
    # and its rows carry a structure_id but no structure_offset -- that column resolves
    # outward to the reaction the element belongs to, which is the offset its ID is
    # meant to be read against. Correct by name resolution rather than by construction,
    # so it is compared against the elements rather than assumed.
    where = {
        "op": "exists",
        "path": "inputs.components",
        "where": {
            "op": "and",
            "clauses": [
                {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
                {"op": "not_null", "path": "smiles"},
            ],
        },
    }
    pivoted = _search(corpus, where)
    monkeypatch.setattr(execute.Corpus, "_pivoted_table", _no_pivoted_table)
    assert _search(corpus, where) == pivoted
    # And the answer is the pyridine-bearing reactions, not everything or nothing.
    assert pivoted == {"ord-aa01", "ord-aa02"}


def test_a_pivoted_quantifier_still_reaches_the_structures_artifact(corpus):
    # Guards the offset arithmetic across datasets: ethanol lives in the second file,
    # so finding it through a pivot proves the outward structure_offset is the right
    # one rather than the first file's.
    found = _search(
        corpus,
        {
            "op": "exists",
            "path": "inputs.components",
            "where": {
                "op": "and",
                "clauses": [
                    {"op": "substructure", "path": "smiles", "smarts": "[OX2H]"},
                    {"op": "not_null", "path": "smiles"},
                ],
            },
        },
    )
    assert found == {"ord-bb01"}


def test_a_pivot_under_a_directory_that_looks_like_a_pattern_is_read(
    wide_root, tmp_path
):
    # read_parquet reads its arguments as globs, so a bracket in a directory name is a
    # character class, and it matches the sibling this decoy provides rather than the
    # artifact itself. The decoy holds no rows, so a query reading it answers nothing
    # and says nothing about having done so.
    pivots = _write_pivots(
        wide_root, ("outcomes.products",), into=tmp_path / "pivots[v2]"
    )
    real = next(iter(pivots.rglob("*.parquet")))
    decoy = tmp_path / "pivotsv" / "outcomes.products" / real.parent.name
    decoy.mkdir(parents=True)
    empty = pq.read_table(real)
    pq.write_table(
        empty.slice(0, 0).replace_schema_metadata(empty.schema.metadata),
        decoy / real.name,
    )
    assert pq.read_table(real).num_rows > 0

    where = {
        "op": "exists",
        "path": "outcomes.products",
        "where": {"op": "not_null", "path": "is_desired_product"},
    }
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
    ) as corpus:
        found = corpus.search(query.Query.model_validate({"where": where}))
    assert found.num_rows > 0


def test_check_pivots_reports_the_levels_held_as_artifacts(corpus_dir, tmp_path):
    # An operator reads these counts to see that a sharded corpus has all its shards, so
    # this runs over the two-dataset corpus: a count of one, or one taken from the
    # levels rather than the files, would agree with the answer over a single shard. A
    # level with no artifacts is absent rather than zero, since most levels have none.
    pivots = _write_pivots(corpus_dir, ("outcomes.products",), into=tmp_path)
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
    ) as corpus:
        assert corpus.check_pivots() == {"outcomes.products": 2}


def test_check_pivots_is_empty_without_a_directory(wide_corpus):
    # Not an error: every level is built in process, which is the phase the executor
    # started in and still supports.
    assert wide_corpus.check_pivots() == {}


def test_check_pivots_refuses_a_stranger_at_startup(wide_root, corpus_dir):
    # The whole point: this refusal would otherwise arrive with whichever query first
    # wanted the level, at whatever hour that was.
    stranger = _write_pivots(corpus_dir, ("outcomes.products",))
    with (
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(stranger),
            warm=False,
        ) as corpus,
        pytest.raises(execute.PairingError, match="another corpus"),
    ):
        corpus.check_pivots()


def test_check_pivots_leaves_the_levels_ready(wide_root, caplog):
    # Publishing is the check, so a level it accepted is one the first query does not
    # have to glob for or build.
    pivots = _write_pivots(wide_root, ("outcomes.products",))
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
    ) as corpus:
        corpus.check_pivots()
        # Cleared, because check_pivots logs the publish it just did and the question
        # here is what the *query* had left to do.
        caplog.clear()
        with caplog.at_level(logging.INFO, logger="ord_schema.search.execute"):
            found = _search(corpus, _white("exists"))
    assert found == {"ord-wd01", "ord-wd02", "ord-wd06"}
    messages = [record.message for record in caplog.records]
    assert not any("read 1 pivot artifacts" in message for message in messages)
    assert not any("building the pivot" in message for message in messages)


def test_a_singular_struct_under_a_level_is_answered_by_that_level_s_pivot(
    deep_corpus, monkeypatch
):
    # An authentic standard is one compound per measurement rather than a list of its
    # own, so no pivot is named for it; the measurements' pivot carries it. This is the
    # one occurrence-indexed path that is not a repeated level, so it is also what an
    # index the pivot was meant to subsume would have quietly stopped covering.
    where = {
        "op": "exists",
        "path": "outcomes.products.measurements.authentic_standard",
        "where": {
            "op": "and",
            "clauses": [
                {"op": "substructure", "path": "smiles", "smarts": "c1ccccc1"},
                {"op": "not_null", "path": "smiles"},
            ],
        },
    }
    pivoted = _search(deep_corpus, where)
    assert pivoted == {"ord-cc01"}
    monkeypatch.setattr(execute.Corpus, "_pivoted_table", _no_pivoted_table)
    assert _search(deep_corpus, where) == pivoted


def test_a_nested_correlation_binds_the_measurement_to_its_own_product(wide_corpus):
    # ord-wd01 has a 90% yield on an undesired product in one outcome and a desired
    # product with 10% in another; ord-wd02 has both on one product. Only the second
    # answers, and a correlation that lost either ordinal would return both.
    found = _search(
        wide_corpus, _DIFFERENTIAL["a desired product with a yield above 50%"]
    )
    # ord-wd07 keeps both in one outcome, so only the product ordinal separates them.
    assert found == {"ord-wd02"}


@pytest.fixture(scope="module")
def drawn_two_ways(tmp_path_factory) -> Iterator[execute.Corpus]:
    """A corpus storing acetate, an enol, and a salt, none spelled the query's way."""
    acetate = reaction_pb2.Reaction(reaction_id="ord-dd01")
    _component(acetate.inputs["in"], "CC(=O)[O-]", _ROLE.REAGENT)
    enol = reaction_pb2.Reaction(reaction_id="ord-dd02")
    _component(enol.inputs["in"], "CC(O)=C", _ROLE.REACTANT)
    salt = reaction_pb2.Reaction(reaction_id="ord-dd03")
    _component(salt.inputs["in"], "CC(=O)[O-].[Na+]", _ROLE.REAGENT)
    other = reaction_pb2.Reaction(reaction_id="ord-dd04")
    _component(other.inputs["in"], "CCO", _ROLE.SOLVENT)
    root = tmp_path_factory.mktemp("drawn")
    source = root / "data" / "ord_dataset-dd.parquet"
    source.parent.mkdir(parents=True, exist_ok=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-dd",
            name="test",
            description="test",
            reactions=[acetate, enol, salt, other],
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
        str(projected),
        str(structured),
        resolver={"acetic_acid": "CC(=O)O"}.__getitem__,
    ) as value:
        yield value


def test_same_compound_matches_another_protonation_state(drawn_two_ways):
    # The corpus stores the acetate; an eq on the spelling answers nothing and says so
    # by returning no rows, which is the worst way for a query to be wrong.
    assert (
        _search(
            drawn_two_ways,
            _exists(
                {"op": "eq", "path": "smiles", "value": {"literal": "CC(=O)O"}},
            ),
        )
        == set()
    )
    assert _search(
        drawn_two_ways,
        _exists({"op": "same_compound", "path": "smiles", "smiles": "CC(=O)O"}),
    ) == {"ord-dd01"}


def test_same_compound_matches_another_tautomer(drawn_two_ways):
    assert _search(
        drawn_two_ways,
        _exists({"op": "same_compound", "path": "smiles", "smiles": "CC(=O)C"}),
    ) == {"ord-dd02"}


def test_same_compound_does_not_reach_a_salt(drawn_two_ways):
    # Sodium acetate is a different reagent from acetic acid rather than a different
    # drawing of it, and dd03 is the reaction that would come back if it were not.
    matched = _search(
        drawn_two_ways,
        _exists({"op": "same_compound", "path": "smiles", "smiles": "CC(=O)O"}),
    )
    assert "ord-dd03" not in matched


def test_same_compound_does_not_match_a_different_molecule(drawn_two_ways):
    assert (
        _search(
            drawn_two_ways,
            _exists({"op": "same_compound", "path": "smiles", "smiles": "c1ccccc1"}),
        )
        == set()
    )


def test_same_parent_reaches_the_salt_that_same_compound_does_not(drawn_two_ways):
    # dd03 is sodium acetate. It is a different reagent from acetic acid and the same
    # parent, which is the whole difference between the two operators.
    assert _search(
        drawn_two_ways,
        _exists({"op": "same_parent", "path": "smiles", "smiles": "CC(=O)O"}),
    ) == {"ord-dd01", "ord-dd03"}


def test_same_parent_still_does_not_match_a_different_compound(drawn_two_ways):
    assert (
        _search(
            drawn_two_ways,
            _exists({"op": "same_parent", "path": "smiles", "smiles": "c1ccccc1"}),
        )
        == set()
    )


def test_a_compound_name_resolves_to_a_same_compound_query(drawn_two_ways):
    assert _search(
        drawn_two_ways,
        _exists({"op": "same_compound", "path": "smiles", "compound": "acetic_acid"}),
    ) == {"ord-dd01"}


def _without_mol_hash(corpus_dir, tmp_path, *, files: int) -> pathlib.Path:
    """Copies a corpus and drops ``mol_hash`` from some of its structures artifacts.

    The stamps are carried over untouched, which is the whole point: they say nothing
    about a file's columns, so a file predating one still reads as current.

    Args:
        corpus_dir: The corpus to copy.
        tmp_path: Where to put the copy.
        files: How many structures artifacts to drop the column from.

    Returns:
        The copy's root.
    """
    root = tmp_path / "aged"
    shutil.copytree(corpus_dir, root)
    for path in sorted((root / "structures").glob("*.parquet"))[:files]:
        table = pq.read_table(path)
        pq.write_table(
            table.drop_columns(["mol_hash"]).replace_schema_metadata(
                table.schema.metadata
            ),
            path,
        )
    return root


def test_a_structures_artifact_without_mol_hash_is_refused(corpus_dir, tmp_path):
    # Stale is not what this is: the stamps still match, so nothing rebuilds the file
    # and nothing else would notice until DuckDB failed to bind the column.
    root = _without_mol_hash(corpus_dir, tmp_path, files=2)
    assert structures.is_current(
        next((root / "structures").glob("*.parquet")),
        base.load_stamps(next((root / "structures").glob("*.parquet"))).source_md5,
    )
    with pytest.raises(execute.PairingError, match="mol_hash"):
        _open(root)


def test_one_aged_structures_artifact_is_refused_with_the_others_current(
    corpus_dir, tmp_path
):
    # The likelier state, and the worse one: DuckDB reads a glob of mismatched schemas
    # as an error on every structure query rather than on the one that wants the column.
    root = _without_mol_hash(corpus_dir, tmp_path, files=1)
    with pytest.raises(execute.PairingError, match="derive it again"):
        _open(root)


def test_a_corpus_reads_each_artifacts_stamps_once(corpus_dir, monkeypatch):
    # Pairing opens a footer per artifact to verify it, and the fingerprint is taken
    # from those stamps. Reading them a second time costs 0.65 seconds on the real
    # corpus for values already in hand.
    opened: list[str] = []
    real = base.load_stamps

    def counting(path):
        opened.append(str(path))
        return real(path)

    monkeypatch.setattr(base, "load_stamps", counting)
    projections = str(corpus_dir / "projections" / "*.parquet")
    structures_glob = str(corpus_dir / "structures" / "*.parquet")
    with execute.Corpus(projections, structures_glob) as value:
        assert value.fingerprint
        assert value.fingerprint  # A second read opens nothing either.
    assert len(opened) == len(set(opened))


def test_a_corpus_fingerprint_names_the_artifacts_it_opened(corpus_dir):
    # The fingerprint travels in a question log so an old record can be reproduced, so
    # it has to be the same for two openings of one corpus and different for a corpus
    # holding different artifacts.
    both = str(corpus_dir / "projections" / "*.parquet")
    one = str(corpus_dir / "projections" / "ord_dataset-aa.parquet")
    structures_glob = str(corpus_dir / "structures" / "*.parquet")
    with (
        execute.Corpus(both, structures_glob) as first,
        execute.Corpus(both, structures_glob) as second,
        execute.Corpus(
            one, str(corpus_dir / "structures" / "ord_dataset-aa.parquet")
        ) as subset,
    ):
        assert first.fingerprint == second.fingerprint
        assert first.fingerprint != subset.fingerprint


def test_a_corpus_fingerprint_moves_when_only_the_structures_are_rebuilt(
    corpus_dir, tmp_path
):
    # Structures answer every structure predicate and are derived separately from the
    # projections they pair with, so a corpus can be rebuilt with its projections
    # untouched and still answer differently. A fingerprint blind to them says two
    # corpora agree when they do not, and a question log promising to reproduce an old
    # query from one would quietly be wrong. The stamp changed here is one currency
    # does not check, so the pair still opens.
    root = tmp_path / "restructured"
    shutil.copytree(corpus_dir, root)
    for path in (root / "structures").glob("*.parquet"):
        table = pq.read_table(path)
        stamps = base.load_stamps(path)
        pq.write_table(
            table.replace_schema_metadata(
                base.to_metadata(
                    dataclasses.replace(stamps, source_dataset_id="ord_dataset-other")
                )
            ),
            path,
        )
    with (
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
        ) as before,
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
        ) as after,
    ):
        assert before.fingerprint != after.fingerprint


# Reading the occurrence index from artifacts


def _write_occurrences(
    pivots: pathlib.Path,
    into: pathlib.Path,
    paths: Sequence[str] = tuple(execute.INDEXED_PATHS),
) -> pathlib.Path:
    """Derives occurrence artifacts for ``paths`` from the pivots under ``pivots``.

    Written under ``<path>/<shard>/`` because that is where ``derive_occurrences`` puts
    them, mirroring the pivots it descends from.

    Args:
        pivots: Holds the pivot artifacts to derive from.
        into: Where to write.
        paths: Which indexed paths to write, every one by default.

    Returns:
        ``into``.
    """
    for path in paths:
        level = occurrences.PATHS[path][0].path
        for pivoted in sorted((pivots / level).glob("*/*.parquet")):
            shard = into / path / pivoted.parent.name
            shard.mkdir(parents=True, exist_ok=True)
            occurrences.write_occurrences(pivoted, shard / pivoted.name, path=path)
    return into


@pytest.fixture
def occurrence_dirs(tmp_path, corpus_dir) -> tuple[pathlib.Path, pathlib.Path]:
    """The pivots the occurrences descend from, and the occurrences themselves."""
    pivots = _write_pivots(corpus_dir, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    return pivots, _write_occurrences(pivots, tmp_path / "occurrences")


@pytest.mark.parametrize("smarts", ["c1ccncc1", "[OX2H]", "c1ccccc1", "[#6]"])
def test_an_index_read_from_artifacts_answers_what_unnesting_answers(
    corpus_dir, occurrence_dirs, smarts
):
    # The artifact is the same occurrences by a third route, so it cannot answer
    # differently from either of the two that build them. Over both datasets, because
    # an occurrence's ID is its element's plus its own dataset's offset and a corpus of
    # one dataset has an offset of zero to get wrong.
    _, occurrence_dir = occurrence_dirs
    where = _exists({"op": "substructure", "path": "smiles", "smarts": smarts})
    with _indexed_corpus(corpus_dir, None) as unnested:
        expected = _search(unnested, where)
    assert expected
    with _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read:
        assert _search(read, where) == expected


def test_an_index_read_from_artifacts_binds_the_role_to_the_element(
    corpus_dir, occurrence_dirs
):
    # The role travels on the occurrence row, which is what keeps "pyridine as the
    # solvent" a condition on one row rather than an intersection of two reaction sets.
    _, occurrence_dir = occurrence_dirs
    where = _role_and_structure("c1ccncc1", "SOLVENT")
    with _indexed_corpus(corpus_dir, None) as unnested:
        expected = _search(unnested, where)
    with _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read:
        assert _search(read, where) == expected


def test_an_index_read_from_artifacts_keeps_a_structure_in_its_own_dataset(
    corpus_dir, occurrence_dirs
):
    # An offset dropped or taken from the wrong dataset still lands inside the corpus's
    # structures, so it returns a reaction rather than raising -- the wrong one.
    _, occurrence_dir = occurrence_dirs
    with (
        _indexed_corpus(corpus_dir, None) as unnested,
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        where = _exists({"op": "substructure", "path": "smiles", "smarts": "[OX2H]"})
        expected = _search(unnested, where)
        assert expected
        assert _search(read, where) == expected


def test_a_covered_corpus_publishes_the_index_as_a_view(
    corpus_dir, occurrence_dirs, caplog
):
    # The whole point of the artifact: the rows stay on disk. Over ORD the table holds
    # 1.19 GB and wants several more to build, and the view holds none of it.
    _, occurrence_dir = occurrence_dirs
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        read.check_index()
    assert any("as a view" in record.message for record in caplog.records)
    assert any(
        f"{len(execute.INDEXED_PATHS)} of {len(execute.INDEXED_PATHS)} paths read "
        "from occurrence artifacts" in record.message
        for record in caplog.records
    )


def test_one_uncovered_path_materializes_the_whole_index(
    corpus_dir, occurrence_dirs, caplog
):
    # A view whose branches unnest the projection would repeat that traversal on every
    # query rather than pay it once, so partial coverage is a table.
    pivots, _ = occurrence_dirs
    partial = _write_occurrences(
        pivots, pathlib.Path(str(pivots) + "-partial"), ("outcomes.products",)
    )
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(partial)) as read,
    ):
        read.check_index()
    assert any("as a table" in record.message for record in caplog.records)
    assert any(
        f"1 of {len(execute.INDEXED_PATHS)} paths read from occurrence artifacts"
        in record.message
        for record in caplog.records
    )


def test_an_uncovered_path_still_answers(corpus_dir, occurrence_dirs):
    # Partial coverage is a partial speedup, never a missing answer: the paths with no
    # artifact are read the way they were before there were any.
    pivots, _ = occurrence_dirs
    partial = _write_occurrences(
        pivots, pathlib.Path(str(pivots) + "-one"), ("outcomes.products",)
    )
    where = _exists({"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"})
    with _indexed_corpus(corpus_dir, None) as unnested:
        expected = _search(unnested, where)
    assert expected
    with _indexed_corpus(corpus_dir, None, occurrences_dir=str(partial)) as read:
        assert _search(read, where) == expected


def test_a_directory_with_no_artifacts_is_not_an_error(corpus_dir, tmp_path):
    # Naming an empty directory reads as a corpus that has not derived them yet, which
    # is the state every corpus starts in.
    empty = tmp_path / "none"
    empty.mkdir()
    with _indexed_corpus(corpus_dir, None, occurrences_dir=str(empty)) as read:
        assert sum(read.check_index().values()) > 0


# What the reader refuses


def test_occurrence_artifacts_of_another_corpus_are_refused(
    tmp_path, corpus_dir, occurrence_dirs
):
    # An occurrence names its reaction by ID and its structure by a dataset-local one,
    # and both resolve against whatever this corpus holds -- so artifacts from another
    # answer a quantifier confidently and wrongly rather than failing.
    _, occurrence_dir = occurrence_dirs
    one = tmp_path / "one"
    for path in execute.INDEXED_PATHS:
        found = occurrences.artifact_paths(occurrence_dir, path)
        target = one / path / found[0].parent.name
        target.mkdir(parents=True)
        shutil.copy(found[0], target / found[0].name)
    with pytest.raises(execute.PairingError, match="another corpus names reactions"):
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(one)).check_index()


def test_an_artifact_stamped_with_another_path_is_refused(
    tmp_path, corpus_dir, occurrence_dirs
):
    # Every indexed path writes the same three columns, so where a file sits says
    # nothing about what it holds and the stamped path is the only thing that does.
    _, occurrence_dir = occurrence_dirs
    misfiled = tmp_path / "misfiled"
    shutil.copytree(occurrence_dir, misfiled)
    stranger = occurrences.artifact_paths(misfiled, "inputs.components")[0]
    into = occurrences.artifact_paths(misfiled, "outcomes.products")[0]
    shutil.copy(stranger, into)
    with pytest.raises(execute.PairingError, match="holds the occurrences at"):
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(misfiled)).check_index()


def test_two_artifacts_of_one_source_dataset_are_refused(
    tmp_path, corpus_dir, occurrence_dirs
):
    # Both are read, so every occurrence at the path is stated twice. A semi-join
    # returns a reaction once however many rows name it, and the reached-structure
    # count is over distinct IDs, so nothing downstream can see it -- but check_index
    # reports the doubled count as the corpus's own.
    _, occurrence_dir = occurrence_dirs
    doubled = tmp_path / "doubled"
    shutil.copytree(occurrence_dir, doubled)
    found = occurrences.artifact_paths(doubled, "inputs.components")[0]
    shutil.copy(found, found.parent / f"copy-{found.name}")
    with pytest.raises(execute.PairingError, match="of the same source dataset"):
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(doubled)).check_index()


def test_a_file_that_is_not_an_occurrences_artifact_is_refused(
    tmp_path, corpus_dir, occurrence_dirs
):
    _, occurrence_dir = occurrence_dirs
    mixed = tmp_path / "mixed"
    shutil.copytree(occurrence_dir, mixed)
    projected = sorted((corpus_dir / "projections").glob("*.parquet"))[0]
    shutil.copy(projected, mixed / "inputs.components" / projected.name)
    with pytest.raises(execute.PairingError, match="is a projection, not an"):
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(mixed)).check_index()


def test_a_stale_occurrence_artifact_is_refused(tmp_path, corpus_dir, occurrence_dirs):
    # The stamps are the only thing saying an artifact was written by the library
    # reading it, and an RDKit that has changed canonicalization changes which
    # structure a row names.
    _, occurrence_dir = occurrence_dirs
    staled = tmp_path / "stale"
    shutil.copytree(occurrence_dir, staled)
    found = occurrences.artifact_paths(staled, "inputs.components")[0]
    table = pq.read_table(found)
    metadata = dict(table.schema.metadata)
    metadata[b"ord.rdkit_version"] = b"0000.00.0"
    pq.write_table(table.replace_schema_metadata(metadata), found)
    with pytest.raises(execute.PairingError, match="is a stale occurrences artifact"):
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(staled)).check_index()
    # Read anyway where the caller has said staleness is theirs to judge, which is the
    # same latitude require_current gives every other artifact.
    with _indexed_corpus(
        corpus_dir, None, occurrences_dir=str(staled), require_current=False
    ) as read:
        assert sum(read.check_index().values()) > 0


def test_requiring_occurrences_refuses_a_partial_directory(
    corpus_dir, occurrence_dirs, tmp_path
):
    # A deployment sized for the view is not sized for the table, so the fallback that
    # keeps a corpus answering is the one that gets the container killed.
    pivots, _ = occurrence_dirs
    partial = _write_occurrences(pivots, tmp_path / "partial", ("outcomes.products",))
    with pytest.raises(execute.PairingError, match="holds no occurrences artifacts"):
        _indexed_corpus(
            corpus_dir, None, occurrences_dir=str(partial), require_occurrences=True
        )


def test_requiring_occurrences_accepts_a_covered_directory(corpus_dir, occurrence_dirs):
    _, occurrence_dir = occurrence_dirs
    with _indexed_corpus(
        corpus_dir,
        None,
        occurrences_dir=str(occurrence_dir),
        require_occurrences=True,
    ) as read:
        assert sum(read.check_index().values()) > 0


def test_requiring_occurrences_without_a_directory_is_refused(corpus_dir):
    with pytest.raises(ValueError, match="require_occurrences was set without"):
        _indexed_corpus(corpus_dir, None, require_occurrences=True)


# What the timeout covers


def test_the_query_is_given_what_the_earlier_phases_left(corpus, monkeypatch):
    # The bound is over the whole call, so a phase that spends half of it leaves the
    # query the other half. Passing the full bound down would let a search take twice
    # what it was allowed and report no error.
    seen: list[float] = []
    original = execute._run_with_timeout

    def recording(cursor, sql, parameters, timeout_seconds):
        seen.append(timeout_seconds)
        return original(cursor, sql, parameters, timeout_seconds)

    monkeypatch.setattr(execute, "_run_with_timeout", recording)
    slow = execute.Corpus._matches

    def dawdle(self, cursor, parameter, resolve):
        time.sleep(0.2)
        return slow(self, cursor, parameter, resolve)

    monkeypatch.setattr(execute.Corpus, "_matches", dawdle)
    request = query.Query.model_validate(
        {
            "where": _exists(
                {"op": "substructure", "path": "smiles", "compound": "pyridine"}
            )
        }
    )
    corpus.search(request, timeout_seconds=10)
    (left,) = seen
    assert 0 < left < 10 - 0.2


def test_a_phase_that_outlasts_the_bound_raises_naming_itself(corpus, monkeypatch):
    # The phases that cannot be interrupted -- an external resolver, RDKit with the GIL
    # released, a build another search is waiting on -- are checked as they finish. The
    # error says which one, because "the search timed out" over a corpus that answers in
    # milliseconds is a question rather than an answer.
    def dawdle(self, cursor, parameter, resolve):
        time.sleep(0.05)
        return self._bitmap([])

    monkeypatch.setattr(execute.Corpus, "_matches", dawdle)
    request = query.Query.model_validate(
        {
            "where": _exists(
                {"op": "substructure", "path": "smiles", "compound": "pyridine"}
            )
        }
    )
    # Named by the pattern the caller asked for, not by the compiler's placeholder.
    with pytest.raises(TimeoutError, match="matching substructure 'pyridine'"):
        corpus.search(request, timeout_seconds=0.01)


def test_no_timeout_leaves_every_phase_unbounded(corpus, monkeypatch):
    # None has to stay a real "no bound" rather than a very large one: a phase that
    # checked a deadline of None would raise on a corpus that is merely slow. The
    # contrast is the assertion -- the same slow phase under a bound does raise, so a
    # passing unbounded search says something about None rather than about the fixture.
    def dawdle(self, cursor, parameter, resolve):
        time.sleep(0.05)
        return self._bitmap([])

    monkeypatch.setattr(execute.Corpus, "_matches", dawdle)
    request = query.Query.model_validate(
        {
            "where": _exists(
                {"op": "substructure", "path": "smiles", "compound": "pyridine"}
            )
        }
    )
    with pytest.raises(TimeoutError):
        corpus.search(request, timeout_seconds=0.01)
    assert corpus.search(request).num_rows == 0


@pytest.mark.parametrize("spent", [0, -0.5])
def test_a_bound_spent_before_the_query_starts_does_not_run_it(corpus, spent):
    # A timer armed for a bound already spent fires while the cursor is idle, and DuckDB
    # clears an interrupt that arrived with nothing running -- so the query would go on
    # to run unbounded and answer after the deadline it was given. Refused instead, from
    # the same value that would have armed the timer, so nothing can expire in between.
    #
    # The SQL names a table that does not exist: if it ran at all, what comes back is a
    # DuckDB binder error rather than the TimeoutError this asserts.
    cursor = corpus._connection.cursor()
    try:
        with pytest.raises(TimeoutError, match="whole bound before the query"):
            execute._run_with_timeout(cursor, "SELECT * FROM no_such_table", {}, spent)
    finally:
        cursor.close()


class _LostTimer:
    """A timer that never fires, standing in for an interrupt DuckDB dropped."""

    def __init__(self, seconds: float, function) -> None:
        """Accepts and discards what ``threading.Timer`` would act on.

        Args:
            seconds: Ignored.
            function: Ignored; never called, which is the point.
        """
        del seconds, function

    def start(self) -> None:
        """Does nothing."""

    def cancel(self) -> None:
        """Does nothing."""


def test_an_answer_that_arrives_past_the_bound_is_a_timeout(corpus, monkeypatch):
    # The interrupt is not guaranteed to land: a timer armed for a very short bound can
    # fire while the cursor is still idle, and DuckDB clears an interrupt that arrived
    # with nothing running -- measured, a bound of a nanosecond loses it about one run
    # in five, after which the query runs unbounded. So the answer is checked against
    # the clock rather than trusted because no interrupt arrived.
    #
    # The lost interrupt is simulated rather than raced for: a timer that never fires
    # is what the caller sees either way, and a test that waited for the real race
    # would pass four times in five while the bug was present.
    monkeypatch.setattr(execute.threading, "Timer", _LostTimer)
    cursor = corpus._connection.cursor()
    try:
        with pytest.raises(TimeoutError, match=r"past the .* it was given"):
            execute._run_with_timeout(
                cursor,
                "SELECT count(*) FROM range(30_000_000) t(i) WHERE i % 7 = 0",
                {},
                0.001,
            )
    finally:
        cursor.close()


def test_a_corpus_bounds_what_it_returns_without_being_asked(corpus_dir):
    # The grammar leaves limit optional, so an unbounded default builds every matching
    # reaction into an Arrow table in whatever process asked -- millions of them at
    # corpus scale, from a query that reads as ordinary.
    assert execute.DEFAULT_MAX_ROWS == 1000
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        assert value._max_rows == execute.DEFAULT_MAX_ROWS


def test_a_caller_computing_over_every_match_can_still_ask_for_one(corpus_dir):
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        max_rows=None,
    ) as value:
        answer = value.search(query.Query.model_validate({}))
    assert answer.num_rows == 3
    assert execute.TRUNCATED.encode() not in (answer.schema.metadata or {})


def test_an_answer_the_bound_may_have_cut_says_so_on_the_table(corpus_dir):
    # On the table rather than logged alone: whoever reads the answer is often not
    # whoever reads the log, and a summary saying "2 reactions" is wrong in a way
    # nothing in the rows shows.
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        max_rows=2,
    ) as value:
        answer = value.search(query.Query.model_validate({}))
    assert answer.num_rows == 2
    assert (answer.schema.metadata or {})[execute.TRUNCATED.encode()] == b"true"


def test_an_answer_the_bound_did_not_reach_carries_no_stamp(corpus_dir):
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        max_rows=100,
    ) as value:
        answer = value.search(query.Query.model_validate({}))
    assert answer.num_rows == 3
    assert execute.TRUNCATED.encode() not in (answer.schema.metadata or {})
