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
import datetime
import gc
import logging
import pathlib
import re
import shutil
import threading
import time
from collections.abc import Iterable, Iterator, Sequence

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary

from ord_schema import message_helpers, parquet
from ord_schema.artifacts import (
    base,
    first_commits,
    history,
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
    """Builds every artifact for both datasets: a whole tree, as a deployment reads.

    Pivots over every level and occurrences at every indexed path sit beside the
    projections and structures, so a corpus opened over the tree routes a quantifier
    as a published build would; a test wanting the projection's lists opens without
    them.
    """
    root = tmp_path_factory.mktemp("corpus")
    for shard, reactions in _REACTIONS.items():
        source = root / "data" / f"ord_dataset-{shard}.parquet"
        source.parent.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            dataset_pb2.Dataset(
                dataset_id=f"ord_dataset-{shard}",
                name=f"dataset {shard}",
                description=f"the {shard} reactions",
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
    _write_pivots(root, sorted(pivot.LEVELS))
    _write_occurrences(root / "pivots", root / "occurrences")
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
    _write_pivots(root, _OCCURRENCE_LEVELS)
    _write_occurrences(root / "pivots", root / "occurrences")
    return root


@pytest.fixture(scope="module")
def deep_corpus(deep_root) -> Iterator[execute.Corpus]:
    """The corpus over ``deep_root``, with pivots over the indexed levels alone."""
    with _open(deep_root, resolver={}.__getitem__) as value:
        yield value


@pytest.fixture(scope="module")
def corpus(corpus_dir) -> Iterator[execute.Corpus]:
    with _open(
        corpus_dir, resolver={"pyridine": "c1ccncc1", "ethanol": "CCO"}.__getitem__
    ) as value:
        yield value


def _open(root, **kwargs) -> execute.Corpus:
    """Opens the corpus under ``root``, with its pivots and occurrences unless told.

    Args:
        root: The tree, laid out as ``corpus_dir`` lays one out.
        **kwargs: Passed to ``Corpus``. ``pivots_dir`` and ``occurrences_dir`` default
            to the tree's own; pass None to open without them.
    """
    kwargs.setdefault("pivots_dir", str(root / "pivots"))
    kwargs.setdefault("occurrences_dir", str(root / "occurrences"))
    return execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        **kwargs,
    )


def _with_sources(root, **kwargs) -> execute.Corpus:
    """Opens the corpus under ``root`` with the sources it was derived from."""
    return _open(
        root,
        sources_pattern=str(root / "data" / "*.parquet"),
        resolver={"pyridine": "c1ccncc1"}.__getitem__,
        **kwargs,
    )


@pytest.fixture(scope="module")
def sourced_corpus(corpus_dir) -> Iterator[execute.Corpus]:
    """The corpus over ``corpus_dir``, opened with its sources."""
    with _with_sources(corpus_dir) as value:
        yield value


def _exists_at(path, where):
    return {"op": "exists", "path": path, "where": where}


def _exists(where):
    return _exists_at("inputs.components", where)


def _reactions(table) -> set[str]:
    return set(table.column("reaction_id").to_pylist())


def _search(corpus, where) -> set[str]:
    return _reactions(corpus.search(query.Query.model_validate({"where": where})))


def test_a_search_narrows_to_the_dataset_a_reaction_came_from(corpus):
    where = {"op": "eq", "path": "dataset_id", "value": {"literal": "ord_dataset-bb"}}
    assert _search(corpus, where) == {"ord-bb01"}


def test_reactions_group_by_the_dataset_they_came_from(corpus):
    table = corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": ["dataset_id"],
                    "measures": [{"fn": "count", "name": "n"}],
                }
            }
        )
    )
    # Named by the path, which for a top-level column is the column's own name.
    assert table.column_names == ["dataset_id", "n"]
    counts = dict(
        zip(
            table.column("dataset_id").to_pylist(),
            table.column("n").to_pylist(),
            strict=True,
        )
    )
    assert counts == {"ord_dataset-aa": 2, "ord_dataset-bb": 1}


@pytest.mark.parametrize("descending", [False, True])
def test_an_aggregate_over_a_level_orders_by_its_group_key(corpus, descending):
    # The group key is read from the level's elements, so the ordering has to name the
    # grouped expression rather than a reaction column of the same path.
    table = corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "over": "inputs.components",
                    "group_by": ["smiles"],
                    "measures": [{"fn": "count", "name": "n"}],
                },
                "order_by": [{"key": "smiles", "descending": descending}],
            }
        )
    )
    # Named by the path as written, relative to the element, not by the expression.
    assert table.column_names == ["smiles", "n"]
    smiles = table.column("smiles").to_pylist()
    assert len(smiles) == 3
    assert smiles == sorted(smiles, reverse=descending)


def test_an_aggregate_over_reactions_orders_by_its_group_key(corpus):
    table = corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": ["dataset_id"],
                    "measures": [{"fn": "count", "name": "n"}],
                },
                "order_by": [{"key": "dataset_id", "descending": True}],
            }
        )
    )
    assert table.column("dataset_id").to_pylist() == [
        "ord_dataset-bb",
        "ord_dataset-aa",
    ]


@pytest.mark.parametrize("name", ["true", "order", "at"])
@pytest.mark.parametrize("descending", [False, True])
def test_a_measure_named_by_a_keyword_orders_by_the_measure(corpus, name, descending):
    # An ORDER BY naming order or at bare fails to parse. aa holds two reactions and bb
    # one, so each direction puts a different dataset first.
    table = corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": ["dataset_id"],
                    "measures": [{"fn": "count", "name": name}],
                },
                "order_by": [{"key": name, "descending": descending}],
            }
        )
    )
    assert table.column_names == ["dataset_id", name]
    expected = [("ord_dataset-bb", 1), ("ord_dataset-aa", 2)]
    if descending:
        expected.reverse()
    assert [(row["dataset_id"], row[name]) for row in table.to_pylist()] == expected


def test_an_aggregate_over_a_level_breaks_count_ties_by_its_group_key(corpus):
    table = corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "over": "inputs.components",
                    "group_by": ["smiles"],
                    "measures": [{"fn": "count", "name": "n"}],
                },
                "order_by": [{"key": "n", "descending": True}, {"key": "smiles"}],
            }
        )
    )
    rows = [(row["n"], row["smiles"]) for row in table.to_pylist()]
    assert rows == sorted(rows, key=lambda row: (-row[0], row[1]))
    assert rows[0][0] == 2


_TEMPERATURE_CONTROL = "conditions.temperature.control.type"
_STIRRING = "conditions.stirring.type"


@pytest.fixture(scope="module")
def conditions_corpus(tmp_path_factory) -> Iterator[execute.Corpus]:
    """A corpus whose temperature control and stirring, both ``type``, differ.

    The two paths share a last segment and hold values from different enums, so a
    column holding the other path's values shows as a name the enum does not have.
    """
    controls = reaction_pb2.TemperatureConditions.TemperatureControl
    stirring = reaction_pb2.StirringConditions
    reactions = []
    for reaction_id, control, method in (
        ("ord-dd01", controls.AMBIENT, stirring.STIR_BAR),
        ("ord-dd02", controls.OIL_BATH, stirring.NONE),
        ("ord-dd03", controls.OIL_BATH, stirring.NONE),
    ):
        reaction = _reaction(reaction_id, components=[("CCO", _ROLE.REACTANT)])
        reaction.conditions.temperature.control.type = control
        reaction.conditions.stirring.type = method
        reactions.append(reaction)
    root = tmp_path_factory.mktemp("conditions")
    source = root / "data" / "ord_dataset-dd.parquet"
    source.parent.mkdir(parents=True, exist_ok=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-dd",
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
    with execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        yield value


def test_group_by_paths_sharing_a_last_segment_name_their_own_columns(
    conditions_corpus,
):
    # DuckDB names an unaliased column after its last segment, which would return two
    # columns both called type and leave a reader by name one of them.
    table = conditions_corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": [_TEMPERATURE_CONTROL, _STIRRING],
                    "measures": [{"fn": "count", "name": "n"}],
                }
            }
        )
    )
    assert table.column_names == [_TEMPERATURE_CONTROL, _STIRRING, "n"]
    rows = {
        (row[_TEMPERATURE_CONTROL], row[_STIRRING], row["n"])
        for row in table.to_pylist()
    }
    assert rows == {("AMBIENT", "STIR_BAR", 1), ("OIL_BATH", "NONE", 2)}


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        # Each ordering puts the other path's values in the opposite order, so sorting
        # by the wrong one shows.
        (_STIRRING, ["STIR_BAR", "NONE"]),
        (_TEMPERATURE_CONTROL, ["NONE", "STIR_BAR"]),
    ],
)
def test_an_aggregate_orders_by_group_by_paths_sharing_a_last_segment(
    conditions_corpus, key, expected
):
    table = conditions_corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": [_TEMPERATURE_CONTROL, _STIRRING],
                    "measures": [{"fn": "count", "name": "n"}],
                },
                "order_by": [{"key": key, "descending": True}],
            }
        )
    )
    assert table.column(_STIRRING).to_pylist() == expected


def test_a_source_without_a_dataset_id_reads_as_null(tmp_path):
    # A source recording no dataset ID gets no stamp; its reactions answer is_null and
    # group under NULL rather than vanishing or borrowing another file's ID.
    for name, dataset_id, reaction_ids in (
        ("named", "ord_dataset-named", ["ord-nm01"]),
        ("unnamed", "", ["ord-un01", "ord-un02"]),
    ):
        source = tmp_path / "data" / f"{name}.parquet"
        source.parent.mkdir(parents=True, exist_ok=True)
        parquet.save_dataset(
            dataset_pb2.Dataset(
                dataset_id=dataset_id,
                name="test",
                description="test",
                reactions=[
                    _reaction(reaction_id, components=[("CCO", _ROLE.REACTANT)])
                    for reaction_id in reaction_ids
                ],
            ),
            str(source),
        )
        projected = tmp_path / "projections" / source.name
        projected.parent.mkdir(parents=True, exist_ok=True)
        projection.write_projection(source, projected)
        structured = tmp_path / "structures" / source.name
        structured.parent.mkdir(parents=True, exist_ok=True)
        structures.write_structures(projected, structured)
    unnamed = tmp_path / "projections" / "unnamed.parquet"
    assert base.load_stamps(unnamed).source_dataset_id is None
    with execute.Corpus(
        str(tmp_path / "projections" / "*.parquet"),
        str(tmp_path / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        where = {"op": "is_null", "path": "dataset_id"}
        assert _search(value, where) == {"ord-un01", "ord-un02"}
        table = value.search(
            query.Query.model_validate(
                {
                    "aggregate": {
                        "group_by": ["dataset_id"],
                        "measures": [{"fn": "count", "name": "n"}],
                    }
                }
            )
        )
    counts = dict(
        zip(
            table.column("dataset_id").to_pylist(),
            table.column("n").to_pylist(),
            strict=True,
        )
    )
    assert counts == {"ord_dataset-named": 1, None: 2}


_ALANINE = "C[C@@H](N)C(=O)O"
_ALANINE_MIRROR = "C[C@H](N)C(=O)O"
_RACEMIC_ALANINE = f"{_ALANINE_MIRROR} |&1:1|"
_UNKNOWN_ALANINE = f"{_ALANINE_MIRROR} |o1:1|"


@pytest.fixture(scope="module")
def chiral_root(tmp_path_factory) -> pathlib.Path:
    """A corpus holding alanine as each enantiomer, without stereo, and in two groups.

    The racemate (``&1``) and the unknown single enantiomer (``o1``) are both drawn as
    the mirror image of ``_ALANINE``, so a pattern drawn as ``_ALANINE`` reaches them
    only through their stereo groups.
    """
    root = tmp_path_factory.mktemp("chiral")
    reactions = [
        _reaction("ord-ch01", components=[(_ALANINE, _ROLE.REACTANT)]),
        _reaction("ord-ch02", components=[(_ALANINE_MIRROR, _ROLE.REACTANT)]),
        _reaction("ord-ch03", components=[("CC(N)C(=O)O", _ROLE.REACTANT)]),
        _reaction("ord-ch04", components=[(_RACEMIC_ALANINE, _ROLE.REACTANT)]),
        _reaction("ord-ch05", components=[(_UNKNOWN_ALANINE, _ROLE.REACTANT)]),
    ]
    source = root / "data" / "ord_dataset-ch.parquet"
    source.parent.mkdir(parents=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-ch",
            name="test",
            description="test",
            reactions=reactions,
        ),
        str(source),
    )
    projected = root / "projections" / source.name
    projected.parent.mkdir(parents=True)
    projection.write_projection(source, projected)
    structured = root / "structures" / source.name
    structured.parent.mkdir(parents=True)
    structures.write_structures(projected, structured)
    return root


@pytest.fixture(scope="module")
def chiral_corpus(chiral_root) -> Iterator[execute.Corpus]:
    """The corpus over ``chiral_root``."""
    with _open(chiral_root, resolver={}.__getitem__) as value:
        yield value


@pytest.fixture(scope="module")
def reaction_root(tmp_path_factory) -> pathlib.Path:
    """Reactions that a reaction SMARTS reading every template keeps apart."""
    root = tmp_path_factory.mktemp("reactions")
    acid, amine, amide = "CC(=O)O", "NCc1ccccc1", "CC(=O)NCc1ccccc1"
    reactions = [
        # An amide coupling: both templates on REACTANT components, the amide a product.
        _reaction(
            "ord-rs01",
            components=[(acid, _ROLE.REACTANT), (amine, _ROLE.REACTANT)],
            product=amide,
        ),
        # The amine recorded as a reagent, so no REACTANT holds the second template.
        _reaction(
            "ord-rs02",
            components=[(acid, _ROLE.REACTANT), (amine, _ROLE.REAGENT)],
            product=amide,
        ),
        # An esterification: the acid template matches, the product template does not.
        _reaction(
            "ord-rs03",
            components=[(acid, _ROLE.REACTANT), ("CO", _ROLE.REACTANT)],
            product="CC(=O)OC",
        ),
        # The acid and an amide product with no amine anywhere: one reactant template of
        # two, which is enough for the cartridge's @>.
        _reaction(
            "ord-rs04",
            components=[(acid, _ROLE.REACTANT), ("Cc1ccccc1", _ROLE.REACTANT)],
            product=amide,
        ),
        # A Suzuki coupling with its palladium recorded as a catalyst.
        _reaction(
            "ord-rs05",
            components=[
                ("Brc1ccccc1", _ROLE.REACTANT),
                ("OB(O)c1ccccc1", _ROLE.REACTANT),
                ("[Pd]", _ROLE.CATALYST),
            ],
            product="c1ccc(-c2ccccc2)cc1",
        ),
        # An esterification of N-acetylglycine, whose one reactant holds both the acid
        # and the nitrogen, with an amine base as a reagent: every template matches,
        # but not on two different REACTANT components.
        _reaction(
            "ord-rs06",
            components=[
                ("CC(=O)NCC(=O)O", _ROLE.REACTANT),
                ("CO", _ROLE.REACTANT),
                ("CCN(CC)CC", _ROLE.REAGENT),
            ],
            product="CC(=O)NCC(=O)OC",
        ),
        # The same esterification with N-acetylglycine recorded twice, as two portions
        # would be: two REACTANT components, one molecule.
        _reaction(
            "ord-rs07",
            components=[
                ("CC(=O)NCC(=O)O", _ROLE.REACTANT),
                ("CC(=O)NCC(=O)O", _ROLE.REACTANT),
                ("CO", _ROLE.REACTANT),
            ],
            product="CC(=O)NCC(=O)OC",
        ),
    ]
    source = root / "data" / "ord_dataset-rs.parquet"
    source.parent.mkdir(parents=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-rs",
            name="test",
            description="test",
            reactions=reactions,
        ),
        str(source),
    )
    projected = root / "projections" / source.name
    projected.parent.mkdir(parents=True)
    projection.write_projection(source, projected)
    structured = root / "structures" / source.name
    structured.parent.mkdir(parents=True)
    structures.write_structures(projected, structured)
    _write_pivots(root, sorted(pivot.LEVELS))
    return root


@pytest.fixture(scope="module", params=["pivoted", "unpivoted"])
def reaction_corpus(request, reaction_root) -> Iterator[execute.Corpus]:
    """The reactions, counted through pivot artifacts and through the lists."""
    with execute.Corpus(
        str(reaction_root / "projections" / "*.parquet"),
        str(reaction_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=_pivots_for(reaction_root, request.param),
    ) as value:
        yield value


def _pivots_for(root: pathlib.Path, route: str) -> str | None:
    """Returns the tree's pivots directory for ``"pivoted"``, and None for the lists."""
    return str(root / "pivots") if route == "pivoted" else None


def test_a_reaction_smarts_needs_every_template_on_its_side(reaction_corpus):
    where = {"op": "reaction_smarts", "smarts": "C(=O)O.N>>C(=O)N"}
    assert _search(reaction_corpus, where) == {"ord-rs01"}


def test_a_grouped_reaction_smarts_matches_both_pieces_in_one_molecule(
    reaction_corpus,
):
    where = {"op": "reaction_smarts", "smarts": "(C(=O)O.N)>>C(=O)N"}
    assert _search(reaction_corpus, where) == {"ord-rs06", "ord-rs07"}


def test_an_agent_template_finds_a_catalyst(reaction_corpus):
    where = {"op": "reaction_smarts", "smarts": "cBr.cB(O)O>[Pd]>c-c"}
    assert _search(reaction_corpus, where) == {"ord-rs05"}


def test_a_substructure_matches_the_drawn_stereoisomer_unless_chirality_is_off(
    chiral_corpus,
):
    # Asked in both orders on one corpus, because the match set is cached: a cache key
    # that left chirality out would answer the second question with the first one's set.
    pattern = {"op": "substructure", "path": "smiles", "smarts": _ALANINE}
    drawn = {"ord-ch01", "ord-ch04", "ord-ch05"}
    assert _search(chiral_corpus, _exists(pattern)) == drawn
    assert _search(chiral_corpus, _exists(pattern | {"chirality": False})) == {
        "ord-ch01",
        "ord-ch02",
        "ord-ch03",
        "ord-ch04",
        "ord-ch05",
    }
    assert _search(chiral_corpus, _exists(pattern)) == drawn


@pytest.mark.parametrize(
    ("smarts", "expected"),
    [
        (_ALANINE, {"ord-ch01", "ord-ch04", "ord-ch05"}),
        (_ALANINE_MIRROR, {"ord-ch02", "ord-ch04", "ord-ch05"}),
    ],
)
def test_a_stereo_group_matches_either_enantiomer_it_allows(
    chiral_corpus, smarts, expected
):
    # The racemate holds both enantiomers and the unknown one may be either, so each is
    # found by a pattern drawn as either, whichever way round its center was written.
    pattern = {"op": "substructure", "path": "smiles", "smarts": smarts}
    assert _search(chiral_corpus, _exists(pattern)) == expected


def test_a_pattern_drawn_as_a_racemate_matches_only_the_racemate(chiral_corpus):
    # A single enantiomer, known or not, is not a racemic mixture.
    pattern = {
        "op": "substructure",
        "path": "smiles",
        "smarts": f"{_ALANINE} |&1:1|",
    }
    assert _search(chiral_corpus, _exists(pattern)) == {"ord-ch04"}


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


def test_a_memory_limit_is_given_to_duckdb(corpus_dir):
    with _open(corpus_dir, memory_limit="512MiB") as corpus:
        setting = corpus._connection.execute(
            "SELECT current_setting('memory_limit')"
        ).fetchone()
    assert setting == ("512.0 MiB",)


@pytest.mark.parametrize(
    "statement",
    [
        "SELECT * FROM read_text('{outside}')",
        "COPY (SELECT 1) TO '{outside}.csv'",
        "ATTACH '{outside}.duckdb' AS elsewhere",
        "SELECT * FROM read_text('http://169.254.170.2/v2/credentials')",
        "INSTALL httpfs",
    ],
)
@pytest.mark.parametrize("opened", ["corpus", "sourced_corpus"])
def test_a_search_cursor_reaches_nothing_outside_the_corpus(
    request, opened, tmp_path, statement
):
    # Opened with sources, a corpus allows their directories as well, and no others.
    corpus = request.getfixturevalue(opened)
    outside = tmp_path / "secret.txt"
    outside.write_text("not the corpus's to read\n")
    cursor = corpus._connection.cursor()
    try:
        with pytest.raises(duckdb.PermissionException):
            cursor.execute(statement.format(outside=outside))
    finally:
        cursor.close()


@pytest.mark.parametrize(
    "statement",
    [
        "SET enable_external_access=true",
        "SET allowed_directories=['/']",
        "SET autoload_known_extensions=true",
    ],
)
def test_a_search_cursor_cannot_lift_the_confinement(corpus, statement):
    cursor = corpus._connection.cursor()
    try:
        with pytest.raises(duckdb.InvalidInputException, match="lock"):
            cursor.execute(statement)
    finally:
        cursor.close()


# A Pd G3 precatalyst from the corpus, whose two dative bonds the projection writes in a
# |C:...| block.
_PRECATALYST = (
    "CC(C)c1cc(C(C)C)c(-c2ccccc2[P](C2CCCCC2)(C2CCCCC2)[Pd]2([NH2]c3ccccc3-c3cccc[c]32)"
    "[O]S(C)(=O)=O)c(C(C)C)c1 |C:16.31,30.32|"
)


def test_a_resolved_name_is_spelled_as_the_projection_spells_it():
    # Plain SMILES writes a dative bond inline, as ->, so a resolver built on RDKit
    # answers in a spelling the projection never writes, and two names for one molecule
    # would key two match sets.
    projected = message_helpers.canonical_smiles(Chem.MolFromSmiles(_PRECATALYST))
    inline = Chem.MolToSmiles(Chem.MolFromSmiles(_PRECATALYST))
    assert "->" in inline
    assert "|C:" in projected
    assert execute._canonical(inline) == projected
    assert execute._canonical(projected) == projected


def test_a_resolved_racemate_keeps_its_stereo_group(chiral_root):
    # Plain SMILES drops the |&1:1| block, which would leave the racemate's query
    # molecule as the enantiomer it was drawn as -- and keyed alongside that
    # enantiomer, so whichever was asked first would answer both.
    resolver = {"racemate": _RACEMIC_ALANINE, "enantiomer": _ALANINE_MIRROR}
    asked = [("enantiomer", {"ord-ch02"}), ("racemate", {"ord-ch04"})]
    with _open(chiral_root, resolver=resolver.__getitem__) as corpus:
        for name, expected in asked:
            where = {"op": "same_compound", "path": "smiles", "compound": name}
            assert _search(corpus, _exists(where)) == expected


def test_a_name_resolving_to_nonsense_is_an_error(corpus_dir):
    # Canonicalizing cannot fix an unparseable answer, and swallowing it here would
    # hide a resolver fault behind an empty result.
    assert execute._canonical("not a molecule") is None
    where = {"op": "same_compound", "path": "smiles", "compound": "nonsense"}
    resolver = {"nonsense": "not a molecule"}
    with (
        _open(corpus_dir, resolver=resolver.__getitem__) as corpus,
        pytest.raises(ValueError, match="unreadable SMILES 'not a molecule'"),
    ):
        _search(corpus, _exists(where))


def _same_compound(name):
    return _exists({"op": "same_compound", "path": "smiles", "compound": name})


def test_a_name_is_resolved_once_per_corpus(corpus_dir):
    # A resolver is an external service, and what a name means does not change between
    # one question and the next.
    asked = []

    def resolver(name):
        asked.append(name)
        return {"pyridine": "c1ccncc1"}[name]

    with _open(corpus_dir, resolver=resolver) as corpus:
        first = _search(corpus, _same_compound("pyridine"))
        second = _search(corpus, _same_compound("pyridine"))
    assert first == second == {"ord-aa01", "ord-aa02"}
    assert asked == ["pyridine"]


def test_the_resolution_cache_drops_the_name_unused_longest(corpus_dir, monkeypatch):
    monkeypatch.setattr(execute, "_CACHED_RESOLUTIONS", 2)
    asked = []
    answers = {"pyridine": "c1ccncc1", "ethanol": "CCO", "benzene": "c1ccccc1"}

    def resolver(name):
        asked.append(name)
        return answers[name]

    # Pyridine is asked again before benzene arrives, so ethanol is what goes, and the
    # resolver is called for it again.
    order = ("pyridine", "ethanol", "pyridine", "benzene", "pyridine", "ethanol")
    with _open(corpus_dir, resolver=resolver) as corpus:
        for name in order:
            _search(corpus, _same_compound(name))
    assert asked == ["pyridine", "ethanol", "benzene", "ethanol"]


def test_a_failed_resolution_is_asked_again(corpus_dir):
    # A service that is down for one question must not answer every later one.
    answers = iter([ValueError("the service is down"), "c1ccncc1"])

    def resolver(name):
        answer = next(answers)
        if isinstance(answer, Exception):
            raise answer
        return answer

    with _open(corpus_dir, resolver=resolver) as corpus:
        with pytest.raises(ValueError, match="the service is down"):
            _search(corpus, _same_compound("pyridine"))
        assert _search(corpus, _same_compound("pyridine")) == {"ord-aa01", "ord-aa02"}


def test_a_name_resolving_to_nonsense_is_asked_again(corpus_dir):
    asked = []

    def resolver(name):
        asked.append(name)
        return "not a molecule"

    with _open(corpus_dir, resolver=resolver) as corpus:
        for _ in range(2):
            with pytest.raises(ValueError, match="unreadable SMILES 'not a molecule'"):
                _search(corpus, _same_compound("nonsense"))
    assert asked == ["nonsense", "nonsense"]


def test_a_search_says_what_each_name_resolved_to(corpus_dir):
    # The query names a compound, and the structure it searched for is what a reader
    # needs to judge the answer. It is stated as the corpus spells it, not as the
    # service answered.
    resolver = {"pyridine": "C1=CC=NC=C1", "ethanol": "OCC"}.__getitem__
    where = {
        "op": "or",
        "clauses": [_same_compound(n) for n in ("pyridine", "ethanol")],
    }
    with _open(corpus_dir, resolver=resolver) as corpus:
        answer = corpus.search(query.Query.model_validate({"where": where}))
    assert execute.resolutions(answer) == {"pyridine": "c1ccncc1", "ethanol": "CCO"}


def test_a_search_naming_no_compound_carries_no_resolutions(corpus_dir):
    where = _exists({"op": "same_compound", "path": "smiles", "smiles": "c1ccncc1"})
    with _open(corpus_dir, resolver={}.__getitem__) as corpus:
        answer = corpus.search(query.Query.model_validate({"where": where}))
    assert answer.num_rows == 2
    assert execute.RESOLVED.encode() not in (answer.schema.metadata or {})
    assert execute.resolutions(answer) == {}


def _copy_corpus(corpus_dir, tmp_path) -> pathlib.Path:
    """Copies the whole tree and its sources so a test can damage one in isolation."""
    for name in ("data", "projections", "structures"):
        (tmp_path / name).mkdir()
        for source in (corpus_dir / name).glob("*.parquet"):
            (tmp_path / name / source.name).write_bytes(source.read_bytes())
    for name in ("pivots", "occurrences"):
        shutil.copytree(corpus_dir / name, tmp_path / name)
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


def test_a_projection_without_statistics_is_refused(corpus_dir, tmp_path):
    # The check above reads its bound from the footer, so a projection recording none
    # would pass it whatever IDs it carries.
    root = _copy_corpus(corpus_dir, tmp_path)
    target = root / "projections" / "ord_dataset-aa.parquet"
    pq.write_table(pq.read_table(target), target, write_statistics=False)
    with pytest.raises(
        execute.PairingError,
        match=re.escape(f"{target} records no statistics for its structure IDs")
        + r".*ord_schema\.artifacts\.scripts\.derive_projection --force",
    ):
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
        )


def _footers_alive() -> int:
    """Returns how many parsed Parquet footers the collector can still see."""
    return sum(isinstance(found, pq.FileMetaData) for found in gc.get_objects())


def test_reading_the_structure_id_bound_frees_the_footer_on_return(corpus_dir):
    # Opening a corpus reads the bound from every projection, so a footer that waits for
    # the cyclic collector rather than going with its file is one more parsed footer per
    # dataset held at once. The collector is held off so that it cannot free one first.
    path = str(corpus_dir / "projections" / "ord_dataset-aa.parquet")
    gc.collect()
    gc.disable()
    try:
        before = _footers_alive()
        assert execute._max_structure_id(path, "") is not None
        after = _footers_alive()
    finally:
        gc.enable()
    assert after == before


def test_the_structure_id_bound_reads_every_row_group_and_nested_column(tmp_path):
    # The columns are found in the first row group and read in every one, so the
    # largest ID sits in the last group, and in a nested column rather than the
    # top-level one.
    path = tmp_path / "split.parquet"
    pq.write_table(
        pa.table(
            {
                "structure_id": [0, 1, 2],
                "products": [[{"structure_id": 1}], [], [{"structure_id": 7}]],
            }
        ),
        path,
        row_group_size=1,
    )
    assert execute._max_structure_id(str(path), "") == 7


def _field_names(fields: Iterable[pa.Field]) -> Iterator[str]:
    """Yields the name of every field, however deeply nested, in ``fields``."""
    for field in fields:
        yield field.name
        kind = field.type
        if pa.types.is_struct(kind):
            yield from _field_names(kind)
        elif pa.types.is_map(kind):
            yield from _field_names([kind.key_field, kind.item_field])
        elif pa.types.is_list(kind) or pa.types.is_large_list(kind):
            yield from _field_names([kind.value_field])


@pytest.mark.parametrize(
    "schema",
    [
        projection.SCHEMA,
        occurrences.SCHEMA,
        *(pivot.schema(level) for level in pivot.LEVELS.values()),
    ],
)
def test_no_artifact_field_name_holds_a_dot(schema):
    # The bound finds its columns by the last component of a dotted chunk path, which
    # names the leaf only while no field name holds a dot itself.
    assert [name for name in _field_names(schema) if "." in name] == []


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


def _derive(root, shard, reactions) -> None:
    """Writes a source under ``root`` and derives its projection and structures."""
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


def _stored(corpus_dir) -> dict[str, bytes]:
    """Returns every record the sources store, keyed by reaction ID."""
    stored = {}
    for source in sorted((corpus_dir / "data").glob("*.parquet")):
        table = pq.read_table(source)
        stored.update(
            zip(
                table.column("reaction_id").to_pylist(),
                table.column("reaction").to_pylist(),
                strict=True,
            )
        )
    return stored


def test_reactions_come_back_in_the_order_asked(sourced_corpus, corpus_dir):
    asked = ["ord-bb01", "ord-aa02", "ord-aa01"]
    table = sourced_corpus.reactions(asked)
    assert table.schema.names == ["reaction_id", "dataset_id", "reaction"]
    assert table.schema.types == [pa.string(), pa.string(), pa.binary()]
    assert table.column("reaction_id").to_pylist() == asked
    assert table.column("dataset_id").to_pylist() == [
        "ord_dataset-bb",
        "ord_dataset-aa",
        "ord_dataset-aa",
    ]
    blobs = table.column("reaction").to_pylist()
    # The record is the source's own bytes rather than a message parsed and written
    # again, and it parses to the reaction it was asked for.
    stored = _stored(corpus_dir)
    assert blobs == [stored[reaction_id] for reaction_id in asked]
    assert [reaction_pb2.Reaction.FromString(blob).reaction_id for blob in blobs] == (
        asked
    )


def test_a_reaction_the_corpus_does_not_hold_has_no_row(sourced_corpus):
    table = sourced_corpus.reactions(["ord-zz99", "ord-aa02"])
    assert table.column("reaction_id").to_pylist() == ["ord-aa02"]


def test_a_request_the_corpus_holds_none_of_returns_no_rows(sourced_corpus):
    table = sourced_corpus.reactions(["ord-zz98", "ord-zz99"])
    assert table.num_rows == 0
    assert table.schema == sourced_corpus.reactions(["ord-aa01"]).schema


def test_a_reaction_asked_twice_comes_back_once(sourced_corpus):
    table = sourced_corpus.reactions(["ord-aa01", "ord-bb01", "ord-aa01"])
    assert table.column("reaction_id").to_pylist() == ["ord-aa01", "ord-bb01"]


def test_an_empty_request_returns_the_columns_and_no_rows(sourced_corpus):
    table = sourced_corpus.reactions([])
    assert table.num_rows == 0
    assert table.schema == sourced_corpus.reactions(["ord-aa01"]).schema


def test_a_request_for_more_reactions_than_max_rows_is_refused(corpus_dir):
    # Refused rather than cut: an ID the corpus does not hold has no row either, so a
    # fetch cut short would read as reactions that are not there.
    with _with_sources(corpus_dir, max_rows=2, warm=False) as value:
        with pytest.raises(ValueError, match="3 reactions, more than the 2"):
            value.reactions(["ord-aa01", "ord-aa02", "ord-bb01"])
        # The bound counts distinct IDs, so a repeat is not charged twice.
        table = value.reactions(["ord-aa01", "ord-aa01", "ord-bb01"])
        assert table.column("reaction_id").to_pylist() == ["ord-aa01", "ord-bb01"]
    with _with_sources(corpus_dir, max_rows=None, warm=False) as value:
        assert value.reactions(["ord-aa01", "ord-aa02", "ord-bb01"]).num_rows == 3


@pytest.mark.parametrize("method", ["reactions", "datasets"])
def test_a_corpus_opened_without_sources_has_none_to_serve(corpus, method):
    arguments = (["ord-aa01"],) if method == "reactions" else ()
    with pytest.raises(ValueError, match="without sources_pattern"):
        getattr(corpus, method)(*arguments)


def test_datasets_describe_each_source(sourced_corpus):
    assert sourced_corpus.datasets().to_pylist() == [
        {
            "dataset_id": "ord_dataset-aa",
            "name": "dataset aa",
            "description": "the aa reactions",
            "reactions": 2,
        },
        {
            "dataset_id": "ord_dataset-bb",
            "name": "dataset bb",
            "description": "the bb reactions",
            "reactions": 1,
        },
    ]


def test_a_search_s_hits_fetch_their_reactions(sourced_corpus):
    found = sourced_corpus.search(
        query.Query.model_validate(
            {
                "where": _exists(
                    {"op": "substructure", "path": "smiles", "compound": "pyridine"}
                ),
                "order_by": [{"key": "reaction_id", "descending": True}],
            }
        )
    )
    hits = found.column("reaction_id").to_pylist()
    assert hits == ["ord-aa02", "ord-aa01"]
    table = sourced_corpus.reactions(hits)
    assert table.column("reaction_id").to_pylist() == hits
    blobs = table.column("reaction").to_pylist()
    for blob, reaction_id in zip(blobs, hits, strict=True):
        reaction = reaction_pb2.Reaction.FromString(blob)
        assert reaction.reaction_id == reaction_id
        assert "c1ccncc1" in {
            identifier.value
            for component in reaction.inputs["in"].components
            for identifier in component.identifiers
        }


def test_a_search_s_column_is_fetched_without_converting_it(sourced_corpus):
    # What a search returns is a pyarrow table, so its column is what a caller holds.
    column = pa.chunked_array([["ord-aa02", "ord-aa01"]])
    for asked in (column, column.combine_chunks()):
        table = sourced_corpus.reactions(asked)
        assert table.column("reaction_id").to_pylist() == ["ord-aa02", "ord-aa01"]


def test_opening_with_sources_logs_what_hashing_them_cost(corpus_dir, caplog):
    # Hashing reads every reaction, which over ORD is most of what opening costs, so a
    # slow start is explained in the log rather than left to look like a hang.
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        _with_sources(corpus_dir, warm=False),
    ):
        pass
    assert any(re.match(r"hashed 2 sources in \d", m) for m in _messages(caplog))


def test_a_source_edited_after_derivation_is_refused(corpus_dir, tmp_path):
    # The projections restate the source as it was hashed when they were derived. An
    # edited source would serve records the search never saw.
    root = _copy_corpus(corpus_dir, tmp_path)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-bb",
            name="dataset bb",
            description="the bb reactions, edited",
            reactions=_REACTIONS["bb"],
        ),
        str(root / "data" / "ord_dataset-bb.parquet"),
    )
    with pytest.raises(
        execute.PairingError, match=r"holds the source these projections restate"
    ) as raised:
        _with_sources(root)
    assert str(root / "projections" / "ord_dataset-bb.parquet") in str(raised.value)


def test_a_missing_source_is_refused(corpus_dir, tmp_path):
    root = _copy_corpus(corpus_dir, tmp_path)
    (root / "data" / "ord_dataset-aa.parquet").unlink()
    with pytest.raises(
        execute.PairingError, match=r"holds the source these projections restate"
    ) as raised:
        _with_sources(root)
    assert str(root / "projections" / "ord_dataset-aa.parquet") in str(raised.value)


def test_a_source_no_projection_restates_is_refused(corpus_dir, tmp_path):
    # A fetch would return its reactions though no search can find them.
    root = _copy_corpus(corpus_dir, tmp_path)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-cc",
            name="test",
            description="test",
            reactions=[_reaction("ord-cc01", components=[("CCN", _ROLE.REACTANT)])],
        ),
        str(root / "data" / "ord_dataset-cc.parquet"),
    )
    with pytest.raises(execute.PairingError, match="no projection restates") as raised:
        _with_sources(root)
    assert str(root / "data" / "ord_dataset-cc.parquet") in str(raised.value)


def test_two_identical_sources_are_refused(corpus_dir, tmp_path):
    # Which of them answers a fetch would be arbitrary, so neither may.
    root = _copy_corpus(corpus_dir, tmp_path)
    source = root / "data" / "ord_dataset-aa.parquet"
    shutil.copy(source, root / "data" / "copy.parquet")
    with pytest.raises(execute.PairingError, match="same source dataset"):
        _with_sources(root)


def test_a_sources_pattern_matching_nothing_is_refused(corpus_dir, tmp_path):
    with pytest.raises(execute.PairingError, match="no sources matched"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            sources_pattern=str(tmp_path / "absent" / "*.parquet"),
        )


def test_a_source_whose_reactions_cannot_be_read_is_refused(corpus_dir, tmp_path):
    # The footer reads, so the file opens; hashing its damaged reactions is what fails.
    root = _copy_corpus(corpus_dir, tmp_path)
    source = root / "data" / "ord_dataset-aa.parquet"
    chunk = pq.ParquetFile(source).metadata.row_group(0).column(1)
    damaged = bytearray(source.read_bytes())
    start = chunk.data_page_offset + 8
    damaged[start : start + 32] = b"\x00" * 32
    source.write_bytes(bytes(damaged))
    with pytest.raises(execute.PairingError, match="not a readable source dataset"):
        _with_sources(root)


def test_a_sources_pattern_matching_an_artifact_is_refused(corpus_dir):
    with pytest.raises(execute.PairingError, match="not a readable source dataset"):
        execute.Corpus(
            str(corpus_dir / "projections" / "*.parquet"),
            str(corpus_dir / "structures" / "*.parquet"),
            sources_pattern=str(corpus_dir / "projections" / "*.parquet"),
        )


def test_sources_that_do_not_join_to_their_dataset_ids_are_refused(
    corpus_dir, tmp_path, monkeypatch
):
    # The same loss _prepare counts its views against: a source DuckDB resolves
    # differently than Python did drops out of the join, and a reaction a search finds
    # is then one a fetch does not return. Staged as that test stages it, by pointing
    # the view at copies filed elsewhere within a tree the corpus reads.
    root = _copy_corpus(corpus_dir, tmp_path)
    elsewhere = root / "data" / "elsewhere"
    elsewhere.mkdir()
    for source in sorted((root / "data").glob("*.parquet")):
        shutil.copy(source, elsewhere / source.name)
    original = execute._sql_paths

    def redirected(paths):
        listed = [str(path) for path in paths]
        if all(pathlib.Path(path).parent.name == "data" for path in listed):
            return original(str(elsewhere / pathlib.Path(path).name) for path in listed)
        return original(listed)

    monkeypatch.setattr(execute, "_sql_paths", redirected)
    with pytest.raises(execute.PairingError, match="the sources hold 3 reactions"):
        _with_sources(root)


def test_a_reaction_two_sources_state_is_refused(tmp_path):
    # Each source pairs with its own projection, so the corpus opens; but a fetch of
    # the reaction they share has two records to choose between.
    for shard, smiles in (("aa", "c1ccncc1"), ("bb", "CCO")):
        _derive(
            tmp_path,
            shard,
            [_reaction("ord-x01", components=[(smiles, _ROLE.SOLVENT)])],
        )
    with (
        _with_sources(tmp_path, warm=False) as value,
        pytest.raises(execute.PairingError, match="ord-x01 is stated 2 times"),
    ):
        value.reactions(["ord-x01"])


def test_a_source_without_a_dataset_id_is_fetched_and_described(tmp_path):
    # The dataset ID comes from the projection's stamps, the value a search's
    # dataset_id column carries, so a source recording none reads as NULL in both.
    _derive(
        tmp_path,
        "named",
        [_reaction("ord-nm01", components=[("CCO", _ROLE.REACTANT)])],
    )
    source = tmp_path / "data" / "unnamed.parquet"
    parquet.save_dataset(
        dataset_pb2.Dataset(
            name="unnamed",
            description="test",
            reactions=[_reaction("ord-un01", components=[("CCO", _ROLE.REACTANT)])],
        ),
        str(source),
    )
    projection.write_projection(source, tmp_path / "projections" / source.name)
    structures.write_structures(
        tmp_path / "projections" / source.name,
        tmp_path / "structures" / source.name,
    )
    with _with_sources(tmp_path, warm=False) as value:
        table = value.reactions(["ord-un01", "ord-nm01"])
        datasets = value.datasets()
    assert table.column("dataset_id").to_pylist() == [None, "ord_dataset-named"]
    assert datasets.column("dataset_id").to_pylist() == ["ord_dataset-named", None]
    assert datasets.column("name").to_pylist() == ["test", "unnamed"]


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
    with _open(corpus_dir, resolver={}.__getitem__) as value:
        # Publishing a view is DDL on the corpus's own connection, so the pivots are
        # published first and what is counted is the search.
        value.check_pivots()
        value._connection = counter = _CursorCounter(value._connection)
        assert value.search(request).num_rows == 1
        # Steady state, with every view published: one search, one cursor.
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


def test_a_search_does_not_wait_on_the_scan_order(corpus):
    # Set globally for the reason the footer cache is: a search runs on its own cursor,
    # which a session-scoped setting would leave holding rows back until every earlier
    # row group had finished.
    cursor = corpus._connection.cursor()
    try:
        setting = cursor.execute(
            "SELECT current_setting('preserve_insertion_order')"
        ).fetchone()
    finally:
        cursor.close()
    assert setting == (False,)


def test_a_rewritten_artifact_is_not_answered_from_its_cached_footer(
    corpus_dir, tmp_path
):
    # Holding the footers means holding something the file can go on to contradict. The
    # rewrite lands in the same second as the read that cached it, which is the case a
    # cache keyed on a coarse timestamp would answer stale. Opened without pivots, so
    # both searches reach the projection rather than an artifact derived from it.
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
    # The library is large at corpus scale. First searches arriving together must not
    # each build their own copy, so the build is serialized and whoever waits takes the
    # finished one.
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
        # together is what makes the race reachable at all. The corpus keeps what a name
        # resolved to, but no thread can store an answer before all of them arrive, so
        # each misses and calls this. Were concurrent misses ever made to wait on one
        # call, the threads that skipped it would leave this barrier one short, so it
        # times out rather than parking the suite.
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

    def index(path, allocate):
        nonlocal spent
        condition = execute._index_condition(path, allocate)
        spent = spent or condition is not None
        return condition

    query.compile_query(query.Query.model_validate(body), index=index)
    return spent


def test_a_role_beside_the_structure_is_left_to_the_pivot():
    # The role is an element field, and an occurrence row carries none, so the
    # quantifier goes to the level's pivot, which holds the whole element, and the role
    # binds there as a parameter like any other literal.
    compiled = query.compile_query(
        query.Query.model_validate(
            {"where": _exists({"op": "and", "clauses": [_SUBSTRUCTURE, _SOLVENT]})}
        ),
        index=execute._index_condition,
        pivot=lambda path: "pivoted" if path == "inputs.components" else None,
    )
    assert "FROM occurrences" not in compiled.sql
    assert "FROM pivoted" in compiled.sql
    assert "SOLVENT" not in compiled.sql
    assert list(compiled.literals.values()) == ["SOLVENT"]


def test_a_reaction_smarts_spends_the_occurrence_index_on_its_products_alone():
    # A product template compiles to an exists over products with a structure alone,
    # which the index answers. A reactant template binds the role beside its structure,
    # which an occurrence row does not carry.
    assert _index_spent(
        {"where": {"op": "reaction_smarts", "smarts": "C(=O)O.N>>C(=O)N"}}
    )
    assert not _index_spent(
        {"where": {"op": "reaction_smarts", "smarts": "C(=O)O.N>>"}}
    )


def _no_index_condition(path, allocate):
    """Stands in for _index_condition so every quantifier compiles over the elements."""
    del path, allocate  # Unused.


# Every shape where the index takes a clause. The index answers one quantifier, not the
# query, so everything around that quantifier -- aggregates, orderings, limits, scalars,
# negation, disjunction, a second quantifier -- composes and belongs here.
_ROUTABLE = {
    "structure alone": {"where": _exists(_SUBSTRUCTURE)},
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
                _exists(_SUBSTRUCTURE),
                _exists({"op": "substructure", "path": "smiles", "smarts": "c1ccccc1"}),
            ],
        }
    },
    # The shape the index answers half of: over the corpus, what the other half costs
    # is set by whether its column was materialized, which is why the budget is what it
    # is.
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
    # The first quantifier binds a role and compiles over the elements; the second is
    # the index's. One statement holds both routes.
    "a role-bound quantifier beside an indexed one": {
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
    # The role is an element field like any other, and an occurrence row carries none.
    # Routed, the condition would be dropped and pyridine as a reactant would match.
    "structure and role": {
        "where": _exists({"op": "and", "clauses": [_SUBSTRUCTURE, _SOLVENT]})
    },
    "role and structure, reversed": {
        "where": _exists({"op": "and", "clauses": [_SOLVENT, _SUBSTRUCTURE]})
    },
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
    # Routed, the inequality would be dropped and the pyridine that is the solvent would
    # match too.
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
    # Two roles on one element is a contradiction the projection answers with no rows.
    # Routed, both would be dropped and every pyridine would match.
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
    # isolated_color is an element field no occurrence row carries, like every other.
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
        # Pyridine as the solvent, which is aa01's and not aa02's. Routed with the role
        # dropped, this would come back as both.
        ("structure and role", {"ord-aa01"}),
        # A pyridine that is not the solvent, which is aa02's reactant. Routed with the
        # role dropped, this would come back as aa01 too, the one reaction it excludes.
        ("a role inequality", {"ord-aa02"}),
        # One component required to be two roles at once, which nothing is. Routed with
        # both roles dropped, this would come back as every pyridine.
        ("two roles at once", set()),
        # An equality on smiles is a condition beside the structure like any other.
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
    # wrong. Asserted for the cases whose routed answer would not merely differ but
    # invert or over-return.
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
    # A role is text a query supplies. Spliced into the SQL unescaped, this literal
    # closes the string and leaves `OR '1'='1` as a condition, which matches every
    # element -- so the assertion that bites is that nothing comes back, not that
    # nothing raises.
    matched = _search(corpus, _role_and_structure("c1ccncc1", "SOLVENT' OR '1'='1"))
    assert matched == set()


def test_the_index_is_published_once_and_only_when_wanted(corpus_dir, caplog):
    with _open(corpus_dir, resolver={}.__getitem__, warm=False) as value:
        # An element predicate with no structure in it cannot use the index, so it must
        # not pay to check one.
        value.search(
            query.Query.model_validate(
                {
                    "where": _exists(
                        {"op": "eq", "path": "smiles", "value": {"literal": "CCO"}}
                    )
                }
            )
        )
        assert value._indexed is None
        caplog.set_level(logging.INFO, logger="ord_schema.search.execute")
        request = query.Query.model_validate(_ROUTABLE["structure alone"])
        assert value.search(request).num_rows == 2
        assert value._indexed is True
        # A second search reuses it rather than publishing again, which the answer
        # alone cannot show -- a second publish answers identically and costs a scan.
        assert value.search(request).num_rows == 2
    published = [
        record for record in caplog.records if record.message.startswith("indexed ")
    ]
    assert len(published) == 1


def test_the_indexed_paths_come_from_the_artifact():
    # One walk, in the artifact, and this pins that it stays one. A path routed here
    # that no artifact carries is a quantifier answered from nothing.
    assert set(execute.INDEXED_PATHS) == set(occurrences.PATHS)


def test_a_corpus_publishes_the_occurrence_index_at_open(corpus_dir):
    # The default, so the first structure query answers from an index that is already
    # checked rather than paying the scan that checks it.
    with _open(corpus_dir, resolver={}.__getitem__) as value:
        assert value._indexed is True


def test_a_cold_corpus_leaves_the_index_to_the_first_query(corpus_dir):
    with _open(corpus_dir, resolver={}.__getitem__, warm=False) as value:
        assert value._indexed is None


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
    root = tmp_path / "corpus"
    shutil.copytree(corpus_dir, root)
    # Filed under a tree the corpus reads, since a copy outside every one is refused by
    # the confinement before a join could drop its rows.
    elsewhere = root / "projections" / "elsewhere"
    elsewhere.mkdir()
    for projected in sorted((root / "projections").glob("*.parquet")):
        shutil.copy(projected, elsewhere / projected.name)
    original = execute._sql_paths

    def redirected(paths):
        listed = [str(path) for path in paths]
        if all(pathlib.Path(path).parent.name == "projections" for path in listed):
            return original(str(elsewhere / pathlib.Path(path).name) for path in listed)
        return original(listed)

    monkeypatch.setattr(execute, "_sql_paths", redirected)
    with pytest.raises(execute.PairingError, match="the projections hold"):
        execute.Corpus(
            str(root / "projections" / "*.parquet"),
            str(root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
        )


def _doubled_tree(root: pathlib.Path) -> pathlib.Path:
    """Writes a whole tree whose two datasets both state the reaction ``ord-x01``.

    Args:
        root: Where to write it.

    Returns:
        ``root``.
    """
    for shard, smiles in (("aa", "c1ccncc1"), ("bb", "CCO")):
        source = root / "data" / f"ord_dataset-{shard}.parquet"
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
        projected = root / "projections" / source.name
        projected.parent.mkdir(parents=True, exist_ok=True)
        projection.write_projection(source, projected)
        structured = root / "structures" / source.name
        structured.parent.mkdir(parents=True, exist_ok=True)
        structures.write_structures(projected, structured)
    _write_pivots(root, _OCCURRENCE_LEVELS)
    _write_occurrences(root / "pivots", root / "occurrences")
    return root


def test_warming_refuses_at_open_what_the_index_would_refuse_on_a_query(tmp_path):
    # Two datasets stating one reaction pair cleanly and only the index refuses them,
    # so this is the refusal warming is able to move. Opening such a corpus and failing
    # later leaves one that looks searchable and is not.
    with pytest.raises(execute.PairingError, match="ord-x01 is stated by 2 rows"):
        _open(_doubled_tree(tmp_path), resolver={}.__getitem__)


def test_checking_the_index_publishes_it_and_counts_every_path(corpus_dir):
    # The point of the check is that it happens at startup rather than inside whichever
    # query first wants the index, so it has to leave the index published. Every
    # indexed path is reported, including the ones this corpus records nothing at: a
    # zero is ordinary, and omitting it would read as a path the walk lost.
    with _open(corpus_dir, resolver={}.__getitem__) as value:
        counts = value.check_index()
        assert value._indexed is True
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
        _open(
            root, resolver={}.__getitem__, require_current=False, warm=False
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
    # Every other role condition here asks about inputs.components. A role that came
    # out right for the inputs and NULL elsewhere would answer "pyridine as the solvent
    # in the workup" with nothing, silently, while the same query without the role
    # still worked.
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
    body = _ROUTABLE["structure alone"]
    assert _index_spent(body), "expected the index to take the clause"
    request = query.Query.model_validate(body)
    assert _reactions(corpus.search(request, timeout_seconds=60)) == {
        "ord-aa01",
        "ord-aa02",
    }
    # What reaches the query is what the earlier phases left, not the whole bound: the
    # screen and the index check come out of the same budget.
    (spent,) = [timeout for _, timeout in seen]
    assert 0 < spent <= 60
    assert "FROM occurrences" in seen[0][0]


def test_a_spent_clause_reads_the_index_and_a_declined_one_does_not(
    corpus_dir, tmp_path
):
    # Everything else here compares the two paths, which agree -- so a search that
    # quietly stopped spending the index would keep passing while the index answered
    # nothing. A row only the index holds separates them: the semi-join filters
    # reactions by membership, so redirecting an existing reaction through a planted
    # occurrence row changes the indexed answer and cannot change the pivot's.
    root = _copy_corpus(corpus_dir, tmp_path)
    with _open(root, resolver={}.__getitem__) as value:
        spent = query.Query.model_validate(_ROUTABLE["structure alone"])
        assert _reactions(value.search(spent)) == {"ord-aa01", "ord-aa02"}
        # Copies aa01's occurrence rows, its pyridine among them, and attributes them to
        # bb01, which holds no pyridine at all. Written into the artifact, which the
        # index is a view over, since a view takes no INSERT.
        artifact = (
            root / "occurrences" / "inputs.components" / "aa" / "ord_dataset-aa.parquet"
        )
        table = pq.read_table(artifact)
        planted = table.filter(
            pa.compute.equal(table.column("reaction_id"), "ord-aa01")
        ).set_column(
            0,
            table.schema.field("reaction_id"),
            pa.array(["ord-bb01"] * 2, type=pa.string()),
        )
        pq.write_table(pa.concat_tables([table, planted]), artifact)
        assert _reactions(value.search(spent)) == {"ord-aa01", "ord-aa02", "ord-bb01"}
        # The same structure bound to the element's role, which no occurrence row
        # carries, so the quantifier compiles over the elements and the planted row is
        # invisible.
        declined_body = _NOT_ROUTABLE["structure and role"]
        assert not _index_spent(declined_body)
        declined = query.Query.model_validate(declined_body)
        assert _reactions(value.search(declined)) == {"ord-aa01"}


def test_concurrent_first_searches_publish_one_index(corpus_dir, caplog, monkeypatch):
    # Publishing the index scans the artifacts to check them. First searches arriving
    # together must not each publish their own; whoever waits takes the checked one.
    # The barrier sits at the door of the publish rather than in the resolver, which a
    # search reaches only afterwards, so the race is arranged rather than hoped for --
    # and the publish's own log line is what counts them, since concurrent CREATEs
    # would also collide in DuckDB's catalog and raise.
    caplog.set_level(logging.INFO, logger="ord_schema.search.execute")
    threads_count = 4
    ready = threading.Barrier(threads_count)
    publishing = execute.Corpus._occurrences

    def synchronized(self):
        ready.wait(timeout=_WAIT)
        return publishing(self)

    monkeypatch.setattr(execute.Corpus, "_occurrences", synchronized)
    with _open(
        corpus_dir, resolver={"pyridine": "c1ccncc1"}.__getitem__, warm=False
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
    published = [
        record for record in caplog.records if record.message.startswith("indexed ")
    ]
    assert failures == []
    assert alive == []
    assert len(published) == 1
    assert results == [2] * threads_count


def test_an_index_that_reaches_nothing_is_refused(corpus_dir, monkeypatch):
    # An index that reaches no structure at all answers every structure query with no
    # matches, which reads like an answer rather than like a corpus nobody can search.
    # The fixture has no authentic standards, so an index over that path alone reaches
    # nothing -- what artifacts written from the wrong schema would look like.
    path = "outcomes.products.measurements.authentic_standard"
    monkeypatch.setattr(execute, "INDEXED_PATHS", (path,))
    request = query.Query.model_validate(
        {
            "where": {
                "op": "exists",
                "path": path,
                "where": {"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"},
            }
        }
    )
    with _open(corpus_dir, resolver={}.__getitem__, warm=False) as value:
        with pytest.raises(execute.PairingError, match="reached 0 of the corpus's 5"):
            value.search(request)
        # The view is in the catalog, reaching nothing; what was not set is the flag,
        # so nobody reads it.
        assert value._indexed is None
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
    with (
        _open(_doubled_tree(tmp_path), resolver={}.__getitem__, warm=False) as value,
        pytest.raises(execute.PairingError, match="ord-x01 is stated by 2 rows"),
    ):
        value.search(query.Query.model_validate(_ROUTABLE["structure alone"]))


def test_a_reaction_two_files_both_state_is_refused_by_a_pivot(tmp_path):
    # A pivot finds reactions by ID as the index does, so with no index at all a
    # question narrowed to aa, about a component only bb's copy holds, would answer
    # with aa's reaction. Refused when the first pivot is published.
    root = _doubled_tree(tmp_path)
    narrowed = {
        "op": "and",
        "clauses": [
            {"op": "eq", "path": "dataset_id", "value": {"literal": "ord_dataset-aa"}},
            _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}}),
        ],
    }
    with (
        _open(root, resolver={}.__getitem__, occurrences_dir=None, warm=False) as value,
        pytest.raises(execute.PairingError, match="ord-x01 is stated by 2 rows"),
    ):
        _search(value, narrowed)


def test_warming_refuses_two_files_stating_a_reaction_without_an_index(tmp_path):
    # Warming publishes no index here, so the check has to run on its own: the
    # corpus a pivot would refuse on its first query is refused at open.
    root = _doubled_tree(tmp_path)
    with pytest.raises(execute.PairingError, match="ord-x01 is stated by 2 rows"):
        _open(root, resolver={}.__getitem__, occurrences_dir=None)


def test_an_index_that_misses_one_path_is_refused(corpus_dir, monkeypatch):
    # The dangerous shape is not an index that reaches nothing -- it is one that reaches
    # almost everything. A path whose artifact finds no element leaves the other paths
    # carrying the total, so a count of occurrences looks healthy while every query at
    # the dead path answers with silence. Counting the structures reached is what tells
    # the two apart.
    kept = tuple(path for path in execute.INDEXED_PATHS if path != "outcomes.products")
    monkeypatch.setattr(execute, "INDEXED_PATHS", kept)
    # Toluene is every reaction's product and nobody's input, so dropping the products
    # path leaves exactly its two structures unreached.
    with (
        _open(corpus_dir, resolver={}.__getitem__, warm=False) as value,
        pytest.raises(
            execute.PairingError,
            match=r"reached 3 of the corpus's 5 structures, .*an occurrence artifact "
            "missing rows",
        ),
    ):
        value.search(query.Query.model_validate(_ROUTABLE["structure alone"]))


def test_a_path_the_index_does_not_cover_is_left_to_the_projection(
    corpus_dir, monkeypatch
):
    # The planner's path check is unfalsifiable against the real schema, since every
    # path a structure can sit at is indexed. It becomes load-bearing the moment that
    # stops being true, and routing a query to a path the index does not carry answers
    # it with the empty set.
    kept = tuple(path for path in execute.INDEXED_PATHS if path != "outcomes.products")
    monkeypatch.setattr(execute, "INDEXED_PATHS", kept)
    body = _ROUTABLE["products rather than inputs"]
    assert not _index_spent(body)
    request = query.Query.model_validate(body)
    with _open(corpus_dir, resolver={}.__getitem__, warm=False) as value:
        # Answered without the index, which is therefore never published.
        assert _reactions(value.search(request)) == {"ord-aa01", "ord-aa02", "ord-bb01"}
        assert value._indexed is None


def test_the_indexed_paths_are_the_ones_a_structure_can_sit_at():
    # Derived from the projection schema, through the artifact's own walk of it.
    assert sorted(execute.INDEXED_PATHS) == [
        "inputs.components",
        "outcomes.products",
        "outcomes.products.measurements.authentic_standard",
        "workups.input.components",
    ]


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


def _tables(value) -> set[str]:
    """Returns the names of the tables the corpus connection holds."""
    cursor = value._connection.cursor()
    try:
        rows = cursor.execute("SELECT table_name FROM duckdb_tables()").fetchall()
    finally:
        cursor.close()
    return {row[0] for row in rows}


def _in_a_thread(work) -> tuple[threading.Thread, list]:
    """Runs ``work`` in a thread, and returns it with the list its result lands in."""
    done: list = []
    thread = threading.Thread(target=lambda: done.append(work()))
    thread.start()
    return thread, done


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


def test_a_pattern_and_a_name_that_read_alike_are_not_one_answer(tmp_path):
    # A stated substructure pattern is SMARTS; a resolved compound is SMILES. The same
    # text is two query molecules across those parsers: [NH4+] read as SMARTS holds its
    # nitrogen to four hydrogens, and read as SMILES leaves the count unconstrained, so
    # it matches tetramethylammonium too. A resolved name is keyed in its canonical
    # spelling, which for ammonium is that same text, so a key holding only the text
    # would answer one with the other.
    _derive(
        tmp_path,
        "nh",
        [
            _reaction("ord-nh01", components=[("[NH4+]", _ROLE.REACTANT)]),
            _reaction("ord-nh02", components=[("C[N+](C)(C)C", _ROLE.REACTANT)]),
        ],
    )
    by_name = query.Query.model_validate(
        {"where": _exists({"op": "substructure", "path": "smiles", "compound": "a"})}
    )

    def answers(first, second):
        with _open(tmp_path, resolver={"a": "[NH4+]"}.__getitem__) as value:
            return _reactions(value.search(first)), _reactions(value.search(second))

    # Whichever runs first, both hold.
    both = {"ord-nh01", "ord-nh02"}
    assert answers(_substructure("[NH4+]"), by_name) == ({"ord-nh01"}, both)
    assert answers(by_name, _substructure("[NH4+]")) == (both, {"ord-nh01"})


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

    Over the whole corpus, no reaction records a NULL or empty ``outcomes`` and few
    record a NULL ``products`` under one, so corpus agreement is weak evidence about
    levels that are absent, empty, or carry a NULL leaf. These say it
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

    # A percentage that is not a yield, beside one that is. Selectivity, purity, and
    # area are recorded as percentages too, so a reduction over the bare path answers a
    # broader question than "yield" asks.
    mixed = reaction_pb2.Reaction(reaction_id="ord-wd08")
    only = mixed.outcomes.add().products.add(is_desired_product=True)
    only.measurements.add(type="YIELD").percentage.value = 40
    only.measurements.add(type="SELECTIVITY").percentage.value = 95
    reactions.append(mixed)

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
    pivots = _write_pivots(root, sorted(pivot.LEVELS))
    with execute.Corpus(
        str(projected),
        str(structured),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
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
    monkeypatch.setattr(execute.Corpus, "_pivot_view", _no_pivot_view, raising=True)
    assert _search(wide_corpus, where) == pivoted


def _no_pivot_view(self, path: str) -> str | None:
    """Stands in for _pivot_view so a quantifier compiles over the elements."""
    del self, path  # Unused.
    return None


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
    name = wide_corpus._pivot_view("outcomes.products")
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
    name = wide_corpus._pivot_view("outcomes.products")
    assert name is not None
    blank = wide_corpus._connection.execute(
        f"SELECT count(*) FROM {name} WHERE reaction_id = 'ord-wd05'"  # noqa: S608
    ).fetchone()
    # A level that is present but empty contributes no rows, and one the source never
    # recorded contributes none either -- which is what makes forall vacuously true of
    # both.
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


def _reduced_per_reaction(corpus, reducer, where=None):
    """Returns the reduction each reaction answers, keyed by reaction."""
    key = {
        "reduce": reducer,
        "path": "outcomes.products.measurements.percentage.value",
    }
    if where is not None:
        key["where"] = where
    table = corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": ["reaction_id"],
                    "measures": [{"fn": "min", "path": key, "name": "reduced"}],
                }
            }
        )
    )
    return {row["reaction_id"]: row["reduced"] for row in table.to_pylist()}


@pytest.mark.parametrize("reducer", ["min", "max", "avg", "sum", "count"])
def test_a_pivoted_reduction_answers_what_the_projection_answers(
    wide_root, tmp_path, reducer
):
    # The pivot reaches the same elements by another route, so it cannot answer
    # differently -- including for the reactions that carry no element at all, which
    # is where the two spellings of "nothing" diverge if the count is left to len().
    pivots = _write_pivots(
        wide_root, ("outcomes.products.measurements",), into=tmp_path / "pivots"
    )
    projected = str(wide_root / "projections" / "*.parquet")
    structured = str(wide_root / "structures" / "*.parquet")
    with execute.Corpus(projected, structured, resolver={}.__getitem__) as lists:
        expected = _reduced_per_reaction(lists, reducer)
    with execute.Corpus(
        projected, structured, resolver={}.__getitem__, pivots_dir=str(pivots)
    ) as pivoted:
        assert _reduced_per_reaction(pivoted, reducer) == expected


@pytest.mark.parametrize("reducer", ["min", "max", "avg", "sum", "count"])
def test_a_filtered_reduction_answers_the_same_on_either_route(
    wide_root, tmp_path, reducer
):
    # The filter narrows the same elements whichever relation holds them, and the
    # wide corpus records measurements that are not yields for exactly this.
    pivots = _write_pivots(
        wide_root, ("outcomes.products.measurements",), into=tmp_path / "pivots"
    )
    projected = str(wide_root / "projections" / "*.parquet")
    structured = str(wide_root / "structures" / "*.parquet")
    where = {"op": "eq", "path": "type", "value": {"literal": "YIELD"}}
    with execute.Corpus(projected, structured, resolver={}.__getitem__) as lists:
        expected = _reduced_per_reaction(lists, reducer, where)
    with execute.Corpus(
        projected, structured, resolver={}.__getitem__, pivots_dir=str(pivots)
    ) as pivoted:
        assert _reduced_per_reaction(pivoted, reducer, where) == expected


def test_a_reduction_filter_may_screen_structures(deep_root, tmp_path):
    # A structure test reads structure_offset, which the reactions carry and a pivot
    # row does not. Inside the correlated subquery the unqualified name resolves to the
    # outer relation, which is the offset that makes the element's dataset-local
    # structure_id corpus-wide; getting that wrong screens against another file's
    # molecules and answers plausibly.
    pivots = _write_pivots(deep_root, ("inputs.components",), into=tmp_path / "pivots")
    counted = {
        "aggregate": {
            "group_by": ["reaction_id"],
            "measures": [
                {
                    "fn": "min",
                    "path": {
                        "reduce": "count",
                        "path": "inputs.components.reaction_role",
                        "where": {
                            "op": "substructure",
                            "path": "smiles",
                            "smarts": "[OX2H]",
                        },
                    },
                    "name": "alcohols",
                }
            ],
        }
    }
    answers = {}
    for pivots_dir in (None, str(pivots)):
        with execute.Corpus(
            str(deep_root / "projections" / "*.parquet"),
            str(deep_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=pivots_dir,
        ) as corpus:
            for smarts in ("[OX2H]", "[n]"):
                counted["aggregate"]["measures"][0]["path"]["where"]["smarts"] = smarts
                table = corpus.search(query.Query.model_validate(counted))
                answers[pivots_dir is not None, smarts] = {
                    row["reaction_id"]: row["alcohols"] for row in table.to_pylist()
                }
    # Every input here is ethanol, and the pyridine is a workup's. So the alcohol
    # screen keeps both -- an agreement on zeros would say nothing about the offset --
    # and the aromatic-nitrogen screen keeps neither, which is what a dropped filter
    # could not produce.
    for pivoted in (False, True):
        assert answers[pivoted, "[OX2H]"] == {"ord-cc01": 1, "ord-cc02": 1}
        assert answers[pivoted, "[n]"] == {"ord-cc01": 0, "ord-cc02": 0}


@pytest.mark.parametrize("pivoted", [False, True])
def test_a_filter_changes_which_elements_are_reduced(wide_root, tmp_path, pivoted):
    # The point of the filter, and the reason the bare path is the wrong question:
    # ord-wd08 records a selectivity percentage beside its yield, so counting every
    # percentage counts two where the yield is one, and the highest is the selectivity.
    pivots = _write_pivots(
        wide_root, ("outcomes.products.measurements",), into=tmp_path / "pivots"
    )
    yields = {"op": "eq", "path": "type", "value": {"literal": "YIELD"}}
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots) if pivoted else None,
    ) as corpus:
        assert _reduced_per_reaction(corpus, "count")["ord-wd08"] == 2
        assert _reduced_per_reaction(corpus, "count", yields)["ord-wd08"] == 1
        assert _reduced_per_reaction(corpus, "max")["ord-wd08"] == 95
        assert _reduced_per_reaction(corpus, "max", yields)["ord-wd08"] == 40


@pytest.mark.parametrize("pivoted", [False, True])
def test_counting_a_level_a_reaction_lacks_is_zero(wide_root, tmp_path, pivoted):
    # ord-wd04 records no outcomes, which the projection writes as NULL rather than as
    # an empty list, and ord-wd03 an outcome carrying no products. Neither reaction
    # holds an unknown number of measurements; both hold none.
    pivots = _write_pivots(
        wide_root, ("outcomes.products.measurements",), into=tmp_path / "pivots"
    )
    with execute.Corpus(
        str(wide_root / "projections" / "*.parquet"),
        str(wide_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots) if pivoted else None,
    ) as corpus:
        counted = _reduced_per_reaction(corpus, "count")
    assert counted["ord-wd04"] == 0
    assert counted["ord-wd03"] == 0


def _integral_reaction(reaction_id, *, rpm, kelvin, inputs):
    """A reaction stirred at ``rpm`` and held at ``kelvin``, adding ``inputs`` in turn.

    ``inputs`` maps each input's name to its addition order and the SMILES of its
    components.
    """
    reaction = reaction_pb2.Reaction(reaction_id=reaction_id)
    reaction.conditions.stirring.rate.rpm = rpm
    setpoint = reaction.conditions.temperature.setpoint
    setpoint.value = kelvin
    setpoint.units = reaction_pb2.Temperature.KELVIN
    for name, (order, smiles) in inputs.items():
        reaction.inputs[name].addition_order = order
        for each in smiles:
            _component(reaction.inputs[name], each, _ROLE.REACTANT)
    return reaction


@pytest.fixture(scope="module")
def integral_root(tmp_path_factory) -> pathlib.Path:
    """Reactions recording integer leaves, which DuckDB sums as HUGEINT."""
    reactions = [
        _integral_reaction(
            "ord-in01",
            rpm=300,
            kelvin=350.5,
            inputs={"a": (1, ["CCO", "CO"]), "b": (2, ["c1ccccc1"])},
        ),
        _integral_reaction(
            "ord-in02", rpm=200, kelvin=300.25, inputs={"a": (5, ["CCO"])}
        ),
    ]
    root = tmp_path_factory.mktemp("integral")
    source = root / "data" / "ord_dataset-in.parquet"
    source.parent.mkdir(parents=True)
    parquet.save_dataset(
        dataset_pb2.Dataset(
            dataset_id="ord_dataset-in",
            name="test",
            description="test",
            reactions=reactions,
        ),
        str(source),
    )
    projected = root / "projections" / source.name
    projected.parent.mkdir(parents=True)
    projection.write_projection(source, projected)
    structured = root / "structures" / source.name
    structured.parent.mkdir(parents=True)
    structures.write_structures(projected, structured)
    _write_pivots(root, sorted(pivot.LEVELS))
    return root


@pytest.fixture(scope="module", params=["pivoted", "unpivoted"])
def integral_corpus(request, integral_root) -> Iterator[execute.Corpus]:
    """The integral reactions, reduced through pivot artifacts and through the lists."""
    with execute.Corpus(
        str(integral_root / "projections" / "*.parquet"),
        str(integral_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=_pivots_for(integral_root, request.param),
    ) as value:
        yield value


def _measured(corpus, measure, group_by=()) -> pa.Table:
    """Returns the table a single measure named ``measured`` comes back in."""
    return corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": list(group_by),
                    "measures": [{**measure, "name": "measured"}],
                }
            }
        )
    )


def test_a_sum_over_an_integer_leaf_is_an_int(integral_corpus):
    # DuckDB sums an INTEGER as HUGEINT, which Arrow carries as decimal128(38, 0) and
    # Python reads as a Decimal, and json.dumps(default=str) then writes as "500".
    table = _measured(
        integral_corpus, {"fn": "sum", "path": "conditions.stirring.rate.rpm"}
    )
    assert table.schema.field("measured").type == pa.int64()
    (total,) = table.column("measured").to_pylist()
    assert total == 500
    assert isinstance(total, int)


def test_a_sum_over_a_count_is_an_int(integral_corpus):
    # A sum of counts beside a count of the same rows, both answering "how many".
    count = {"reduce": "count", "path": "inputs.components"}
    table = _measured(integral_corpus, {"fn": "sum", "path": count})
    assert table.schema.field("measured").type == pa.int64()
    (total,) = table.column("measured").to_pylist()
    assert total == 4
    assert isinstance(total, int)


def test_a_summed_reduction_is_an_int_where_it_surfaces(integral_corpus):
    # Read through a min rather than a sum, which would cast the reduction again and
    # hide whether the reduction itself comes back an integer.
    summed = {"reduce": "sum", "path": "inputs.addition_order"}
    table = _measured(integral_corpus, {"fn": "min", "path": summed}, ["reaction_id"])
    assert table.schema.field("measured").type == pa.int64()
    per_reaction = {row["reaction_id"]: row["measured"] for row in table.to_pylist()}
    assert per_reaction == {"ord-in01": 3, "ord-in02": 5}
    assert all(isinstance(value, int) for value in per_reaction.values())


def test_a_sum_over_a_floating_point_leaf_stays_a_float(integral_corpus):
    table = _measured(
        integral_corpus,
        {"fn": "sum", "path": "conditions.temperature.setpoint_kelvin"},
    )
    assert table.schema.field("measured").type == pa.float64()
    (total,) = table.column("measured").to_pylist()
    assert total == 650.75
    assert isinstance(total, float)


_MODIFIED = "provenance.record_modified.time.timestamp"


def _modified_reaction(reaction_id, *stamps) -> reaction_pb2.Reaction:
    """A reaction whose record was modified at each of ``stamps``."""
    reaction = _reaction(reaction_id, components=[("CCO", _ROLE.REACTANT)])
    for stamp in stamps:
        reaction.provenance.record_modified.add().time.value = stamp
    return reaction


@pytest.fixture(scope="module")
def modified_root(tmp_path_factory) -> pathlib.Path:
    """Three datasets whose earliest and latest modifications disagree on an order.

    aa was submitted first and edited last, so ranking datasets by their latest
    modification and by their earliest give different answers. cc records none.
    """
    root = tmp_path_factory.mktemp("modified")
    _derive(
        root,
        "aa",
        [
            _modified_reaction("ord-md01", "2021-03-01", "2024-06-01"),
            _modified_reaction("ord-md02", "2021-03-02"),
        ],
    )
    _derive(
        root,
        "bb",
        [
            _modified_reaction("ord-md03", "2023-01-14 09:30:00", "2023-02-01"),
            _modified_reaction("ord-md04", "2023-01-15 10:00:00"),
        ],
    )
    _derive(root, "cc", [_modified_reaction("ord-md05")])
    _write_pivots(root, sorted(pivot.LEVELS))
    return root


@pytest.fixture(scope="module", params=["pivoted", "unpivoted"])
def modified_corpus(request, modified_root) -> Iterator[execute.Corpus]:
    """The modified reactions, reduced through pivot artifacts and through the lists."""
    with execute.Corpus(
        str(modified_root / "projections" / "*.parquet"),
        str(modified_root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=_pivots_for(modified_root, request.param),
    ) as value:
        yield value


@pytest.mark.parametrize(
    ("fn", "reduce", "expected"),
    [
        ("min", "min", {"ord_dataset-bb": datetime.date(2023, 1, 14)}),
        ("max", "max", {"ord_dataset-aa": datetime.date(2024, 6, 1)}),
        ("max", "min", {"ord_dataset-bb": datetime.date(2023, 1, 15)}),
        ("min", "max", {"ord_dataset-bb": datetime.date(2023, 1, 15)}),
    ],
)
def test_datasets_rank_by_a_day_reduced_from_their_reactions(
    modified_corpus, fn, reduce, expected
):
    # min of each reaction's min is when the dataset was first recorded, and max of
    # each reaction's max is when it was last edited; aa is the newest by one and not
    # the other. The mixed pairs answer differently from both, so neither the measure
    # nor the reduction can be ignored. cc records no day, so its NULL sorts last. Each
    # answer is the day alone, whatever time of it the record wrote.
    table = modified_corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": ["dataset_id"],
                    "measures": [
                        {
                            "fn": fn,
                            "path": {"reduce": reduce, "path": _MODIFIED},
                            "name": "modified",
                        }
                    ],
                },
                "order_by": [{"key": "modified", "descending": True}],
                "limit": 1,
            }
        )
    )
    assert table.schema.field("modified").type == pa.date32()
    assert {row["dataset_id"]: row["modified"] for row in table.to_pylist()} == expected


def test_reactions_order_by_their_latest_modification(modified_corpus):
    table = modified_corpus.search(
        query.Query.model_validate(
            {
                "order_by": [
                    {
                        "key": {"reduce": "max", "path": _MODIFIED},
                        "descending": True,
                    }
                ]
            }
        )
    )
    assert table.column("reaction_id").to_pylist() == [
        "ord-md01",
        "ord-md03",
        "ord-md04",
        "ord-md02",
        "ord-md05",
    ]


def test_an_edit_on_a_day_matches_whatever_time_it_was_made(modified_corpus):
    # md04 was edited at 10:00 that day, so only a comparison by day finds it; an
    # instant compared against the day's midnight would find nothing. Read through the
    # pivot and the lists alike.
    where = {
        "op": "exists",
        "path": "provenance.record_modified",
        "where": {
            "op": "eq",
            "path": "time.timestamp",
            "value": {"literal": "2023-01-15"},
        },
    }
    assert _search(modified_corpus, where) == {"ord-md04"}


_CREATED = "provenance.record_created.time.timestamp"


def _created_reaction(reaction_id, stamp) -> reaction_pb2.Reaction:
    """A reaction whose record was created at ``stamp``."""
    reaction = _reaction(reaction_id, components=[("CCO", _ROLE.REACTANT)])
    reaction.provenance.record_created.time.value = stamp
    return reaction


@pytest.fixture(scope="module")
def created_corpus(tmp_path_factory) -> Iterator[execute.Corpus]:
    """Reactions on three days: three on March 1 and one a second either side.

    A comparison finer than a day tells dy02, dy03, and dy04 apart; one by day cannot.
    """
    root = tmp_path_factory.mktemp("created")
    _derive(
        root,
        "aa",
        [
            _created_reaction("ord-dy01", "2024-02-29 23:59:59"),
            _created_reaction("ord-dy02", "2024-03-01 00:00:00"),
            _created_reaction("ord-dy03", "2024-03-01 12:30:00"),
            _created_reaction("ord-dy04", "2024-03-01T23:59:59"),
            _created_reaction("ord-dy05", "2024-03-02 00:00:00"),
        ],
    )
    with execute.Corpus(
        str(root / "projections" / "*.parquet"),
        str(root / "structures" / "*.parquet"),
        resolver={}.__getitem__,
    ) as value:
        yield value


@pytest.mark.parametrize(
    ("op", "expected"),
    [
        ("eq", {"ord-dy02", "ord-dy03", "ord-dy04"}),
        ("ne", {"ord-dy01", "ord-dy05"}),
        ("lt", {"ord-dy01"}),
        ("le", {"ord-dy01", "ord-dy02", "ord-dy03", "ord-dy04"}),
        ("gt", {"ord-dy05"}),
        ("ge", {"ord-dy02", "ord-dy03", "ord-dy04", "ord-dy05"}),
    ],
)
def test_a_timestamp_compares_by_day(created_corpus, op, expected):
    # eq takes every time on the day, le takes the whole day in, and gt leaves it all
    # out. Compared as instants against the day's midnight, eq would find dy02 alone and
    # le would stop there, with gt taking in the rest of the day.
    where = {"op": op, "path": _CREATED, "value": {"literal": "2024-03-01"}}
    assert _search(created_corpus, where) == expected


def test_a_search_refuses_a_literal_carrying_a_time(created_corpus):
    where = {"op": "le", "path": _CREATED, "value": {"literal": "2024-03-01T12:00:00"}}
    with pytest.raises(query.QueryError, match="timestamps compare by day"):
        _search(created_corpus, where)


def test_reactions_group_by_the_day_they_were_created(created_corpus):
    table = created_corpus.search(
        query.Query.model_validate(
            {
                "aggregate": {
                    "group_by": [_CREATED],
                    "measures": [{"fn": "count", "name": "reactions"}],
                }
            }
        )
    )
    assert table.schema.field(_CREATED).type == pa.date32()
    assert {row[_CREATED]: row["reactions"] for row in table.to_pylist()} == {
        datetime.date(2024, 2, 29): 1,
        datetime.date(2024, 3, 1): 3,
        datetime.date(2024, 3, 2): 1,
    }


@pytest.mark.parametrize(
    ("fn", "expected"),
    [
        ("min", datetime.date(2024, 2, 29)),
        ("max", datetime.date(2024, 3, 2)),
        ("count_distinct", 3),
    ],
)
def test_measures_over_a_timestamp_read_its_day(created_corpus, fn, expected):
    # Five instants, three days: a distinct count counts the days.
    table = _measured(created_corpus, {"fn": fn, "path": _CREATED})
    if fn != "count_distinct":
        assert table.schema.field("measured").type == pa.date32()
    assert table.column("measured").to_pylist() == [expected]


def test_reactions_order_by_the_day_they_were_created(created_corpus):
    # Ties within a day fall to the second key, which an ordering by instant would never
    # reach: by instant, descending, dy04 would lead March 1 rather than trail it.
    table = created_corpus.search(
        query.Query.model_validate(
            {
                "order_by": [
                    {"key": _CREATED, "descending": True},
                    {"key": "reaction_id"},
                ]
            }
        )
    )
    assert table.column("reaction_id").to_pylist() == [
        "ord-dy05",
        "ord-dy02",
        "ord-dy03",
        "ord-dy04",
        "ord-dy01",
    ]


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


def test_a_role_beside_the_structure_is_answered_by_the_pivot(tmp_path, corpus_dir):
    # The role binds a query to one element: aa01 holds benzene as a reactant and
    # pyridine as its solvent. No occurrence row carries it, so the pivot answers, and
    # the index is never built for the asking. A role read wrong from the pivot answers
    # "pyridine as the solvent" with silence and "benzene as the solvent" with aa01,
    # neither of which raises.
    pivots = _write_pivots(corpus_dir, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    with _indexed_corpus(corpus_dir, pivots, warm=False) as corpus:
        assert _search(corpus, _role_and_structure("c1ccncc1", "SOLVENT")) == {
            "ord-aa01"
        }
        assert _search(corpus, _role_and_structure("c1ccccc1", "SOLVENT")) == set()
        assert "inputs.components" in corpus._pivot_views
        assert corpus._indexed is None


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


def test_requiring_pivots_refuses_a_level_with_no_artifacts(wide_root):
    # The silent case: a missing level is otherwise answered from the projection's
    # lists on every query over it, slower and with no error anywhere.
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


def test_a_level_without_artifacts_is_answered_from_the_projection(wide_root, caplog):
    # A partial set of artifacts is a partial speedup, not a missing answer: the level
    # with no artifact is read from the projection's lists, and the log says so.
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
        assert "outcomes.products.measurements" not in corpus._pivot_views
    assert found == {"ord-wd01", "ord-wd02", "ord-wd07", "ord-wd08"}
    assert "this query reads the projection alone" in _messages(caplog)


@pytest.mark.parametrize(
    ("levels", "logged"),
    [
        # The inner level's pivot is looked up before the outer level declines, and the
        # whole quantifier is then answered from the lists.
        (("outcomes.products.measurements",), "the projection alone"),
        (
            ("outcomes.products", "outcomes.products.measurements"),
            "the pivots over ['outcomes.products', 'outcomes.products.measurements']",
        ),
    ],
)
def test_the_log_names_the_pivots_the_query_reads(
    wide_root, tmp_path, caplog, levels, logged
):
    # Someone asking why a query was slow reads this line, so a pivot named there and
    # absent from the SQL sends them to the wrong relation.
    pivots = _write_pivots(wide_root, levels, into=tmp_path / "pivots")
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
            corpus, _DIFFERENTIAL["a desired product with a yield above 50%"]
        )
    assert found
    assert f"this query reads {logged}" in _messages(caplog)


def _messages(caplog) -> list[str]:
    """Returns what the executor logged."""
    return [
        record.message
        for record in caplog.records
        if record.name == "ord_schema.search.execute"
    ]


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


def test_a_pivot_row_naming_no_reaction_does_not_empty_a_negated_exists(
    corpus_dir, tmp_path
):
    # An exists over a pivot is spelled IN, and NOT IN against a subquery holding a
    # NULL answers NULL for every reaction, returning nothing. The subquery leaves such
    # rows out, so the negation still answers every reaction without a solvent.
    pivots = _write_pivots(corpus_dir, ("inputs.components",), into=tmp_path / "pivots")
    # The aa shard, which holds the one solvent.
    target = pivots / "inputs.components" / "aa" / "ord_dataset-aa.parquet"
    table = pq.read_table(target)
    roles = pa.compute.struct_field(table.column("element"), "reaction_role")
    solvent = table.filter(pa.compute.equal(roles, "SOLVENT")).slice(0, 1)
    assert solvent.num_rows == 1
    unnamed = solvent.set_column(
        0, table.schema.field("reaction_id"), pa.array([None], type=pa.string())
    )
    pq.write_table(
        pa.concat_tables([table, unnamed]).replace_schema_metadata(
            table.schema.metadata
        ),
        target,
    )
    with execute.Corpus(
        str(corpus_dir / "projections" / "*.parquet"),
        str(corpus_dir / "structures" / "*.parquet"),
        resolver={}.__getitem__,
        pivots_dir=str(pivots),
        warm=False,
    ) as corpus:
        assert _search(corpus, {"op": "not", "clause": _exists(_SOLVENT)}) == {
            "ord-aa02",
            "ord-bb01",
        }


def test_a_stranger_s_pivot_is_refused_by_the_check(wide_root, corpus_dir):
    # check_pivots publishes every level at startup, so a pivots_dir belonging to
    # another corpus fails the deployment there rather than whichever query first
    # quantified over the level.
    stranger = _write_pivots(corpus_dir, ("inputs.components",))
    with (
        execute.Corpus(
            str(wide_root / "projections" / "*.parquet"),
            str(wide_root / "structures" / "*.parquet"),
            resolver={}.__getitem__,
            pivots_dir=str(stranger),
        ) as corpus,
        pytest.raises(execute.PairingError, match="another corpus"),
    ):
        corpus.check_pivots()


def test_two_pivots_of_one_dataset_under_a_level_are_refused(corpus_dir, tmp_path):
    # A stray copy under the level satisfies the set comparison over source hashes, and
    # both artifacts are then read. What that costs is in _check_pivots; what makes it
    # worth a refusal here is that the corpus opens and answers, so nothing else says.
    pivots = _write_pivots(corpus_dir, ("inputs.components",), into=tmp_path / "pivots")
    stray = pivots / "inputs.components" / "again"
    stray.mkdir(parents=True)
    for artifact in sorted((pivots / "inputs.components").rglob("*.parquet")):
        shutil.copy(artifact, stray / artifact.name)
    with (
        _open(corpus_dir, resolver={}.__getitem__, pivots_dir=str(pivots)) as corpus,
        pytest.raises(execute.PairingError, match="the same source dataset"),
    ):
        corpus.check_pivots()


def test_a_stale_pivot_is_refused(corpus_dir, tmp_path):
    # The stamps are the only thing saying an artifact was written by the library
    # reading it; a pivot derived by an older one can hold columns this schema walk
    # does not find.
    pivots = _write_pivots(corpus_dir, ("inputs.components",), into=tmp_path / "pivots")
    target = next((pivots / "inputs.components").rglob("*.parquet"))
    with pq.ParquetFile(target) as artifact:
        table = artifact.read()
        schema = artifact.schema_arrow
    metadata = dict(schema.metadata)
    metadata[b"ord.rdkit_version"] = b"0000.00.0"
    pq.write_table(table.cast(schema.with_metadata(metadata)), target)
    with (
        _open(corpus_dir, resolver={}.__getitem__, pivots_dir=str(pivots)) as corpus,
        pytest.raises(execute.PairingError, match="is a stale"),
    ):
        corpus.check_pivots()


def test_a_pivot_missing_a_column_the_schema_declares_is_refused(corpus_dir, tmp_path):
    # The element struct grows whenever the projection's does, and the stamps say
    # nothing about columns -- so an artifact derived before the column is current by
    # every other measure here. Read as a corpus member it fails deep in DuckDB, as a
    # binder error naming a struct key rather than the file that predates it, which
    # sends whoever reads it looking at the query instead of the tree.
    pivots = _write_pivots(corpus_dir, ("inputs.components",), into=tmp_path / "pivots")
    target = next((pivots / "inputs.components").rglob("*.parquet"))
    table = pq.read_table(target)
    metadata = table.schema.metadata
    element = table.column("element").combine_chunks()
    dropped = element.field(element.type.field(0).name)
    table = table.drop_columns(["element"]).append_column(
        pa.field("element", pa.struct([element.type.field(0)])),
        pa.StructArray.from_arrays([dropped], fields=[element.type.field(0)]),
    )
    pq.write_table(table.replace_schema_metadata(metadata), target)
    with (
        _open(corpus_dir, resolver={}.__getitem__, pivots_dir=str(pivots)) as corpus,
        pytest.raises(execute.PairingError, match="without"),
    ):
        corpus.check_pivots()


def test_a_pivot_past_its_dataset_s_structures_is_refused(corpus_dir, tmp_path):
    # The occurrence artifacts' blind spot, one derivation earlier: a quantifier
    # answered from a pivot tests the bit at the element's ID plus its reaction's
    # offset, and an ID past the dataset lands in the next one's range.
    pivots = _write_pivots(corpus_dir, _OCCURRENCE_LEVELS, into=tmp_path / "pivots")
    target = pivot.artifact_paths(pivots, "outcomes.products")[0]
    _shift_structure_ids(target, 100)
    with (
        pytest.raises(
            execute.PairingError,
            match=re.escape(
                f"{target} carries structure_id 102 but its dataset holds only 3 "
                "structures"
            )
            + r".*ord_schema\.artifacts\.scripts\.derive_pivots --force",
        ),
        _indexed_corpus(corpus_dir, pivots, warm=False) as read,
    ):
        read.check_pivots()


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
        pytest.raises(execute.PairingError, match="holds the pivot over"),
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
    monkeypatch.setattr(execute.Corpus, "_pivot_view", _no_pivot_view)
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


def test_check_pivots_is_empty_without_a_directory(corpus_dir):
    # Not an error: every level is then answered from the projection's lists.
    with _open(
        corpus_dir, resolver={}.__getitem__, pivots_dir=None, occurrences_dir=None
    ) as value:
        assert value.check_pivots() == {}


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
    monkeypatch.setattr(execute.Corpus, "_pivot_view", _no_pivot_view)
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
        resolver={"acetic acid": "CC(=O)O"}.__getitem__,
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
        _exists({"op": "same_compound", "path": "smiles", "compound": "acetic acid"}),
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


def test_an_artifact_stamped_with_an_empty_name_is_refused_by_name(
    corpus_dir, tmp_path
):
    # load_stamps requires the artifact key to be present, not to be non-empty, so a
    # footer written by hand or truncated mid-write reaches the check with an empty
    # name. Naming the article off it must not be what fails: an IndexError raised
    # while composing the message would bury the file this is refusing.
    root = tmp_path / "unnamed"
    shutil.copytree(corpus_dir, root)
    target = next((root / "structures").glob("*.parquet"))
    table = pq.read_table(target)
    metadata = dict(table.schema.metadata)
    metadata[b"ord.artifact"] = b""
    pq.write_table(table.replace_schema_metadata(metadata), target)
    with pytest.raises(execute.PairingError, match=re.escape(str(target))):
        _open(root)


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
    # from those stamps. Reading them a second time costs a footer read per artifact
    # for values already in hand.
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


def test_an_occurrences_artifact_missing_a_column_is_refused(
    corpus_dir, occurrence_dirs
):
    # The same check _check_pivots makes, for an artifact whose columns cannot drift
    # today: the two are the same at every path. It is here
    # because the stamps say nothing about columns whatever the artifact, so the day
    # this schema grows one, every file written before it would read as current and
    # fail deep in DuckDB rather than here.
    _, occurrence_dir = occurrence_dirs
    target = occurrences.artifact_paths(occurrence_dir, "inputs.components")[0]
    table = pq.read_table(target)
    metadata = table.schema.metadata
    pq.write_table(
        table.drop_columns(["structure_id"]).replace_schema_metadata(metadata), target
    )
    with (
        pytest.raises(execute.PairingError, match="without"),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        _search(read, _white("exists"))


@pytest.mark.parametrize("smarts", ["c1ccncc1", "[OX2H]", "c1ccccc1", "[#6]"])
def test_an_index_read_from_artifacts_answers_what_the_pivots_answer(
    corpus_dir, occurrence_dirs, smarts
):
    # The artifact is the same occurrences by another route, so it cannot answer
    # differently from the pivots they were derived from. Over both datasets, because
    # an occurrence's ID is its element's plus its own dataset's offset and a corpus of
    # one dataset has an offset of zero to get wrong.
    pivots, occurrence_dir = occurrence_dirs
    where = _exists({"op": "substructure", "path": "smiles", "smarts": smarts})
    with _indexed_corpus(corpus_dir, pivots) as pivoted:
        expected = _search(pivoted, where)
    assert expected
    with _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read:
        assert _search(read, where) == expected


def test_an_index_read_from_artifacts_keeps_a_structure_in_its_own_dataset(
    corpus_dir, occurrence_dirs
):
    # An offset dropped or taken from the wrong dataset still lands inside the corpus's
    # structures, so it returns a reaction rather than raising -- the wrong one.
    pivots, occurrence_dir = occurrence_dirs
    with (
        _indexed_corpus(corpus_dir, pivots) as pivoted,
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        where = _exists({"op": "substructure", "path": "smiles", "smarts": "[OX2H]"})
        expected = _search(pivoted, where)
        assert expected
        assert _search(read, where) == expected


def test_a_covered_corpus_publishes_the_index_as_a_view(
    corpus_dir, occurrence_dirs, caplog
):
    # The rows stay on disk: the corpus holds a view over them and nothing else.
    _, occurrence_dir = occurrence_dirs
    with (
        caplog.at_level(logging.INFO, logger="ord_schema.search.execute"),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        assert sum(read.check_index().values()) > 0
        assert read._indexed is True
    assert any(
        f"as a view over {occurrence_dir}" in record.message
        for record in caplog.records
    )


def test_one_uncovered_path_leaves_the_corpus_without_an_index(
    corpus_dir, occurrence_dirs, caplog
):
    # The check that makes the index trustworthy counts the structures every path
    # reaches, so an index over some paths is no index at all, and a warning says which
    # paths are missing.
    pivots, _ = occurrence_dirs
    partial = _write_occurrences(
        pivots, pathlib.Path(str(pivots) + "-partial"), ("outcomes.products",)
    )
    with (
        caplog.at_level(logging.WARNING, logger="ord_schema.search.execute"),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(partial)) as read,
    ):
        assert read.check_index() == {}
        assert read._indexed is False
    missing = sorted(set(execute.INDEXED_PATHS) - {"outcomes.products"})
    assert any(
        f"holds no occurrences artifacts for {missing}" in message
        for message in _messages(caplog)
    )


def test_an_uncovered_path_still_answers(corpus_dir, occurrence_dirs):
    # Partial coverage is a missing index, never a missing answer: every structure
    # predicate is read from the pivots or the projection, as a corpus without any
    # artifacts reads it.
    pivots, _ = occurrence_dirs
    partial = _write_occurrences(
        pivots, pathlib.Path(str(pivots) + "-one"), ("outcomes.products",)
    )
    where = _exists({"op": "substructure", "path": "smiles", "smarts": "c1ccncc1"})
    with _indexed_corpus(corpus_dir, pivots) as pivoted:
        expected = _search(pivoted, where)
    assert expected
    with _indexed_corpus(corpus_dir, None, occurrences_dir=str(partial)) as read:
        assert _search(read, where) == expected


def test_a_directory_with_no_artifacts_is_not_an_error(corpus_dir, tmp_path):
    # Naming an empty directory reads as a corpus that has not derived them yet, which
    # is the state every corpus starts in: no index, and every answer still given.
    empty = tmp_path / "none"
    empty.mkdir()
    with _indexed_corpus(corpus_dir, None, occurrences_dir=str(empty)) as read:
        assert read.check_index() == {}
        assert _search(read, _exists(_SUBSTRUCTURE)) == {"ord-aa01", "ord-aa02"}


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
    # Every indexed path writes the same two columns, so where a file sits says
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


def _empty_products(occurrence_dir) -> None:
    """Empties one dataset's products artifact, keeping its schema and stamps."""
    # Toluene is every reaction's product and nobody's input, so this leaves exactly
    # that dataset's toluene unreached.
    target = occurrences.artifact_paths(occurrence_dir, "outcomes.products")[0]
    pq.write_table(pq.read_table(target).slice(0, 0), target)


def test_an_occurrence_artifact_short_of_rows_is_named_as_the_cause(
    corpus_dir, occurrence_dirs
):
    # A truncated file keeps the stamps and columns the reader checks, so the only sign
    # is a structure the index never reaches. Blaming the projections there sends the
    # reader to the one input that is fine.
    _, occurrence_dir = occurrence_dirs
    _empty_products(occurrence_dir)
    with (
        pytest.raises(
            execute.PairingError,
            match=(
                r"reached 4 of the corpus's 5 structures.*an occurrence artifact "
                r"missing rows.*re-derive .* with ord_schema\.artifacts\.scripts\."
                r"derive_occurrences --force"
            ),
        ),
        _indexed_corpus(
            corpus_dir,
            None,
            occurrences_dir=str(occurrence_dir),
            require_occurrences=True,
        ) as read,
    ):
        read.check_index()


def _shift_structure_ids(target: pathlib.Path, by: int) -> None:
    """Adds ``by`` to every structure ID one artifact states, keeping its stamps.

    Args:
        target: An occurrence artifact, which holds its IDs in a column of their own, or
            a pivot, which holds them in its element struct.
        by: What to add.
    """
    table = pq.read_table(target)
    if "structure_id" in table.column_names:
        name = "structure_id"
        values = table.column(name)
        shifted = pa.compute.add(values, by).cast(values.type)
    else:
        name = pivot.ELEMENT
        element = table.column(name).combine_chunks()
        shifted = pa.StructArray.from_arrays(
            [
                (
                    pa.compute.add(element.field(field.name), by).cast(field.type)
                    if field.name == "structure_id"
                    else element.field(field.name)
                )
                for field in element.type
            ],
            fields=list(element.type),
            mask=element.is_null(),
        )
    index = table.schema.get_field_index(name)
    pq.write_table(table.set_column(index, table.schema.field(index), shifted), target)


def _first_dataset(corpus) -> tuple[str, int]:
    """Returns the basename of the dataset whose IDs start at zero, and its size.

    Args:
        corpus: Any corpus over ``corpus_dir``.

    Returns:
        The basename its artifacts share, and how many structures it holds.
    """
    structured = corpus._connection.execute(
        "SELECT structures_filename FROM structure_offsets WHERE structure_offset = 0"
    ).fetchone()[0]
    return pathlib.Path(structured).name, pq.read_metadata(structured).num_rows


def _occurrence_artifact(occurrence_dir, path: str, name: str) -> pathlib.Path:
    """Returns the occurrence artifact at ``path`` with basename ``name``."""
    return next(
        found
        for found in occurrences.artifact_paths(occurrence_dir, path)
        if found.name == name
    )


def test_an_occurrence_past_the_corpus_is_refused(corpus, corpus_dir, occurrence_dirs):
    # Toluene is each dataset's last structure and every reaction's product, so this
    # moves exactly one dataset's toluene. The ID it lost and the one it gained cancel
    # in a count of distinct IDs, which accepts the index, and every structure query at
    # the path then fails on a bit past the end of the corpus.
    name, count = _first_dataset(corpus)
    _, occurrence_dir = occurrence_dirs
    target = _occurrence_artifact(occurrence_dir, "outcomes.products", name)
    _shift_structure_ids(target, 100)
    with (
        pytest.raises(
            execute.PairingError,
            match=re.escape(
                f"{target} carries structure_id {count + 99} but its dataset holds "
                f"only {count} structures"
            )
            + r".*ord_schema\.artifacts\.scripts\.derive_occurrences --force",
        ),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        read.check_index()


def test_an_occurrence_moved_into_the_next_dataset_is_refused(
    corpus, corpus_dir, occurrence_dirs
):
    # Moved one along, the first dataset's toluene lands on the next dataset's first
    # structure, inside the corpus. Stated once more at another path, the toluene stays
    # reached, so the count matches and the reaction answers as holding the other
    # dataset's molecule: a wrong answer rather than an error, which only a bound per
    # dataset sees.
    name, count = _first_dataset(corpus)
    _, occurrence_dir = occurrence_dirs
    products = _occurrence_artifact(occurrence_dir, "outcomes.products", name)
    inputs = _occurrence_artifact(occurrence_dir, "inputs.components", name)
    held = pq.read_table(inputs)
    restated = pq.read_table(products).slice(0, 1).cast(held.schema)
    pq.write_table(pa.concat_tables([held, restated]), inputs)
    _shift_structure_ids(products, 1)
    with (
        pytest.raises(
            execute.PairingError,
            match=re.escape(
                f"{products} carries structure_id {count} but its dataset holds only "
                f"{count} structures"
            ),
        ),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        read.check_index()


def test_an_occurrence_artifact_without_statistics_is_refused(
    corpus_dir, occurrence_dirs
):
    # The bound is read from the footer, so a file recording none states IDs nothing
    # checked, and the count after it accepts one moved out of its dataset.
    _, occurrence_dir = occurrence_dirs
    target = occurrences.artifact_paths(occurrence_dir, "outcomes.products")[0]
    pq.write_table(pq.read_table(target), target, write_statistics=False)
    with (
        pytest.raises(
            execute.PairingError,
            match=re.escape(f"{target} records no statistics for its structure IDs")
            + r".*ord_schema\.artifacts\.scripts\.derive_occurrences --force",
        ),
        _indexed_corpus(corpus_dir, None, occurrences_dir=str(occurrence_dir)) as read,
    ):
        read.check_index()


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
    # with nothing running, after which the query runs unbounded. So the answer is
    # checked against the clock rather than trusted because no interrupt arrived.
    #
    # The lost interrupt is simulated rather than raced for: a timer that never fires
    # is what the caller sees either way, and a test that waited for the real race
    # would mostly pass while the bug was present.
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
    # reaction into an Arrow table in whatever process asked -- potentially the whole
    # corpus, from a query that reads as ordinary.
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


# When each reaction of corpus_dir first reached ord-data: one in an early commit under
# a dataset ID it no longer holds, the others in a later one.
_ARRIVALS = {
    "ord-aa01": {
        "sha": "a" * 40,
        "timestamp": datetime.datetime(2021, 3, 3, 13, 33, 15),
        "utc_offset_minutes": -300,
        "dataset_id": "ord_dataset-old",
    },
    "ord-aa02": {
        "sha": "b" * 40,
        "timestamp": datetime.datetime(2022, 12, 6, 13, 51, 13),
        "utc_offset_minutes": -300,
        "dataset_id": "ord_dataset-aa",
    },
    "ord-bb01": {
        "sha": "b" * 40,
        "timestamp": datetime.datetime(2022, 12, 6, 13, 51, 13),
        "utc_offset_minutes": -300,
        "dataset_id": "ord_dataset-bb",
    },
}


def _write_first_commits(root, output, arrivals=None) -> str:
    """Derives the first_commits artifacts of ``root``'s projections into ``output``.

    Returns:
        A glob matching them, for ``Corpus``.
    """
    rows = [
        {"reaction_id": reaction_id, **row}
        for reaction_id, row in (arrivals or _ARRIVALS).items()
    ]
    index = history.Index(
        pa.Table.from_pylist(rows, schema=history.INDEX_SCHEMA), ["a" * 40, "b" * 40]
    )
    for projected in sorted((root / "projections").glob("*.parquet")):
        output.mkdir(parents=True, exist_ok=True)
        first_commits.write_first_commits(
            projected, output / projected.name, index=index
        )
    return str(output / "*.parquet")


@pytest.fixture(scope="module")
def arrivals_pattern(corpus_dir, tmp_path_factory) -> str:
    """The first_commits artifacts of ``corpus_dir``, in a directory of their own."""
    return _write_first_commits(corpus_dir, tmp_path_factory.mktemp("first_commits"))


@pytest.fixture(scope="module")
def history_corpus(corpus_dir, arrivals_pattern) -> Iterator[execute.Corpus]:
    with _open(
        corpus_dir, resolver={}.__getitem__, first_commits_pattern=arrivals_pattern
    ) as value:
        yield value


_ADDED_IN_2022 = [
    {"op": "ge", "path": "first_commit.timestamp", "value": {"literal": "2022-01-01"}},
    {"op": "lt", "path": "first_commit.timestamp", "value": {"literal": "2023-01-01"}},
]


class TestFirstCommits:
    def test_a_filter_reads_the_first_commits(self, history_corpus, caplog):
        with caplog.at_level(logging.INFO, logger=execute.logger.name):
            found = _search(history_corpus, {"op": "and", "clauses": _ADDED_IN_2022})
        assert found == {"ord-aa02", "ord-bb01"}
        assert "this query reads the first commits" in caplog.text

    def test_a_negated_filter_answers_the_rest(self, history_corpus):
        where = {"op": "not", "clause": {"op": "and", "clauses": _ADDED_IN_2022}}
        assert _search(history_corpus, where) == {"ord-aa01"}

    def test_a_filter_combines_with_one_on_the_reaction(self, history_corpus):
        old = {
            "op": "eq",
            "path": "first_commit.dataset_id",
            "value": {"literal": "ord_dataset-old"},
        }
        ethanol = _exists({"op": "eq", "path": "smiles", "value": {"literal": "CCO"}})
        assert _search(history_corpus, {"op": "or", "clauses": [old, ethanol]}) == {
            "ord-aa01",
            "ord-bb01",
        }
        assert (
            _search(history_corpus, {"op": "and", "clauses": [old, ethanol]}) == set()
        )

    def test_a_query_that_does_not_name_them_reads_the_projection_alone(
        self, history_corpus, caplog
    ):
        with caplog.at_level(logging.INFO, logger=execute.logger.name):
            _search(
                history_corpus,
                {"op": "eq", "path": "reaction_id", "value": {"literal": "ord-aa01"}},
            )
        assert "this query reads the projection alone" in caplog.text

    def test_groups_by_where_a_reaction_first_appeared(self, history_corpus):
        answer = history_corpus.search(
            query.Query.model_validate(
                {
                    "aggregate": {
                        "group_by": ["first_commit.dataset_id"],
                        "measures": [
                            {
                                "fn": "count_distinct",
                                "path": "reaction_id",
                                "name": "reactions",
                            }
                        ],
                    }
                }
            )
        )
        assert dict(
            zip(
                answer.column("first_commit.dataset_id").to_pylist(),
                answer.column("reactions").to_pylist(),
                strict=True,
            )
        ) == {"ord_dataset-old": 1, "ord_dataset-aa": 1, "ord_dataset-bb": 1}

    def test_measures_the_earliest_arrival(self, history_corpus):
        answer = history_corpus.search(
            query.Query.model_validate(
                {
                    "aggregate": {
                        "measures": [
                            {
                                "fn": "min",
                                "path": "first_commit.timestamp",
                                "name": "earliest",
                            }
                        ]
                    }
                }
            )
        )
        assert answer.column("earliest").to_pylist() == [datetime.date(2021, 3, 3)]

    def test_orders_by_arrival(self, history_corpus):
        answer = history_corpus.search(
            query.Query.model_validate(
                {"order_by": [{"key": "first_commit.timestamp"}], "limit": 1}
            )
        )
        assert answer.column("reaction_id").to_pylist() == ["ord-aa01"]

    def test_a_corpus_without_them_refuses_a_query_naming_them(self, corpus):
        with pytest.raises(query.QueryError, match="holds no first commits"):
            _search(corpus, {"op": "and", "clauses": _ADDED_IN_2022})

    def test_their_stamps_reach_the_fingerprint(self, history_corpus, corpus):
        assert history_corpus.fingerprint != corpus.fingerprint

    def test_the_history_they_read_reaches_the_fingerprint(self, corpus_dir, tmp_path):
        # Rebuilt from a rewritten history, the artifacts can answer differently while
        # every stamp the other artifacts carry stays the same.
        fingerprints = []
        for head in ("c" * 40, "d" * 40):
            output = tmp_path / head[0]
            rows = [
                {"reaction_id": reaction_id, **row}
                for reaction_id, row in _ARRIVALS.items()
            ]
            index = history.Index(
                pa.Table.from_pylist(rows, schema=history.INDEX_SCHEMA), [head]
            )
            output.mkdir()
            for projected in sorted((corpus_dir / "projections").glob("*.parquet")):
                first_commits.write_first_commits(
                    projected, output / projected.name, index=index
                )
            with _open(
                corpus_dir,
                resolver={}.__getitem__,
                first_commits_pattern=str(output / "*.parquet"),
            ) as value:
                fingerprints.append(value.fingerprint)
        assert fingerprints[0] != fingerprints[1]


class TestFirstCommitsPairing:
    def test_refuses_a_dataset_without_them(self, corpus_dir, tmp_path):
        pattern = _write_first_commits(corpus_dir, tmp_path / "first_commits")
        (tmp_path / "first_commits" / "ord_dataset-bb.parquet").unlink()
        with pytest.raises(
            execute.PairingError, match="projections and first_commits artifacts do not"
        ):
            _open(corpus_dir, resolver={}.__getitem__, first_commits_pattern=pattern)

    def test_refuses_them_from_another_walk_when_current_is_required(
        self, corpus_dir, arrivals_pattern, monkeypatch
    ):
        monkeypatch.setattr(history, "VERSION", "2")
        with pytest.raises(
            execute.PairingError, match="not derived from a version 2 walk"
        ):
            _open(
                corpus_dir,
                resolver={}.__getitem__,
                first_commits_pattern=arrivals_pattern,
            )
        with _open(
            corpus_dir,
            resolver={}.__getitem__,
            first_commits_pattern=arrivals_pattern,
            require_current=False,
            warm=False,
        ) as value:
            assert _search(value, {"op": "and", "clauses": _ADDED_IN_2022})

    def test_refuses_a_pattern_matching_nothing(self, corpus_dir, tmp_path):
        with pytest.raises(
            execute.PairingError, match="no first_commits artifacts matched"
        ):
            _open(
                corpus_dir,
                resolver={}.__getitem__,
                first_commits_pattern=str(tmp_path / "*.parquet"),
            )

    def test_refuses_one_short_of_its_projection(self, corpus_dir, tmp_path):
        pattern = _write_first_commits(corpus_dir, tmp_path / "first_commits")
        short = tmp_path / "first_commits" / "ord_dataset-aa.parquet"
        table = pq.read_table(short)
        pq.write_table(table.slice(0, 1), short)
        with pytest.raises(
            execute.PairingError, match=r"holds 1 first commits but .* holds 2"
        ):
            _open(corpus_dir, resolver={}.__getitem__, first_commits_pattern=pattern)

    def test_refuses_a_doubled_reaction_before_reading_them(self, tmp_path):
        root = _doubled_tree(tmp_path / "doubled")
        pattern = _write_first_commits(
            root,
            tmp_path / "first_commits",
            arrivals={"ord-x01": _ARRIVALS["ord-aa01"]},
        )
        with (
            _open(
                root,
                resolver={}.__getitem__,
                warm=False,
                first_commits_pattern=pattern,
            ) as value,
            pytest.raises(execute.PairingError, match="stated by 2 rows"),
        ):
            _search(value, {"op": "and", "clauses": _ADDED_IN_2022})

    @pytest.mark.parametrize(
        "replacement",
        ["ord-zz01", "ord-aa01"],
        ids=["one the corpus lacks", "one stated twice"],
    )
    def test_refuses_them_naming_another_reaction(
        self, corpus_dir, tmp_path, replacement
    ):
        # As many rows as the projection and stamps still current, but one names
        # another reaction: read alone, a count would include it.
        pattern = _write_first_commits(corpus_dir, tmp_path / "first_commits")
        edited = tmp_path / "first_commits" / "ord_dataset-aa.parquet"
        table = pq.read_table(edited)
        ids = table.column("reaction_id").to_pylist()
        ids[ids.index("ord-aa02")] = replacement
        pq.write_table(
            table.set_column(0, table.schema.field(0), pa.array(ids, pa.string())),
            edited,
        )
        refused = "reaction ord-aa02 has no first commit"
        with pytest.raises(execute.PairingError, match=refused):
            _open(corpus_dir, resolver={}.__getitem__, first_commits_pattern=pattern)
        with _open(
            corpus_dir,
            resolver={}.__getitem__,
            first_commits_pattern=pattern,
            warm=False,
        ) as value:
            # Cold, the check waits for a query reading them, and only those refuse.
            assert _search(value, {"op": "not_null", "path": "reaction_id"})
            for _ in range(2):
                with pytest.raises(execute.PairingError, match=refused):
                    _search(value, {"op": "and", "clauses": _ADDED_IN_2022})

    @pytest.mark.parametrize(
        "left_out",
        [*(field.name for field in history.TYPE), None],
        ids=[*(field.name for field in history.TYPE), "the whole commit"],
    )
    def test_refuses_a_first_commit_missing_a_field(
        self, corpus_dir, tmp_path, left_out
    ):
        # The derivation writes every field, and a semi-join and a row carrying the
        # field would answer a negation over a NULL one differently.
        pattern = _write_first_commits(corpus_dir, tmp_path / "first_commits")
        edited = tmp_path / "first_commits" / "ord_dataset-aa.parquet"
        table = pq.read_table(edited)
        rows = table.to_pylist()
        for row in rows:
            if row["reaction_id"] == "ord-aa02":
                if left_out is None:
                    row["first_commit"] = None
                else:
                    row["first_commit"][left_out] = None
        pq.write_table(pa.Table.from_pylist(rows, schema=table.schema), edited)
        absent = (
            ", ".join(field.name for field in history.TYPE)
            if left_out is None
            else left_out
        )
        with pytest.raises(
            execute.PairingError,
            match=re.escape(
                f"reaction ord-aa02 has a first commit without its {absent},"
            ),
        ):
            _open(corpus_dir, resolver={}.__getitem__, first_commits_pattern=pattern)


# Questions about arrivals and reaction IDs and nothing else, which the first commits
# answer without the projection.
_ARRIVAL_QUESTIONS = {
    "added in 2022": {"where": {"op": "and", "clauses": _ADDED_IN_2022}},
    "not added in 2022": {
        "where": {"op": "not", "clause": {"op": "and", "clauses": _ADDED_IN_2022}}
    },
    "how many arrived in 2022": {
        "where": {"op": "and", "clauses": _ADDED_IN_2022},
        "aggregate": {"measures": [{"fn": "count", "name": "reactions"}]},
    },
    "how many first appeared in each dataset": {
        "aggregate": {
            "group_by": ["first_commit.dataset_id"],
            "measures": [
                {"fn": "count_distinct", "path": "reaction_id", "name": "reactions"}
            ],
        }
    },
    "the earliest and latest arrivals": {
        "aggregate": {
            "measures": [
                {"fn": "min", "path": "first_commit.timestamp", "name": "earliest"},
                {"fn": "max", "path": "first_commit.timestamp", "name": "latest"},
            ]
        }
    },
    "commits by how many reactions they brought": {
        "aggregate": {
            "group_by": ["first_commit.sha"],
            "measures": [{"fn": "count", "name": "first_commit"}],
        },
        "order_by": [
            {"key": "first_commit", "descending": True},
            {"key": "first_commit.sha"},
        ],
    },
    "the latest two arrivals": {
        "order_by": [
            {"key": "first_commit.timestamp", "descending": True},
            {"key": "reaction_id"},
        ],
        "limit": 2,
    },
    "first held elsewhere, or one reaction": {
        "where": {
            "op": "or",
            "clauses": [
                {
                    "op": "eq",
                    "path": "first_commit.dataset_id",
                    "value": {"literal": "ord_dataset-old"},
                },
                {"op": "eq", "path": "reaction_id", "value": {"literal": "ord-bb01"}},
            ],
        }
    },
    "a recorded commit": {"where": {"op": "not_null", "path": "first_commit.sha"}},
}


class TestFirstCommitsAlone:
    @pytest.mark.parametrize(
        "payload", _ARRIVAL_QUESTIONS.values(), ids=_ARRIVAL_QUESTIONS.keys()
    )
    def test_answers_as_the_reactions_they_join_to_do(
        self, history_corpus, monkeypatch, caplog, payload
    ):
        request = query.Query.model_validate(payload)
        with caplog.at_level(logging.INFO, logger=execute.logger.name):
            alone = history_corpus.search(request)
        assert "this query reads the first commits alone" in _messages(caplog)
        monkeypatch.setattr(query, "_reads_history_alone", lambda _: False)
        caplog.clear()
        with caplog.at_level(logging.INFO, logger=execute.logger.name):
            joined = history_corpus.search(request)
        assert "this query reads the first commits" in _messages(caplog)
        assert alone.num_rows > 0
        assert alone.schema.remove_metadata() == joined.schema.remove_metadata()
        if request.order_by:
            assert alone.to_pylist() == joined.to_pylist()
        else:
            assert sorted(alone.to_pylist(), key=repr) == sorted(
                joined.to_pylist(), key=repr
            )

    def test_orders_by_a_measure_named_like_them(self, history_corpus):
        # The measure's alias shadows the first_commit column, and the ordering has to
        # bind the alias: by the column, it would not bind at all beside a GROUP BY.
        answer = history_corpus.search(
            query.Query.model_validate(
                _ARRIVAL_QUESTIONS["commits by how many reactions they brought"]
            )
        )
        assert answer.to_pylist() == [
            {"first_commit.sha": "b" * 40, "first_commit": 2},
            {"first_commit.sha": "a" * 40, "first_commit": 1},
        ]

    def test_a_question_naming_the_projection_too_reads_it(
        self, history_corpus, caplog
    ):
        hot = {
            "op": "gt",
            "path": "conditions.temperature.setpoint_kelvin",
            "value": {"literal": 0},
        }
        with caplog.at_level(logging.INFO, logger=execute.logger.name):
            _search(history_corpus, {"op": "and", "clauses": [hot, *_ADDED_IN_2022]})
        assert "this query reads the first commits" in _messages(caplog)
        assert "this query reads the first commits alone" not in _messages(caplog)
