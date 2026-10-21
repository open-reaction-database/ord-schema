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

"""Running a compiled query against a corpus of projections and structures base.

The compiler (:mod:`ord_schema.search.query`) leaves two kinds of work undone on
purpose, and this module is where both happen:

* **Names become structures.** Compound names bind as resolved SMILES, through
  :mod:`ord_schema.resolvers` by default or any injected resolver. The corpus keeps
  what each name resolved to, and a result states it under ``RESOLVED``.
* **Structure predicates become bitmaps.** Each compiles to a bitmap test on the
  element's ``structure_id``, and the chemistry runs here. Substructure goes through an
  RDKit ``SubstructLibrary`` built over the corpus: it screens on the pattern
  fingerprints the artifact already stores and matches the survivors exactly, across
  its own threads with the GIL released, so a server can share one corpus between
  concurrent requests without forking. Similarity stays in SQL and has no verification
  step at all -- Tanimoto is defined on the Morgan fingerprint, so the screen is the
  whole answer.

Structure IDs are dataset-local, so a corpus-wide bitmap needs an offset per file.
``Corpus`` pairs every projection with its structures artifact, refuses a pair whose
stamps disagree about the source dataset -- the two files are two statements of one
derivation, and a mismatched pair would join IDs to the wrong molecules silently --
and publishes a ``reactions`` relation carrying each row's offset in the column the
compiled SQL expects.

A ``Corpus`` is safe to share between concurrent searches. A ``DuckDBPyConnection``
holds the pending result of its last ``execute``, so two threads stepping on one read
each other's rows rather than raising; each search therefore takes its own cursor. A
cursor in turn sees only the catalog, not the objects ``register`` attaches to the
connection that registered them, so every relation here is a catalog table or view.

Finding the reactions that hold a match is a scan of every reaction, which costs more
than the chemistry for all but the narrowest queries. An **occurrence index** -- one row
per structure occurrence, carrying the corpus-wide ID and the path it sat at -- turns
that scan into a semi-join against a narrow table. It is spent through the compiler's
``ElementIndex`` hook, one quantifier at a time: an ``exists`` whose body is one
structure predicate on the element's own ``smiles`` and nothing else becomes a condition
on the row, and everything around it compiles exactly as it would have. An aggregate
over the quantifier, a scalar beside it, an ordering, a limit, a negation, a
disjunction, or a second quantifier each compose with the condition rather than
disqualifying it.

What still compiles over the elements, and so reads the level's pivot where there is
one and the projection where there is not: a quantifier whose body binds any element
field beside the structure -- ``reaction_role`` as much as any other -- reaches a nested
one, holds no structure predicate or more than one, or is not a plain conjunction, and
every ``forall``, since the index shows the elements that match, never that every
element does. Either one reads the whole element, so "pyridine as the solvent" stays a
condition on one element rather than an intersection of two reaction sets, which
over-returns.

An occurrence and a pivot row each name their reaction by ID, so neither is read over a
corpus that states a reaction twice; two files carrying one reaction would have each
copy's elements answering for the other.

Reading the projection at all costs a re-parse of every file's footer, and over a schema
this wide that parse is most of what a query spends however few leaves it goes on to
read. The corpus keeps the parsed footers instead; see ``_cache_footers``.

The index is a view over occurrence artifacts covering every indexed path, and the
pivots are views over pivot artifacts: the corpus builds neither, and a level or a path
with no artifact is answered from the projection's own lists, which is slower and
nothing else. The substructure library is the one thing built, at open unless the corpus
was given ``warm=False``, which leaves it to the first substructure predicate.

The grammar bounds what a query can cost (one pass and a sort), and ``search``
takes a wall-clock timeout over the whole call, bounding when the caller is answered
rather than when work stops; see ``_Deadline`` for which phases can be interrupted and
which are only checked.

Given the source datasets, a corpus also serves what a search found. A search returns
reaction IDs, and ``Corpus.reactions`` reads each one's serialized ``Reaction`` out of
the source file it came from, so the search and the records are one snapshot: every
source is paired with the projection whose ``ord.source_md5`` stamp names its hash, and
sources and projections that do not pair one to one are refused at open.
"""

import array
import collections
import dataclasses
import glob
import hashlib
import json
import math
import pathlib
import re
import threading
import time
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Self

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, DataStructs
from rdkit.Chem import rdSubstructLibrary

from ord_schema import parquet, resolvers
from ord_schema.artifacts import (
    base,
    first_commits,
    history,
    occurrences,
    pivot,
    projection,
    structures,
)
from ord_schema.logging import get_logger
from ord_schema.search import query

logger = get_logger(__name__)

# Structures read per batch while building the library. Bounds peak memory during the
# build without making the per-batch Python overhead matter.
_BUILD_BATCH = 50_000

# Match sets kept for reuse. Each is a bitmap with one bit per corpus structure ID, so
# the cache is small beside the chemistry whose repeat it saves.
_CACHED_MATCHES = 16

# Compound names kept resolved for reuse. An entry is a name and its SMILES, and a miss
# is a call to an external service, so the bound is generous.
_CACHED_RESOLUTIONS = 1024

# The modules that write each derived file, which a refusal names as the remedy.
_DERIVE_PROJECTION = "ord_schema.artifacts.scripts.derive_projection"
_DERIVE_PIVOTS = "ord_schema.artifacts.scripts.derive_pivots"
_DERIVE_OCCURRENCES = "ord_schema.artifacts.scripts.derive_occurrences"

# Most rows a search returns where the caller states no bound. The grammar leaves
# ``limit`` optional and a question like "reactions above 350 K" matches a large share
# of the corpus, so an unbounded default builds an Arrow table of them in whatever
# process asked -- from a query that reads as ordinary. A thousand is a page of
# reactions to look at rather than a dataset to compute over; a caller that wants the
# whole match set says so with max_rows=None, and one that wants a different page says
# so with a number.
DEFAULT_MAX_ROWS = 1000

# Stamped on a result the corpus's own bound may have cut short. On the table rather
# than logged alone, because whoever reads the answer is often not whoever reads the
# log: a summary saying "1000 reactions" is wrong in a way nothing in the rows shows.
TRUNCATED = "ord.truncated"

# Stamped on a result whose query named a compound: a JSON object mapping each name to
# the SMILES it was searched as. The query keeps the name, and the structure is what a
# reader needs to judge the answer; read it with ``resolutions``.
RESOLVED = "ord.resolved"

# What Corpus.reactions returns. ``reaction`` is the serialized Reaction as its source
# stores it, so serving one parses nothing.
_RECORD_SCHEMA = pa.schema(
    [
        pa.field("reaction_id", pa.string()),
        pa.field(query.DATASET_ID, pa.string()),
        pa.field("reaction", pa.binary()),
    ]
)

# What Corpus.datasets returns: the dataset scalars a source's footer carries, and the
# row count it states.
_DATASETS_SCHEMA = pa.schema(
    [
        pa.field(query.DATASET_ID, pa.string()),
        pa.field("name", pa.string()),
        pa.field("description", pa.string()),
        pa.field("reactions", pa.int64()),
    ]
)


# What a structure predicate's answer depends on: the operation, whether the string
# is read as SMARTS or as SMILES, the string, the similarity threshold, and whether a
# substructure match honors stereocenters.
_MatchKey = tuple[str, bool, str, float | None, bool]


def _reads_as_smarts(parameter: query.StructureParameter) -> bool:
    """Returns whether the predicate's string is read as SMARTS rather than SMILES.

    Which parser reads a string is part of what the string means -- ``C1=CC=CC=C1`` is
    benzene as SMILES and six aliphatic carbons as SMARTS -- so the query molecule and
    the key a match set is cached under have to make this call the same way.

    Args:
        parameter: The structure predicate.

    Returns:
        True for a substructure predicate stating its own pattern, which is the only
        form a query writes in SMARTS; a resolved compound is always SMILES.
    """
    return parameter.pattern is not None and parameter.op == "substructure"


# The predicates the artifact answers with one equality: the column each compares, and
# the function that puts a query molecule into that column's terms.
_HASH_MATCHES: dict[str, tuple[str, Callable[[Chem.Mol], str | None]]] = {
    "same_compound": ("mol_hash", structures.mol_hash),
    "same_parent": ("parent_hash", structures.parent_hash),
}


class PairingError(ValueError):
    """The artifacts, or the sources beside them, do not form a consistent corpus."""


def _resolve_with_resolvers(name: str) -> str:
    """Resolves a compound name to SMILES through ord_schema.resolvers."""
    smiles, _ = resolvers.resolve_name("name", name)
    return smiles


def _canonical(smiles: str) -> str | None:
    """Returns the SMILES as the corpus spells it, or None if RDKit cannot parse it.

    A resolver answers in whatever form its service uses -- PubChem returns pyridine
    as ``C1=CC=NC=C1`` -- and a structure predicate's match set is cached under the
    resolved string, so two names for one molecule share an entry only once both are
    spelled one way. That spelling is the projection's, which keeps the stereo groups
    plain SMILES drops: without them a racemate's query molecule is one enantiomer, and
    shares that enantiomer's cache entry.

    Args:
        smiles: What the resolver returned.

    Returns:
        The canonical form, or None if RDKit cannot parse it.
    """
    try:
        return resolvers.canonicalize_smiles(smiles)
    except ValueError:
        return None


def resolutions(table: pa.Table) -> dict[str, str]:
    """Returns each compound name a search resolved, mapped to the SMILES it matched.

    Args:
        table: What ``Corpus.search`` returned.

    Returns:
        Each name mapped to the canonical SMILES the search matched it as, read from
        the ``RESOLVED`` stamp; empty for a query that named no compound.
    """
    stamp = (table.schema.metadata or {}).get(RESOLVED.encode())
    return json.loads(stamp) if stamp is not None else {}


def sql_string(value: str) -> str:
    """Returns ``value`` as a SQL string literal, single quotes escaped."""
    escaped = value.replace("'", "''")
    return f"'{escaped}'"


def _sql_paths(paths: Iterable[str]) -> str:
    """Returns ``paths`` as a SQL list literal naming each file and no other.

    read_parquet reads its arguments as globs, so a path is not simply itself: a
    directory named with a bracket becomes a character class that matches a sibling
    instead, and one with a star or a question mark matches its neighbors as well as
    itself. Either way the query answers off files nobody asked for, which is the kind
    of wrong that arrives as data rather than as an error.

    Args:
        paths: Paths to files, as this process found them.

    Returns:
        A SQL list literal, single quotes escaped and glob characters spelled as the
        one-character classes that match them literally.
    """
    return "[" + ", ".join(sql_string(glob.escape(path)) for path in paths) + "]"


def _max_structure_id(path: str, derive: str) -> int | None:
    """Returns the largest ``structure_id`` a projection or an artifact of one carries.

    Read from Parquet footer statistics, so it decodes no column data however large the
    file is. Every column named ``structure_id`` counts, however deeply it is nested.

    The columns are found by the paths their chunks record rather than through
    ``metadata.schema``. pyarrow's ``ParquetSchema`` and the ``FileMetaData`` it belongs
    to refer to each other, so a footer whose schema was read stays parsed after the
    file closes, until the cyclic collector runs. Reading every file of a corpus that
    way holds many footers at once.

    Args:
        path: Path to a projection, or to a pivot or occurrence artifact.
        derive: The module that writes ``path``, which a refusal names.

    Returns:
        The maximum across every ``structure_id`` column, or None when the file holds no
        compound with a structure.

    Raises:
        PairingError: If a row group's footer records no statistics for a column of
            structure IDs, which leaves them unbounded without reading every row.
    """
    with pq.ParquetFile(path) as parquet_file:
        metadata = parquet_file.metadata
        columns: list[int] = []
        largest = None
        for group in range(metadata.num_row_groups):
            row_group = metadata.row_group(group)
            if group == 0:
                # Every row group holds the same leaves in the same order. A path is
                # dotted, and no field name in an artifact schema holds a dot.
                columns = [
                    index
                    for index in range(row_group.num_columns)
                    if row_group.column(index).path_in_schema.rsplit(".", 1)[-1]
                    == "structure_id"
                ]
            for index in columns:
                chunk = row_group.column(index)
                statistics = chunk.statistics
                if statistics is not None and statistics.has_min_max:
                    value = statistics.max
                    largest = value if largest is None else max(largest, value)
                elif statistics is None and chunk.num_values:
                    # Statistics without a bound are a column holding no ID; no
                    # statistics at all say nothing, and this library always writes
                    # them, so the file was rewritten by something else.
                    raise PairingError(
                        f"{path} records no statistics for its structure IDs, so they "
                        "cannot be checked against its dataset's structures without "
                        f"reading every row; re-derive it with {derive} --force, since "
                        "a run without it skips a file whose stamps are current"
                    )
    return largest


def _article(word: str) -> str:
    """Returns the indefinite article ``word`` takes, or "a" for an empty name.

    ``load_stamps`` requires the artifact key to be present rather than non-empty, so
    an empty name reaches here from a footer written by hand or truncated mid-write.
    Tested with ``startswith`` rather than by indexing, which raises on it, or against
    a string of vowels, which contains the empty one and would answer "an".

    Args:
        word: The name to precede.

    Returns:
        "an" before a vowel, "a" otherwise.
    """
    return "an" if word.startswith(("a", "e", "i", "o", "u")) else "a"


def _refuse_unreadable(
    path: str,
    stamps: base.Stamps,
    *,
    artifact: str,
    schema: pa.Schema,
    require_current: bool,
    held: Callable[[str], str | None] | None = None,
    expected: str | None = None,
    reading: str = "",
) -> None:
    """Refuses a file a corpus cannot read as the artifact it is filed as.

    The one definition of what makes an artifact readable, for every artifact a corpus
    opens. One definition because the checks have to agree and nothing else makes them:
    a copy per artifact still runs when it falls behind the others, so the artifact
    whose copy is short a check is admitted and fails somewhere far from here.

    Order matters. What the file *is* comes first, since the later checks read it as
    the artifact it claims to be, and a level or path stamped on something else says
    nothing.

    Args:
        path: The file to check, named in every message.
        stamps: Its stamps, already read.
        artifact: Artifact name it is expected to hold.
        schema: The columns this library's version of that artifact carries.
        require_current: Refuse an artifact not written by the current versions. A
            column is refused whatever this says, because a reader binding one that is
            absent fails however old the file is; stamps alone are the caller's policy.
        held: Reads the level or path a file is stamped with, for an artifact split
            across a tree. None where the artifact is one file per dataset and where it
            sits says nothing.
        expected: The level or path ``path`` sits under, which ``held`` must match.
        reading: How a message names the relation to ``expected`` -- "over" for a pivot
            level, "at" for an indexed path.

    Raises:
        PairingError: If the file holds another artifact, is stamped with another level
            or path, lacks a column this library reads, or is stale under
            ``require_current``.
    """
    if stamps.artifact != artifact:
        raise PairingError(
            f"{path} is {_article(stamps.artifact)} {stamps.artifact}, not "
            f"{_article(artifact)} {artifact}"
        )
    at = f" {reading} {expected}" if expected is not None else ""
    if held is not None:
        stamped = held(path)
        if stamped != expected:
            raise PairingError(
                f"{path} holds the {artifact} {reading} {stamped}, but sits where "
                f"{expected} is read from, so a quantifier would be answered by the "
                "wrong elements"
            )
    missing = base.missing_columns(path, schema)
    if missing:
        # The stamps say nothing about a file's columns, so this is the only check that
        # sees a schema change. Without it a file predating a column is read as a corpus
        # member and fails deep in DuckDB -- as a binder error on the query that wants
        # the column, or, where some files have it and some do not, as a schema mismatch
        # that takes down every structure query.
        raise PairingError(
            f"{path} is {_article(artifact)} {artifact} artifact{at} without "
            f"{missing}, which this library reads; derive it again first"
        )
    if require_current and not base.stamps_are_current(stamps, artifact):
        raise PairingError(
            f"{path} is a stale {artifact} artifact; derive it again first"
        )


def _index(
    pattern: str, artifact: str, schema: pa.Schema, require_current: bool
) -> dict[str, tuple[str, base.Stamps]]:
    """Returns the artifacts matching ``pattern``, keyed by source dataset.

    Keyed by the source hash rather than the filename because that is what an
    artifact *is*: two files pair when they restate the same dataset, whatever
    they are called or wherever they sit. Keying on the basename would let two
    files in different directories collapse onto one another silently.

    Args:
        pattern: Glob matching the artifact files.
        artifact: Artifact name every match must hold.
        schema: The columns this library's version of the artifact carries.
        require_current: Refuse artifacts not written by the current versions.

    Returns:
        A mapping from source dataset hash to the path and the stamps read off it. The
        stamps come back rather than being dropped because reading them opens a footer,
        and the corpus fingerprint wants exactly the same ones.

    Raises:
        PairingError: If a match holds another artifact, lacks a column this library
            reads, is stale under ``require_current``, or restates a dataset another
            match already did.
    """
    index: dict[str, tuple[str, base.Stamps]] = {}
    for path in sorted(glob.glob(pattern, recursive=True)):
        stamps = base.load_stamps(path)
        _refuse_unreadable(
            path,
            stamps,
            artifact=artifact,
            schema=schema,
            require_current=require_current,
        )
        if stamps.source_md5 in index:
            raise PairingError(
                f"{path} and {index[stamps.source_md5][0]} are both {artifact} "
                "artifacts of the same source dataset; which one answers a query "
                "would be arbitrary"
            )
        index[stamps.source_md5] = (path, stamps)
    return index


def _fingerprint(records: Iterable[tuple[Any, ...]]) -> str:
    """Returns a short digest naming the artifacts a corpus opened.

    Args:
        records: One per artifact in the corpus, in any order: its stamps as a tuple,
            and for a first_commits artifact the walk of ord-data's history it read as
            well. Each artifact paired with a projection belongs here: the structures
            artifact answers every structure predicate and the first_commits artifact
            every question about when a reaction arrived, and each is derived
            separately, so a corpus whose partners were rebuilt can answer differently
            with its projections untouched.

    Returns:
        Sixteen hex characters, sorted so the same corpus opened through a different
        glob digests the same. The whole stamp goes in rather than the source hash
        alone, because two artifacts built from identical sources by different library
        versions can disagree too.
    """
    digest = hashlib.sha256()
    # Sorted as text rather than as tuples: a stamp field is optional, and comparing
    # None against a string to order two artifacts raises.
    for value in sorted(repr(record) for record in records):
        digest.update(value.encode())
    return digest.hexdigest()[:16]


def _pair(
    projection_pattern: str,
    structures_pattern: str,
    require_current: bool,
    first_commits_pattern: str | None = None,
) -> tuple[list[tuple[str, str, str, str | None]], str, dict[str, str]]:
    """Returns (projection, structures, source, dataset ID) tuples, verified by stamps.

    Args:
        projection_pattern: Glob matching the projection base.
        structures_pattern: Glob matching the structures base.
        require_current: Refuse artifacts not written by the current versions.
        first_commits_pattern: Glob matching the first_commits artifacts, or None for
            a corpus that holds none.

    Returns:
        One tuple per source dataset, ordered by the source hash, carrying that hash
        beside the paths -- every artifact derived from the dataset names it, so it is
        what pairs a pivot with the offset this corpus gave its projection -- and the
        dataset ID the projection's stamps record, or None where its source recorded
        none. Also the fingerprint over the stamps this already read to verify them,
        and each source's first_commits artifact, empty where the corpus holds none.

    Raises:
        PairingError: If no projection matches, if any side lacks a column this library
            reads, or if any side states a dataset another does not, since a
            projection's IDs index its own partner's molecules and nobody else's, and a
            first commit asked of a dataset missing one would read as never having
            arrived. Given ``first_commits_pattern``, also if it matches nothing, or
            under ``require_current`` if an artifact was derived from another walk of
            ord-data's history.
    """
    projections = _index(
        projection_pattern, projection.ARTIFACT, projection.SCHEMA, require_current
    )
    structure_files = _index(
        structures_pattern, structures.ARTIFACT, structures.SCHEMA, require_current
    )
    if not projections:
        raise PairingError(f"no projections matched: {projection_pattern}")
    unpaired = projections.keys() ^ structure_files.keys()
    if unpaired:
        orphans = sorted((projections | structure_files)[key][0] for key in unpaired)
        raise PairingError(
            "projections and structures artifacts do not pair up; these have no "
            f"counterpart derived from the same source dataset: {orphans}"
        )
    arrivals: dict[str, tuple[str, base.Stamps]] = {}
    if first_commits_pattern is not None:
        arrivals = _index(
            first_commits_pattern,
            first_commits.ARTIFACT,
            first_commits.SCHEMA,
            require_current,
        )
        if not arrivals:
            raise PairingError(
                f"no first_commits artifacts matched: {first_commits_pattern}"
            )
        unpaired = projections.keys() ^ arrivals.keys()
        if unpaired:
            # Every dataset or none: a semi-join on a dataset with no first commits
            # answers a negated question as though its reactions matched.
            orphans = sorted((projections | arrivals)[key][0] for key in unpaired)
            raise PairingError(
                "projections and first_commits artifacts do not pair up; these have no "
                f"counterpart derived from the same source dataset: {orphans}"
            )
    # The walk of ord-data's history each first_commits artifact read, which its other
    # stamps do not record: rebuilt from another, it answers differently under them.
    walks = {key: history.read_stamp(path) for key, (path, _) in arrivals.items()}
    if require_current:
        stale = sorted(
            arrivals[key][0]
            for key, walk in walks.items()
            if walk is None or walk[0] != history.VERSION
        )
        if stale:
            raise PairingError(
                "these first_commits artifacts were not derived from a version "
                f"{history.VERSION} walk of ord-data's history; derive them again "
                f"first: {stale}"
            )
    pairs = [
        (
            projections[key][0],
            structure_files[key][0],
            key,
            projections[key][1].source_dataset_id,
        )
        for key in sorted(projections)
    ]
    fingerprint = _fingerprint(
        [
            *(
                dataclasses.astuple(stamps)
                for index in (projections, structure_files)
                for _, stamps in index.values()
            ),
            *(
                (*dataclasses.astuple(stamps), walks[key])
                for key, (_, stamps) in arrivals.items()
            ),
        ]
    )
    return pairs, fingerprint, {key: path for key, (path, _) in arrivals.items()}


@dataclasses.dataclass(frozen=True)
class _Source:
    """A source dataset, paired with the projection that restates it.

    Attributes:
        path: The source file.
        dataset_id: The dataset ID the projection's stamps record, which is the value
            a search's ``dataset_id`` column carries; None where the source recorded
            none.
        name: The dataset name, from the source's footer.
        description: The dataset description, from the source's footer.
        reactions: Rows the source's footer states.
    """

    path: str
    dataset_id: str | None
    name: str
    description: str
    reactions: int


def _pair_sources(
    pattern: str, pairs: list[tuple[str, str, str, str | None]]
) -> list[_Source]:
    """Returns the sources matching ``pattern``, each paired with its projection.

    Paired by content hash, the one every artifact derived from a source records as
    ``ord.source_md5``, for the reason ``_index`` keys on it: a file is the dataset it
    holds, whatever it is called or wherever it sits. Hashing reads every reaction, so
    this is most of what opening a corpus with sources costs.

    Args:
        pattern: Glob matching the source dataset files.
        pairs: The corpus's pairs, from ``_pair``.

    Returns:
        One record per source, in the order of ``pairs``.

    Raises:
        PairingError: If no file matches, a match cannot be read as a source dataset,
            two matches hold the same content, a projection's source is not among the
            matches, or a match restates a dataset no projection does.
    """
    start = time.perf_counter()
    views: dict[str, parquet.DatasetView] = {}
    for path in sorted(glob.glob(pattern, recursive=True)):
        # The footer is read on opening and every reaction on hashing, so a file that is
        # not a source and one whose reactions cannot be read both surface here.
        try:
            view = parquet.DatasetView(path)
            digest = view.md5()
        except (OSError, ValueError) as error:
            raise PairingError(
                f"{path} is not a readable source dataset: {error}"
            ) from error
        if digest in views:
            raise PairingError(
                f"{path} and {views[digest].path} are the same source dataset; which "
                "one answers a fetch would be arbitrary"
            )
        views[digest] = view
    if not views:
        raise PairingError(f"no sources matched: {pattern}")
    logger.info(
        "hashed %d sources in %.1fs to pair them with their projections",
        len(views),
        time.perf_counter() - start,
    )
    projected = {source: path for path, _, source, _ in pairs}
    missing = sorted(projected[key] for key in projected.keys() - views.keys())
    if missing:
        raise PairingError(
            f"no file matching {pattern} holds the source these projections restate, "
            f"so the reactions a search finds in them could not be fetched: {missing}"
        )
    orphans = sorted(views[key].path for key in views.keys() - projected.keys())
    if orphans:
        raise PairingError(
            "no projection restates these sources, so a fetch could return a reaction "
            f"no search can find: {orphans}"
        )
    return [
        _Source(
            path=views[source].path,
            dataset_id=dataset_id,
            name=views[source].name,
            description=views[source].description,
            reactions=len(views[source].reactions),
        )
        for _, _, source, dataset_id in pairs
    ]


# A structure that did not parse still gets a library entry, so every structure ID
# resolves to one and the mapping needs no special case. An empty molecule fingerprints
# to no bits, and a query with no atoms is refused, so it can never match.
_UNPARSEABLE = Chem.Mol().ToBinary()
_NO_BITS = DataStructs.ExplicitBitVect(structures.PATTERN_FP_SIZE)


# The paths the occurrence index covers, in the order the artifact lists them. They are
# ``occurrences.PATHS`` rather than a walk of this module's own, so the compiler routes
# a quantifier to the index exactly where an artifact carries the rows to answer it: two
# walks of the schema would be two answers waiting to disagree, and the disagreement is
# silent, since a path routed here that no artifact holds is a quantifier answered from
# nothing.
INDEXED_PATHS = tuple(occurrences.PATHS)


def _warn_when_the_bound_was_reached(
    request: query.Query, compiled: query.Compiled, returned: int
) -> bool:
    """Says when this corpus's bound, not the query's, may have cut the answer short.

    Said after the query rather than while compiling it, because a bound applied is not
    an answer cut short: a query the bound never reached is the ordinary case, and
    warning about it on every generous limit is how a reader learns to skip the line
    that matters. A result that comes back exactly full is the most a row count can say
    -- it may hold every match and end there -- and nothing in the table says which.

    An aggregated query is the louder case. A truncated list of reactions is a sample of
    the matching ones; a truncated list of groups is part of a distribution read as the
    whole of it, and where the query ordered by nothing it is an arbitrary part.

    Args:
        request: The query as asked, whose own ``limit`` is not this corpus's doing.
        compiled: What it compiled to, carrying the limit that reached the SQL.
        returned: Rows the query returned.

    Returns:
        Whether the answer may be short of the matches, for stamping on the table.
    """
    if compiled.limit is None or compiled.limit == request.limit:
        return False
    if returned < compiled.limit:
        return False
    if request.aggregate is not None:
        logger.warning(
            "this query grouped into at least the %d rows the corpus bounds it at, so "
            "the groups it returned are %s of the distribution rather than all of it; "
            "ask for fewer groups, or raise max_rows",
            compiled.limit,
            "an arbitrary part" if not request.order_by else "the leading part",
        )
    else:
        logger.warning(
            "this query returned the %d rows the corpus bounds it at, so there may be "
            "matches it did not return; raise max_rows to see them",
            compiled.limit,
        )
    return True


# How an index condition names the index, which is also how the log recognizes one.
_INDEX_RELATION = "occurrences AS occurrence"


def _index_condition(path: str, allocate: Callable[[], str]) -> str | None:
    """Returns a row condition standing for an element quantifier, or None.

    The compiler's ``ElementIndex``: it offers only a quantifier whose body is a lone
    structure predicate on the element's own ``smiles``, since a row holds the path
    and the corpus-wide structure ID and no element field. This declines a path the
    index does not cover, which is left to the elements.

    What comes back is a semi-join rather than a whole query, which is what lets the
    rest of the query stay exactly what it was: the aggregate, the ordering, the
    limit, and every clause beside this one compile unchanged, and this condition
    narrows the same rows they are about. A reaction is one row of ``reactions``, so
    the join cannot multiply rows and needs no DISTINCT.

    Args:
        path: The path the quantifier ranges over.
        allocate: Names the structure parameter the match set binds under.

    Returns:
        The condition, or None to leave the quantifier compiled over the elements.
    """
    if path not in INDEXED_PATHS:
        return None
    conditions = [
        f"occurrence.path = {sql_string(path)}",
        f"get_bit(CAST(${allocate()} AS BITSTRING), occurrence.global_id::INTEGER) = 1",
    ]
    # S608: every fragment is a path this module indexes or a compiler-issued
    # parameter name. The inner relation is aliased so that the outer reaction_id and
    # the inner one cannot be confused for each other.
    return (
        f"reaction_id IN (SELECT occurrence.reaction_id FROM {_INDEX_RELATION} "  # noqa: S608
        f"WHERE {' AND '.join(conditions)})"
    )


def _relations_read(sql: str, looked_up: dict[str, str]) -> str:
    """Returns which relations a compiled query reads, for a person reading the log.

    Read off the SQL rather than recorded as the compiler asks for views: compiling a
    nested quantifier looks up the inner level's pivot before the outer level can
    decline, and the whole quantifier is then answered from the projection's lists.

    Args:
        sql: The compiled query.
        looked_up: The view each pivot lookup returned, keyed by level.

    Returns:
        The occurrence index, the pivots over the levels the SQL names, the first
        commits, or the projection alone.
    """
    reads = []
    if _INDEX_RELATION in sql:
        reads.append("the occurrence index")
    if _reads_history(sql):
        reads.append("the first commits")
    # Matched as whole words, since one level's view name can begin another's:
    # pivot_outcomes_products is a prefix of pivot_outcomes_products_measurements.
    pivoted = sorted(
        path
        for path, name in looked_up.items()
        if re.search(rf"\b{re.escape(name)}\b", sql)
    )
    if pivoted:
        reads.append(f"the pivots over {pivoted}")
    return " and ".join(reads) or "the projection alone"


def _reads_history(sql: str) -> bool:
    """Returns whether compiled SQL reads the first commits, matched as a whole word."""
    return re.search(rf"\b{query.HISTORY_TABLE}\b", sql) is not None


def _group(entry_of: array.array, count: int) -> tuple[array.array, array.array]:
    """Inverts structure-to-entry into entry-to-structures.

    Args:
        entry_of: The library entry each structure ID resolves to, indexed by ID.
        count: How many entries the library holds.

    Returns:
        ``(members, starts)``, where the structure IDs sharing entry ``i`` are
        ``members[starts[i]:starts[i + 1]]``. Counting sort, so one pass each way and
        no per-entry list object.
    """
    starts = array.array("I", bytes(4 * (count + 1)))
    for entry in entry_of:
        starts[entry + 1] += 1
    for index in range(1, count + 1):
        starts[index] += starts[index - 1]
    members = array.array("I", bytes(4 * len(entry_of)))
    cursor = starts[:-1]
    for global_id, entry in enumerate(entry_of):
        members[cursor[entry]] = global_id
        cursor[entry] += 1
    return members, starts


def _cache_footers(connection: duckdb.DuckDBPyConnection) -> None:
    """Keeps the Parquet footers parsed between queries.

    A scan re-reads and re-parses the footer of every file it touches, and the
    projection's footers are large: every leaf of the schema in every file, with
    statistics per row group. Over ORD that parse is most of what a query costs, and it
    is charged again on every query however little column data the query goes on to
    read; the search README has the measurements.

    The footers are bounded by the files rather than by what is asked of them. The
    corpus parses every projection footer at open, counting the reactions view to check
    its pairing. DuckDB's own accounting books only part of what they hold, so most of
    it sits outside ``memory_limit``. DuckDB re-reads a file that has been rewritten, so
    an artifact replaced under an open corpus is read as it stands.

    Set globally rather than on this connection: a search runs on its own cursor, which
    is its own session, and a session-scoped setting would leave every search paying the
    parse.

    Args:
        connection: The corpus connection, on the database every cursor shares.
    """
    connection.execute("SET GLOBAL parquet_metadata_cache=true")


def _free_result_order(connection: duckdb.DuckDBPyConnection) -> None:
    """Lets rows come back in whatever order the scan's threads finish them.

    A query without ``ORDER BY`` promises no order, and keeping the scan's order anyway
    holds every row back until the row groups before it are done: over the projection's
    row groups, a wide scan returning many rows waits longer on that than on the scan,
    where a count of the same rows does not wait at all. ``ORDER BY`` is honored, and
    ``reactions`` orders its own rows.

    Set globally, as ``_cache_footers`` is, so that each search's cursor is covered.

    Args:
        connection: The corpus connection, on the database every cursor shares.
    """
    connection.execute("SET GLOBAL preserve_insertion_order=false")


def _confine(connection: duckdb.DuckDBPyConnection, directories: Iterable[str]) -> None:
    """Leaves the database able to reach the files under ``directories`` and no others.

    The compiler is what keeps a query to one relation, and this is what holds if the
    compiler is ever wrong: a statement that reached DuckDB with text it should not
    carry would otherwise run with this process's filesystem and network. Confined, it
    cannot read a file outside the corpus's own trees -- ``/proc/self/environ`` among
    them -- fetch a URL, attach a database, or load an extension, and the configuration
    is locked, so no statement can lift any of that. Spilling to DuckDB's temporary
    directory works under it.

    Disabling external access alone would refuse the lazy Parquet views too; a list of
    allowed directories beside it is the one form that leaves them readable. Writes
    under those directories are allowed, since DuckDB's list has no read-only form, and
    mounting the trees read-only is what closes that.

    Each directory is listed as the reads spell it, which DuckDB resolves the way it
    resolves the reads themselves, relative paths and ``..`` included. One whose name
    holds a glob metacharacter is listed twice. ``_sql_paths`` escapes it for
    ``read_parquet``, and DuckDB checks that pattern's text against the list as well as
    every file the pattern expands to, so either form alone refuses the read. The
    escaped form opens nothing else: as a pattern it matches only the directory it
    escapes.

    Set on the database rather than the connection, so every search's cursor is held to
    it, and after ``_cache_footers``, which the lock would otherwise refuse.

    Args:
        connection: The corpus connection, on the database every cursor shares.
        directories: The trees the corpus reads: the projections', the structures', the
            sources' where it was given them, and any pivot or occurrence directory.
    """
    spelled = set(directories)
    allowed = sorted(spelled | {glob.escape(directory) for directory in spelled})
    connection.execute("SET allowed_directories = $allowed", {"allowed": allowed})
    connection.execute("SET enable_external_access = false")
    connection.execute("SET lock_configuration = true")


class _Deadline:
    """When a search has to be finished, and what it may still spend.

    A search is several phases -- compiling, which publishes the views it routes to on
    their first use; each structure predicate, its compound resolved and its matches
    screened and verified; and the query itself -- and only the last is a DuckDB
    statement a timer can interrupt. The rest is RDKit with the GIL released, an
    external resolver, or a build another search waits on. So this bounds when the
    caller gets an answer or an error rather than when work stops: each phase is
    checked as it finishes, and the query gets whatever is left.

    A shared build is deliberately not interrupted. The substructure library is built
    once for the corpus under a lock, and the occurrence index checked once, so a timer
    that killed either would fail every search queued behind it -- callers that asked
    for no bound included.
    """

    def __init__(self, seconds: float | None) -> None:
        """Starts the clock.

        Args:
            seconds: Wall-clock bound, or None for no deadline.
        """
        self._at = None if seconds is None else time.monotonic() + seconds
        self._seconds = seconds

    def remaining(self) -> float | None:
        """Returns the seconds left, or None where there is no deadline."""
        if self._at is None:
            return None
        return self._at - time.monotonic()

    def check(self, phase: str) -> None:
        """Raises if the deadline has passed.

        Args:
            phase: What has just finished, named in the error so the caller learns
                which part of the search was the slow one.

        Raises:
            TimeoutError: If the deadline has passed.
        """
        left = self.remaining()
        if left is not None and left <= 0:
            raise TimeoutError(
                f"the search exceeded {self._seconds} seconds while {phase}; it ran to "
                "completion because that phase cannot be interrupted"
            )


def _run_with_timeout(
    cursor: duckdb.DuckDBPyConnection,
    sql: str,
    parameters: dict[str, Any],
    timeout_seconds: float,
) -> pa.Table:
    """Runs ``sql``, interrupting it if it outlasts ``timeout_seconds``.

    The one phase of a search that is interruptible, and the only one whose work is this
    caller's alone: the cursor is the search's own, so interrupting it fails nothing
    else.

    ``Timer.cancel`` only sets a flag, so a timer that has already passed its own
    check fires anyway. The lock makes the interrupt and the teardown exclusive, so
    it either lands while this query is running or not at all -- and in particular
    never after the caller closes the cursor.

    Args:
        cursor: The cursor to run on, and the one an expired timer interrupts.
        sql: The compiled query.
        parameters: Values to bind.
        timeout_seconds: Wall-clock bound.

    Returns:
        The result as an Arrow table.

    Raises:
        TimeoutError: If the query is interrupted by the timer, if ``timeout_seconds``
            has already run out before it starts, or if the query answers after the
            bound because its interrupt was lost.
    """
    # Refused rather than armed. A timer set for a bound already spent fires while the
    # cursor is still idle, and DuckDB clears an interrupt nothing was running when it
    # arrived -- so the query would then run unbounded and answer after the deadline it
    # was given. Decided from the same value that would arm the timer, so no time passes
    # between the test and the arming for the bound to expire in.
    if timeout_seconds <= 0:
        raise TimeoutError(
            "the search spent its whole bound before the query could start"
        )
    lock = threading.Lock()
    running = True

    def interrupt() -> None:
        with lock:
            if running:
                cursor.interrupt()

    timer = threading.Timer(timeout_seconds, interrupt)
    started = time.monotonic()
    timer.start()
    try:
        answered = cursor.execute(sql, parameters).to_arrow_table()
    except duckdb.InterruptException as error:
        raise TimeoutError(f"query exceeded {timeout_seconds} seconds") from error
    finally:
        timer.cancel()
        with lock:
            running = False
    # The interrupt is not guaranteed to land. A timer armed for a very short bound can
    # fire while the cursor is still idle, and DuckDB clears an interrupt that arrived
    # with nothing running, after which the query goes on unbounded. So the answer is
    # checked against the clock rather than trusted because no interrupt arrived: the
    # bound is a promise about when the caller is answered, and one kept only when a
    # race goes the right way is not a promise.
    elapsed = time.monotonic() - started
    if elapsed > timeout_seconds:
        raise TimeoutError(
            f"the query answered after {elapsed:.3f} seconds, past the "
            f"{timeout_seconds} it was given; its interrupt did not land"
        )
    return answered


class Corpus:
    """A searchable pairing of projections with their structures base.

    Opens one DuckDB connection over both artifact sets and publishes the relation a
    compiled query runs against. Given the source datasets, it also serves the reactions
    a search finds; see ``reactions``. The connection reaches the corpus's own
    directories and nothing else; see ``_confine``. Use as a context manager, or call
    ``close``.
    """

    def __init__(
        self,
        projection_pattern: str,
        structures_pattern: str,
        *,
        resolver: Callable[[str], str] | None = None,
        threads: int = -1,
        require_current: bool = True,
        pivots_dir: str | None = None,
        require_pivots: bool = False,
        occurrences_dir: str | None = None,
        require_occurrences: bool = False,
        first_commits_pattern: str | None = None,
        memory_limit: str | None = None,
        warm: bool = True,
        max_rows: int | None = DEFAULT_MAX_ROWS,
        sources_pattern: str | None = None,
    ) -> None:
        """Pairs the artifacts, checks their stamps, and prepares the relations.

        Args:
            projection_pattern: Glob matching the projection files.
            structures_pattern: Glob matching the structures base.
            resolver: Maps a compound name to SMILES. Defaults to
                ``ord_schema.resolvers``, which calls external services; inject
                something local for tests or offline use. An answer RDKit can read is
                kept for later searches, for the ``_CACHED_RESOLUTIONS`` names used most
                recently. A failure is not kept, so the next search naming that compound
                asks again.
            threads: Threads RDKit uses to match a substructure query; -1 means one
                per core. These are real threads with the GIL released, so concurrent
                searches overlap instead of serializing. The count is per search, not
                per corpus, so a server expecting several at once should divide the
                core count by that concurrency rather than leaving this at -1.
            require_current: Refuse artifacts not written by the current library,
                artifact, and RDKit versions. The screen's completeness guarantee
                assumes the query and the artifact fingerprint identically, so turning
                this off is only for reading a corpus known to match anyway.
            pivots_dir: Directory holding derived pivot artifacts, one subdirectory per
                repeated level. Given one, a quantifier over a level is a semi-join
                against a view over the level's artifacts, which answers in a fraction
                of the time a filter over the projection's lists takes and costs nothing
                to hold. A level with no subdirectory is answered from the lists, so a
                partial set of artifacts is a partial speedup rather than a missing
                answer. Derived by ``ord_schema.artifacts.scripts.derive_pivots``.
            require_pivots: Refuse a corpus whose ``pivots_dir`` does not hold every
                level, rather than answering the missing ones from the projection's
                lists. Needs a ``pivots_dir``. What it costs is a footer read per
                artifact at startup; what it buys is that no query is slow for a level
                nobody meant to leave out.
            occurrences_dir: Directory holding derived occurrence artifacts, one
                subdirectory per indexed path. Given one covering every path, the
                occurrence index is a view over those files, and an ``exists`` whose
                body is one structure predicate is a semi-join against it rather than a
                filter over the elements. A directory covering only some paths leaves
                the corpus with no index, said in a warning: every structure predicate
                is then answered from the pivots or the projection, which is slower and
                nothing else. ``require_occurrences`` refuses it instead. Derived by
                ``ord_schema.artifacts.scripts.derive_occurrences``.
            require_occurrences: Refuse a corpus whose ``occurrences_dir`` does not
                cover every indexed path, rather than running without an index. Needs
                an ``occurrences_dir``. What it costs is a footer read per artifact at
                startup.
            first_commits_pattern: Glob matching the first_commits artifacts, which hold
                when each reaction first reached ord-data. Given one, every projection
                has to pair with one, and a query naming ``first_commit`` reads them;
                no other query does. Without one, such a query is refused. Derived by
                ``ord_schema.artifacts.scripts.derive_first_commits``.
            memory_limit: What DuckDB may hold, spelled as it spells sizes: ``6500MiB``,
                ``8GB``. Left unset, DuckDB takes about 80% of the machine, or of the
                container's cap, which it does read. That default suits a process that
                is mostly DuckDB, and this one is not: RDKit, the interpreter, the
                substructure library, and most of what DuckDB's own cache of parsed
                footers holds sit outside that accounting. Exceeding a cap is a kill
                rather than an exception, so a container wants this set low enough to
                leave room for them.
            warm: Build the substructure library at open and publish the occurrence
                index, on by default. Both otherwise wait for the first query that
                wants them, so a corpus the library or the index refuses fails that
                query rather than the deployment, and a container short of the
                library's memory finds out at whoever runs the first substructure
                search. A corpus asked only scalar or similarity queries needs neither
                and can stay cold; ``check_index`` publishes the index alone.
            max_rows: Most rows a search may return, applied whether or not the query
                asked for a limit and clamping one that asks for more. The grammar
                leaves ``limit`` optional, so a query without one returns every matching
                reaction -- potentially the whole corpus, built into an Arrow table in
                this process, from a query that reads as ordinary. Defaults to
                ``DEFAULT_MAX_ROWS``; None leaves every query's own limit, or none,
                exactly as stated, which is what a caller computing over the whole match
                set wants. A result this bound may have cut short carries
                ``TRUNCATED`` in its schema metadata, since nothing in the rows tells a
                short answer from a complete one. It also bounds ``reactions``, which
                refuses a request naming more distinct IDs rather than cutting it: an ID
                the corpus does not hold has no row either, so a cut fetch would read
                as reactions that are not there.
            sources_pattern: Glob matching the source dataset files the projections
                were derived from, in ``ord_schema.parquet``'s format. Given one,
                ``reactions`` serves the records a search finds and ``datasets``
                describes each dataset, both read from the files the projections
                restate. Each source is hashed at open to pair it with its projection,
                which reads every source once, and is read on demand rather than held.

        Raises:
            PairingError: If any pattern matches nothing, a file on one side has no
                partner on another, a pair's stamps disagree about the source dataset,
                a first_commits artifact holds a different count of reactions from its
                projection,
                a projection carries a structure ID its structures artifact does not
                hold or records no statistics bounding them, or -- with
                ``require_current`` -- any artifact is stale. Given
                ``sources_pattern``, also if it matches nothing or a file it cannot read
                as a source dataset, two of its files hold the same content, a
                projection's source is not among them, or one restates a dataset no
                projection does. With ``warm``, also whatever publishing the occurrence
                index and building the substructure library refuse: a corpus stating a
                reaction twice, one whose occurrence artifacts do not reach every
                structure it holds, one whose IDs are not every integer from zero, or
                pivot or occurrence artifacts that belong to another corpus or carry a
                structure ID their dataset does not hold, all of which otherwise
                surface on the first query that wants what refused them rather than at
                open.
            ValueError: If ``max_rows`` is less than one, which is a bound no query can
                satisfy rather than a small one; if ``require_pivots`` is set without a
                ``pivots_dir``; or if ``require_occurrences`` is set without an
                ``occurrences_dir``.
        """
        self._resolver = resolver if resolver is not None else _resolve_with_resolvers
        # Names the resolver answered, canonical, most recently used last; see _resolve.
        self._resolutions: collections.OrderedDict[str, str] = collections.OrderedDict()
        self._resolutions_lock = threading.Lock()
        if max_rows is not None and max_rows <= 0:
            raise ValueError(
                f"max_rows is {max_rows}, which no query can return; leave it unset to "
                "bound nothing"
            )
        self._max_rows = max_rows
        self._threads = threads
        # Built on first substructure query; see _library. The library holds one entry
        # per distinct molecule, and _members with _starts map an entry back to the
        # structure IDs sharing it.
        self._substructure_library: rdSubstructLibrary.SubstructLibrary | None = None
        self._members = array.array("I")
        self._starts = array.array("I")
        self._library_lock = threading.Lock()
        # Whether the corpus has an occurrence index, settled by the first query that
        # could spend one; see _occurrences. A check that finds the corpus unindexable
        # keeps its reason here, since the corpus cannot change while this object is
        # open and checking again costs a scan of the artifacts to reach it.
        self._indexed: bool | None = None
        self._refused: str | None = None
        self._occurrences_lock = threading.Lock()
        # Whether the corpus states each reaction once, settled by the first route that
        # finds reactions by ID; see _require_unique_reactions. A refusal is kept for
        # the reason _refused is.
        self._unique_checked = False
        self._duplicated: str | None = None
        self._unique_lock = threading.Lock()
        # Recent match sets, most recently used last; see _matches.
        self._matched: collections.OrderedDict[_MatchKey, str] = (
            collections.OrderedDict()
        )
        # Match sets a thread is computing, so callers asking the same question while it
        # runs wait for that answer instead of repeating it; see _matches.
        self._matching: dict[_MatchKey, threading.Event] = {}
        self._matches_lock = threading.Lock()
        if require_pivots and pivots_dir is None:
            raise ValueError(
                "require_pivots was set without a pivots_dir, so there is nowhere to "
                "look for what it would require"
            )
        if require_occurrences and occurrences_dir is None:
            raise ValueError(
                "require_occurrences was set without an occurrences_dir, so there is "
                "nowhere to look for what it would require"
            )
        pairs, self._fingerprint, arrivals = _pair(
            projection_pattern,
            structures_pattern,
            require_current,
            first_commits_pattern,
        )
        self._history = bool(arrivals)
        # The fingerprint digests only the artifacts' stamps, which name the sources
        # too: a source pairs only with the projection whose stamps record its hash.
        self._sources = (
            _pair_sources(sources_pattern, pairs)
            if sources_pattern is not None
            else None
        )
        self._pivots_dir = pivots_dir
        self._occurrences_dir = occurrences_dir
        self._require_current = require_current
        # The projections a pivot artifact has to have been derived from, and the views
        # already published over the artifacts that were; see _pivot_view.
        self._projections = [pair[0] for pair in pairs]
        self._pivot_views: dict[str, str] = {}
        self._views_lock = threading.Lock()
        # The artifacts found for each level and for each indexed path, with the source
        # each names, under their own lock: finding them holds a footer read per
        # artifact, which is not something to hold the lock a published view is looked
        # up under.
        self._pivot_artifacts_found: dict[str, dict[str, str] | None] = {}
        self._occurrence_artifacts_found: dict[str, dict[str, str] | None] = {}
        self._artifacts_lock = threading.Lock()
        # What a dataset's structure IDs add to reach a corpus-wide one, and how many it
        # holds, keyed by the source every artifact derived from it names; filled by
        # _prepare.
        self._offset_of_source: dict[str, int] = {}
        self._structures_of_source: dict[str, int] = {}
        self._connection = duckdb.connect(
            config={"memory_limit": memory_limit} if memory_limit is not None else {}
        )
        try:
            _cache_footers(self._connection)
            _free_result_order(self._connection)
            _confine(
                self._connection,
                [
                    *(
                        str(pathlib.Path(path).parent)
                        for pair in pairs
                        for path in pair[:2]
                    ),
                    *(
                        str(pathlib.Path(source.path).parent)
                        for source in self._sources or ()
                    ),
                    *(str(pathlib.Path(path).parent) for path in arrivals.values()),
                    *(
                        directory
                        for directory in (pivots_dir, occurrences_dir)
                        if directory is not None
                    ),
                ],
            )
            self._total, self._searchable = self._prepare(pairs)
            if arrivals:
                self._prepare_history(pairs, arrivals)
            if self._sources is not None:
                self._prepare_sources(self._sources)
            if require_pivots:
                self._require_every_pivot()
            if require_occurrences:
                self._require_every_occurrence()
            if warm:
                self._warm()
        except Exception:
            # Nothing else holds the connection yet, and the caller has no object to
            # close one from if __init__ does not return.
            self._connection.close()
            raise

    def _warm(self) -> None:
        """Publishes the index and builds the library at open rather than on a query.

        Checks that the corpus states each reaction once whether or not it holds any
        artifact, so a corpus a pivot or the index would refuse is refused at open.

        Raises:
            PairingError: Whatever the duplicate check, the index, or the library
                refuses; see ``_require_unique_reactions``, ``_occurrences``, and
                ``_library``.
        """
        self._require_unique_reactions()
        self._occurrences()
        self._library()

    def _require_unique_reactions(self) -> None:
        """Refuses a corpus stating a reaction more than once, checked once.

        A pivot, the occurrence index, and the first commits all find reactions by ID:
        a semi-join on ``reaction_id`` answers for every row carrying that ID, and a
        join repeats each row once per copy. Two files stating the
        same reaction pair cleanly -- what pairing checks is that each projection has
        its structures, not that the corpus states a reaction once -- and each copy's
        elements would then answer the other's queries, so a question narrowed to one
        dataset is answered by the other's components. Costs a scan of one column.

        Raises:
            PairingError: If any reaction is stated by more than one row. The reason
                is kept and re-raised, since a corpus does not change under an open
                ``Corpus``.
        """
        with self._unique_lock:
            if self._duplicated is not None:
                raise PairingError(self._duplicated)
            if self._unique_checked:
                return
            # Its own cursor: this runs while searches are in flight, and the shared
            # connection holds their results.
            cursor = self._connection.cursor()
            try:
                duplicated = cursor.execute(
                    f"SELECT reaction_id, count(*) FROM {query.TABLE} "  # noqa: S608
                    "GROUP BY reaction_id HAVING count(*) > 1 "
                    "ORDER BY reaction_id LIMIT 1"
                ).fetchone()
            finally:
                cursor.close()
            if duplicated is not None:
                self._duplicated = (
                    f"reaction {duplicated[0]} is stated by {duplicated[1]} rows of "
                    "the corpus, and the pivots, the occurrence index, and the first "
                    "commits find reactions by ID, so one copy would answer for the "
                    "other; glob one artifact per reaction"
                )
                raise PairingError(self._duplicated)
            self._unique_checked = True

    def _prepare(
        self, pairs: list[tuple[str, str, str, str | None]]
    ) -> tuple[int, int]:
        """Publishes the relations, and returns the total and searchable row counts."""
        offsets = []
        total = 0
        stated = 0
        for projected, structured, source, dataset_id in pairs:
            with pq.ParquetFile(structured) as artifact:
                count = artifact.metadata.num_rows
            with pq.ParquetFile(projected) as artifact:
                stated += artifact.metadata.num_rows
            # The IDs the projection carries have to land inside the partner's rows.
            # Stamps cannot see this: an artifact rederived from a rewritten
            # projection keeps the same source hash, so a short one pairs cleanly and
            # then aliases its neighbor's molecules -- in range for get_bit, wrong
            # about the chemistry, and silent. Read from footer statistics, so it
            # costs no column data.
            largest = _max_structure_id(projected, _DERIVE_PROJECTION)
            if largest is not None and largest >= count:
                raise PairingError(
                    f"{projected} carries structure_id {largest} but {structured} "
                    f"holds only {count} structures, so its IDs would join to another "
                    "dataset's molecules; derive the structures artifact from this "
                    "projection again"
                )
            offsets.append((projected, structured, total, dataset_id))
            # Keyed by the source rather than by the file, so an artifact derived from
            # this projection finds the offset wherever it is filed; see
            # ``_artifact_offsets``.
            self._offset_of_source[source] = total
            self._structures_of_source[source] = count
            total += count
        registered = pa.table(
            {
                "projection_filename": [offset[0] for offset in offsets],
                "structures_filename": [offset[1] for offset in offsets],
                query.STRUCTURE_OFFSET: pa.array(
                    [offset[2] for offset in offsets], type=pa.int64()
                ),
                query.DATASET_ID: pa.array(
                    [offset[3] for offset in offsets], type=pa.string()
                ),
            }
        )
        # A registered object belongs to the connection that registered it, so a
        # per-search cursor raises a catalog error on any view built over one. Copying
        # it into the catalog puts the views below within reach of every cursor.
        self._connection.register("registered_offsets", registered)
        self._connection.execute(
            "CREATE TABLE structure_offsets AS SELECT * FROM registered_offsets"
        )
        self._connection.unregister("registered_offsets")
        # Inlined rather than bound because DuckDB refuses parameters in DDL. The
        # paths came from this process's own glob, and the quoting still escapes them.
        projection_files = _sql_paths(pair[0] for pair in pairs)
        structure_files = _sql_paths(pair[1] for pair in pairs)
        # S608: the fragments are this module's constants and this process's own
        # glob results, quoted with their single quotes escaped.
        self._connection.execute(
            f"""
            CREATE VIEW {query.TABLE} AS
            SELECT p.* EXCLUDE (filename), o.{query.STRUCTURE_OFFSET},
                   o.{query.DATASET_ID}
            FROM read_parquet({projection_files}, filename=true) p
            JOIN structure_offsets o ON p.filename = o.projection_filename
            """  # noqa: S608
        )
        self._connection.execute(
            f"""
            CREATE VIEW corpus_structures AS
            SELECT s.* EXCLUDE (filename),
                   (s.structure_id + o.{query.STRUCTURE_OFFSET})::BIGINT AS global_id
            FROM read_parquet({structure_files}, filename=true) s
            JOIN structure_offsets o ON s.filename = o.structures_filename
            """  # noqa: S608
        )
        # Both views join on a path that Python read and DuckDB re-resolved. Anything
        # that makes the two disagree -- a glob metacharacter in a directory name is
        # the reachable one, since read_parquet globs each element it is handed --
        # drops rows from the join rather than failing, leaving a corpus that answers
        # every query with silence. Each side is counted against the total its own
        # footers state, which is where the loss is visible: a reaction the view lost is
        # one no query can find, and the occurrence index cannot stand in for the count,
        # since a path read from a pivot artifact reaches the structures whatever the
        # view holds.
        counts = self._connection.execute(
            "SELECT count(*), count(pattern_fp) FROM corpus_structures"
        ).fetchone()
        assert counts is not None  # An aggregate over any relation returns one row.
        joined, searchable = counts
        if joined != total:
            raise PairingError(
                f"the structures artifacts hold {total} rows but only {joined} joined "
                "to their offsets; a path did not survive read_parquet"
            )
        reactions = self._connection.execute(
            f"SELECT count(*) FROM {query.TABLE}"  # noqa: S608
        ).fetchone()
        assert reactions is not None  # An aggregate over any relation returns one row.
        if reactions[0] != stated:
            raise PairingError(
                f"the projections hold {stated} reactions but only {reactions[0]} "
                "joined to their offsets; a path did not survive read_parquet"
            )
        return total, searchable

    def _prepare_history(
        self, pairs: list[tuple[str, str, str, str | None]], arrivals: dict[str, str]
    ) -> None:
        """Publishes ``query.HISTORY_TABLE``: each reaction's first commit to ord-data.

        A view over the artifacts rather than a table, so the corpus holds none of it,
        and nothing reads it but a query naming ``first_commit``.

        Args:
            pairs: The corpus's projections, paired by source.
            arrivals: Each source's first_commits artifact.

        Raises:
            PairingError: If an artifact holds a different count of reactions from its
                projection, or the view loses rows to a path read_parquet re-resolved.
        """
        stated = 0
        for projected, _, source, _ in pairs:
            # Both counts are footer reads. Pairing by source hash cannot see an
            # artifact edited after it was written, which keeps its stamps, and one
            # short of its projection leaves reactions with no first commit.
            with pq.ParquetFile(projected) as artifact:
                expected = artifact.metadata.num_rows
            with pq.ParquetFile(arrivals[source]) as artifact:
                count = artifact.metadata.num_rows
            if count != expected:
                raise PairingError(
                    f"{arrivals[source]} holds {count} first commits but {projected} "
                    f"holds {expected} reactions; derive it from this projection again"
                )
            stated += count
        files = _sql_paths(arrivals[pair[2]] for pair in pairs)
        # S608: the relation name is this module's constant and the paths this
        # process's own glob results, quoted with their single quotes escaped.
        self._connection.execute(
            f"""
            CREATE VIEW {query.HISTORY_TABLE} AS
            SELECT reaction_id, {first_commits.COLUMN}
            FROM read_parquet({files})
            """  # noqa: S608
        )
        held = self._connection.execute(
            f"SELECT count(*) FROM {query.HISTORY_TABLE}"  # noqa: S608
        ).fetchone()
        assert held is not None  # An aggregate over any relation returns one row.
        if held[0] != stated:
            raise PairingError(
                f"the first_commits artifacts hold {stated} rows but only {held[0]} "
                "were read; a path did not survive read_parquet"
            )

    def _prepare_sources(self, sources: list[_Source]) -> None:
        """Publishes ``corpus_sources``: every source's records beside its dataset ID.

        A view over the files rather than a table, so the corpus holds none of the
        reactions. The ID filter a fetch applies reaches the Parquet scan through the
        join as a range on ``reaction_id``, the same scan a read of the files alone
        gets; see ``reactions`` for what that range costs.

        Args:
            sources: The paired sources, from ``_pair_sources``.

        Raises:
            PairingError: If the view holds fewer reactions than the sources' footers
                state.
        """
        registered = pa.table(
            {
                "source_filename": pa.array(
                    [source.path for source in sources], type=pa.string()
                ),
                query.DATASET_ID: pa.array(
                    [source.dataset_id for source in sources], type=pa.string()
                ),
            }
        )
        # Copied into the catalog for the reason structure_offsets is; see _prepare.
        self._connection.register("registered_sources", registered)
        self._connection.execute(
            "CREATE TABLE source_catalog AS SELECT * FROM registered_sources"
        )
        self._connection.unregister("registered_sources")
        source_files = _sql_paths(source.path for source in sources)
        # S608: the fragments are this module's constants and this process's own glob
        # results, quoted with their single quotes escaped.
        self._connection.execute(
            f"""
            CREATE VIEW corpus_sources AS
            SELECT s.reaction_id, c.{query.DATASET_ID}, s.reaction
            FROM read_parquet({source_files}, filename=true) s
            JOIN source_catalog c ON s.filename = c.source_filename
            """  # noqa: S608
        )
        # The join loses a file whose path DuckDB resolves differently than Python did,
        # as _prepare's do, and a reaction lost there is one a search finds and a fetch
        # does not return.
        stated = sum(source.reactions for source in sources)
        joined = self._connection.execute(
            "SELECT count(*) FROM corpus_sources"
        ).fetchone()
        assert joined is not None  # An aggregate over any relation returns one row.
        if joined[0] != stated:
            raise PairingError(
                f"the sources hold {stated} reactions but only {joined[0]} joined to "
                "their dataset IDs; a path did not survive read_parquet"
            )

    @property
    def fingerprint(self) -> str:
        """Returns a short digest naming the artifacts this corpus opened.

        A question log records this beside the query rather than the rows the query
        returned, so an old record stays reproducible: the same query against the same
        fingerprint answers the same way, and a different fingerprint says why it does
        not.
        """
        return self._fingerprint

    def _require_every_pivot(self) -> None:
        """Refuses a corpus whose pivots_dir does not hold every level, current.

        ``check_pivots`` reports what is held; this insists. A level with no artifacts
        is otherwise answered from the projection's lists on every query over it, which
        is slower and no error anywhere -- the failure a caller naming a pivots_dir is
        least likely to be expecting, since naming one says the artifacts are meant to
        answer.

        Raises:
            PairingError: If any level has no artifacts. Wrong provenance, a level
                filed under another, a structure ID outside its dataset, and staleness
                are raised by the publishing itself.
        """
        missing = [
            path for path in sorted(pivot.LEVELS) if self._pivot_view(path) is None
        ]
        if missing:
            raise PairingError(
                f"{self._pivots_dir} holds no pivot artifacts for {missing}, which "
                "this corpus was opened requiring; derive them first with "
                f"{_DERIVE_PIVOTS}"
            )

    def check_pivots(self) -> dict[str, int]:
        """Publishes every pivot artifact this corpus can read, and returns the counts.

        Meant for startup. Artifacts are otherwise found by the first query that wants
        a level, so a tree derived from another corpus, filed under the wrong level, or
        left stale is a ``PairingError`` raised at whatever hour that query arrives.
        Calling this makes the same refusal happen while a server is still starting,
        where it is a failed deployment rather than a failed request.

        Publishing is the check: a level whose artifacts are wrong cannot be published,
        and one that is fine is left ready, so the first real query does not pay the
        glob. Cheap either way -- a view over Parquet reads footers, not rows.

        Returns:
            The number of artifacts found per level, for the levels that have any.
            Empty when the corpus was opened without ``pivots_dir``, which is not an
            error: every level is then answered from the projection's lists.

        Raises:
            PairingError: If any level's artifacts were not derived from this corpus's
                projections, hold a different level than the one they are filed under,
                carry a structure ID their dataset does not hold, or -- with
                ``require_current`` -- are stale.
        """
        if self._pivots_dir is None:
            return {}
        found = {}
        for path in sorted(pivot.LEVELS):
            if self._pivot_view(path) is not None:
                found[path] = len(pivot.artifact_paths(self._pivots_dir, path))
        logger.info(
            "%d of %d levels are held as artifacts: %s",
            len(found),
            len(pivot.LEVELS),
            ", ".join(f"{path} {count}" for path, count in found.items()) or "none",
        )
        return found

    def check_index(self) -> dict[str, int]:
        """Publishes the occurrence index and returns what it reached, per path.

        The sibling of ``check_pivots``, and for the same reason: a corpus opened
        without ``warm`` publishes the index on the first query that can spend it, so a
        corpus it refuses -- a reaction stated twice, a structure no artifact reaches,
        an artifact of another corpus -- fails that query rather than the deployment.
        A corpus opened with ``warm`` has published it already, and the call is then
        its counts. Leaving a corpus cold suits one asked only for scalar or similarity
        queries, which never spend the library ``warm`` would also build.

        Returns:
            The number of occurrences found per indexed path, including the paths that
            hold none -- most corpora record no authentic standards, and a zero there is
            ordinary rather than a sign the walk missed something. Empty when the corpus
            has no index: it was opened without ``occurrences_dir``, or with one that
            does not cover every indexed path.

        Raises:
            PairingError: If the corpus states a reaction more than once, the index
                does not reach every structure the corpus holds, or an occurrence
                artifact belongs to another corpus or carries a structure ID its
                dataset does not hold.
        """
        if not self._occurrences():
            return {}
        cursor = self._connection.cursor()
        try:
            counts = dict(
                cursor.execute(
                    "SELECT path, count(*) FROM occurrences GROUP BY path"
                ).fetchall()
            )
        finally:
            cursor.close()
        return {path: counts.get(path, 0) for path in INDEXED_PATHS}

    def __enter__(self) -> Self:
        """Returns the corpus itself; closing on exit is the whole protocol."""
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Closes the connection."""
        self.close()

    def close(self) -> None:
        """Closes the connection."""
        self._connection.close()

    def _query_molecule(
        self, parameter: query.StructureParameter, resolve: Callable[[str], str]
    ) -> Chem.Mol:
        """Returns the query molecule for a structure predicate.

        Args:
            parameter: The predicate to build a molecule for.
            resolve: Maps a compound name to SMILES.

        Returns:
            A pattern from SMARTS for a substructure predicate stated as one; a
            molecule from SMILES otherwise (including every resolved compound).

        Raises:
            ValueError: If a resolved compound's SMILES does not parse -- the
                resolver's contract is RDKit-readable SMILES, so this is a resolver
                bug, not a query error.
        """
        if _reads_as_smarts(parameter):
            assert parameter.pattern is not None  # What reading it as SMARTS means.
            return Chem.MolFromSmarts(parameter.pattern)  # Validated at model time.
        if parameter.pattern is not None:
            smiles = parameter.pattern
        else:
            assert parameter.compound is not None  # One of the two is always set.
            smiles = resolve(parameter.compound)
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(
                f"resolver returned unreadable SMILES {smiles!r} "
                f"for {parameter.compound!r}"
            )
        return molecule

    def _library(self) -> rdSubstructLibrary.SubstructLibrary:
        """Returns the substructure library, building it on first use.

        Built at open under ``warm``, and otherwise on the first substructure query: it
        costs a build and the memory to hold it, and a corpus asked only for similarity
        or scalar queries never needs it.

        Structures are deduplicated per dataset, so a molecule used by several datasets
        occupies a row in each and would otherwise be matched once per copy. The library
        holds one entry per distinct SMILES -- over ORD, markedly fewer than the rows,
        and the matching shrinks with them -- and ``_members`` maps an entry back to
        every structure ID sharing its molecule. Two structures with one SMILES have one
        molecule and one fingerprint, both derived from that SMILES by the RDKit the
        stamps pin, so which row is read does not matter.

        Concurrent first searches serialize on the build rather than each building a
        copy, which would hold a library's memory apiece. Whoever waits needs the
        finished library anyway.

        Returns:
            A library over every distinct molecule in the corpus, screened by the
            pattern fingerprints the artifact already stores.

        Raises:
            PairingError: If the corpus IDs are not exactly ``0`` to ``total - 1``, so
                that the position an ID is assumed to occupy holds another structure; if
                the library and the entry table disagree about how many distinct
                molecules there are; or if a structure holds one of its derived columns
                without the other.
        """
        with self._library_lock:
            if self._substructure_library is not None:
                return self._substructure_library
            start = time.perf_counter()
            molecules = rdSubstructLibrary.CachedMolHolder()
            patterns = rdSubstructLibrary.PatternHolder(structures.PATTERN_FP_SIZE)
            cursor = self._connection.cursor()
            # Which library entry each structure ID resolves to, in ID order. An array
            # of unsigned ints rather than a list, which holds four bytes per ID where a
            # list holds a pointer to an int object.
            entry_of = array.array("I")
            entries: dict[str, int] = {}
            position = 0
            try:
                reader = cursor.execute(
                    "SELECT global_id, smiles, mol_binary, pattern_fp "
                    "FROM corpus_structures ORDER BY global_id"
                ).to_arrow_reader(_BUILD_BATCH)
                for batch in reader:
                    identifiers = batch.column("global_id").to_pylist()
                    smiles_values = batch.column("smiles").to_pylist()
                    blobs = batch.column("mol_binary").to_pylist()
                    fingerprints = batch.column("pattern_fp").to_pylist()
                    for global_id, smiles, blob, fingerprint in zip(
                        identifiers, smiles_values, blobs, fingerprints, strict=True
                    ):
                        # An ID is read as the position it occupies in this pass, which
                        # it is only while the IDs are every integer from zero, each
                        # once. A duplicate or a gap slides every later structure one
                        # position along, which a row count cannot see.
                        if global_id != position:
                            raise PairingError(
                                f"the corpus states structure ID {global_id} where "
                                f"{position} was expected, so its IDs are not one "
                                "unbroken run and every later structure sits at a "
                                "position other than its own ID; derive the structures "
                                "artifacts again"
                            )
                        # A structures artifact writes the derived columns together or
                        # not at all. One without the other means the row came from
                        # somewhere else, and taking either branch would be a guess:
                        # dropping a live fingerprint makes the structure unmatchable
                        # while it still counts as searchable, and a missing one
                        # reaches RDKit as a type error naming no row.
                        if (blob is None) != (fingerprint is None):
                            raise PairingError(
                                f"structure {global_id} has "
                                f"{'no ' if blob is None else 'a '}molecule but "
                                f"{'no ' if fingerprint is None else 'a '}fingerprint; "
                                "a structures artifact writes both or neither, so this "
                                "one was not written by this library"
                            )
                        entry = entries.get(smiles)
                        if entry is None:
                            entry = entries[smiles] = len(entries)
                            if blob is None:
                                molecules.AddBinary(_UNPARSEABLE)
                                patterns.AddFingerprint(_NO_BITS)
                            else:
                                molecules.AddBinary(blob)
                                patterns.AddFingerprint(
                                    DataStructs.CreateFromBinaryText(fingerprint)
                                )
                        entry_of.append(entry)
                        position += 1
            finally:
                cursor.close()
            library = rdSubstructLibrary.SubstructLibrary(molecules, patterns)
            if position != self._total:
                raise PairingError(
                    f"the library covers {position} structures but the corpus has "
                    f"{self._total}, so a structure ID would fall outside the run this "
                    "read and resolve to no entry"
                )
            # The holders and the entry table are filled in the same branch, so a
            # divergence means one of them missed a row -- after which every entry above
            # the divergence maps to a neighbor's structures: in range, wrong molecule,
            # and invisible to a count of either one alone.
            if len(library) != len(entries):
                raise PairingError(
                    f"the library holds {len(library)} molecules against "
                    f"{len(entries)} distinct SMILES, so an entry index does not name "
                    "the molecule it was built from"
                )
            self._members, self._starts = _group(entry_of, len(entries))
            logger.info(
                "built a substructure library over %d distinct molecules covering "
                "%d structures in %.1fs",
                len(library),
                position,
                time.perf_counter() - start,
            )
            self._substructure_library = library
            return self._substructure_library

    def _artifact_offsets(self, published: str, sources: dict[str, str]) -> str:
        """Publishes what each artifact's rows add to reach a corpus-wide ID.

        An artifact and the projection it descends from restate one source dataset, so
        they carry the same ``source_md5`` however many derivations apart they sit, and
        that hash is what pairs an artifact with the offset its projection was given.
        Pairing on the filename instead would rest on the two trees being laid out
        alike, which nothing enforces.

        Args:
            published: The relation to create, named for the level or path it serves so
                that a retry after a refusal replaces these rows rather than leaving a
                relation behind per attempt.
            sources: The source each artifact names, from the check that accepted them.

        Returns:
            ``published``.
        """
        rows = []
        for name, source in sources.items():
            # The checks require these artifacts' sources to be exactly the
            # projections', and the offsets are keyed by those same hashes.
            offset = self._offset_of_source.get(source)
            assert offset is not None, f"{name} names a source the corpus does not hold"
            rows.append((name, offset))
        registered = pa.table(
            {
                "artifact_filename": [row[0] for row in rows],
                query.STRUCTURE_OFFSET: pa.array(
                    [row[1] for row in rows], type=pa.int64()
                ),
            }
        )
        # Called under _occurrences_lock, so two cursors never register one name at
        # once. S608: the name comes from this module's walk of the schema.
        self._connection.register("registered_artifact_offsets", registered)
        self._connection.execute(
            f"CREATE OR REPLACE TABLE {published} AS "  # noqa: S608
            "SELECT * FROM registered_artifact_offsets"
        )
        self._connection.unregister("registered_artifact_offsets")
        return published

    def _occurrence_artifacts(self, path: str) -> dict[str, str] | None:
        """Returns the artifacts holding ``path``'s occurrences, checked, or None.

        Globbed and checked once per path and remembered, including the answer that
        there are none, so a corpus without artifacts reads the directory once rather
        than on every query.

        Args:
            path: The indexed path, as the query grammar names it.

        Returns:
            The source dataset each artifact was derived from, keyed by its path, or
            None where the corpus has no directory or the directory holds none for this
            path.

        Raises:
            PairingError: If the artifacts do not belong to this corpus; see
                ``_check_split``.
        """
        if self._occurrences_dir is None:
            return None
        with self._artifacts_lock:
            if path in self._occurrence_artifacts_found:
                return self._occurrence_artifacts_found[path]
            files = [
                str(found)
                for found in occurrences.artifact_paths(self._occurrences_dir, path)
            ]
            if not files:
                self._occurrence_artifacts_found[path] = None
                return None
            sources = self._check_split(
                path,
                files,
                artifact=occurrences.ARTIFACT,
                schema=occurrences.SCHEMA,
                held=occurrences.occurrence_path,
                reading="at",
                derive=_DERIVE_OCCURRENCES,
            )
            self._occurrence_artifacts_found[path] = sources
            return sources

    def _occurrences_from_artifact(self, path: str) -> str | None:
        """Returns the SELECT reading ``path``'s occurrences from artifacts, or None.

        The artifact is already one row per occurrence, so the read adds only the
        corpus-wide ID: the element's own plus its dataset's offset, which is a running
        total over the corpus's pairs and so belongs to no artifact. It is joined here
        from the file each row came from, keyed by the source that file and its
        projection both name.

        Args:
            path: The indexed path.

        Returns:
            The SELECT, or None where the corpus has no artifacts for the path.

        Raises:
            PairingError: If the artifacts do not belong to this corpus; see
                ``_check_split``.
        """
        sources = self._occurrence_artifacts(path)
        if sources is None:
            return None
        offsets = self._artifact_offsets(
            f"occurrences_{path.replace('.', '_')}_offsets", sources
        )
        # The artifact holds no row for an element carrying no structure, so nothing
        # here filters: what a pivot needs a WHERE for, the derivation already did.
        # S608: the path is a key of INDEXED_PATHS and the filenames come from this
        # process's own glob, quoted with their single quotes escaped.
        return f"""
            SELECT x.reaction_id, {sql_string(path)} AS path,
                   (x.structure_id + o.{query.STRUCTURE_OFFSET})::UINTEGER AS global_id
            FROM read_parquet({_sql_paths(sources)}, filename=true) x
            JOIN {offsets} o ON x.filename = o.artifact_filename
            """  # noqa: S608

    def _require_every_occurrence(self) -> None:
        """Refuses a corpus whose occurrences_dir does not hold every indexed path.

        A path with no artifact otherwise leaves the corpus with no index at all, and
        every structure predicate answered from the pivots or the projection, which is
        slower on every query and no error anywhere -- the failure a caller naming an
        occurrences_dir is least likely to be expecting.

        Raises:
            PairingError: If any indexed path has no artifacts. Wrong provenance, a
                path filed under another, a structure ID outside its dataset, and
                staleness are raised by the check itself.
        """
        missing = [
            path
            for path in sorted(INDEXED_PATHS)
            if self._occurrence_artifacts(path) is None
        ]
        if missing:
            raise PairingError(
                f"{self._occurrences_dir} holds no {occurrences.ARTIFACT} artifacts "
                f"for {missing}, which this corpus was opened requiring; derive them "
                f"first with {_DERIVE_OCCURRENCES}"
            )

    def _occurrences(self) -> bool:
        """Publishes the occurrence index as a view over its artifacts, once.

        One row per structure occurrence, carrying the corpus-wide ID and the path it
        sat at. Published only over an ``occurrences_dir`` covering every indexed path:
        a corpus with no directory, or one covering some paths and not others, has no
        index, and the compiler answers every structure quantifier from the pivots or
        the projection. Decided for the corpus rather than per path because the check
        below is what makes the index trustworthy, and it needs every path to count.

        Costs a scan of the artifacts to check that they reach every structure the
        corpus holds, charged at open under ``warm`` and otherwise to the first query
        that could spend the index. Concurrent first queries serialize here and share
        the answer.

        Returns:
            Whether the corpus has an index.

        Raises:
            PairingError: If the corpus states a reaction more than once, or the
                artifacts do not reach every structure the corpus holds. An unreached
                structure is one whose reactions no routed query can find, and the
                answer comes back empty rather than wrong -- which reads like an answer
                rather than like a corpus that cannot be searched. The reason is kept
                and re-raised, since a corpus does not change under an open ``Corpus``
                and checking again is a scan to reach the same refusal.
        """
        with self._occurrences_lock:
            if self._refused is not None:
                raise PairingError(self._refused)
            if self._indexed is not None:
                return self._indexed
            start = time.perf_counter()
            reads = []
            missing = []
            for path in INDEXED_PATHS:
                reading = self._occurrences_from_artifact(path)
                if reading is None:
                    missing.append(path)
                else:
                    reads.append(reading)
            if missing:
                if self._occurrences_dir is None:
                    logger.info(
                        "opened without occurrences_dir, so structure predicates are "
                        "answered from the pivots or the projection"
                    )
                else:
                    logger.warning(
                        "%s holds no %s artifacts for %s, so this corpus has no "
                        "occurrence index and answers every structure predicate from "
                        "the pivots or the projection; derive them with %s",
                        self._occurrences_dir,
                        occurrences.ARTIFACT,
                        sorted(missing),
                        _DERIVE_OCCURRENCES,
                    )
                self._indexed = False
                return False
            self._require_unique_reactions()
            selects = "\nUNION ALL\n".join(reads)
            # Its own cursor: this runs while searches are in flight, and the shared
            # connection holds their results.
            cursor = self._connection.cursor()
            try:
                # OR REPLACE, so a publish interrupted after the view exists is one the
                # next query repeats rather than one that collides with itself forever.
                # S608: the fragments are this module's own reads of the artifacts, not
                # anything a query supplies.
                cursor.execute(f"CREATE OR REPLACE VIEW occurrences AS {selects}")
                counts = dict(
                    cursor.execute(
                        "SELECT path, count(*) FROM occurrences GROUP BY path"
                    ).fetchall()
                )
                reached = cursor.execute(
                    "SELECT count(DISTINCT global_id) FROM occurrences"
                ).fetchone()
                assert reached is not None  # An aggregate returns one row.
                # A structures artifact holds one row per distinct structure its
                # projection uses, so every ID it states is one some element carries and
                # the index has to find all of them. Counting them is what catches an
                # artifact short of rows -- truncated or edited -- which keeps the
                # stamps and columns its reader checks; a total over every path cannot,
                # since the other paths carry the total and the short one answers its
                # queries with silence. A count equal to the total is every structure
                # only because every ID read lies inside its own dataset's range, which
                # _prepare holds the projections to and _check_split the artifacts.
                if reached[0] != self._total:
                    self._refused = (
                        f"the occurrence index reached {reached[0]} of the corpus's "
                        f"{self._total} structures, so the reactions holding the rest "
                        "cannot be found; the likeliest cause is an occurrence "
                        f"artifact missing rows; re-derive {self._occurrences_dir} "
                        f"with {_DERIVE_OCCURRENCES} --force, since a run without it "
                        "skips a file whose stamps are current"
                    )
                    raise PairingError(self._refused)
                self._indexed = True
            finally:
                cursor.close()
            # Per path, so how the occurrences are distributed is visible rather than
            # only their total. A path holding none is normal -- most corpora have no
            # authentic standards -- which is why the refusal above counts structures
            # rather than paths.
            logger.info(
                "indexed %d structure occurrences in %.1fs as a view over %s: %s",
                sum(counts.values()),
                time.perf_counter() - start,
                self._occurrences_dir,
                ", ".join(f"{path} {counts.get(path, 0)}" for path in INDEXED_PATHS),
            )
            return True

    def _substructure_ids(
        self, parameter: query.StructureParameter, resolve: Callable[[str], str]
    ) -> list[int]:
        """Screens and verifies a substructure predicate; returns global IDs.

        RDKit screens and matches in one call, across its own threads with the GIL
        released, so concurrent searches overlap here instead of serializing.

        Args:
            parameter: The predicate to match.
            resolve: Maps a compound name to SMILES.

        Returns:
            The corpus-wide structure IDs whose molecules contain the query.
        """
        molecule = self._query_molecule(parameter, resolve)
        library = self._library()
        # Stereo groups count wherever stereocenters do. A racemate recorded as
        # "C[C@H](N)C(=O)O |&1:1|" holds both enantiomers, and matching its center as
        # drawn would make the answer depend on which way round the record wrote it.
        parameters = Chem.SubstructMatchParameters()
        parameters.useChirality = parameter.chirality
        parameters.useEnhancedStereo = parameter.chirality
        # maxResults defaults to 1000, which would silently truncate: a broad pattern
        # matches far more of ORD's distinct molecules than that.
        matched = library.GetMatches(
            molecule,
            parameters,
            numThreads=self._threads,
            maxResults=len(library) or 1,
        )
        # A library entry is a molecule, not a structure: every ID sharing that molecule
        # matched too, and the answer is stated in IDs.
        identifiers: list[int] = []
        for entry in matched:
            identifiers.extend(
                self._members[self._starts[entry] : self._starts[entry + 1]]
            )
        return identifiers

    def _similarity_ids(
        self,
        cursor: duckdb.DuckDBPyConnection,
        parameter: query.StructureParameter,
        resolve: Callable[[str], str],
    ) -> list[int]:
        """Screens a similarity predicate; the fingerprint is the whole answer."""
        molecule = self._query_molecule(parameter, resolve)
        fingerprint = structures.morgan_fingerprint(molecule)
        blob = DataStructs.BitVectToBinaryText(fingerprint)
        popcount = fingerprint.GetNumOnBits()
        threshold = parameter.threshold
        assert threshold is not None  # A similarity predicate always carries one.
        # Tanimoto >= t bounds the candidate popcount to [t * q, q / t], which the
        # artifact's morgan_popcount column answers without touching the fingerprints.
        rows = cursor.execute(
            """
            SELECT global_id FROM corpus_structures
            WHERE morgan_fp IS NOT NULL
              AND morgan_popcount BETWEEN $lo AND $hi
              AND bit_count(CAST(morgan_fp AS BITSTRING)
                            & CAST($q AS BITSTRING))::DOUBLE
                  / ($n + morgan_popcount
                     - bit_count(CAST(morgan_fp AS BITSTRING) & CAST($q AS BITSTRING)))
                  >= $t
            """,
            {
                "q": blob,
                "n": popcount,
                "lo": math.ceil(threshold * popcount),
                "hi": math.floor(popcount / threshold),
                "t": threshold,
            },
        ).fetchall()
        return [row[0] for row in rows]

    def _hashed_ids(
        self,
        cursor: duckdb.DuckDBPyConnection,
        parameter: query.StructureParameter,
        resolve: Callable[[str], str],
    ) -> list[int]:
        """Matches on one of the artifact's hash columns; returns global IDs.

        The artifact already holds every structure's hashes, so the whole answer is one
        equality against a column -- no screen, no verification, and no chemistry
        beyond hashing the query molecule the way the artifact hashed its own.

        Args:
            cursor: The cursor the match runs on.
            parameter: The predicate to match, whose ``op`` chooses the column.
            resolve: Maps a compound name to SMILES.

        Returns:
            The corpus-wide structure IDs sharing the query molecule's hash.

        Raises:
            ValueError: If a resolved compound's SMILES does not parse, or if RDKit
                refuses to hash the query molecule -- either leaves the question
                unaskable rather than quietly unanswerable.
        """
        column, hasher = _HASH_MATCHES[parameter.op]
        molecule = self._query_molecule(parameter, resolve)
        hashed = hasher(molecule)
        if hashed is None:
            raise ValueError(
                f"could not hash the query molecule for {parameter.name!r}; the corpus "
                "compares compounds by hash, so there is nothing to compare it to"
            )
        # S608: the column comes from this module's own table, keyed by an op the
        # grammar holds to a Literal; the molecule reaches the query as a parameter.
        rows = cursor.execute(
            f"SELECT global_id FROM corpus_structures WHERE {column} = $h",  # noqa: S608
            {"h": hashed},
        ).fetchall()
        return [row[0] for row in rows]

    def _bitmap(self, matched: Sequence[int]) -> str:
        """Returns the match set as a bitmap over the corpus-wide ID space."""
        bits = bytearray(b"0" * max(self._total, 1))
        for global_id in matched:
            bits[global_id] = ord("1")
        return bits.decode()

    def _pivot_view(self, path: str) -> str | None:
        """Returns a view over the pivot artifacts for ``path``, or None if none exist.

        Published once per level, over whichever artifacts ``_pivot_artifacts`` found
        and checked, and remembered under its own lock so a query over a published level
        waits on nothing. The compiler's ``PivotIndex``: a quantifier over a level with
        a view is a semi-join against it, and one over a level without is a filter over
        the projection's lists.

        Args:
            path: The repeated level, as the query grammar names it.

        Returns:
            The view's name, or None to leave the level to the projection's lists --
            which includes a path that names no level at all.

        Raises:
            PairingError: If the corpus states a reaction more than once, which a
                semi-join on ``reaction_id`` cannot tell apart; if the artifacts for
                this level were not derived from the projections this corpus reads -- a
                pivot of another corpus answers confidently and wrongly, since its
                reaction IDs resolve against the wrong rows -- if any lacks a column
                this library reads, if any carries a structure ID its dataset does not
                hold, or, with ``require_current``, if any of them is stale.
        """
        if path not in pivot.LEVELS:
            return None
        sources = self._pivot_artifacts(path)
        if sources is None:
            return None
        self._require_unique_reactions()
        files = list(sources)
        with self._views_lock:
            published = self._pivot_views.get(path)
            if published is not None:
                return published
            name = pivot.table_name(path)
            # S608: the name comes from this module's walk of the schema, and the paths
            # from this process's own glob, quoted with their single quotes escaped.
            self._connection.execute(
                f"CREATE OR REPLACE VIEW {name} AS "  # noqa: S608
                f"SELECT * FROM read_parquet({_sql_paths(files)})"
            )
            logger.info("read %d pivot artifacts for %s", len(files), path)
            self._pivot_views[path] = name
            return name

    def _pivot_artifacts(self, path: str) -> dict[str, str] | None:
        """Returns the artifacts holding the pivot over ``path``, checked, or None.

        Globbed and checked once per level and remembered, including the answer that
        there are none, so a corpus without artifacts reads the directory once rather
        than on every query. The stamps the check reads are kept with the paths, since
        the source an artifact names is what pairs it with an offset and reading the
        footers again would be one open per artifact per level.

        Args:
            path: The repeated level, as the query grammar names it.

        Returns:
            The source dataset each artifact was derived from, keyed by its path, or
            None where the directory holds none for this level.

        Raises:
            PairingError: If the artifacts were not derived from the projections this
                corpus reads, if any lacks a column this library reads, if any carries
                a structure ID its dataset does not hold, or -- with
                ``require_current`` -- if any is stale.
        """
        if self._pivots_dir is None:
            return None
        with self._artifacts_lock:
            if path in self._pivot_artifacts_found:
                return self._pivot_artifacts_found[path]
            files = [
                str(found) for found in pivot.artifact_paths(self._pivots_dir, path)
            ]
            if not files:
                self._pivot_artifacts_found[path] = None
                return None
            sources = self._check_split(
                path,
                files,
                artifact=pivot.ARTIFACT,
                schema=pivot.schema(pivot.LEVELS[path]),
                held=pivot.pivot_path,
                reading="over",
                derive=_DERIVE_PIVOTS,
            )
            self._pivot_artifacts_found[path] = sources
            return sources

    def _check_split(
        self,
        path: str,
        files: list[str],
        *,
        artifact: str,
        schema: pa.Schema,
        held: Callable[[str], str | None],
        reading: str,
        derive: str,
    ) -> dict[str, str]:
        """Refuses artifacts of a tree-split artifact that do not belong to this corpus.

        These artifacts name their reactions by ID, and an ID resolves against whatever
        rows the projections hold. Artifacts derived from another corpus therefore
        answer a quantifier confidently and wrongly rather than failing, which is why
        the check is on the source datasets rather than on the file count.

        Args:
            path: The level or indexed path the files claim to hold.
            files: The artifacts found for it.
            artifact: Artifact name every file must hold.
            schema: The columns this library's version of that artifact carries at
                ``path``.
            held: Reads the level or path a file is stamped with.
            reading: How a message names the relation to ``path`` -- "over" for a pivot
                level, "at" for an indexed path.
            derive: The module that writes these artifacts, which a refusal names.

        Returns:
            The source dataset each artifact was derived from, keyed by its path, which
            is what pairs an artifact with the offset its projection was given; see
            ``_artifact_offsets``.

        Raises:
            PairingError: If a file is not this artifact over ``path``, if it lacks a
                column this library reads, if two artifacts restate one source dataset,
                if the set of source datasets differs from the projections', if one
                carries a structure ID its dataset does not hold or records no
                statistics bounding them, or -- with ``require_current`` -- if any
                artifact is stale.
        """
        wanted = {base.load_stamps(name).source_md5 for name in self._projections}
        sources: dict[str, str] = {}
        derived_from: dict[str, str] = {}
        for name in files:
            stamps = base.load_stamps(name)
            _refuse_unreadable(
                name,
                stamps,
                artifact=artifact,
                schema=schema,
                require_current=self._require_current,
                held=held,
                expected=path,
                reading=reading,
            )
            # Two artifacts of one dataset pass the set comparison below and are both
            # read, so every element is stated twice. A quantifier cannot see it -- a
            # semi-join returns a reaction once however many rows name it -- and
            # neither can the structure count that decides whether the index reaches
            # the corpus, which is taken over distinct IDs. What does see it is
            # check_index, which reports the doubled count as the corpus's own.
            if stamps.source_md5 in derived_from:
                raise PairingError(
                    f"{name} and {derived_from[stamps.source_md5]} are both {artifact} "
                    f"artifacts {reading} {path} of the same source dataset, so its "
                    "elements would each be read twice"
                )
            derived_from[stamps.source_md5] = name
            sources[name] = stamps.source_md5
        if set(sources.values()) != wanted:
            found = set(sources.values())
            raise PairingError(
                f"the {artifact} artifacts {reading} {path} were derived from "
                f"{len(found)} source datasets and the projections from {len(wanted)}, "
                f"and {len(wanted - found)} of the projections' are missing; an "
                "artifact of another corpus names reactions this one does not hold"
            )
        for name, source in sources.items():
            # Stamps and columns say nothing about the IDs a file states, so an artifact
            # edited after it was written passes everything above. Offset onto the
            # corpus-wide range, an ID its dataset does not hold names another
            # dataset's structure or none, and the occurrence index's count of distinct
            # IDs cannot see it: the ID a row lost and the one it gained cancel. Bounded
            # per dataset rather than by the corpus total, since an ID just past its own
            # dataset lands inside the next one's range.
            largest = _max_structure_id(name, derive)
            count = self._structures_of_source[source]
            if largest is not None and largest >= count:
                raise PairingError(
                    f"{name} carries structure_id {largest} but its dataset holds only "
                    f"{count} structures, so its IDs would join to another dataset's "
                    f"molecules or to none; re-derive it with {derive} --force, since "
                    "a run without it skips a file whose stamps are current"
                )
        return sources

    def _resolve(self, name: str) -> str:
        """Returns a compound name's SMILES, calling the resolver only on a miss.

        The resolver is called outside the lock, so a slow service holds up only the
        search that asked. Two searches missing on one name at once both call it, and
        whichever stores last is kept.

        Args:
            name: The compound name a structure predicate gives.

        Returns:
            The canonical SMILES, or the resolver's answer unchanged if RDKit cannot
            parse it. That answer is left for the query molecule to refuse by name, and
            is not kept, so the next search naming the compound asks again.

        Raises:
            ValueError: If the default resolver cannot resolve the name; an injected
                resolver's own errors pass through unchanged. Nothing is kept.
        """
        with self._resolutions_lock:
            smiles = self._resolutions.get(name)
            if smiles is not None:
                self._resolutions.move_to_end(name)
                return smiles
        answered = self._resolver(name)
        smiles = _canonical(answered)
        if smiles is None:
            return answered
        logger.info("resolved %r to %r", name, smiles)
        with self._resolutions_lock:
            self._resolutions[name] = smiles
            self._resolutions.move_to_end(name)
            while len(self._resolutions) > _CACHED_RESOLUTIONS:
                self._resolutions.popitem(last=False)
        return smiles

    def _matches(
        self,
        cursor: duckdb.DuckDBPyConnection,
        parameter: query.StructureParameter,
        resolve: Callable[[str], str],
    ) -> str:
        """Returns a structure predicate's match set as a bitmap, reusing a recent one.

        The chemistry depends on the query molecule, the operation, and the threshold,
        so the same predicate asked twice has the same answer -- and asking it costs a
        pass over every molecule in the corpus. A compound is keyed by what it resolved
        to rather than by its name, so two names for one molecule share an entry and a
        resolver that starts answering differently is a different key.

        One question is answered once even while it is still being answered: a caller
        arriving during another thread's pass waits for that pass rather than starting
        a second, so a burst of identical requests costs one. The wait is per predicate,
        so unrelated searches still overlap. A pass that raises publishes nothing and a
        waiter behind it takes the matching over, which leaves the error with the caller
        that hit it rather than with everyone queued behind.

        Args:
            cursor: The cursor a similarity screen or a hash match runs on.
            parameter: The predicate to evaluate.
            resolve: Maps a compound name to SMILES.

        Returns:
            The bitmap over corpus-wide structure IDs.

        Raises:
            ValueError: If the predicate names neither a pattern nor a compound, if a
                resolved compound's SMILES does not parse, or if RDKit refuses to hash
                the query molecule of a predicate that matches on a hash.
            PairingError: If the library does not come out one entry per distinct
                molecule over an unbroken run of IDs.
        """
        if parameter.pattern is not None:
            pattern = parameter.pattern
        elif parameter.compound is not None:
            pattern = resolve(parameter.compound)
        else:
            raise ValueError(
                f"structure parameter {parameter.name!r} names neither a pattern nor a "
                "compound, which the grammar does not allow"
            )
        # Which parser reads the string is part of what the string means, so it is part
        # of the key -- and it is asked of the same function the query molecule comes
        # from, so the two cannot come apart. A name is keyed in its canonical spelling,
        # which a pattern can share: [NH4+] is ammonium's, and as SMARTS it excludes the
        # quaternary ammoniums the SMILES reading matches.
        key = (
            parameter.op,
            _reads_as_smarts(parameter),
            pattern,
            parameter.threshold,
            parameter.chirality,
        )
        while True:
            with self._matches_lock:
                cached = self._matched.get(key)
                if cached is not None:
                    # Reinserted so the eviction below drops what has gone longest
                    # unused rather than what was computed longest ago.
                    self._matched.move_to_end(key)
                    logger.info(
                        "%s %r reused a cached match set", parameter.op, pattern
                    )
                    return cached
                computing = self._matching.get(key)
                if computing is None:
                    self._matching[key] = threading.Event()
                    break
            logger.info("%s %r is already running; waiting", parameter.op, pattern)
            computing.wait()
        try:
            if parameter.op == "substructure":
                matched = self._substructure_ids(parameter, resolve)
            elif parameter.op in _HASH_MATCHES:
                matched = self._hashed_ids(cursor, parameter, resolve)
            else:
                matched = self._similarity_ids(cursor, parameter, resolve)
            logger.info(
                "%s %r matched %d of %d structures (%d unsearchable)",
                parameter.op,
                pattern,
                len(matched),
                self._searchable,
                self._total - self._searchable,
            )
            bitmap = self._bitmap(matched)
            with self._matches_lock:
                self._matched[key] = bitmap
                while len(self._matched) > _CACHED_MATCHES:
                    self._matched.popitem(last=False)
        finally:
            # Woken after the result is published, so a waiter finds it on the next turn
            # of the loop; a failed pass leaves nothing there and the waiter recomputes.
            with self._matches_lock:
                waiters = self._matching.pop(key, None)
            if waiters is not None:
                waiters.set()
        return bitmap

    def search(
        self, request: query.Query, *, timeout_seconds: float | None = None
    ) -> pa.Table:
        """Compiles and runs a query, returning the result as an Arrow table.

        Runs on its own cursor, so concurrent searches sharing this corpus do not read
        each other's results. Two searches wait on each other at the one-time library
        build and index check, and at a structure predicate one of them is already
        matching.

        A quantifier the occurrence index can answer becomes a condition on the row
        rather than a filter over the elements, and the rest of the query -- the clauses
        beside it, the aggregate, the ordering, the limit -- compiles exactly as it
        would have. The screening and verification are the compiler's either way, so an
        indexed query and a projection one differ only in how they reach the reactions
        holding a match. Which happened is logged, because it is the one thing about a
        search that the result does not show.

        Args:
            request: The query to run.
            timeout_seconds: Wall-clock bound on the whole call, not on the final
                query alone. The query is interrupted; a phase already running is not,
                and raises when it finishes instead -- so a search can take longer than
                this and still report the overrun. ``_Deadline`` says which phases those
                are and why. None runs unbounded.

        Returns:
            The selected columns: ``reaction_id`` for a plain query, and for an
            aggregated one a column per ``group_by`` path, named by the whole path,
            then one per measure, named by its ``name``. At most this corpus's
            ``max_rows``, where it has one: an answer that reached the bound is logged
            as possibly cut and carries ``TRUNCATED`` in its schema metadata.
            ``Compiled.limit`` is what actually reached the SQL, for a caller that
            wants the number rather than the flag. A query naming a compound carries
            ``RESOLVED`` as well, which ``resolutions`` reads.

        Raises:
            query.QueryError: If the query does not compile.
            ValueError: If a compound name cannot be resolved.
            PairingError: If the corpus IDs do not come out one unbroken run, or the
                occurrence index refuses the corpus -- a reaction stated by more than
                one row, or a structure no indexed path reaches.
            TimeoutError: If the search exceeds ``timeout_seconds``, naming the phase
                that was running when the bound passed.
        """

        def index(path: str, allocate: Callable[[], str]) -> str | None:
            # Asked of the corpus before the condition is written, since the condition
            # names the index's view and a corpus without one has nothing for it to
            # read -- and only for a path the index covers, so a query over any other
            # never pays the check that publishes it.
            if path not in INDEXED_PATHS or not self._occurrences():
                return None
            return _index_condition(path, allocate)

        # Kept per search as well as on the corpus, whose cache another search can evict
        # from meanwhile: a compound named in two structure predicates stays one
        # molecule, and the result's RESOLVED stamp states what this search used.
        resolved: dict[str, str] = {}

        def resolve(name: str) -> str:
            if name not in resolved:
                resolved[name] = self._resolve(name)
            return resolved[name]

        # Every view the compiler looked up, by level, so the log can name the ones the
        # SQL kept.
        looked_up: dict[str, str] = {}

        def pivot_view(path: str) -> str | None:
            name = self._pivot_view(path)
            if name is not None:
                looked_up[path] = name
            return name

        deadline = _Deadline(timeout_seconds)
        compiled = query.compile_query(
            request,
            index=index,
            pivot=pivot_view,
            history=query.HISTORY_TABLE if self._history else None,
            max_rows=self._max_rows,
        )
        # Compiling publishes the views it routes to on their first use, which reads
        # a footer per artifact and, for the index, scans the artifacts once.
        deadline.check("compiling the query")
        if compiled.limit != request.limit:
            logger.info("the corpus bounds this query at %d rows", compiled.limit)
        logger.info("this query reads %s", _relations_read(compiled.sql, looked_up))
        if _reads_history(compiled.sql):
            # The first commits are found by reaction ID, as a pivot's elements are.
            self._require_unique_reactions()
        cursor = self._connection.cursor()
        try:
            parameters: dict[str, Any] = dict(compiled.literals)
            for parameter in compiled.structures:
                parameters[parameter.name] = self._matches(cursor, parameter, resolve)
                # Per predicate rather than after all of them, so a query with two
                # does not pay for the second once the first has spent the budget.
                # A compound's resolution is checked here too: a resolver is an
                # external service, and one that hangs is the likeliest way a
                # search overruns without a single slow query in it.
                # Named by what the caller wrote rather than by parameter.name,
                # which is the compiler's own placeholder and identifies nothing
                # outside the SQL it was allocated for.
                asked = parameter.pattern or parameter.compound
                deadline.check(f"matching {parameter.op} {asked!r}")
            left = deadline.remaining()
            if left is None:
                answered = cursor.execute(compiled.sql, parameters).to_arrow_table()
            else:
                answered = _run_with_timeout(cursor, compiled.sql, parameters, left)
            stamps: dict[str, bytes] = {}
            if _warn_when_the_bound_was_reached(request, compiled, answered.num_rows):
                stamps[TRUNCATED] = b"true"
            if resolved:
                stamps[RESOLVED] = json.dumps(resolved).encode()
            if stamps:
                answered = answered.replace_schema_metadata(
                    (answered.schema.metadata or {}) | stamps
                )
            return answered
        finally:
            cursor.close()

    def _paired_sources(self) -> list[_Source]:
        """Returns the sources this corpus was opened with.

        Raises:
            ValueError: If it was opened without ``sources_pattern``.
        """
        if self._sources is None:
            raise ValueError(
                "this corpus was opened without sources_pattern, so it holds no "
                "reactions to fetch or datasets to describe; open it with a glob "
                "matching the source dataset files"
            )
        return self._sources

    def reactions(
        self, reaction_ids: Iterable[str] | pa.Array | pa.ChunkedArray
    ) -> pa.Table:
        """Returns the source records of the given reactions.

        What a search returns is reaction IDs; this is how a caller reaches the
        reactions themselves, read from the files the searched projections restate.
        One ID is cheap, since only its record is decoded. A batch reaches the Parquet
        scan as the range its IDs span, which for IDs spread across the corpus decodes
        nearly every record, and costs about the same whether it names ten IDs or a
        thousand.

        Runs on its own cursor, as ``search`` does, so concurrent fetches and searches
        sharing this corpus do not read each other's rows.

        Args:
            reaction_ids: The reactions wanted, typically a search's ``reaction_id``
                column, as the pyarrow column or as a list.

        Returns:
            Columns ``reaction_id``, ``dataset_id``, and ``reaction``, the serialized
            ``Reaction`` exactly as its source stores it. One row per distinct ID the
            corpus holds, in the order first asked; an ID it does not hold has no row,
            and an empty request returns an empty table with these columns.

        Raises:
            ValueError: If the corpus was opened without ``sources_pattern``, or the
                request names more distinct IDs than ``max_rows``.
            PairingError: If the sources state one of the reactions more than once, so
                which record answers would be arbitrary.
        """
        self._paired_sources()
        if isinstance(reaction_ids, (pa.Array, pa.ChunkedArray)):
            # DuckDB binds Python strings, not the pyarrow scalars iterating one yields.
            reaction_ids = reaction_ids.to_pylist()
        asked = list(dict.fromkeys(reaction_ids))
        if self._max_rows is not None and len(asked) > self._max_rows:
            raise ValueError(
                f"asked for {len(asked)} reactions, more than the {self._max_rows} "
                "this corpus returns at once; ask in batches of at most max_rows, or "
                "open the corpus with a larger one"
            )
        if not asked:
            return _RECORD_SCHEMA.empty_table()
        cursor = self._connection.cursor()
        try:
            found = cursor.execute(
                f"SELECT reaction_id, {query.DATASET_ID}, reaction "  # noqa: S608
                "FROM corpus_sources WHERE reaction_id IN (SELECT unnest($asked))",
                {"asked": asked},
            ).to_arrow_table()
        finally:
            cursor.close()
        found = found.cast(_RECORD_SCHEMA)
        returned = found.column("reaction_id").to_pylist()
        if len(set(returned)) != len(returned):
            repeated, copies = collections.Counter(returned).most_common(1)[0]
            raise PairingError(
                f"reaction {repeated} is stated {copies} times across the sources, so "
                "which record answers for it would be arbitrary"
            )
        position = {reaction_id: index for index, reaction_id in enumerate(asked)}
        # Typed, because a plain empty list (nothing asked for is held) reaches pyarrow
        # as a null array, which take has no kernel for.
        order = pa.array(
            sorted(range(len(returned)), key=lambda row: position[returned[row]]),
            type=pa.int64(),
        )
        return found.take(order)

    def datasets(self) -> pa.Table:
        """Returns what each source dataset says about itself.

        Read from the footers at open, so this reads no file.

        Returns:
            One row per source: ``dataset_id``, ``name``, ``description``, and
            ``reactions``, the row count its footer states. Ordered by
            ``dataset_id``, with any source that recorded none last.

        Raises:
            ValueError: If the corpus was opened without ``sources_pattern``.
        """
        ordered = sorted(
            self._paired_sources(),
            key=lambda source: (
                source.dataset_id is None,
                source.dataset_id or "",
                source.path,
            ),
        )
        return pa.Table.from_pylist(
            [
                {
                    query.DATASET_ID: source.dataset_id,
                    "name": source.name,
                    "description": source.description,
                    "reactions": source.reactions,
                }
                for source in ordered
            ],
            schema=_DATASETS_SCHEMA,
        )
