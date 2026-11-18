# Derived artifacts

A source dataset in [ord-data](https://github.com/open-reaction-database/ord-data) stores
reactions as opaque wire-format bytes, so any question about them costs a full
deserialization. A **derived artifact** restates those same reactions in a shape some
consumer can read cheaply. It adds no information, which is the whole discipline here: an
artifact is *wrong* if it disagrees with its source, and *stale* if its source has moved
on. Both have to be detectable without reading a single column of data. One artifact
reads something besides its source: `arrivals`, from ord-data's git history, which
carries a stamp of that history too; see [Arrivals](#arrivals).

Five artifacts, derived in a chain:

```mermaid
flowchart LR
    D["source dataset<br/>ord-data, wire-format bytes"] --> P["<b>projection</b><br/>every field, as nested Parquet"]
    P --> S["<b>structures</b><br/>one row per distinct molecule<br/>fingerprints + serialized mol"]
    P --> F["<b>arrivals</b><br/>one row per reaction<br/>its first commit to ord-data"]
    H["ord-data's git history<br/>index of first commits"] --> F
    P --> V["<b>pivot</b><br/>one row per element of<br/>one repeated level"]
    V --> O["<b>occurrences</b><br/>one row per structure occurrence<br/>at one indexed path"]
    S -.->|"paired by source dataset"| C["ord_schema.search.Corpus"]
    P -.-> C
    F -.->|"read only by a query<br/>naming first_commit"| C
    V -.-> C
    O -.-> C
```

| artifact | one file per | holds | read by |
| --- | --- | --- | --- |
| [`projection`](projection.py) | source dataset | every field of every message reachable from `Reaction`, as `STRUCT` / `LIST` / `MAP` | any query engine; the compiler in [`ord_schema.search`](../search/) |
| [`structures`](structures.py) | projection | one row per distinct structure: `pattern_fp`, `morgan_fp`, `morgan_popcount`, `mol_binary`, keyed by `structure_id` | the screen and verify steps of substructure search |
| [`arrivals`](arrivals.py) | projection | one row per reaction: its first commit to ord-data's `main`, that commit's time and offset, and the dataset it arrived in | a query naming `first_commit`, and no other |
| [`pivot`](pivot.py) | projection × repeated level | one row per element, carrying `reaction_id`, the ordinal of every enclosing level, and the element's own non-repeated fields | `exists` and `forall` over the level, nested ones joined on the ordinals, rather than a filter over the projection's lists; per-reaction reductions such as a ranking by yield; `aggregate.over`, which reads nothing else; and the occurrences derivation |
| [`occurrences`](occurrences.py) | pivot × indexed path | one row per structure occurrence: `reaction_id`, `structure_id` | the semi-join a structure query with no other element condition is routed to, as a view over the artifacts |

The projection leaves **no field out** — a field left out is a question nobody can ask.
Its schema is generated from the proto descriptors, so a field added upstream appears
without anyone deciding it is worth carrying. Two normalizations are applied, both
because the un-normalized form is a trap: united `{value, precision, units}` triples
become columns named for their unit (`setpoint_kelvin`), so `WHERE temperature > 350`
cannot silently miss the rows recorded in Celsius; and the structural identifiers
collapse to one `smiles`, so "what is this molecule" has one answer. Everything else
stays in the source `reaction` column, which remains authoritative for byte-exact
round-tripping.

## Stamps

Every artifact carries four footer keys in the `ord.` namespace, and
[`base.is_current`](base.py) requires all of them to match before a file counts as
current:

| key | why it can invalidate |
| --- | --- |
| `ord.artifact` | tells one artifact from another without inspecting the schema |
| `ord.artifact_schema_version` | the version of the artifact's schema: its columns and what each holds; one for every artifact, since they are published and fetched as one build |
| `ord.source_md5` | the content of the source dataset this restates |
| `ord.rdkit_version` | RDKit releases have changed both canonicalization and pattern fingerprints, either of which silently breaks a file whose consumers run the newer one |

No key names the ord-schema release. A release need not change any artifact, and a
development build reports the last release's number, so a release stamp says nothing
about currency while making a rebuild of an unchanged artifact write different bytes.

A change a reader of the current schema would misread — a column renamed, removed, or
given a different meaning — bumps `ARTIFACT_SCHEMA_VERSION`, which makes every artifact stale. A
change a reader of either schema reads correctly does not: a column added is seen by
`missing_columns` and derives again on its own, and values improved by a fixed
normalization are derived again with `--force` and published in place. `base_test.py`
pins what each artifact holds when derived from two fixtures, so a change to an artifact
fails there until someone decides which kind it is.

A fifth, `ord.source_dataset_id`, is written only when the source records one, and is the
only key not required to read stamps back. An `arrivals` artifact also carries
`ord.arrivals_version` and `ord.arrivals_head`, naming the index of ord-data's
history it read; see [Arrivals](#arrivals).

Rebuilding an artifact from the same source writes the same bytes, under any ord-schema
release that leaves its definition alone, until pyarrow or DuckDB changes: every Parquet
footer names the pyarrow that wrote it, and DuckDB's batches are a pivot's row groups.
Pivots are sorted on `reaction_id` and their ordinals, and occurrences on every column,
since the unnest that produces them runs without holding any order.

`ord.artifact` carries a second job the other keys do not: nothing outside this library
writes it, so its presence alone says a file is ours to replace. That is what
`derive_tree` asks before writing over anything, rather than asking for the whole stamp
set — a file written before a key was added or renamed reads as no artifact at all, and a
guard wanting all five would refuse the very re-derive that brings such a tree up to date.

An artifact derived from another artifact **passes `source_md5` through** rather than
hashing its parent, so every artifact names the dataset it reflects however many
derivations away it sits, and one comparison answers "is this current for that dataset?"
A chain is invisible to a consumer checking currency.

## Deriving a tree

[`base.derive_tree`](base.py) is the driver every script uses. It walks a glob and, for
each match:

- mirrors the input's directory layout beneath `output_dir`, rooted at the leading
  components of the pattern that hold no wildcard;
- **skips** anything whose footer already records the current source, schema, and RDKit, so
  re-running is cheap and an interrupted run resumes rather than restarts;
- **ignores** a match that is not the kind of parent this derivation reads — a source
  dataset where a projection was wanted, or an artifact from an earlier run — so
  `output_dir` may sit inside a recursive pattern's reach;
- **refuses** to write over any of its inputs, checked across every destination against
  every source rather than each against its own, because an `output_dir` nested in the
  source tree maps one dataset onto a *different* one, which a per-source check passes
  while destroying a source the run has not reached yet.

Five scripts wrap it, in dependency order:

```bash
python -m ord_schema.artifacts.scripts.derive_projection \
    --input_pattern='data/*/*.parquet' --output_dir=projections

python -m ord_schema.artifacts.scripts.derive_structures \
    --input_pattern='projections/*/*.parquet' --output_dir=structures

python -m ord_schema.artifacts.scripts.derive_arrivals \
    --input_pattern='projections/*/*.parquet' --output_dir=arrivals \
    --index=arrivals_index.parquet

python -m ord_schema.artifacts.scripts.derive_pivots \
    --input_pattern='projections/*/*.parquet' --output_dir=pivots \
    --levels $(python -m ord_schema.artifacts.scripts.derive_occurrences --print_levels)

python -m ord_schema.artifacts.scripts.derive_occurrences \
    --pivots_dir=pivots --output_dir=occurrences
```

Each takes `--force` to rewrite what is already current. `derive_pivots` defaults to all
39 repeated levels, which is far more than a deployment answers questions over; naming
the levels it does is the cheaper run — and `derive_occurrences --print_levels` names the
four it reads, which is why the command above asks it rather than spelling them.
`derive_occurrences` covers every path a structure can sit at, which is not a choice: a
path left out is one whose structures no query can find. `derive_arrivals` reads the
index `index_arrivals` walks out of ord-data's history; see
[Arrivals](#arrivals).

## Arrivals

A record says when its depositor created it, not when it entered ORD's public corpus. The
submission pipeline's own `record_modified` event can precede the merge by days, so the
answer comes from ord-data's git history instead. It is not in the record, so it is not in
the projection either, which restates the record: the
[`arrivals`](arrivals.py) artifact carries it, one row per reaction, beside the
projection it derives from. A query names it as `first_commit`:

| field | holds |
| --- | --- |
| `sha` | the first commit on ord-data's `main`, followed by first parent, whose version of a dataset file under `data/` held the reaction's ID |
| `timestamp` | that commit's committer date as git recorded its wall clock, naive, like every projected time |
| `utc_offset_minutes` | the offset git recorded with it, in minutes east of UTC |
| `dataset_id` | the dataset whose file held it then, which need not be the one holding it now |

First parent, because early submissions merged in two steps and a reaction entered the
corpus when the second reached `main`. Every version of every file, because a file's own
first commit is not the answer: every dataset has been rewritten into another format at
least once, and USPTO has been resharded more than once. Per reaction, because a dataset
updated later holds reactions that arrived at different commits.

A reaction whose ID a migration reassigned first appears under its current ID at that
migration, since nothing in either file links the two IDs. ord-data holds two:

- The Nanomole flow dataset, merged in ord-data#41, was deleted and re-added under new
  IDs as `ord_dataset-68cb8b4b2b384e3d85b5b1efae58b203` in #73.
- The USPTO grants merged in #66 through #69 were removed in #149 and re-added under new
  IDs by the reshard in #151, #153, and #155. Their first commits are December 2022.

[`index_arrivals`](scripts/index_arrivals.py) walks the history
([`history`](history.py)) into an index holding every reaction ID `main` has held,
current or not, with the chain of commits it walked. History only grows, so an existing
index is extended from its head and only the commits since are read; `--rebuild` starts
from the root, which an index walked under another `history.VERSION`, or along a chain
the branch no longer holds, needs. Each version is read for its reaction and dataset IDs
alone, through a two-field message that parses any version of a `Dataset`, so text
written under a retired schema still reads.

The walk needs the whole history as committed: it refuses a shallow clone, whose oldest
commit would read as where every reaction it holds arrived, and a grafts file, and it
ignores replace refs. It reads only LFS objects already fetched. The Hugging Face mirror
ord-data's `.lfsconfig` points at holds no object older than October 2024, nor every one
since, while GitHub's LFS store holds them all, so the first walk on a machine fetches
from there once:

```bash
git -C ord-data -c lfs.url=https://github.com/open-reaction-database/ord-data.git/info/lfs \
    lfs fetch origin $(git -C ord-data rev-list --first-parent main)
```

[`derive_arrivals`](scripts/derive_arrivals.py) looks up each projection's
reactions in the index and stamps the artifact with the index's walk version and head. An
artifact stays current while its stamped head is on the chain the index walked: an index
only appends, so the first commit it gives a reaction never changes, and a new ord-data
commit leaves every artifact current. One stamped with a head the chain does not hold
derives again, in seconds for the whole corpus, and no projection is touched. A projection
holding a reaction the index lacks refuses, since a missing first commit would read as a
reaction ord-data never held.

A `Corpus` given `arrivals_pattern` refuses unless one artifact pairs with every
projection and the artifacts name only reactions the corpus holds, with no field left
NULL. It reads them only for a query naming `first_commit`; every other query reads what
it would without them. A query reading nothing but `first_commit` and `reaction_id` reads
the artifacts alone. Otherwise the compiler turns the filters a conjunction places on
`first_commit` into one semi-join, and grouping, measuring, or ordering by it into a join
on `reaction_id`; see [the search README](../search/README.md#arrivals).

Measured on 2026-10-20 over ord-data at `7048706`, 135 first-parent commits: the full
walk read 1,264 LFS objects (7.0 GB) and 285 MB of blobs stored in git in 46 seconds, at a
5.5 GB peak resident set, the largest single `.pb.gz` alone reaching 2.6 GB. The index
holds 4,205,068 reaction IDs, 2,428,291 of them current, in 71 MB. The artifacts for the
93-dataset corpus hold 41 MiB and derive in 10 seconds.

## Published builds

Derived trees are published to the Hugging Face dataset
[ord-artifacts](https://huggingface.co/datasets/open-reaction-database/ord-artifacts),
apart from the ord-data mirror so that rebuilding them leaves that history alone. The
repository holds one build at a time, laid out exactly as a local tree is:

```text
projections/aa/ord_dataset-<id>.parquet
pivots/inputs.components/aa/ord_dataset-<id>.parquet
```

An earlier build, of an older artifact schema version or another RDKit, stays in the
repository's history, and each upload's commit description names the pair it holds.
[`huggingface.fetch_artifacts`](../huggingface.py) reads the build's stamps in place
before downloading anything and refuses a schema version or RDKit this library does not
read, so a local tree holding a build it reads is not overwritten with one it cannot; a
reader left on the old pair passes the revision of a commit holding it.
A fetch deletes local files the build has since dropped and refreshes any pivot level
already held. The RDKit version has to match exactly: pattern fingerprints and
registration hashes both drift between releases, patch releases included, and a
structure query computes its side of each comparison with the RDKit it runs.

[`upload_artifacts`](../scripts/upload_artifacts.py) publishes a local tree in one
commit, replacing the build the repository held, and reports the schema version and RDKit of the
one it replaces: a change switches every reader at once. It accepts only a tree a
`Corpus(require_pivots=True, require_occurrences=True)` opens: one artifact schema version
and RDKit version, every pivot level and indexed path present and filed where its footer
says, and each restating the projections' source datasets once. Unchanged files drop out
of the commit, so a new dataset goes up by deriving it into the full tree and uploading
that. Published files the tree lacks are deleted only with `--allow_deletes`; otherwise
the upload is refused, because a tree holding only the new dataset is a complete build of
it, and uploading it would delete the rest. Its docstring lists the refusals; it needs
the `search` and `huggingface` extras.

```bash
python -m ord_schema.scripts.upload_artifacts \
    --root=artifacts --message="Derive from ord-data <commit>" --dry_run
```

A value-only change is published in place, and a reader's next fetch brings its copy up
to date. A change a reader of the current schema would misread needs the
`ARTIFACT_SCHEMA_VERSION` bump, so that readers still on the old schema are refused rather than
served files they would misread.

## What pairs with what

`structure_id` is assigned by `write_projection` in first-seen order and is **local to
one dataset**. The projection and its structures artifact are two statements of one
derivation and are meaningful only together:

- IDs are stable for a fixed source and a fixed library — assignment is first-seen order
  over the reactions, so the same bytes under the same ord-schema and RDKit give the same
  numbers, which is what lets an artifact pair with a projection by source hash rather
  than by which build produced it. They are *not* stable across a library upgrade that
  changes canonicalization, which is what the version stamps are for. Nothing outside the
  artifacts should record them either way; the column is marked internal and stays out of
  what a model is told it may query.
- A projection rewritten in place needs its structures artifact rederived with it. The
  stamps name the *source dataset* rather than the projection file, so the skip check
  cannot see that the pairing changed — which is why
  [`Corpus`](../search/execute.py) additionally refuses a structures artifact too short
  to hold every `structure_id` its projection carries.
- Pairing is by source dataset rather than by filename, so the two trees need not share a
  layout.

A pivot names its reactions by ID, so one derived from another corpus resolves against
rows this one does not hold — confidently and wrongly. `Corpus` refuses pivot artifacts
whose source datasets are not its projections'. It also refuses a pivot or occurrences
artifact carrying a `structure_id` its own dataset's structures artifact does not hold:
an artifact edited after it was written keeps its stamps, and an ID moved out of its
dataset names another dataset's molecule or none, which the occurrence index's count of
distinct IDs cannot see. The bound is read from footer statistics, so a projection or
artifact recording none for its `structure_id` columns is refused too.

**No artifact carries a corpus-wide ID.** A `structure_id` is local to its dataset; the
corpus-wide one a bitmap is indexed by is that plus an offset, and the offset is a running
total over the datasets the corpus opened. It therefore belongs to no artifact: the same
rows read beside a different set of datasets need a different offset, and one baked into a
file would stay *in range* while naming another dataset's molecule — passing on the whole
corpus and failing only on a subset of it, which is the shape of bug that reaches
production. The corpus assigns offsets at open and joins them on, keyed by the source
dataset every artifact names. That rule is what lets the occurrence index be an artifact
at all rather than a relation rebuilt at every startup: given an `occurrences_dir`
covering every indexed path, [`Corpus`](../search/execute.py) publishes it as a view over
those files and builds nothing. Measured over the full corpus that is 0.18 s and 0.40 GiB
resident at open. A directory covering some paths and not others leaves the corpus with
no index, since the check that makes one trustworthy counts the structures every path
reaches; the structure predicates are then answered from the pivots or the projection.

## Querying a projection

Use list lambdas over repeated levels, never `UNNEST` in a `FROM` clause. Measured in
DuckDB it is 27–200× slower where it finishes at all, and over the full corpus 0.9
seconds against no result in four minutes. It also emits one row per matching identifier
rather than one per reaction, so the two spellings need a `DISTINCT` to even count alike:

```sql
-- fast
SELECT reaction_id FROM p
WHERE len(list_filter(
        flatten(list_transform(map_values(inputs), i -> i.components)),
        c -> len(list_filter(c.identifiers,
                             x -> x.type = 'NAME' AND x.value = 'THF')) > 0)) > 0

-- same answer, does not finish
SELECT DISTINCT reaction_id
FROM p, UNNEST(map_values(inputs)) t(i), UNNEST(i.components) u(c),
        UNNEST(c.identifiers) v(x)
WHERE x.type = 'NAME' AND x.value = 'THF'
```

Writing that by hand is what [`ord_schema.search`](../search/) exists to avoid: it
compiles a validated `Query` to the fast spelling, and routes quantifiers to the pivots
and to a structure index rather than to the projection where it can.
