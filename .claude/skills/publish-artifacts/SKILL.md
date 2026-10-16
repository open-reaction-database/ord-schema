---
name: publish-artifacts
description: Use when publishing derived artifacts (projections, structures, pivots, occurrences) to the ord-artifacts Hugging Face repository, adding a new ord-data dataset to the published artifacts, re-deriving after an RDKit upgrade or an ARTIFACT_SCHEMA bump, running ord_schema.scripts.upload_artifacts, or updating the ord-artifacts dataset card.
---

# Publishing derived artifacts

`upload_artifacts` publishes whatever complete build it is given at the repository root,
replacing the build there; earlier builds stay in its history. It refuses a tree that is
not one consistent build, and refuses to delete published files without
`--allow_deletes`, so a tree missing datasets the repository holds is refused. A tree
holding only some datasets is still a consistent build, though, and the deletes cannot
see a dataset the repository does not hold yet, nor anything on the first publish:
`fetch_artifacts` would serve that partial build to every reader.
`scripts/check_tree.py` closes that gap: **run it before every upload.**

## Setup

Run everything with the repository's `.venv`, synced with the `search` and `huggingface`
extras. Its RDKit is the one the tree is stamped with; under another (the conda `ord`
env, for one) every file reads stale and a derive rebuilds the whole corpus.

```sh
ORD_DATA=~/ord/ord-data
ARTIFACTS=~/ord/artifacts
```

## Steps

1. **Update the source.** `git -C $ORD_DATA pull`, then `git -C $ORD_DATA lfs pull`:
   `data/` is in LFS, and a pointer file is not a dataset.

2. **Derive into the full tree,** in this order. Each step skips files whose stamps are
   current, so a new dataset costs only its own files. The USPTO projection alone takes
   about 35 minutes, and an RDKit upgrade or a schema bump makes every file stale.

   ```sh
   python -m ord_schema.artifacts.scripts.derive_projection \
       --input_pattern="$ORD_DATA/data/*/*.parquet" --output_dir=$ARTIFACTS/projections
   python -m ord_schema.artifacts.scripts.derive_structures \
       --input_pattern="$ARTIFACTS/projections/*/*.parquet" --output_dir=$ARTIFACTS/structures
   python -m ord_schema.artifacts.scripts.derive_pivots \
       --input_pattern="$ARTIFACTS/projections/*/*.parquet" --output_dir=$ARTIFACTS/pivots
   python -m ord_schema.artifacts.scripts.derive_occurrences \
       --pivots_dir=$ARTIFACTS/pivots --output_dir=$ARTIFACTS/occurrences
   ```

3. **Check the tree against ord-data.** Needs the network for the mirror listing:

   ```sh
   python .claude/skills/publish-artifacts/scripts/check_tree.py \
       --root=$ARTIFACTS --ord_data=$ORD_DATA
   ```

   It requires every Parquet dataset in the checkout to have exactly one projection
   stamped with its current content hash, and the checkout to hold the same datasets as
   the ord-data mirror on `main`, byte for byte. A mismatch is fixed by deriving or by
   updating the checkout, never by uploading anyway.

4. **Dry run.**

   ```sh
   python -m ord_schema.scripts.upload_artifacts --root=$ARTIFACTS --dry_run \
       --message="Derive from ord-data $(git -C $ORD_DATA rev-parse --short HEAD)"
   ```

   Read the artifact schema and RDKit it reports against the build it replaces. The same
   pair updates that build; a new one switches every reader at once, and readers left on
   the old pair have to pass an earlier revision, so confirm the switch with the user.
   The count of files to upload is the whole tree; unchanged files drop out only inside
   the real commit.

5. **Upload,** with the same command minus `--dry_run`, once the user has confirmed: it
   publishes. It needs `hf auth login` or `$HF_TOKEN` with write access, and logs the
   commit URL.

6. **Push the dataset card** on the first publish, and whenever
   [`dataset_card.md`](dataset_card.md) changes. It is the repository's `README.md`,
   which `upload_artifacts` never touches. Change it here, through a PR, so the Hub's
   README is always a copy of this file; a schema bump that changes the columns the
   `projections` config names needs its `dataset_info` changed too. Push it once the
   user has confirmed: it publishes too.

   ```sh
   hf upload open-reaction-database/ord-artifacts \
       .claude/skills/publish-artifacts/dataset_card.md README.md --repo-type=dataset \
       --commit-message="Update the dataset card"
   ```

   Then check that the [dataset viewer](https://huggingface.co/datasets/open-reaction-database/ord-artifacts/viewer)
   shows rows for both configs. From the shell, each config should print a row count
   rather than an error:

   ```sh
   for config in projections structures; do
       curl -s "https://datasets-server.huggingface.co/first-rows?dataset=open-reaction-database/ord-artifacts&config=$config&split=train" \
           | jq '.error // (.rows | length)'
   done
   ```

   The viewer reprocesses the repository after each commit, so give it a few minutes
   before reading an error as final. `datasets` has no type for the projection's MAP
   columns, so the `projections` config selects three columns and `dataset_info`
   declares their features. Keep both when editing the card: `columns` alone fails,
   because features are inferred from the full schema first. If `projections` still
   errors, cut the card to the `structures` config alone: drop the `projections` config
   and its `dataset_info`, make `structures` the default, rewrite the viewer sentence,
   land that here, and push again.

## When the upload refuses to delete

The refusal names a count per artifact directory. Rerun the dry run with `--allow_deletes` to list
the files, and pass `--allow_deletes` only when each one is meant to go: a dataset
removed from ord-data, or a pivot level or indexed path the library no longer derives.
Files of datasets that ord-data still holds mean the tree is short, which step 3 should
already have caught.

## Common mistakes

| Mistake | What happens |
| --- | --- |
| Deriving a new dataset into an empty directory and uploading that | A one-dataset build: refused once the repository holds other datasets; on the first publish it goes up, and every reader fetches one dataset. Step 3 catches both. |
| Deriving under the conda `ord` env | Every file reads stale; the corpus re-derives under the older RDKit, and uploading it replaces the published build with one that readers on the newer RDKit refuse. |
| Passing `--allow_deletes` to get past a refusal | Deletes published datasets from the build readers fetch. |
| Reading the dry run's upload count as the size of the change | It is every file in the tree; a new dataset still reads as all of them. |
| Editing the README on the Hub | The next push of `dataset_card.md` overwrites the edit. |
