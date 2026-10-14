---
name: artifact-publish
description: Use when publishing derived artifacts (projections, structures, pivots, occurrences) to the ord-artifacts Hugging Face repository, adding a new ord-data dataset to the published artifacts, re-deriving after an RDKit upgrade or an ARTIFACT_VERSIONS bump, or running ord_schema.scripts.upload_artifacts.
---

# Publishing derived artifacts

`upload_artifacts` publishes whatever complete build it is given. It refuses a tree that
is not one consistent build, and refuses to delete published files without
`--allow_deletes`. A tree holding only some datasets is still a consistent build, though.
Under a new lineage or RDKit prefix it deletes nothing, so neither check sees it, and
`fetch_artifacts` would serve that partial build to every reader on that version.
`scripts/check_tree.py` closes that gap: **run it before every upload.**

## Setup

Run everything with the repository's `.venv`, synced with the `search` and `huggingface`
extras. Its RDKit is the one the tree is stamped with; under another (the conda `ord`
env, for one) every file reads stale and a derive rebuilds the whole corpus under a new
prefix.

```sh
ORD_DATA=~/ord/ord-data
ARTIFACTS=~/ord/artifacts
```

## Steps

1. **Update the source.** `git -C $ORD_DATA pull`, then `git -C $ORD_DATA lfs pull`:
   `data/` is in LFS, and a pointer file is not a dataset.

2. **Derive into the full tree,** in this order. Each step skips files whose stamps are
   current, so a new dataset costs only its own files. The USPTO projection alone takes
   about 35 minutes, and an RDKit upgrade or a version bump makes every file stale.

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
   python .claude/skills/artifact-publish/scripts/check_tree.py \
       --root=$ARTIFACTS --ord_data=$ORD_DATA
   ```

   It requires every Parquet dataset in the checkout to have exactly one projection
   stamped with its current content hash, and the checkout to hold the same dataset IDs
   as the ord-data mirror on `main`. A mismatch is fixed by deriving or by updating the
   checkout, never by uploading anyway.

4. **Dry run.**

   ```sh
   python -m ord_schema.scripts.upload_artifacts --root=$ARTIFACTS --dry_run \
       --message="Derive from ord-data $(git -C $ORD_DATA rev-parse --short HEAD)"
   ```

   Read the four prefixes: the same as the last publish means the upload updates that
   build in place, a new one means a new build. The count of files to upload is the
   whole tree; unchanged files drop out only inside the real commit.

5. **Upload,** with the same command minus `--dry_run`, once the user has confirmed: it
   publishes. It needs `hf auth login` or `$HF_TOKEN` with write access, and logs the
   commit URL.

## When the upload refuses to delete

The refusal names a count per prefix. Rerun the dry run with `--allow_deletes` to list
the files, and pass `--allow_deletes` only when each one is meant to go: a dataset
removed from ord-data, or a pivot level or indexed path the library no longer derives.
Files of datasets that ord-data still holds mean the tree is short, which step 3 should
already have caught.

## Common mistakes

| Mistake | What happens |
| --- | --- |
| Deriving a new dataset into an empty directory and uploading that | A one-dataset build: refused under the published prefix, published as the whole corpus under a new one. Step 3 catches both. |
| Deriving under the conda `ord` env | Every file reads stale; the corpus re-derives under the older RDKit and a new prefix. |
| Passing `--allow_deletes` to get past a refusal | Deletes published datasets from the build readers fetch. |
| Reading the dry run's upload count as the size of the change | It is every file in the tree; a new dataset still reads as all of them. |
