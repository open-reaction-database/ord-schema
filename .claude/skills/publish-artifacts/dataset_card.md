---
license: cc-by-sa-4.0
pretty_name: Open Reaction Database derived artifacts
tags:
- chemistry
- reactions
- cheminformatics
configs:
- config_name: projections
  default: true
  data_files: "projections/projection-1/rdkit-2026.03.5/**/*.parquet"
  columns:
  - reaction_id
  - smiles
  - identifiers
- config_name: structures
  data_files: "structures/structures-1_projection-1/rdkit-2026.03.5/**/*.parquet"
dataset_info:
- config_name: projections
  features:
  - name: reaction_id
    dtype: string
  - name: smiles
    dtype: string
  - name: identifiers
    list:
    - name: type
      dtype: string
    - name: details
      dtype: string
    - name: value
      dtype: string
    - name: is_mapped
      dtype: bool
---

# Open Reaction Database derived artifacts

Each dataset in the [Open Reaction Database](https://huggingface.co/datasets/open-reaction-database/ord-data) restated in shapes a query engine reads cheaply. The artifacts are derived by [`ord_schema.artifacts`](https://github.com/open-reaction-database/ord-schema/tree/main/ord_schema/artifacts) and searched by [`ord_schema.search`](https://github.com/open-reaction-database/ord-schema/tree/main/ord_schema/search). They add nothing to the source datasets, which stay the authoritative record.

## Artifacts

| artifact | one file per | holds |
| --- | --- | --- |
| `projections` | source dataset | every field of every reaction as nested Parquet, with units folded into column names (`setpoint_kelvin`) and each molecule's structural identifiers resolved to one `smiles` |
| `structures` | source dataset | one row per distinct molecule: SMILES, registration hashes, pattern and Morgan fingerprints, and a serialized RDKit mol, keyed by `structure_id` |
| `pivots` | source dataset and repeated level | one row per element of that level, with its `reaction_id`, the ordinal of every enclosing level, and the element's own non-repeated fields |
| `occurrences` | source dataset and indexed path | one row per molecule occurrence: `reaction_id`, `structure_id`, `reaction_role` |

A `structure_id` is local to its dataset: it joins a projection, pivot, or occurrences file only to the structures file derived from the same source dataset.

## Layout

```text
projections/projection-1/rdkit-2026.03.5/<shard>/ord_dataset-<id>.parquet
structures/structures-1_projection-1/rdkit-2026.03.5/<shard>/ord_dataset-<id>.parquet
pivots/pivot-1_projection-1/rdkit-2026.03.5/<level>/<shard>/ord_dataset-<id>.parquet
occurrences/occurrences-1_pivot-1_projection-1/rdkit-2026.03.5/<path>/<shard>/ord_dataset-<id>.parquet
```

- The shard and file name are the source dataset's own, as in ord-data's `data/<shard>/ord_dataset-<id>.parquet`.
- Each artifact sits under a build prefix, `<artifact>/<lineage>/rdkit-<version>`. The lineage names the version of the artifact's definition and of each artifact it derives from; the RDKit version is the release it was derived under.
- A new definition or a new RDKit release publishes a new build beside the existing ones. A reader needs the build matching its own RDKit exactly, since fingerprints and registration hashes change between releases.
- Every file's Parquet footer names the artifact (`ord.artifact`), its lineage (`ord.artifact_lineage`), the RDKit version (`ord.rdkit_version`), and the MD5 of the source dataset it restates (`ord.source_md5`).

The dataset viewer shows the reaction ID, reaction SMILES, and identifiers of the projections, and the structures in full, for the `rdkit-2026.03.5` build.

## Usage

`fetch_artifacts` downloads the build that the installed ord-schema and RDKit read, and a repeat call fetches only what changed:

```bash
pip install 'ord-schema[search,huggingface]'
```

```python
from ord_schema import huggingface
from ord_schema.search import execute

artifacts = huggingface.fetch_artifacts("ord-artifacts")
corpus = execute.Corpus(
    artifacts.projection_pattern,
    artifacts.structures_pattern,
    pivots_dir=str(artifacts.pivots),
    occurrences_dir=str(artifacts.occurrences),
)
```

Pass `revision=` a commit to pin one. Each upload's commit message names the ord-data revision it was derived from. See the [search README](https://github.com/open-reaction-database/ord-schema/blob/main/ord_schema/search/README.md) for queries, and the [artifacts README](https://github.com/open-reaction-database/ord-schema/blob/main/ord_schema/artifacts/README.md) for how each artifact is derived and stamped.

## Citation

Please cite the Open Reaction Database as given on the [ord-data card](https://huggingface.co/datasets/open-reaction-database/ord-data).
