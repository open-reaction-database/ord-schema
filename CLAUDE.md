# ord-schema

Protocol buffer schema for the Open Reaction Database and the `ord_schema` package that builds, validates, serializes, and queries it. Module map: [ord_schema/README.md](ord_schema/README.md); `artifacts/`, `search/`, and `orm/` each have a README.

## Commands

```bash
uv sync --extra tests  # development environment in .venv
uv run pytest -n auto  # tests are foo_test.py beside foo.py
uv run ruff check --fix . && uv run ruff format .
uv run ty check ord_schema  # same scope as CI Lint
./compile_proto_wrappers.sh  # after editing proto/ord-schema/proto/*.proto; needs Node
uv lock  # after changing dependencies in pyproject.toml
```

`pre-commit install` runs ruff, ty, addlicense, clang-format on `.proto`, markdownlint-cli2, and nbstripout on each commit, and refuses commits to `main`. It does not run the tests.

## Gotchas

- ORM tests start PostgreSQL through `testing.postgresql`, so `initdb` must be on `PATH`; without it they fail with `command not found: initdb`. It must come from a Postgres with the RDKit cartridge, since `rdkit_mappers_test.py` has no skip for its absence. CI installs conda-forge `rdkit-postgresql` (`.github/conda-postgres.yml`).
- Live PubChem/OPSIN resolver tests run by default; set `ORD_LIVE_RESOLVERS=false` to skip them offline.
- `ord_schema/proto/*_pb2.py` and `*_pb2.pyi` are generated: regenerate them, never hand-edit, and commit them with the license header the addlicense hook adds. CI fails on drift and runs `buf breaking` against `main` and the last release.
- Three submission notebooks (listed in `.github/workflows/run_tests.yml`) compare stored output under nbval; after changing what one prints, run `uv run jupyter execute --inplace <notebook>`. Running notebooks locally rewrites the tracked `examples/tools/set_dative_bonds/output.csv`; revert it unless updating it is the point.
- Changing what an artifact holds, or the RDKit it is derived under, fails `_GOLDEN_OUTPUT` in `ord_schema/artifacts/base_test.py`, and the failure prints the line to re-pin. `ARTIFACT_SCHEMA_VERSION` moves only when a reader of the old schema would misread new files.
- The package version moves only through `publish.yml` (bump-my-version), never by hand.
- ord-interface, ord-app, and ord-data import from `ord_schema`; grep their `origin/main` before removing a public symbol.
