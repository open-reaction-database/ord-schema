# Review report demo

The review report bundle for the
[Nielsen deoxyfluorination screen](../../submissions/2_Nielsen_Deoxyfluorination_Screen/),
as `review_dataset` writes it:

- [`Deoxyfluorination_screen_review.html`](Deoxyfluorination_screen_review.html): the
  page a reviewer reads. Download it and open it in a browser. It is self-contained and
  follows the browser's light or dark theme.
- [`tables/`](tables/): every table behind the page, complete:
  - `findings.csv`: dataset-level and validation findings;
  - `compounds.csv`: every compound;
  - `flat.csv`: the varying fields, one row per reaction;
  - `constants.csv`: the fields every reaction holds the same way;
  - `column_stats.csv`: statistics for each varying field;
  - `exemplars.txtpb`: the reactions picked for the reviewer to read, as a dataset
    whose description says why each was picked. It can be imported into the ORD
    reaction editor with *Create Dataset from File*; the page also has a download
    for each exemplar, to import with *Create Reaction from File*.

The screen is recorded consistently, so the dataset-level checks find nothing. The one
validation warning, raised by all 80 reactions, is a record time without a UTC offset.
The exemplars are one reaction for its single shape, and the reactions with the lowest
and highest yield.

## Rebuilding it

From the repository root:

```bash
uv run python -m ord_schema.scripts.review_dataset \
    --input=examples/submissions/2_Nielsen_Deoxyfluorination_Screen/nielsen_fig1_dataset.pbtxt \
    --output_dir=examples/demo/review_report
```

The page records when it was generated and with which ord-schema, so rebuilding it
changes those lines even when nothing else has changed.
