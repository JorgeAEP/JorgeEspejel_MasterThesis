# RAVEN curation audit guide

The authoritative full curation audit is `curation_tracker.xlsx`.

Its main evidence tables are:

- `Removed Reactions`: reaction deletions, with pre-removal equations,
  bounds, gene associations, evidence and validation state.
- `Removed Metabolites`: metabolite removals and associated target decisions.
- `Directionality Changes`: bound edits with original and revised values.
- `Retained and Protected`: reviewed branches not removed, with rationale.
- `Validation History` and `Checkpoint Comparison`: model-level checks for
  each curation transition.

The terminal Step-88 change is additionally supplied in machine-readable CSV
form: `removed_objects.csv`, `component_summary.csv`,
`biomass_precursor_checks.csv`, and `validation_summary.csv`.
