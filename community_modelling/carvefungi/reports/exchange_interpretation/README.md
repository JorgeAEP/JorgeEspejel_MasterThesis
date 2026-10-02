# Exchange interpretation tables

This folder reorganizes the large MICOM exchange-flux output into smaller,
question-oriented tables. No MICOM simulations were run by this step.

## Generated tables

- `01_exchange_flux_long_nonzero_combined.csv`
  - Combined non-zero exchange fluxes from Script 96 and Script 100.
  - One row = one organism × exchange × run.

- `02_exchange_summary_by_run_exchange.csv`
  - One row = one run × exchange.
  - Aggregates medium flux, fungal uptake/secretion, bacterial uptake/secretion,
    producer/consumer counts, and candidate interaction flags.

- `03_medium_uptake_by_run.csv`
  - Medium components consumed from the community medium.

- `04_fungal_exchange_fluxes_by_run.csv`
  - Fungal uptake and secretion flux rows.

- `05_bacterial_exchange_fluxes_by_run.csv`
  - Bacterial uptake and secretion flux rows.

- `06_candidate_bacteria_to_fungus_crossfeeding.csv`
  - Exchanges where at least one bacterium secretes and at least one fungus
    consumes the same exchange.

- `07_candidate_fungus_to_bacteria_crossfeeding.csv`
  - Exchanges where at least one fungus secretes and at least one bacterium
    consumes the same exchange.

- `08_resource_competition_candidates.csv`
  - Exchanges where fungi and bacteria both consume the same exchange.

- `09_hvsc1_vs_hvsc2_exchange_differences.csv`
  - Exchange-level HvSC1 minus HvSC2 differences for matched contexts.

- `10_top_exchange_summary_by_run.csv`
  - Top 50 largest exchange signals per run.

## Row counts

```text
combined_long_rows: 61837
summary_run_exchange_rows: 9435
medium_uptake_rows: 1052
fungal_exchange_rows: 30
bacterial_exchange_rows: 22721
bacteria_to_fungus_candidate_rows: 1
fungus_to_bacteria_candidate_rows: 4
resource_competition_rows: 1
hvsc1_vs_hvsc2_rows: 3925
top_exchange_rows: 649
```

## How to interpret flux signs

For organism rows:

- positive flux = secretion/export;
- negative flux = uptake/import.

For medium rows:

- negative medium flux indicates medium uptake by the community;
- positive medium flux indicates net secretion to the medium.

Candidate cross-feeding and competition rows are hypotheses from flux patterns,
not proof of direct transfer. They should be used to prioritize follow-up
inspection or targeted blocking/supplementation tests.
