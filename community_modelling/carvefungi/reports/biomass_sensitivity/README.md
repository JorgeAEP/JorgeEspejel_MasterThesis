# Biomass × tradeoff 2×2 summary

This folder compares native-open fungal exchange simulations under two biomass formulations and two MICOM community-growth fractions.

- Biomass formulations: shared 52-component biomass vs species-specific biomass.
- Fractions: 1.0 strict maximum-community-growth solution vs 0.9 relaxed tradeoff.
- Medium: proline-relaxed final shared barley medium (`final_shared_proline_1.csv`).

The purpose is to separate biomass-formulation effects from optimization/tradeoff effects.

Files:
- `biomass_tradeoff_2x2_growth_long.csv`: one row per manifest/fraction.
- `biomass_tradeoff_2x2_growth_pivot.csv`: compact growth table.
- `biomass_tradeoff_2x2_delta_table.csv`: same compact table with 0.9 minus 1.0 deltas.
