# Fungal-species × SynCom factorial MICOM panel

This folder contains the main community modelling panel for the CarveMe fungal models.

Main ecological medium:

- `original_proline_0p1`: original shared barley/community medium with `EX_proL_m = 0.1`.

Targeted control:

- `proline_relaxed_1`: Fusarium-specific monoculture control with `EX_proL_m = 1.0`.

The main design compares:

- Bipolaris alone
- Fusarium alone
- Bipolaris + Fusarium
- HvSC1 alone
- HvSC2 alone
- each single fungus + each SynCom
- both fungi together + each SynCom

Each fungal context is tested with the shared 52-component biomass and the species-specific biomass.

Runs summarized: 9

Important summary tables are in `summaries/`.
