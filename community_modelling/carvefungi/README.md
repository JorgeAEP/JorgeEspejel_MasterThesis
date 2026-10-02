# CarveFungi fungal--SynCom community analysis

This directory supports the main fungal--bacterial MICOM comparisons and the
biomass/medium/optimisation sensitivity analyses.

`models/native/` contains the curated fungal model before interface
harmonisation. `models/community_ready/` contains the derived fungal copies
used in MICOM. They expose exchange identifiers compatible with the
dual-namespace medium but retain the fungal internal networks. The 27 HvSC1
and 26 HvSC2 bacterial SBML models are under `models/bacterial_models/`.

`manifests/` contains the nine shared-biomass layouts and the corresponding
species-specific-biomass layouts. `media/` contains the common dual-namespace
apoplastic medium, its source audit, and the additional support media used by
the sensitivity analyses. `biomass/` retains the shared baseline and the
species-specific biomass construction scripts.

`harmonisation/` contains the fungal interface changes, external-identifier
crosswalks, medium-mapping provenance and exchange-coverage audit. The
scripts cover community preparation, MICOM runs, growth/precursor analyses,
niche overlap and exchange/cross-feeding screens. The former resilience script
for Section 3.3.3 is intentionally absent.
