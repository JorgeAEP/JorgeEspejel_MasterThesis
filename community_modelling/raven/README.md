# RAVEN fungal--SynCom community analysis

This directory contains the inputs and script for the RAVEN
*B. sorokiniana* phosphate-gradient analysis. The final curated fungal model
is duplicated here as an analysis input; its authoritative model lineage and
curation history remain under `models/raven/`.

`scripts/116_run_raven_phosphate_gradient_micom.py` is the MICOM workflow.
The manifests make the fungal monoculture, HvSC1 + fungus and HvSC2 + fungus
layouts explicit. `media/sensitivity_panel/` contains the phosphate-bound
conditions, while `media/phaseA_kegg_full_medium.csv` is the relaxed medium.

`harmonisation/` documents the bridge between the fungal KEGG-oriented
namespace and the bacterial community interfaces. Archived result records are
in `reports/`. No script in this directory is run during archive preparation.
