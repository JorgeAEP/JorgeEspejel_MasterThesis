# RAVEN *Bipolaris sorokiniana* reconstruction

`original/` contains the imported KEGG/RAVEN-track SBML files. The
COBRApy-compatible file is the starting reconstruction retained for this
archive; the gap-filled file is preserved as the initial working variant.

`curated/88_post_remove_disconnected_bh4_pair.xml` is the final independent
RAVEN-track model used for the phosphate-gradient community analysis. The
terminal curation event removed the disconnected residual BH4 reaction
`R11762` and metabolites `C04244` and `C21640`. Its implementation, checks
and removed-object audit are retained in `scripts/` and `reports/`.

The preceding multi-step manual curation history is recorded in
`reports/curation_tracker.xlsx`; this is the authoritative reaction-level
tracker for RAVEN edits. It includes removals, repairs, evidence decisions and
their validation states. `reports/validation_summary.csv` records the final
Step-88 checks.

The model has cytosolic and extracellular compartments and uses a
KEGG-oriented namespace. Its different compartmental structure and exchange
representation from the CarveFungi models are important limits to a direct
one-to-one reconstruction comparison.
