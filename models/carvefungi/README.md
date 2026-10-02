# CarveFungi-derived fungal reconstructions

This branch contains separately reconstructed *Bipolaris sorokiniana* and
*Fusarium graminearum* models. Each organism directory distinguishes the
protein and annotation inputs, unmodified native reconstruction, corrected
derived model, construction scripts, and curation validation records.

The reconstruction scripts use the CarveFungi workflow with eggNOG-based
functional annotation and predicted subcellular localisation. They preserve
the original local execution setup and require third-party CarveFungi,
eggNOG-mapper and PSIPRED installations. Those large external tools and their
databases are not included in this archive.

The corrected model files are not overwritten native reconstructions. They are
derived copies created after the corresponding curation validation. The
reaction-level before/after audit and final validation records are included
with each organism.
