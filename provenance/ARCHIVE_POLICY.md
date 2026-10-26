# Archive policy

This repository is a frozen reproducibility package, not a live analysis
workspace. Archive assembly did not rerun reconstruction, curation, FBA or
MICOM analyses. It copied the final model artifacts, source inputs, scripts,
configuration files, mappings, validation records and selected output tables
from the existing project workspaces.

The archive intentionally excludes external CarveFungi, eggNOG-mapper and
PSIPRED software trees and reference databases. Their size and third-party
licensing make them unsuitable for inclusion. The required interfaces and
recorded environment are documented under `environment/` and in the archived
construction scripts.

`MANIFEST.tsv` records each archived file and its size. The SHA256 checksum
file covers every archived file except the manifest and checksum file
themselves, avoiding self-referential checksums.
