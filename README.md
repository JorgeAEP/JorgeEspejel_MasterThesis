# JorgeEspejel_MasterThesis

This private repository is the reproducibility archive for the computational
work supporting My Thesis on fungal--
bacterial metabolic interactions in barley root synthetic communities.

It preserves the terminal model lineage, the scripts used to construct and
curate the fungal models, the media and identifier mappings used for joint
community modelling, and the scripts and selected outputs supporting the
MICOM analyses. No analysis is run by this archive; the included scripts are
provided to document and, where dependencies are available, reproduce the
workflow.

## Repository map

- `models/raven/`: original and final curated *Bipolaris sorokiniana* RAVEN
  reconstruction, its final curation script, and curation/validation records.
- `models/carvefungi/`: input proteomes, native reconstructions, corrected
  derived models, construction scripts, curation scripts, and validation
  records for *B. sorokiniana* and *Fusarium graminearum*.
- `community_modelling/raven/`: the curated RAVEN model, phosphate-gradient
  medium, manifests, MICOM script, mappings, and archived result records.
- `community_modelling/carvefungi/`: bacterial models, native and
  community-ready fungal copies, biomass definitions, dual-namespace media,
  manifests, MICOM scripts, identifier mappings, and selected result tables.
- `environment/`: the recorded Python dependency specifications and software
  versions.
- `provenance/`: model and analysis registries, a complete archive manifest,
  checksums, and the verification script.

## Model lineage

The archive distinguishes three model states:

1. **Original/native models** are the first exported reconstruction files.
2. **Curated models** are derived SBML/JSON files after the documented
   curation or stoichiometric corrections.
3. **Community-ready models** are derived copies that expose an integration
   interface required for MICOM. They do not replace the native or curated
   fungal models and do not modify the internal fungal metabolic networks.

See `provenance/model_lineage.tsv` for the source-to-archive mapping.

## Running order

The scripts are retained in the order used in the original workspaces. A
future reproduction should proceed as follows:

1. Install the software described in `environment/` and the external
   CarveFungi/eggNOG/PSIPRED resources described in the per-model READMEs.
2. Reconstruct each fungal model from its input proteome using the archived
   CarveFungi workflow script.
3. Apply the documented final curation script for the relevant fungal model.
4. Use the final curated model, biomass definition, medium, manifest and
   identifier crosswalk from the appropriate community-modelling folder.
5. Run the corresponding MICOM script only in a separate output directory.

The construction scripts are preserved as executed, including their original
workspace assumptions. Their required source files, command logic and output
paths are documented in the model README files. They should be adapted to a
new local installation rather than run against the archived model files in
place.

## Analyses and outputs

The archive contains selected final tables and figures plus run records for
the thesis analyses. `provenance/analysis_registry.tsv` states which script,
models, medium, biomass and manifests belong to each analysis, and whether its
outputs are archived here.

## Verification

Run the structural archive check without rerunning any modelling analysis:

```bash
python3 provenance/verify_archive.py
```

`provenance/checksums.sha256` can be checked with:

```bash
sha256sum -c provenance/checksums.sha256
```

## Scope and interpretation

The archive documents model-derived hypotheses and reproducibility inputs. It
does not make the MICOM exchange candidates direct evidence of metabolite
transfer, cross-feeding, pathogen suppression, or in-planta community
stability.
