#!/usr/bin/env python3
"""Verify archive structure without running model construction or MICOM."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]

REQUIRED = [
    "README.md",
    "environment/requirements.txt",
    "environment/software_versions.txt",
    "models/raven/original/Bipolaris_sorokiniana_COBRApy_compatible.xml",
    "models/raven/curated/88_post_remove_disconnected_bh4_pair.xml",
    "models/raven/scripts/308_apply_step88_remove_disconnected_bh4_pair.py",
    "models/carvefungi/bipolaris/original/Bipolaris_sorokiniana_carvefungi.xml",
    "models/carvefungi/bipolaris/curated/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.xml",
    "models/carvefungi/fusarium/original/Fusarium_graminearum_carvefungi.xml",
    "models/carvefungi/fusarium/curated/Fusarium_graminearum_carvefungi_after_stoichiometry_correction.xml",
    "community_modelling/raven/scripts/116_run_raven_phosphate_gradient_micom.py",
    "community_modelling/carvefungi/scripts/96_run_fungal_species_syncom_factorial_panel.py",
    "community_modelling/carvefungi/media/final_shared_barley_community_medium_dual_namespace.csv",
    "community_modelling/carvefungi/harmonisation/fungal_curated76_interface_changes.csv",
    "provenance/model_lineage.tsv",
    "provenance/script_registry.tsv",
    "provenance/analysis_registry.tsv",
]

def main() -> int:
    missing = [path for path in REQUIRED if not (ROOT / path).is_file()]

    if missing:
        print("Missing required files:")
        print("\n".join(f"- {path}" for path in missing))
    if missing:
        return 1

    bacterial_models = list((ROOT / "community_modelling/carvefungi/models/bacterial_models").glob("HvSC*/*.xml"))
    if len(bacterial_models) != 53:
        print(f"Expected 53 bacterial SBML models, found {len(bacterial_models)}.")
        return 1

    print("Archive structure verified.")
    print("Bacterial SBML models: 53 (HvSC1: 27; HvSC2: 26).")
    print("No modelling analyses were run by this verification.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
