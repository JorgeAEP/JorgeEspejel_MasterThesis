#!/usr/bin/env python3
"""Prepare dual-fungi biomass-variant MICOM panel.

Creates derived model files and manifests for community modelling with:

* Bipolaris CarveFungi with current 52-component biomass
* Bipolaris CarveFungi with Alternaria-derived biomass
* Fusarium CarveFungi with current 52-component biomass
* Fusarium CarveFungi with finalized Fusarium fujikuroi-derived biomass

The panel includes:

* bacterial-only HvSC1/HvSC2 controls
* individual fungal controls
* Bipolaris + Fusarium fungal-pair controls
* Bipolaris + Fusarium + HvSC1
* Bipolaris + Fusarium + HvSC2

Two media are prepared:

* final_shared_proline_1: official final shared barley/community medium
* final_shared_proline_0p1: same composition, proline-limited control

No source model is overwritten.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import shutil
from pathlib import Path

import cobra
import pandas as pd


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
WORKSPACE = PACKAGE.parent

OUT = PACKAGE / "reports/94_prepare_dual_fungi_biomass_variant_panel"
MODEL_OUT = PACKAGE / "fungal_models/biomass_variant_panel"
MEDIA_OUT = PACKAGE / "media/biomass_variant_panel"
MANIFEST_DIR = PACKAGE / "manifests"

FINAL_MEDIUM = (
    PACKAGE
    / "media/final_shared_barley_medium/final_shared_barley_community_medium_dual_namespace.csv"
)

BIPOLARIS_BASE = (
    PACKAGE
    / "fungal_models/carvefungi_shared_interface/"
    "Bipolaris_sorokiniana_carvefungi_shared_medium_interface.xml"
)
FUSARIUM_BASE = (
    WORKSPACE
    / "Fusarium_graminearum_carveme/04_models/"
    "Fusarium_graminearum_carvefungi_after_stoichiometry_correction.xml"
)

BIPOLARIS_SPECIFIC_BIOMASS = (
    PACKAGE
    / "biomass/bipolaris_alternaria_harmonized/"
    "biomass_bipolaris_alternaria_harmonized.py"
)
FUSARIUM_SPECIFIC_BIOMASS = (
    PACKAGE
    / "biomass/fusarium_graminearum_f_fujikuroi_harmonized/"
    "biomass_fusarium_fujikuroi_harmonized.py"
)

HVSC1_BACTERIA_ONLY = MANIFEST_DIR / "hvsc1_bacteria_only.csv"
HVSC2_BACTERIA_ONLY = MANIFEST_DIR / "hvsc2_bacteria_only.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolve_metabolite(model: cobra.Model, met_id: str) -> cobra.Metabolite:
    candidates = [met_id]
    if met_id.startswith("M_"):
        candidates.append(met_id[2:])
    else:
        candidates.append("M_" + met_id)
    for candidate in candidates:
        try:
            return model.metabolites.get_by_id(candidate)
        except KeyError:
            pass
    raise KeyError(f"Missing metabolite: {met_id}")


def add_specific_biomass(model: cobra.Model, module_path: Path) -> str:
    module = load_module(module_path, f"biomass_{module_path.parent.name}")
    biomass_id = module.BIOMASS_ID
    old_biomass_like = [
        reaction
        for reaction in model.reactions
        if (
            "BIOMASS" in reaction.id.upper()
            or "BIOMASS" in str(reaction.name).upper()
            or "GROWTH" in str(reaction.name).upper()
        )
    ]
    if old_biomass_like:
        model.remove_reactions(old_biomass_like)
    reaction = cobra.Reaction(biomass_id)
    reaction.name = getattr(module, "BIOMASS_NAME", biomass_id)
    reaction.lower_bound = 0.0
    reaction.upper_bound = 1000.0
    reaction.add_metabolites(
        {
            resolve_metabolite(model, met_id): coeff
            for met_id, coeff in module.STOICHIOMETRY.items()
        }
    )
    reaction.annotation["sbo"] = "SBO:0000629"
    model.add_reactions([reaction])
    model.objective = reaction
    return biomass_id


def save_model(model: cobra.Model, xml_path: Path, json_path: Path) -> None:
    xml_path.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(xml_path))
    cobra.io.save_json_model(model, str(json_path))


def create_model_variants() -> dict[str, dict[str, object]]:
    MODEL_OUT.mkdir(parents=True, exist_ok=True)

    variants = {}

    # Bipolaris current 52-component biomass.
    bip_52 = cobra.io.read_sbml_model(str(BIPOLARIS_BASE))
    variants["bipolaris_carvefungi_52"] = {
        "model": bip_52,
        "xml": MODEL_OUT / "Bipolaris_sorokiniana_carvefungi_52biomass.xml",
        "json": MODEL_OUT / "Bipolaris_sorokiniana_carvefungi_52biomass.json",
        "model_id": "Bipolaris_sorokiniana_carvefungi_52",
        "species": "Bipolaris_sorokiniana",
        "source": "Bipolaris_carvefungi_52biomass",
        "biomass": "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF",
    }

    # Bipolaris Alternaria-derived biomass.
    bip_specific = cobra.io.read_sbml_model(str(BIPOLARIS_BASE))
    bip_specific_biomass = add_specific_biomass(bip_specific, BIPOLARIS_SPECIFIC_BIOMASS)
    variants["bipolaris_carvefungi_specific"] = {
        "model": bip_specific,
        "xml": MODEL_OUT / "Bipolaris_sorokiniana_carvefungi_alternaria_biomass.xml",
        "json": MODEL_OUT / "Bipolaris_sorokiniana_carvefungi_alternaria_biomass.json",
        "model_id": "Bipolaris_sorokiniana_carvefungi_alternaria_specific",
        "species": "Bipolaris_sorokiniana",
        "source": "Bipolaris_carvefungi_Alternaria_biomass",
        "biomass": bip_specific_biomass,
    }

    # Fusarium current 52-component biomass.
    fus_52 = cobra.io.read_sbml_model(str(FUSARIUM_BASE))
    variants["fusarium_carvefungi_52"] = {
        "model": fus_52,
        "xml": MODEL_OUT / "Fusarium_graminearum_carvefungi_52biomass.xml",
        "json": MODEL_OUT / "Fusarium_graminearum_carvefungi_52biomass.json",
        "model_id": "Fusarium_graminearum_carvefungi_52",
        "species": "Fusarium_graminearum",
        "source": "Fusarium_carvefungi_52biomass",
        "biomass": "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF",
    }

    # Fusarium fujikuroi-derived biomass.
    fus_specific = cobra.io.read_sbml_model(str(FUSARIUM_BASE))
    fus_specific_biomass = add_specific_biomass(fus_specific, FUSARIUM_SPECIFIC_BIOMASS)
    variants["fusarium_carvefungi_specific"] = {
        "model": fus_specific,
        "xml": MODEL_OUT / "Fusarium_graminearum_carvefungi_fujikuroi_biomass.xml",
        "json": MODEL_OUT / "Fusarium_graminearum_carvefungi_fujikuroi_biomass.json",
        "model_id": "Fusarium_graminearum_carvefungi_fujikuroi_specific",
        "species": "Fusarium_graminearum",
        "source": "Fusarium_carvefungi_Fujikuroi_biomass",
        "biomass": fus_specific_biomass,
    }

    for variant in variants.values():
        save_model(variant["model"], variant["xml"], variant["json"])
        del variant["model"]
        variant["xml_sha256"] = sha256(variant["xml"])
        variant["json_sha256"] = sha256(variant["json"])

    return variants


def read_rows(path: Path) -> list[dict[str, object]]:
    return pd.read_csv(path).to_dict(orient="records")


def fungal_row(sample_id: str, variant: dict[str, object]) -> dict[str, object]:
    return {
        "sample_id": sample_id,
        "id": variant["model_id"],
        "file": str(variant["xml"]),
        "abundance": 1.0,
        "source": variant["source"],
        "kind": "fungus",
    }


def with_sample(rows: list[dict[str, object]], sample_id: str) -> list[dict[str, object]]:
    out = []
    for row in rows:
        new = dict(row)
        new["sample_id"] = sample_id
        out.append(new)
    return out


def write_manifest(name: str, rows: list[dict[str, object]]) -> Path:
    path = MANIFEST_DIR / f"{name}.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def create_manifests(variants: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    hvsc1_rows = read_rows(HVSC1_BACTERIA_ONLY)
    hvsc2_rows = read_rows(HVSC2_BACTERIA_ONLY)

    panel: dict[str, list[dict[str, object]]] = {}

    # Individual fungal controls.
    for key in [
        "bipolaris_carvefungi_52",
        "bipolaris_carvefungi_specific",
        "fusarium_carvefungi_52",
        "fusarium_carvefungi_specific",
    ]:
        name = f"fungus_alone__{key}"
        panel[name] = [fungal_row(name, variants[key])]

    # Dual fungal controls.
    pair_52 = "fungi_pair__carvefungi_52"
    panel[pair_52] = [
        fungal_row(pair_52, variants["bipolaris_carvefungi_52"]),
        fungal_row(pair_52, variants["fusarium_carvefungi_52"]),
    ]
    pair_specific = "fungi_pair__carvefungi_specific"
    panel[pair_specific] = [
        fungal_row(pair_specific, variants["bipolaris_carvefungi_specific"]),
        fungal_row(pair_specific, variants["fusarium_carvefungi_specific"]),
    ]

    # Single fungus + one SynCom at a time.
    for syncom, bacteria_rows in [("hvsc1", hvsc1_rows), ("hvsc2", hvsc2_rows)]:
        for fungus_key in [
            "bipolaris_carvefungi_52",
            "bipolaris_carvefungi_specific",
            "fusarium_carvefungi_52",
            "fusarium_carvefungi_specific",
        ]:
            name = f"{syncom}_plus_{fungus_key}"
            panel[name] = with_sample(bacteria_rows, name) + [fungal_row(name, variants[fungus_key])]

    # Dual fungi + one SynCom at a time.
    for syncom, bacteria_rows in [("hvsc1", hvsc1_rows), ("hvsc2", hvsc2_rows)]:
        name_52 = f"{syncom}_plus_fungi_pair__carvefungi_52"
        panel[name_52] = (
            with_sample(bacteria_rows, name_52)
            + [
                fungal_row(name_52, variants["bipolaris_carvefungi_52"]),
                fungal_row(name_52, variants["fusarium_carvefungi_52"]),
            ]
        )

        name_specific = f"{syncom}_plus_fungi_pair__carvefungi_specific"
        panel[name_specific] = (
            with_sample(bacteria_rows, name_specific)
            + [
                fungal_row(name_specific, variants["bipolaris_carvefungi_specific"]),
                fungal_row(name_specific, variants["fusarium_carvefungi_specific"]),
            ]
        )

    # Bacteria-only controls under the same media.
    panel["hvsc1_bacteria_only"] = with_sample(hvsc1_rows, "hvsc1_bacteria_only")
    panel["hvsc2_bacteria_only"] = with_sample(hvsc2_rows, "hvsc2_bacteria_only")

    audit = []
    for name, rows in panel.items():
        path = write_manifest(name, rows)
        table = pd.DataFrame(rows)
        audit.append(
            {
                "manifest": name,
                "path": str(path),
                "rows": int(table.shape[0]),
                "fungus_rows": int((table["kind"].astype(str).str.lower() == "fungus").sum()),
                "bacterium_rows": int((table["kind"].astype(str).str.lower() == "bacterium").sum()),
                "sha256": sha256(path),
            }
        )
    return audit


def update_medium_proline(path: Path, out_path: Path, proline_flux: float) -> None:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["reaction"] == "EX_proL_m":
                row["flux"] = str(proline_flux)
            rows.append(row)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reaction", "flux"])
        writer.writeheader()
        writer.writerows({"reaction": row["reaction"], "flux": row["flux"]} for row in rows)


def create_media() -> list[dict[str, object]]:
    MEDIA_OUT.mkdir(parents=True, exist_ok=True)
    if not FINAL_MEDIUM.exists():
        raise FileNotFoundError(FINAL_MEDIUM)

    proline_1 = MEDIA_OUT / "final_shared_proline_1.csv"
    proline_0p1 = MEDIA_OUT / "final_shared_proline_0p1.csv"
    shutil.copy2(FINAL_MEDIUM, proline_1)
    update_medium_proline(FINAL_MEDIUM, proline_0p1, 0.1)

    readme = MEDIA_OUT / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Dual-fungi biomass-variant media",
                "",
                "`final_shared_proline_1.csv` is the official final shared barley/community medium.",
                "",
                "`final_shared_proline_0p1.csv` is the proline-limited control. It has the same components but sets:",
                "",
                "```text",
                "EX_proL_m = 0.1",
                "```",
                "",
                "Use the proline-limited control to test whether bacteria can rescue Fusarium-specific biomass growth by cross-feeding proline or related metabolites.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    return [
        {
            "scenario": "final_shared_proline_1",
            "path": str(proline_1),
            "rows": sum(1 for _ in proline_1.open()) - 1,
            "proline_flux": 1.0,
        },
        {
            "scenario": "final_shared_proline_0p1",
            "path": str(proline_0p1),
            "rows": sum(1 for _ in proline_0p1.open()) - 1,
            "proline_flux": 0.1,
        },
    ]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    required = [
        FINAL_MEDIUM,
        BIPOLARIS_BASE,
        FUSARIUM_BASE,
        BIPOLARIS_SPECIFIC_BIOMASS,
        FUSARIUM_SPECIFIC_BIOMASS,
        HVSC1_BACTERIA_ONLY,
        HVSC2_BACTERIA_ONLY,
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

    source_checksums_before = {str(path): sha256(path) for path in [BIPOLARIS_BASE, FUSARIUM_BASE]}
    variants = create_model_variants()
    manifest_audit = create_manifests(variants)
    media_audit = create_media()
    source_checksums_after = {str(path): sha256(path) for path in [BIPOLARIS_BASE, FUSARIUM_BASE]}

    pd.DataFrame(variants).T.reset_index(names="variant").to_csv(OUT / "model_variant_audit.csv", index=False)
    pd.DataFrame(manifest_audit).to_csv(OUT / "manifest_audit.csv", index=False)
    pd.DataFrame(media_audit).to_csv(OUT / "media_audit.csv", index=False)

    manifest_names = [row["manifest"] for row in manifest_audit]
    summary = {
        "model_variants_folder": str(MODEL_OUT),
        "media_folder": str(MEDIA_OUT),
        "manifest_count": len(manifest_audit),
        "manifests": sorted(manifest_names),
        "media": media_audit,
        "source_models_unchanged": source_checksums_before == source_checksums_after,
        "outputs": str(OUT),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 94 compact summary")
    print(f"model_variants: {len(variants)}")
    print(f"manifest_count: {len(manifest_audit)}")
    print(f"media_scenarios: {';'.join(row['scenario'] for row in media_audit)}")
    print(f"source_models_unchanged: {summary['source_models_unchanged']}")
    print(f"outputs: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
