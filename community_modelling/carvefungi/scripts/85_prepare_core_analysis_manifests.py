#!/usr/bin/env python3
"""Script 85: prepare manifests for the core community-comparison panel.

Creates missing manifests for:

- fungal monocultures:
    fungus_alone__raven_independent
    fungus_alone__carvefungi
- bacterial SynCom baselines:
    hvsc1_bacteria_only
    hvsc2_bacteria_only

The four fungus+SynCom manifests are expected to already exist:

- hvsc1_plus_bipolaris__raven_independent
- hvsc2_plus_bipolaris__raven_independent
- hvsc1_plus_bipolaris__carvefungi
- hvsc2_plus_bipolaris__carvefungi
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
MANIFEST_DIR = PACKAGE / "manifests"
OUT = PACKAGE / "reports/85_prepare_core_analysis_manifests"

RAVEN_MODEL = PACKAGE / "fungal_models/raven_independent/88_post_remove_disconnected_bh4_pair.xml"
CARVEFUNGI_MODEL = PACKAGE / "fungal_models/carvefungi_shared_interface/Bipolaris_sorokiniana_carvefungi_shared_medium_interface.xml"

HVSC1_COMMUNITY_MANIFEST = MANIFEST_DIR / "hvsc1_plus_bipolaris__raven_independent.csv"
HVSC2_COMMUNITY_MANIFEST = MANIFEST_DIR / "hvsc2_plus_bipolaris__raven_independent.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def fungal_manifest(sample_id: str, model_id: str, model_path: Path, source: str) -> list[dict[str, object]]:
    return [
        {
            "sample_id": sample_id,
            "id": model_id,
            "file": str(model_path),
            "abundance": 1.0,
            "source": source,
            "kind": "fungus",
        }
    ]


def bacterial_only_manifest(source_manifest: Path, sample_id: str) -> list[dict[str, object]]:
    table = pd.read_csv(source_manifest)
    bacteria = table.loc[table["kind"].astype(str).str.lower() == "bacterium"].copy()
    if bacteria.empty:
        raise ValueError(f"No bacterial rows found in {source_manifest}")
    bacteria["sample_id"] = sample_id
    return bacteria.to_dict(orient="records")


def inspect_manifest(path: Path, role: str) -> dict[str, object]:
    table = pd.read_csv(path)
    return {
        "manifest": str(path),
        "role": role,
        "rows": int(table.shape[0]),
        "fungus_rows": int((table["kind"].astype(str).str.lower() == "fungus").sum()),
        "bacterium_rows": int((table["kind"].astype(str).str.lower() == "bacterium").sum()),
        "sha256": sha256(path),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for required in [RAVEN_MODEL, CARVEFUNGI_MODEL, HVSC1_COMMUNITY_MANIFEST, HVSC2_COMMUNITY_MANIFEST]:
        if not required.exists():
            raise FileNotFoundError(required)

    generated = {
        "fungus_alone__raven_independent": fungal_manifest(
            "fungus_alone__raven_independent",
            "Bipolaris_sorokiniana_raven_independent",
            RAVEN_MODEL,
            "Bipolaris_raven_independent",
        ),
        "fungus_alone__carvefungi": fungal_manifest(
            "fungus_alone__carvefungi",
            "Bipolaris_sorokiniana_carvefungi",
            CARVEFUNGI_MODEL,
            "Bipolaris_carvefungi_shared_medium_interface",
        ),
        "hvsc1_bacteria_only": bacterial_only_manifest(HVSC1_COMMUNITY_MANIFEST, "hvsc1_bacteria_only"),
        "hvsc2_bacteria_only": bacterial_only_manifest(HVSC2_COMMUNITY_MANIFEST, "hvsc2_bacteria_only"),
    }

    audit_rows = []
    for name, rows in generated.items():
        path = MANIFEST_DIR / f"{name}.csv"
        write_manifest(rows, path)
        audit_rows.append(inspect_manifest(path, "generated_core_manifest"))

    for name in [
        "hvsc1_plus_bipolaris__raven_independent",
        "hvsc2_plus_bipolaris__raven_independent",
        "hvsc1_plus_bipolaris__carvefungi",
        "hvsc2_plus_bipolaris__carvefungi",
    ]:
        path = MANIFEST_DIR / f"{name}.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        audit_rows.append(inspect_manifest(path, "existing_core_manifest"))

    audit = pd.DataFrame(audit_rows)
    audit.to_csv(OUT / "core_manifest_audit.csv", index=False)

    summary = {
        "package": str(PACKAGE),
        "manifest_dir": str(MANIFEST_DIR),
        "generated_manifests": sorted(generated),
        "core_manifest_count": int(len(audit_rows)),
        "outputs": str(OUT),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 85 compact summary")
    print(f"generated_manifests: {';'.join(sorted(generated))}")
    print(f"core_manifest_count: {len(audit_rows)}")
    print(f"outputs: {OUT}")


if __name__ == "__main__":
    main()
