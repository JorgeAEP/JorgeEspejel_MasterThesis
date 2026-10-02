#!/usr/bin/env python3
"""Summarise fungal--SynCom niche overlap and exchange-candidate metrics.

This script uses the all-exchange identifier audit rather than the earlier
EX-only screen. It reports fungal--bacterial co-uptake (niche overlap) and
non-medium shared-pool-compatible bacteria-to-fungus candidates separately.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
AUDIT = PACKAGE / "reports/133_fungal_bacterial_exchange_identifier_audit"
CLASSIFICATION = PACKAGE / "reports/134_exchange_candidate_classification"
OUT = PACKAGE / "reports/135_fungal_syncom_interface_metrics"
THRESHOLD = 1e-4


def label_fungus(organism: str) -> str:
    if "Bipolaris" in organism:
        return "B. sorokiniana"
    if "Fusarium" in organism:
        return "F. graminearum"
    return organism


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    active = pd.read_csv(AUDIT / "01_active_exchange_fluxes_mapped.csv")
    active = active[active["abs_flux"] >= THRESHOLD].copy()
    uptake = active[active["direction"] == "uptake"].copy()
    fungal = uptake[uptake["role"] == "fungus"].copy()
    bacterial = uptake[uptake["role"] == "bacteria"].copy()

    candidates = pd.read_csv(CLASSIFICATION / "05_shared_pool_compatible_nonmedium_candidates.csv")
    candidates = candidates.drop_duplicates(["run_id", "organism_donor", "organism_recipient", "chemical_key"])

    rows: list[dict[str, object]] = []
    metabolite_rows: list[dict[str, object]] = []
    for (run_id, fungal_member), fungal_data in fungal.groupby(["run_id", "organism"]):
        metadata = fungal_data.iloc[0]
        if metadata["syncom"] == "none":
            continue
        fungal_profile = fungal_data.groupby("chemical_key")["abs_flux"].sum()
        bacterial_data = bacterial[bacterial["run_id"] == run_id]
        bacterial_profile = bacterial_data.groupby("chemical_key")["abs_flux"].sum()
        keys = sorted(set(fungal_profile.index) | set(bacterial_profile.index))
        f = fungal_profile.reindex(keys, fill_value=0.0).to_numpy(float)
        b = bacterial_profile.reindex(keys, fill_value=0.0).to_numpy(float)
        shared = (f > 0) & (b > 0)
        denominator = float(np.linalg.norm(f) * np.linalg.norm(b))
        shared_min = float(np.minimum(f[shared], b[shared]).sum())
        shared_keys = [key for key, is_shared in zip(keys, shared) if is_shared]
        bacterial_co_consumers = bacterial_data[bacterial_data["chemical_key"].isin(shared_keys)]["organism"].nunique()

        candidate_data = candidates[(candidates["run_id"] == run_id) & (candidates["organism_recipient"] == fungal_member)]
        candidate_keys = candidate_data["chemical_key"].unique().tolist()
        supply_coverage = 0.0
        if candidate_keys:
            for key in candidate_keys:
                bacterial_secretion = active[(active["run_id"] == run_id) & (active["role"] == "bacteria") & (active["direction"] == "secretion") & (active["chemical_key"] == key)]["abs_flux"].sum()
                fungal_uptake = fungal_profile.get(key, 0.0)
                supply_coverage += min(bacterial_secretion, fungal_uptake)

        rows.append(
            {
                "run_id": run_id,
                "syncom": metadata["syncom"],
                "fungus_context": metadata["fungus_context"],
                "fungal_recipient": label_fungus(fungal_member),
                "fungal_uptake_metabolites": int((f > 0).sum()),
                "bacterial_uptake_metabolites": int((b > 0).sum()),
                "shared_uptake_metabolites": int(shared.sum()),
                "fungal_bacterial_cosine": float(np.dot(f, b) / denominator) if denominator else 0.0,
                "shared_min_uptake_flux": shared_min,
                "shared_fraction_of_fungal_uptake": shared_min / float(f.sum()) if f.sum() else 0.0,
                "bacterial_members_co_consuming": int(bacterial_co_consumers),
                "shared_pool_compatible_nonmedium_pairs": int(len(candidate_data)),
                "shared_pool_compatible_nonmedium_metabolites": int(len(candidate_keys)),
                "bacterial_candidate_donors": int(candidate_data["organism_donor"].nunique()),
                "potential_bacterial_supply_coverage_of_fungal_uptake": supply_coverage / float(f.sum()) if f.sum() else 0.0,
            }
        )
        for key in shared_keys:
            metabolite_rows.append(
                {
                    "run_id": run_id,
                    "syncom": metadata["syncom"],
                    "fungus_context": metadata["fungus_context"],
                    "fungal_recipient": label_fungus(fungal_member),
                    "chemical_key": key,
                    "fungal_uptake": fungal_profile[key],
                    "bacterial_uptake": bacterial_profile[key],
                    "shared_min_uptake": min(fungal_profile[key], bacterial_profile[key]),
                }
            )

    summary = pd.DataFrame(rows).sort_values(["fungal_recipient", "syncom", "fungus_context"])
    summary.to_csv(OUT / "01_fungal_syncom_interface_metric_summary.csv", index=False)
    pd.DataFrame(metabolite_rows).to_csv(OUT / "02_shared_fungal_bacterial_uptake_metabolites.csv", index=False)
    (OUT / "summary.json").write_text(json.dumps({
        "contexts": int(len(summary)),
        "flux_threshold": THRESHOLD,
        "note": "Co-uptake metrics describe resource-use overlap. Supply-coverage metrics describe the fraction of fungal uptake on non-medium, shared-pool-compatible candidate metabolites that could be covered by total bacterial secretion in the same solution; they do not prove transfer or competition.",
    }, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
