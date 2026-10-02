#!/usr/bin/env python3
"""Extract proline-related quantities from the existing shared-52 MICOM runs.

This script does not rerun MICOM. It combines fungal growth, biomass
coefficients and fungal exchange summaries for the nine shared-52 contexts.

Important limitation: the exchange summaries contain organism-level exchange
fluxes, not intracellular reaction fluxes. Intracellular proline synthesis is
therefore reported as not identifiable here rather than inferred from growth.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
REPORTS = PACKAGE / "reports"
SOURCE = REPORTS / "96_fungal_species_syncom_factorial_panel/summaries"
OUT = REPORTS / "120_shared52_proline_community_fluxes"

# Coefficients are negative in the biomass reactions; positive values below
# represent the biomass-associated proline requirement per unit fungal growth.
PROLINE_COEFFICIENT = {"Bs": 0.0866064142634991, "Fg": 0.0866064142634991}
STANDALONE_CAPACITY = {"Bs": 952.2298245614036, "Fg": 1.0}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    growth = pd.read_csv(SOURCE / "growth_summary_all_runs.csv")
    exchanges = pd.read_csv(REPORTS / "103_exchange_interpretation_tables/04_fungal_exchange_fluxes_by_run.csv")

    growth = growth[(growth["biomass_mode"] == "shared_52") & (growth["scenario"] == "proline_relaxed_1")].copy()
    fungal = exchanges[
        (exchanges["biomass_mode"] == "shared_52")
        & (exchanges["scenario"] == "proline_relaxed_1")
        & (exchanges["role"] == "fungus")
    ].copy()

    # Growth summaries may contain multiple organism rows per run. Keep the
    # fungal member corresponding to the context and retain its growth value.
    fungal_growth = growth[(growth["fungus"].isin(["Bs", "Fg"])) & (growth["role"] == "fungus")].copy()
    fungal_growth = fungal_growth[["run_id", "syncom", "fungus_context", "community_type", "fungus", "growth_rate"]]
    fungal_growth = fungal_growth.rename(columns={"growth_rate": "fungal_growth"})

    rows = []
    for _, row in fungal_growth.iterrows():
        fungus = row["fungus"]
        run_id = row["run_id"]
        matching = fungal[(fungal["run_id"] == run_id) & (fungal["fungus"] == fungus)]
        # The harmonized exchange namespace is expected to be EX_proL_e. The
        # medium namespace may also appear as EX_proL_m in raw outputs.
        proline = matching[matching["exchange"].isin(["EX_proL_e", "EX_proL_m", "EX_proL"])]
        uptake = float(-proline.loc[proline["flux"] < 0, "flux"].sum()) if not proline.empty else 0.0
        secretion = float(proline.loc[proline["flux"] > 0, "flux"].sum()) if not proline.empty else 0.0
        growth_value = float(row["fungal_growth"])
        requirement = growth_value * PROLINE_COEFFICIENT[fungus]
        rows.append({
            "run_id": run_id,
            "context": f"{fungus}+{row['syncom']}" if row["syncom"] != "none" else str(fungus),
            "syncom": row["syncom"],
            "fungus": fungus,
            "community_type": row["community_type"],
            "fungal_growth": growth_value,
            "biomass_associated_proline_requirement": requirement,
            "external_proline_uptake": uptake,
            "external_proline_secretion": secretion,
            "intracellular_proline_synthesis": None,
            "intracellular_synthesis_status": "not identifiable from exchange summaries",
            "standalone_max_proline_capacity": STANDALONE_CAPACITY[fungus],
            "requirement_to_standalone_capacity": requirement / STANDALONE_CAPACITY[fungus],
        })

    result = pd.DataFrame(rows).sort_values(["fungus", "syncom", "community_type"])
    result.to_csv(OUT / "01_shared52_proline_community_fluxes.csv", index=False)
    (OUT / "summary.json").write_text(json.dumps({
        "output_dir": str(OUT),
        "rows": len(result),
        "source_exchange_table": str(REPORTS / "103_exchange_interpretation_tables/04_fungal_exchange_fluxes_by_run.csv"),
        "note": "External exchange is observable; intracellular synthesis requires reaction-level MICOM fluxes or a dedicated constrained rerun.",
    }, indent=2), encoding="utf-8")
    print(json.dumps({"output_dir": str(OUT), "rows": len(result)}, indent=2))


if __name__ == "__main__":
    main()
