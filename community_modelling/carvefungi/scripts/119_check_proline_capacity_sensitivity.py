#!/usr/bin/env python3
"""Test whether proline capacity is associated with fungal growth limitation.

The script performs standalone FBA for the two shared-biomass fungal models
under the final shared medium. It reports:
  1. maximum proline demand flux;
  2. proline flux at the biomass-maximizing solution;
  3. capacity utilization; and
  4. growth sensitivity when proline production is progressively constrained.

This is a diagnostic analysis only; model files are not modified.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cobra
import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
DEFAULT_MEDIUM = PACKAGE / "media/biomass_variant_panel/final_shared_proline_1.csv"
OUT = PACKAGE / "reports/119_proline_capacity_sensitivity"

MODELS = [
    {
        "fungus": "Bs",
        "organism": "Bipolaris sorokiniana",
        "path": PACKAGE / "fungal_models/biomass_variant_panel/Bipolaris_sorokiniana_carvefungi_52biomass.xml",
        "biomass": "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF",
    },
    {
        "fungus": "Fg",
        "organism": "Fusarium graminearum",
        "path": PACKAGE / "fungal_models/biomass_variant_panel/Fusarium_graminearum_carvefungi_52biomass.xml",
        "biomass": "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF",
    },
]


def medium_reaction_to_external(reaction_id: str) -> str | None:
    if reaction_id.startswith("EX_") and reaction_id.endswith("_m"):
        return reaction_id[3:-2] + "_e"
    return None


def map_medium(model: cobra.Model, medium_path: Path) -> dict[str, float]:
    medium_df = pd.read_csv(medium_path)
    by_metabolite: dict[str, list[str]] = {}
    for reaction in model.exchanges:
        for metabolite in reaction.metabolites:
            by_metabolite.setdefault(metabolite.id, []).append(reaction.id)

    mapped: dict[str, float] = {}
    for _, row in medium_df.iterrows():
        rid = str(row["reaction"])
        flux = float(row["flux"])
        matches = []
        if rid in model.reactions:
            matches.append(rid)
        external = medium_reaction_to_external(rid)
        if external:
            matches.extend(by_metabolite.get(external, []))
        for match in set(matches):
            mapped[match] = max(mapped.get(match, 0.0), flux)
    return mapped


def add_proline_demand(model: cobra.Model) -> cobra.Reaction:
    reaction = cobra.Reaction("DM_TEST_proL_c")
    reaction.lower_bound = 0.0
    reaction.upper_bound = 1000.0
    reaction.add_metabolites({model.metabolites.get_by_id("proL_c"): -1.0})
    model.add_reactions([reaction])
    return reaction


def solve_growth(model: cobra.Model, biomass: cobra.Reaction) -> float | None:
    model.objective = biomass
    solution = model.optimize()
    if solution.status != "optimal":
        return None
    return float(solution.objective_value)


def analyze_model(config: dict[str, object], medium_path: Path, fractions: list[float]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    model = cobra.io.read_sbml_model(str(config["path"]))
    model.medium = map_medium(model, medium_path)
    biomass = model.reactions.get_by_id(str(config["biomass"]))

    with model:
        baseline_growth = solve_growth(model, biomass)
        demand = add_proline_demand(model)
        model.objective = demand
        maximum_solution = model.optimize()
        maximum_capacity = float(maximum_solution.objective_value) if maximum_solution.status == "optimal" else None

        model.objective = biomass
        biomass_solution = model.optimize()
        optimal_growth = float(biomass_solution.objective_value) if biomass_solution.status == "optimal" else None

        production_at_optimum = None
        if optimal_growth is not None:
            with model:
                model.reactions.get_by_id(str(biomass.id)).lower_bound = optimal_growth * (1.0 - 1e-7)
                model.objective = demand
                production_solution = model.optimize()
                if production_solution.status == "optimal":
                    production_at_optimum = float(production_solution.fluxes[demand.id])

        utilization = None
        if production_at_optimum is not None and maximum_capacity and maximum_capacity > 0:
            utilization = production_at_optimum / maximum_capacity

        summary = [{
            "fungus": config["fungus"],
            "organism": config["organism"],
            "baseline_growth": baseline_growth,
            "maximum_proline_capacity": maximum_capacity,
            "proline_at_biomass_optimum": production_at_optimum,
            "capacity_utilization": utilization,
        }]

        sensitivity = []
        for fraction in fractions:
            with model:
                demand.upper_bound = maximum_capacity * fraction if maximum_capacity is not None else 0.0
                model.objective = biomass
                solution = model.optimize()
                sensitivity.append({
                    "fungus": config["fungus"],
                    "proline_capacity_fraction": fraction,
                    "proline_upper_bound": demand.upper_bound,
                    "growth": float(solution.objective_value) if solution.status == "optimal" else None,
                    "status": solution.status,
                })
    return summary, sensitivity


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--medium", default=str(DEFAULT_MEDIUM))
    parser.add_argument("--output-dir", default=str(OUT))
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fractions = [1.0, 0.75, 0.5, 0.25, 0.1, 0.01]

    summaries: list[dict[str, object]] = []
    sensitivity: list[dict[str, object]] = []
    for config in MODELS:
        current_summary, current_sensitivity = analyze_model(config, Path(args.medium), fractions)
        summaries.extend(current_summary)
        sensitivity.extend(current_sensitivity)

    pd.DataFrame(summaries).to_csv(output_dir / "01_proline_capacity_summary.csv", index=False)
    pd.DataFrame(sensitivity).to_csv(output_dir / "02_proline_growth_sensitivity.csv", index=False)
    metadata = {
        "medium": str(Path(args.medium).resolve()),
        "models": [str(config["path"]) for config in MODELS],
        "fractions": fractions,
        "note": "Standalone FBA diagnostic; does not establish community-level causality.",
    }
    (output_dir / "summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({"output_dir": str(output_dir), "summary_rows": len(summaries), "sensitivity_rows": len(sensitivity)}, indent=2))


if __name__ == "__main__":
    main()
