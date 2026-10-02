#!/usr/bin/env python3
"""Check producibility of shared 52-component biomass precursors.

This is an optional analysis for Results section 3.3. It loads the two
CarveFungi shared-52 fungal models, applies the final shared barley medium,
and tests whether each biomass reactant can be produced by adding a temporary
demand reaction and maximizing that demand.

It does not edit model files.

Outputs:
    reports/118_shared52_biomass_precursor_producibility/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cobra
import pandas as pd


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
OUT = PACKAGE / "reports/118_shared52_biomass_precursor_producibility"
MEDIUM = PACKAGE / "media/biomass_variant_panel/final_shared_proline_1.csv"

MODEL_CONFIGS = [
    {
        "fungus": "Bs",
        "organism": "Bipolaris sorokiniana",
        "model_path": PACKAGE / "fungal_models/biomass_variant_panel/Bipolaris_sorokiniana_carvefungi_52biomass.xml",
        "biomass_id": "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF",
    },
    {
        "fungus": "Fg",
        "organism": "Fusarium graminearum",
        "model_path": PACKAGE / "fungal_models/biomass_variant_panel/Fusarium_graminearum_carvefungi_52biomass.xml",
        "biomass_id": "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF",
    },
]


def medium_reaction_to_external_metabolite_id(reaction_id: str) -> str | None:
    if not reaction_id.startswith("EX_") or not reaction_id.endswith("_m"):
        return None
    return reaction_id[3:-2] + "_e"


def map_medium_to_model(model: cobra.Model, medium_df: pd.DataFrame) -> tuple[dict[str, float], pd.DataFrame]:
    exchange_by_met: dict[str, list[str]] = {}
    for reaction in model.exchanges:
        for metabolite in reaction.metabolites:
            exchange_by_met.setdefault(metabolite.id, []).append(reaction.id)

    medium: dict[str, float] = {}
    rows: list[dict[str, object]] = []
    for _, row in medium_df.iterrows():
        medium_rxn = str(row["reaction"])
        flux = float(row["flux"])
        external_met = medium_reaction_to_external_metabolite_id(medium_rxn)
        matched = []
        if medium_rxn in model.reactions:
            matched.append(medium_rxn)
        if external_met:
            matched.extend(exchange_by_met.get(external_met, []))
        matched = sorted(set(matched))
        for rid in matched:
            medium[rid] = max(medium.get(rid, 0.0), flux)
        rows.append(
            {
                "medium_reaction": medium_rxn,
                "external_metabolite_id": external_met or "",
                "matched_model_exchange_count": len(matched),
                "matched_model_exchanges": ";".join(matched),
            }
        )
    return medium, pd.DataFrame(rows)


def classify_precursor(met_id: str, met_name: str) -> str:
    text = f"{met_id} {met_name}".lower()
    amino = [
        "ala",
        "arg",
        "asn",
        "asp",
        "cys",
        "gln",
        "glu",
        "gly",
        "his",
        "ile",
        "leu",
        "lys",
        "met",
        "phe",
        "pro",
        "ser",
        "thr",
        "trp",
        "tyr",
        "val",
    ]
    nucleotide = ["amp", "cmp", "gmp", "ump", "datp", "dctp", "dgtp", "dttp", "atp", "ctp", "gtp", "utp"]
    lipid = ["fa", "lipid", "coa", "lnl", "hdca", "ocdca", "stear", "ole", "ergost"]
    carbohydrate = ["gluc", "glycogen", "udp", "uacgam", "chitin", "trehalose", "mannan"]
    cofactor = ["thm", "nad", "fad", "coa", "heme", "pydx", "riboflavin", "folate"]
    if any(token in text for token in amino):
        return "amino acid"
    if any(token in text for token in nucleotide):
        return "nucleotide"
    if any(token in text for token in lipid):
        return "lipid/cofactor"
    if any(token in text for token in carbohydrate):
        return "carbohydrate/cell wall"
    if any(token in text for token in cofactor):
        return "cofactor"
    return "other"


def objective_value(solution) -> float | None:
    if solution.status != "optimal" or solution.objective_value is None:
        return None
    return float(solution.objective_value)


def test_precursor(model: cobra.Model, medium: dict[str, float], metabolite: cobra.Metabolite) -> dict[str, object]:
    demand_id = f"DM_TEST_{metabolite.id}"
    with model:
        model.medium = {rid: value for rid, value in medium.items() if rid in model.reactions}
        demand = cobra.Reaction(demand_id)
        demand.lower_bound = 0.0
        demand.upper_bound = 1000.0
        demand.add_metabolites({metabolite: -1.0})
        model.add_reactions([demand])
        model.objective = demand
        solution = model.optimize()
        value = objective_value(solution)
        return {
            "demand_reaction": demand_id,
            "status": solution.status,
            "max_demand_flux": value,
            "producible": bool(value is not None and value > 1e-9),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--medium", default=str(MEDIUM))
    parser.add_argument("--max-terms", type=int, default=0, help="Debug limit; 0 tests all terms.")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    medium_df = pd.read_csv(args.medium)
    all_rows = []
    mapping_rows = []

    for cfg in MODEL_CONFIGS:
        model = cobra.io.read_sbml_model(str(cfg["model_path"]))
        mapped_medium, mapping = map_medium_to_model(model, medium_df)
        mapping.insert(0, "fungus", cfg["fungus"])
        mapping_rows.append(mapping)

        biomass = model.reactions.get_by_id(cfg["biomass_id"])
        reactants = [(met, coeff) for met, coeff in biomass.metabolites.items() if coeff < 0]
        reactants = sorted(reactants, key=lambda item: item[0].id)
        if args.max_terms > 0:
            reactants = reactants[: args.max_terms]

        with model:
            model.medium = {rid: value for rid, value in mapped_medium.items() if rid in model.reactions}
            model.objective = biomass
            growth_sol = model.optimize()
            baseline_growth = objective_value(growth_sol)
            baseline_status = growth_sol.status

        for met, coeff in reactants:
            result = test_precursor(model, mapped_medium, met)
            all_rows.append(
                {
                    "fungus": cfg["fungus"],
                    "organism": cfg["organism"],
                    "model_path": str(cfg["model_path"]),
                    "biomass_id": cfg["biomass_id"],
                    "baseline_growth_status": baseline_status,
                    "baseline_growth": baseline_growth,
                    "metabolite_id": met.id,
                    "metabolite_name": met.name,
                    "biomass_coefficient": float(coeff),
                    "precursor_class": classify_precursor(met.id, met.name),
                    **result,
                }
            )

    results = pd.DataFrame(all_rows)
    results.to_csv(OUT / "01_shared52_precursor_producibility_long.csv", index=False)
    pd.concat(mapping_rows, ignore_index=True).to_csv(OUT / "02_medium_mapping_audit.csv", index=False)

    matrix = results.pivot_table(
        index=["metabolite_id", "metabolite_name", "precursor_class"],
        columns="fungus",
        values="max_demand_flux",
        aggfunc="first",
    ).reset_index()
    matrix["Bs_to_Fg_ratio"] = matrix["Bs"] / matrix["Fg"].replace({0: pd.NA})
    matrix.to_csv(OUT / "03_shared52_precursor_producibility_matrix.csv", index=False)

    class_summary = (
        results.groupby(["fungus", "precursor_class"])
        .agg(
            term_count=("metabolite_id", "count"),
            producible_count=("producible", "sum"),
            median_max_demand=("max_demand_flux", "median"),
            min_max_demand=("max_demand_flux", "min"),
            max_max_demand=("max_demand_flux", "max"),
        )
        .reset_index()
    )
    class_summary["producible_fraction"] = class_summary["producible_count"] / class_summary["term_count"]
    class_summary.to_csv(OUT / "04_shared52_precursor_class_summary.csv", index=False)

    summary = {
        "output_dir": str(OUT),
        "medium": str(args.medium),
        "rows": len(results),
        "models": [cfg["fungus"] for cfg in MODEL_CONFIGS],
        "nonproducible_terms": int((~results["producible"]).sum()) if not results.empty else 0,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
