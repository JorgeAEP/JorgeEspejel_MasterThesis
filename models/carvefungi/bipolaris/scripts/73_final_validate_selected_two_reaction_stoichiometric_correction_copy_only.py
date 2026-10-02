#!/usr/bin/env python3
"""Script 73: final validation of selected two-reaction stoichiometric correction.

Diagnostic/validation only. The mother SBML is never modified and no SBML is
saved.

This script compares:

    current viable harmonized baseline
    vs
    selected two-reaction stoichiometric representation correction

Selected correction from Script 69/70:

    UF00430_CR:
      h_r coefficient -1 -> -2

    UF00435_C:
      remove 2 h_c from reactants
      h2o_c coefficient +2 -> +1

This was the minimal among-tested passing correction with MEMOTE pass, full
growth retention, and carbon-free growth = 0.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import cobra
import pandas as pd


ROOT = Path("/home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme")
SCRIPT67 = ROOT / "03_scripts/67_find_minimum_b5_inconsistency_hitting_set.py"
OUT = ROOT / "06_reports/73_final_selected_two_reaction_stoichiometric_correction_validation"
RELAXED_BIOMASS_MAPPING = ROOT / "06_reports/relaxed_curated_biomass_variant/relaxed_biomass_mapping.csv"

DIRECT_BIOMASS_ID = "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF"
SELECTED_REACTIONS = ["UF00430_CR", "UF00435_C"]
ZERO = 1e-9


def load_script67_helpers():
    spec = importlib.util.spec_from_file_location("script67_helpers", SCRIPT67)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not load Script 67 helpers from {SCRIPT67}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def met(model: cobra.Model, mid: str):
    if mid not in model.metabolites:
        raise KeyError(f"Missing metabolite {mid}")
    return model.metabolites.get_by_id(mid)


def reaction_snapshot(model: cobra.Model, condition: str, rids: list[str]) -> list[dict]:
    rows = []
    for rid in rids:
        rxn = model.reactions.get_by_id(rid)
        rows.append({
            "condition": condition,
            "reaction_id": rid,
            "reaction_name": rxn.name,
            "reaction_equation": rxn.reaction,
            "mass_balance": json.dumps(rxn.check_mass_balance(), sort_keys=True),
            "bounds": str(rxn.bounds),
            "metabolites": "; ".join(
                f"{m.id}:{coeff:g}" for m, coeff in sorted(rxn.metabolites.items(), key=lambda item: item[0].id)
            ),
        })
    return rows


def apply_selected_correction(model: cobra.Model) -> list[dict]:
    rows = []
    for rid in SELECTED_REACTIONS:
        rows.extend(reaction_snapshot(model, f"before_selected_correction_{rid}", [rid]))

    # UF00430_CR: h_r coefficient -1 -> -2.
    model.reactions.get_by_id("UF00430_CR").add_metabolites({met(model, "h_r"): -1.0})

    # UF00435_C: remove 2 h_c reactants and reduce h2o_c product from 2 to 1.
    model.reactions.get_by_id("UF00435_C").add_metabolites({
        met(model, "h_c"): 2.0,
        met(model, "h2o_c"): -1.0,
    })

    for rid in SELECTED_REACTIONS:
        rows.extend(reaction_snapshot(model, f"after_selected_correction_{rid}", [rid]))
    return rows


def evaluate_model(helpers, model: cobra.Model, condition: str, baseline_growth: float | None) -> dict:
    state = helpers.memote_state(model)
    medium = helpers.load_medium(model)
    growth_status, growth = helpers.optimize_value(model, medium)
    cf_status, cf_growth = helpers.optimize_value(model, helpers.carbon_free_medium(model, medium))
    objective_expr = str(model.objective.expression)
    objective_direction = model.objective.direction
    biomass_equation = model.reactions.get_by_id(DIRECT_BIOMASS_ID).reaction if DIRECT_BIOMASS_ID in model.reactions else ""
    return {
        "condition": condition,
        **state,
        "unconserved_metabolites": ";".join(state["unconserved_metabolites"]),
        "growth_status": growth_status,
        "growth": growth,
        "growth_delta_vs_baseline": growth - baseline_growth if growth is not None and baseline_growth is not None else None,
        "growth_ratio_vs_baseline": growth / baseline_growth if growth is not None and baseline_growth and abs(baseline_growth) > ZERO else None,
        "carbon_free_status": cf_status,
        "carbon_free_growth": cf_growth,
        "objective_expression": objective_expr,
        "objective_direction": objective_direction,
        "biomass_reaction_id": DIRECT_BIOMASS_ID if DIRECT_BIOMASS_ID in model.reactions else "",
        "biomass_equation": biomass_equation,
        "reaction_count": len(model.reactions),
        "metabolite_count": len(model.metabolites),
        "gene_count": len(model.genes),
    }


def biomass_precursor_producibility(helpers, model: cobra.Model, condition: str) -> pd.DataFrame:
    mapping = pd.read_csv(RELAXED_BIOMASS_MAPPING)
    included = mapping.loc[mapping["status"].astype(str) == "included"].copy()
    medium = helpers.load_medium(model)
    rows = []
    for _, row in included.iterrows():
        mid = str(row["carve_id"])
        if mid not in model.metabolites:
            rows.append({
                "condition": condition,
                "metabolite_id": mid,
                "metabolite_name": row.get("metabolite_name", ""),
                "status": "missing_metabolite",
                "max_production": None,
                "producible": False,
            })
            continue
        with model:
            model.medium = {rid: value for rid, value in medium.items() if rid in model.reactions}
            demand = model.add_boundary(model.metabolites.get_by_id(mid), type="demand")
            model.objective = demand
            solution = model.optimize()
            value = float(solution.objective_value or 0.0) if solution.status == "optimal" else None
            rows.append({
                "condition": condition,
                "metabolite_id": mid,
                "metabolite_name": row.get("metabolite_name", ""),
                "status": solution.status,
                "max_production": value,
                "producible": bool(value is not None and value > 1e-8),
            })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    helpers = load_script67_helpers()
    before_checksum = helpers.sha256(helpers.MODEL)

    mother = cobra.io.read_sbml_model(str(helpers.MODEL))
    baseline = helpers.build_harmonized_copy(mother)
    old_C = helpers.load_old_C(baseline)
    inventory = helpers.load_old_C_inventory(old_C)
    old_C_edits_baseline = helpers.apply_old_C_explicit_H_neutralization(baseline, inventory, old_C)

    baseline_growth_status, baseline_growth = helpers.optimize_value(baseline, helpers.load_medium(baseline))
    baseline_eval = evaluate_model(helpers, baseline, "current_viable_harmonized_baseline", baseline_growth)
    baseline_precursors = biomass_precursor_producibility(helpers, baseline, "current_viable_harmonized_baseline")

    corrected = baseline.copy()
    correction_audit = apply_selected_correction(corrected)
    corrected_eval = evaluate_model(helpers, corrected, "selected_two_reaction_stoichiometric_correction", baseline_growth)
    corrected_precursors = biomass_precursor_producibility(helpers, corrected, "selected_two_reaction_stoichiometric_correction")

    evals = pd.DataFrame([baseline_eval, corrected_eval])
    precursor_tests = pd.concat([baseline_precursors, corrected_precursors], ignore_index=True)
    correction_audit_df = pd.DataFrame(correction_audit)
    old_C_edit_table = pd.DataFrame(old_C_edits_baseline)

    biomass_unchanged = baseline_eval["biomass_equation"] == corrected_eval["biomass_equation"]
    objective_unchanged = (
        baseline_eval["objective_expression"] == corrected_eval["objective_expression"]
        and baseline_eval["objective_direction"] == corrected_eval["objective_direction"]
    )
    precursor_summary = precursor_tests.groupby("condition").agg(
        precursors_tested=("metabolite_id", "count"),
        precursors_producible=("producible", "sum"),
    ).reset_index()

    evals.to_csv(OUT / "baseline_vs_selected_correction_validation.csv", index=False)
    precursor_tests.to_csv(OUT / "biomass_precursor_producibility.csv", index=False)
    precursor_summary.to_csv(OUT / "biomass_precursor_producibility_summary.csv", index=False)
    correction_audit_df.to_csv(OUT / "selected_reaction_before_after_audit.csv", index=False)
    old_C_edit_table.to_csv(OUT / "old_C_explicit_H_neutralization_edit_audit.csv", index=False)

    with pd.ExcelWriter(OUT / "final_selected_correction_validation.xlsx", engine="xlsxwriter") as writer:
        evals.to_excel(writer, sheet_name="validation", index=False)
        precursor_summary.to_excel(writer, sheet_name="precursor_summary", index=False)
        precursor_tests.to_excel(writer, sheet_name="precursor_tests", index=False)
        correction_audit_df.to_excel(writer, sheet_name="before_after_reactions", index=False)
        old_C_edit_table.to_excel(writer, sheet_name="H_edit_audit", index=False)
        wrap = writer.book.add_format({"text_wrap": True, "valign": "top"})
        for sheet in writer.sheets.values():
            sheet.freeze_panes(1, 0)
            sheet.set_column(0, 10, 34)
            sheet.set_column(11, 80, 48, wrap)

    after_checksum = helpers.sha256(helpers.MODEL)
    corrected_prec_summary = precursor_summary.loc[
        precursor_summary["condition"] == "selected_two_reaction_stoichiometric_correction"
    ].iloc[0].to_dict()
    pass_gate = {
        "memote_unconserved_zero": bool(corrected_eval["unconserved_count"] == 0),
        "growth_retained_99_999_percent": bool(
            corrected_eval["growth_ratio_vs_baseline"] is not None
            and corrected_eval["growth_ratio_vs_baseline"] >= 0.99999
        ),
        "carbon_free_growth_zero": bool(
            corrected_eval["carbon_free_growth"] is not None
            and abs(corrected_eval["carbon_free_growth"]) <= 1e-9
        ),
        "all_biomass_precursors_producible": bool(
            corrected_prec_summary["precursors_tested"] == corrected_prec_summary["precursors_producible"]
        ),
        "biomass_reaction_unchanged": bool(biomass_unchanged),
        "objective_unchanged": bool(objective_unchanged),
        "mother_model_unchanged": bool(before_checksum == after_checksum),
    }
    summary = {
        "mother_model": str(helpers.MODEL),
        "mother_checksum_before": before_checksum,
        "mother_checksum_after": after_checksum,
        "mother_model_unchanged": before_checksum == after_checksum,
        "diagnostic_only_no_sbml_saved": True,
        "selected_correction": {
            "UF00430_CR": "h_r coefficient -1 -> -2",
            "UF00435_C": "remove 2 h_c reactants; h2o_c coefficient +2 -> +1",
        },
        "baseline": baseline_eval,
        "selected_correction_validation": corrected_eval,
        "biomass_reaction_unchanged": biomass_unchanged,
        "objective_unchanged": objective_unchanged,
        "biomass_precursor_summary": precursor_summary.to_dict(orient="records"),
        "pass_gate": pass_gate,
        "all_pass_gate_checks_true": all(pass_gate.values()),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    md = "\n".join([
        "# Script 73 final selected stoichiometric correction validation", "",
        "Copy-only final validation. No SBML was saved.", "",
        "## Validation", "",
        "```csv",
        evals[[
            "condition", "memote_pass", "unconserved_count", "unconserved_metabolites",
            "growth", "growth_ratio_vs_baseline", "carbon_free_growth",
            "reaction_count", "metabolite_count", "gene_count",
        ]].to_csv(index=False).strip(),
        "```", "",
        "## Pass gate", "",
        "```json",
        json.dumps(pass_gate, indent=2),
        "```", "",
        "## Biomass precursor producibility", "",
        "```csv",
        precursor_summary.to_csv(index=False).strip(),
        "```", "",
        "## Selected reaction before/after", "",
        "```csv",
        correction_audit_df.to_csv(index=False).strip(),
        "```", "",
        "## Integrity", "",
        f"- mother checksum before: `{before_checksum}`",
        f"- mother checksum after: `{after_checksum}`",
        f"- mother model unchanged: `{before_checksum == after_checksum}`", "",
    ])
    (OUT / "summary.md").write_text(md, encoding="utf-8")

    print("Script 73 compact summary")
    print(f"baseline_unconserved_count: {baseline_eval['unconserved_count']}")
    print(f"baseline_unconserved_metabolites: {baseline_eval['unconserved_metabolites']}")
    print(f"baseline_growth: {baseline_eval['growth']}")
    print(f"baseline_carbon_free_growth: {baseline_eval['carbon_free_growth']}")
    print(f"corrected_memote_pass: {corrected_eval['memote_pass']}")
    print(f"corrected_unconserved_count: {corrected_eval['unconserved_count']}")
    print(f"corrected_unconserved_metabolites: {corrected_eval['unconserved_metabolites']}")
    print(f"corrected_growth: {corrected_eval['growth']}")
    print(f"corrected_growth_ratio_vs_baseline: {corrected_eval['growth_ratio_vs_baseline']}")
    print(f"corrected_carbon_free_growth: {corrected_eval['carbon_free_growth']}")
    print(f"biomass_reaction_unchanged: {biomass_unchanged}")
    print(f"objective_unchanged: {objective_unchanged}")
    print(
        "corrected_biomass_precursors_producible: "
        f"{int(corrected_prec_summary['precursors_producible'])}/{int(corrected_prec_summary['precursors_tested'])}"
    )
    print(f"reaction_count_baseline_corrected: {baseline_eval['reaction_count']};{corrected_eval['reaction_count']}")
    print(f"metabolite_count_baseline_corrected: {baseline_eval['metabolite_count']};{corrected_eval['metabolite_count']}")
    print(f"gene_count_baseline_corrected: {baseline_eval['gene_count']};{corrected_eval['gene_count']}")
    print(f"mother_model_unchanged: {before_checksum == after_checksum}")
    print(f"all_pass_gate_checks_true: {all(pass_gate.values())}")
    print(f"outputs: {OUT}")


if __name__ == "__main__":
    main()
