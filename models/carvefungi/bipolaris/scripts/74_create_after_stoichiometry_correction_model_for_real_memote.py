#!/usr/bin/env python3
"""Script 74: create the derived after-stoichiometry-correction model.

This is the first script in this sequence that intentionally writes a corrected
derived SBML model, following user approval. The mother CarveFungi SBML is never
modified.

Derived model contents reproduce the exact Script 73 validated copy:

1. UF01339_C/m/x correction.
2. Direct harmonized biomass with Batch E artificial pool context removed.
3. old_C explicit-H neutralization background.
4. Selected two-reaction stoichiometric correction:

   UF00430_CR:
     h_r coefficient -1 -> -2

   UF00435_C:
     remove 2 h_c from reactants
     h2o_c coefficient +2 -> +1

The script writes the derived SBML and JSON files, records checksums, validates
the same gate checks as Script 73, and writes the real MEMOTE command to run on
the saved SBML.
"""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import cobra
import pandas as pd


ROOT = Path("/home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme")
SCRIPT67 = ROOT / "03_scripts/67_find_minimum_b5_inconsistency_hitting_set.py"
SCRIPT73 = ROOT / "03_scripts/73_final_validate_selected_two_reaction_stoichiometric_correction_copy_only.py"
OUT = ROOT / "06_reports/74_after_stoichiometry_correction_model_creation"

DERIVED_XML = ROOT / "04_models/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.xml"
DERIVED_JSON = ROOT / "04_models/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.json"

DIRECT_BIOMASS_ID = "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF"
SELECTED_REACTIONS = ["UF00430_CR", "UF00435_C"]
ZERO = 1e-9


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not load helper module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
                f"{m.id}:{coeff:g}"
                for m, coeff in sorted(rxn.metabolites.items(), key=lambda item: item[0].id)
            ),
        })
    return rows


def apply_selected_correction(model: cobra.Model) -> list[dict]:
    rows = []
    rows.extend(reaction_snapshot(model, "before_selected_correction", SELECTED_REACTIONS))

    # UF00430_CR: h_r coefficient -1 -> -2.
    model.reactions.get_by_id("UF00430_CR").add_metabolites({
        model.metabolites.get_by_id("h_r"): -1.0,
    })

    # UF00435_C: remove 2 h_c reactants and reduce h2o_c product from 2 to 1.
    model.reactions.get_by_id("UF00435_C").add_metabolites({
        model.metabolites.get_by_id("h_c"): 2.0,
        model.metabolites.get_by_id("h2o_c"): -1.0,
    })

    rows.extend(reaction_snapshot(model, "after_selected_correction", SELECTED_REACTIONS))
    return rows


def evaluate_model(helpers, model: cobra.Model, condition: str, baseline_growth: float | None) -> dict:
    state = helpers.memote_state(model)
    medium = helpers.load_medium(model)
    growth_status, growth = helpers.optimize_value(model, medium)
    cf_status, cf_growth = helpers.optimize_value(model, helpers.carbon_free_medium(model, medium))
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
        "objective_expression": str(model.objective.expression),
        "objective_direction": model.objective.direction,
        "biomass_reaction_id": DIRECT_BIOMASS_ID if DIRECT_BIOMASS_ID in model.reactions else "",
        "biomass_equation": biomass_equation,
        "reaction_count": len(model.reactions),
        "metabolite_count": len(model.metabolites),
        "gene_count": len(model.genes),
    }


def biomass_precursor_producibility(script73, helpers, model: cobra.Model, condition: str) -> pd.DataFrame:
    return script73.biomass_precursor_producibility(helpers, model, condition)


def ensure_safe_output_targets() -> None:
    overwrite = os.environ.get("SCRIPT74_OVERWRITE", "").strip().lower() in {"1", "true", "yes"}
    existing = [path for path in [DERIVED_XML, DERIVED_JSON] if path.exists()]
    if existing and not overwrite:
        paths = "\n".join(f"  - {path}" for path in existing)
        raise SystemExit(
            "Refusing to overwrite existing derived model file(s):\n"
            f"{paths}\n"
            "Set SCRIPT74_OVERWRITE=1 only if you intentionally want to regenerate them."
        )


def real_memote_command() -> str:
    result_json = OUT / "real_memote_after_stoichiometry_correction_results.json.gz"
    report_html = OUT / "real_memote_after_stoichiometry_correction.html"
    result_json_q = f'"{result_json}"'
    report_html_q = f'"{report_html}"'
    derived_xml_q = f'"{DERIVED_XML}"'
    return "\n".join([
        "# Full MEMOTE test collection as JSON:",
        "02_environment/.venv-carveme/bin/memote run \\",
        "  --ignore-git \\",
        "  --solver gurobi \\",
        f"  --filename {result_json_q} \\",
        f"  {derived_xml_q}",
        "",
        "# Full MEMOTE HTML snapshot report:",
        "02_environment/.venv-carveme/bin/memote report snapshot \\",
        "  --solver gurobi \\",
        f"  --filename {report_html_q} \\",
        f"  {derived_xml_q}",
    ])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ensure_safe_output_targets()

    helpers = load_module(SCRIPT67, "script67_helpers")
    script73 = load_module(SCRIPT73, "script73_helpers")

    mother_checksum_before = helpers.sha256(helpers.MODEL)
    mother = cobra.io.read_sbml_model(str(helpers.MODEL))

    baseline = helpers.build_harmonized_copy(mother)
    old_C = helpers.load_old_C(baseline)
    inventory = helpers.load_old_C_inventory(old_C)
    old_C_edits = helpers.apply_old_C_explicit_H_neutralization(baseline, inventory, old_C)

    baseline_growth_status, baseline_growth = helpers.optimize_value(baseline, helpers.load_medium(baseline))
    baseline_eval = evaluate_model(helpers, baseline, "current_viable_harmonized_baseline", baseline_growth)

    corrected = baseline.copy()
    corrected.id = "Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction"
    corrected.name = "Bipolaris sorokiniana CarveFungi after stoichiometry correction"
    correction_audit = apply_selected_correction(corrected)
    corrected_eval = evaluate_model(helpers, corrected, "after_stoichiometry_correction_model", baseline_growth)
    corrected_precursors = biomass_precursor_producibility(
        script73, helpers, corrected, "after_stoichiometry_correction_model"
    )
    precursor_summary = corrected_precursors.groupby("condition").agg(
        precursors_tested=("metabolite_id", "count"),
        precursors_producible=("producible", "sum"),
    ).reset_index()
    corrected_prec_summary = precursor_summary.iloc[0].to_dict()

    cobra.io.write_sbml_model(corrected, str(DERIVED_XML))
    cobra.io.save_json_model(corrected, str(DERIVED_JSON))

    mother_checksum_after = helpers.sha256(helpers.MODEL)
    derived_xml_checksum = helpers.sha256(DERIVED_XML)
    derived_json_checksum = helpers.sha256(DERIVED_JSON)

    biomass_unchanged = baseline_eval["biomass_equation"] == corrected_eval["biomass_equation"]
    objective_unchanged = (
        baseline_eval["objective_expression"] == corrected_eval["objective_expression"]
        and baseline_eval["objective_direction"] == corrected_eval["objective_direction"]
    )

    pass_gate = {
        "saved_derived_xml_exists": DERIVED_XML.exists(),
        "saved_derived_json_exists": DERIVED_JSON.exists(),
        "internal_memote_unconserved_zero": bool(corrected_eval["unconserved_count"] == 0),
        "growth_retained_99_999_percent": bool(
            corrected_eval["growth_ratio_vs_baseline"] is not None
            and corrected_eval["growth_ratio_vs_baseline"] >= 0.99999
        ),
        "carbon_free_growth_zero": bool(
            corrected_eval["carbon_free_growth"] is not None
            and abs(corrected_eval["carbon_free_growth"]) <= ZERO
        ),
        "all_biomass_precursors_producible": bool(
            corrected_prec_summary["precursors_tested"] == corrected_prec_summary["precursors_producible"]
        ),
        "biomass_reaction_unchanged": bool(biomass_unchanged),
        "objective_unchanged": bool(objective_unchanged),
        "mother_model_unchanged": bool(mother_checksum_before == mother_checksum_after),
    }

    evals = pd.DataFrame([baseline_eval, corrected_eval])
    evals.to_csv(OUT / "model_creation_validation.csv", index=False)
    pd.DataFrame(correction_audit).to_csv(OUT / "selected_reaction_before_after_audit.csv", index=False)
    pd.DataFrame(old_C_edits).to_csv(OUT / "old_C_explicit_H_neutralization_edit_audit.csv", index=False)
    corrected_precursors.to_csv(OUT / "biomass_precursor_producibility.csv", index=False)
    precursor_summary.to_csv(OUT / "biomass_precursor_producibility_summary.csv", index=False)
    pd.DataFrame([
        {
            "file_role": "mother_model",
            "path": str(helpers.MODEL),
            "sha256_before": mother_checksum_before,
            "sha256_after": mother_checksum_after,
            "unchanged": mother_checksum_before == mother_checksum_after,
        },
        {
            "file_role": "derived_after_stoichiometry_correction_xml",
            "path": str(DERIVED_XML),
            "sha256_before": "",
            "sha256_after": derived_xml_checksum,
            "unchanged": "",
        },
        {
            "file_role": "derived_after_stoichiometry_correction_json",
            "path": str(DERIVED_JSON),
            "sha256_before": "",
            "sha256_after": derived_json_checksum,
            "unchanged": "",
        },
    ]).to_csv(OUT / "model_checksums.csv", index=False)

    command = real_memote_command()
    (OUT / "run_real_memote_after_stoichiometry_correction.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n\n" + command + "\n",
        encoding="utf-8",
    )

    summary = {
        "mother_model": str(helpers.MODEL),
        "mother_checksum_before": mother_checksum_before,
        "mother_checksum_after": mother_checksum_after,
        "mother_model_unchanged": mother_checksum_before == mother_checksum_after,
        "derived_xml": str(DERIVED_XML),
        "derived_xml_sha256": derived_xml_checksum,
        "derived_json": str(DERIVED_JSON),
        "derived_json_sha256": derived_json_checksum,
        "selected_correction": {
            "UF00430_CR": "h_r coefficient -1 -> -2",
            "UF00435_C": "remove 2 h_c reactants; h2o_c coefficient +2 -> +1",
        },
        "baseline": baseline_eval,
        "after_stoichiometry_correction": corrected_eval,
        "old_C_explicit_H_neutralization_edits_changed": sum(
            1 for row in old_C_edits if row.get("status") == "changed"
        ),
        "biomass_reaction_unchanged": biomass_unchanged,
        "objective_unchanged": objective_unchanged,
        "biomass_precursor_summary": precursor_summary.to_dict(orient="records"),
        "pass_gate": pass_gate,
        "all_pass_gate_checks_true": all(pass_gate.values()),
        "real_memote_command": command,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    (OUT / "summary.md").write_text(
        "\n".join([
            "# Script 74 after-stoichiometry-correction model creation",
            "",
            "This script writes a derived curated model for real MEMOTE testing. The mother model is not modified.",
            "",
            "## Derived files",
            "",
            f"- SBML: `{DERIVED_XML}`",
            f"- SBML SHA256: `{derived_xml_checksum}`",
            f"- JSON: `{DERIVED_JSON}`",
            f"- JSON SHA256: `{derived_json_checksum}`",
            "",
            "## Internal validation before real MEMOTE",
            "",
            "```csv",
            evals[[
                "condition", "memote_pass", "unconserved_count", "unconserved_metabolites",
                "growth", "growth_ratio_vs_baseline", "carbon_free_growth",
                "reaction_count", "metabolite_count", "gene_count",
            ]].to_csv(index=False).strip(),
            "```",
            "",
            "## Pass gate",
            "",
            "```json",
            json.dumps(pass_gate, indent=2),
            "```",
            "",
            "## Real MEMOTE command",
            "",
            "```bash",
            command,
            "```",
            "",
            "## Integrity",
            "",
            f"- mother checksum before: `{mother_checksum_before}`",
            f"- mother checksum after: `{mother_checksum_after}`",
            f"- mother model unchanged: `{mother_checksum_before == mother_checksum_after}`",
            "",
        ]),
        encoding="utf-8",
    )

    print("Script 74 compact summary")
    print(f"derived_xml: {DERIVED_XML}")
    print(f"derived_xml_sha256: {derived_xml_checksum}")
    print(f"derived_json: {DERIVED_JSON}")
    print(f"derived_json_sha256: {derived_json_checksum}")
    print(f"baseline_unconserved_count: {baseline_eval['unconserved_count']}")
    print(f"corrected_internal_memote_pass: {corrected_eval['memote_pass']}")
    print(f"corrected_unconserved_count: {corrected_eval['unconserved_count']}")
    print(f"corrected_growth: {corrected_eval['growth']}")
    print(f"corrected_growth_ratio_vs_baseline: {corrected_eval['growth_ratio_vs_baseline']}")
    print(f"corrected_carbon_free_growth: {corrected_eval['carbon_free_growth']}")
    print(
        "corrected_biomass_precursors_producible: "
        f"{int(corrected_prec_summary['precursors_producible'])}/{int(corrected_prec_summary['precursors_tested'])}"
    )
    print(f"reaction_count_baseline_corrected: {baseline_eval['reaction_count']};{corrected_eval['reaction_count']}")
    print(f"metabolite_count_baseline_corrected: {baseline_eval['metabolite_count']};{corrected_eval['metabolite_count']}")
    print(f"gene_count_baseline_corrected: {baseline_eval['gene_count']};{corrected_eval['gene_count']}")
    print(f"mother_model_unchanged: {mother_checksum_before == mother_checksum_after}")
    print(f"all_pass_gate_checks_true: {all(pass_gate.values())}")
    print(f"outputs: {OUT}")
    print("real_memote_command:")
    print(command)


if __name__ == "__main__":
    main()
