#!/usr/bin/env python3
"""Script 21: final copy-only validation of selected Fusarium correction.

No SBML is saved. The mother model is never modified.

Selected Fusarium stoichiometric-correction candidate:

1. Correct UF01339_C:

       datp_c + 2 h2o_c -> damp_c + 2 h_c + 2 ppi_c
       becomes
       datp_c + h2o_c -> damp_c + h_c + ppi_c

2. Replace artificial biomass/pool context with direct 52-component biomass:

       remove the 15 Batch-E-overlap reactions
       add BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF

3. Remove final 5 trigger/interface reactions:

       UF00283_CM
       UF00284_CM
       UF00806_M
       UF00874_C
       UF01036_C

Script 20 showed that UF00794_c, UF00969_C, and UF00970_C can be retained while
maintaining MEMOTE consistency and growth. Therefore those are not removed.

This script validates:

    MEMOTE support stoichiometric consistency
    unconserved metabolites
    native-medium growth
    selected shared+37-medium growth
    carbon-free growth
    52/52 biomass precursor producibility
    model size
    mother checksum
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import cobra
from cobra import Reaction


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "04_models/Fusarium_graminearum_carvefungi.xml"
SCRIPT_15 = ROOT / "03_scripts/15_test_direct_52_biomass_after_UF01339_and_BatchE_removal_copy_only.py"
OUT = ROOT / "06_reports/21_final_selected_fusarium_stoichiometric_correction_validation"

FINAL_TRIGGER_REMOVALS = [
    "UF00283_CM",
    "UF00284_CM",
    "UF00806_M",
    "UF00874_C",
    "UF01036_C",
]
RETAINED_FROM_TRIGGER_SET = ["UF00794_c", "UF00969_C", "UF00970_C"]
HARMONIZED_BIOMASS_ID = "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF"
ZERO = 1e-9


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S15 = load_module(SCRIPT_15, "script15_for_script21")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row}) if rows else ["empty"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def existing(model: cobra.Model, reaction_ids: Iterable[str]) -> list[str]:
    return sorted(rid for rid in set(reaction_ids) if rid in model.reactions)


def reaction_snapshot(model: cobra.Model, reaction_id: str, condition: str) -> dict[str, object]:
    reaction = model.reactions.get_by_id(reaction_id)
    return {
        "condition": condition,
        "reaction_id": reaction.id,
        "reaction_name": reaction.name,
        "equation": reaction.reaction,
        "bounds": str(reaction.bounds),
        "mass_balance": json.dumps(reaction.check_mass_balance(), sort_keys=True),
        "metabolites": "; ".join(
            f"{met.id}:{coeff:g}"
            for met, coeff in sorted(reaction.metabolites.items(), key=lambda item: item[0].id)
        ),
        "gene_reaction_rule": reaction.gene_reaction_rule,
    }


def build_selected_corrected_copy() -> tuple[cobra.Model, dict[str, float], dict[str, float], list[dict[str, object]]]:
    mother = cobra.io.read_sbml_model(str(MODEL))
    native_medium = dict(mother.medium)

    bipolaris = cobra.io.read_sbml_model(str(S15.S06.BIPOLARIS_FINAL))
    bipolaris_biomass = S15.S06.find_bipolaris_harmonized_biomass(bipolaris)

    work = mother.copy()
    audit_rows: list[dict[str, object]] = []
    audit_rows.extend(S15.apply_uf01339_correction(work))
    removed_batch_e = S15.remove_batch_e_overlap(work)
    for rid in removed_batch_e:
        if rid in mother.reactions:
            audit_rows.append(reaction_snapshot(mother, rid, "removed_Batch_E_overlap_reaction"))

    biomass_info = S15.S06.add_harmonized_biomass(work, bipolaris_biomass)
    biomass_id = str(biomass_info["harmonized_biomass_id"])
    work.reactions.get_by_id(biomass_id).annotation["sbo"] = "SBO:0000629"
    work.reactions.get_by_id(biomass_id).name = (
        "biomass Fusarium direct 52-component harmonized CarveFungi biomass"
    )
    audit_rows.append(reaction_snapshot(work, biomass_id, "added_direct_52_component_biomass"))

    final_removed = existing(work, FINAL_TRIGGER_REMOVALS)
    for rid in final_removed:
        audit_rows.append(reaction_snapshot(work, rid, "removed_final_trigger_reaction"))
    work.remove_reactions(final_removed, remove_orphans=False)

    for rid in RETAINED_FROM_TRIGGER_SET:
        if rid in work.reactions:
            audit_rows.append(reaction_snapshot(work, rid, "retained_from_8_reaction_trigger_set"))

    selected_medium, medium_rows = S15.build_shared_plus_37_selected_medium(work.copy())
    for row in medium_rows:
        row["condition"] = "selected_shared_plus_37_medium"
        audit_rows.append(row)
    return work, native_medium, selected_medium, audit_rows


def optimize_with_medium(model: cobra.Model, medium: dict[str, float]) -> dict[str, object]:
    test = model.copy()
    test.medium = {
        reaction_id: float(value)
        for reaction_id, value in medium.items()
        if test.reactions.has_id(reaction_id) and abs(float(value)) > ZERO
    }
    solution = test.optimize()
    value = float(solution.objective_value) if solution.objective_value is not None else 0.0
    return {
        "status": solution.status,
        "growth": value,
        "medium_entries_applied": len(test.medium),
    }


def carbon_free_medium(model: cobra.Model, medium: dict[str, float]) -> dict[str, float]:
    return S15.carbon_free_medium(model, medium)


def biomass_precursor_producibility(
    model: cobra.Model,
    medium: dict[str, float],
    biomass_id: str,
    medium_condition: str,
) -> list[dict[str, object]]:
    test = model.copy()
    test.medium = {
        reaction_id: float(value)
        for reaction_id, value in medium.items()
        if test.reactions.has_id(reaction_id) and abs(float(value)) > ZERO
    }
    biomass = test.reactions.get_by_id(biomass_id)
    rows: list[dict[str, object]] = []
    for metabolite, coefficient in sorted(biomass.metabolites.items(), key=lambda item: item[0].id):
        if coefficient >= 0:
            continue
        with test:
            demand = Reaction(f"TMP_DM_{metabolite.id}")
            demand.bounds = (0.0, 1000.0)
            demand.add_metabolites({metabolite: -1.0})
            test.add_reactions([demand])
            test.objective = demand
            solution = test.optimize()
            value = float(solution.objective_value) if solution.objective_value is not None else 0.0
        rows.append({
            "medium_condition": medium_condition,
            "metabolite_id": metabolite.id,
            "metabolite_name": metabolite.name,
            "coefficient": float(coefficient),
            "formula": metabolite.formula or "",
            "compartment": metabolite.compartment,
            "maximum_production": value,
            "status": solution.status,
            "producible": value > ZERO,
        })
    return rows


def evaluate_condition(
    model: cobra.Model,
    condition: str,
    medium: dict[str, float],
) -> dict[str, object]:
    state = S15.quiet_memote_state(model)
    growth = optimize_with_medium(model, medium)
    carbon_free = optimize_with_medium(model, carbon_free_medium(model, medium))
    return {
        "condition": condition,
        **state,
        "growth_status": growth["status"],
        "growth": growth["growth"],
        "medium_entries_applied": growth["medium_entries_applied"],
        "carbon_free_status": carbon_free["status"],
        "carbon_free_growth": carbon_free["growth"],
        "carbon_free_medium_entries_applied": carbon_free["medium_entries_applied"],
        "reaction_count": len(model.reactions),
        "metabolite_count": len(model.metabolites),
        "gene_count": len(model.genes),
        "objective": str(model.objective.expression),
        "biomass_reaction": model.reactions.get_by_id(HARMONIZED_BIOMASS_ID).reaction,
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mother_before = sha256(MODEL)
    bipolaris_before = sha256(S15.S06.BIPOLARIS_FINAL)

    corrected, native_medium, selected_medium, audit_rows = build_selected_corrected_copy()
    condition_rows = [
        evaluate_condition(corrected, "selected_correction_native_medium", native_medium),
        evaluate_condition(corrected, "selected_correction_selected_shared_plus_37_medium", selected_medium),
    ]

    precursor_rows = []
    precursor_rows.extend(
        biomass_precursor_producibility(corrected, native_medium, HARMONIZED_BIOMASS_ID, "native_medium")
    )
    precursor_rows.extend(
        biomass_precursor_producibility(
            corrected, selected_medium, HARMONIZED_BIOMASS_ID, "selected_shared_plus_37_medium"
        )
    )
    precursor_summary = []
    for medium_condition in ["native_medium", "selected_shared_plus_37_medium"]:
        rows = [row for row in precursor_rows if row["medium_condition"] == medium_condition]
        precursor_summary.append({
            "medium_condition": medium_condition,
            "precursors_tested": len(rows),
            "precursors_producible": sum(bool(row["producible"]) for row in rows),
        })

    write_csv(OUT / "final_validation_condition_results.csv", condition_rows)
    write_csv(OUT / "selected_correction_reaction_and_medium_audit.csv", audit_rows)
    write_csv(OUT / "biomass_precursor_producibility.csv", precursor_rows)
    write_csv(OUT / "biomass_precursor_producibility_summary.csv", precursor_summary)

    native = condition_rows[0]
    selected = condition_rows[1]
    native_prec = precursor_summary[0]
    selected_prec = precursor_summary[1]

    mother_after = sha256(MODEL)
    bipolaris_after = sha256(S15.S06.BIPOLARIS_FINAL)

    pass_gate = {
        "native_memote_pass": bool(native["memote_pass"]),
        "selected_memote_pass": bool(selected["memote_pass"]),
        "native_unconserved_zero": int(native["unconserved_count"]) == 0,
        "selected_unconserved_zero": int(selected["unconserved_count"]) == 0,
        "native_growth_positive": float(native["growth"]) > ZERO,
        "selected_growth_positive": float(selected["growth"]) > ZERO,
        "native_carbon_free_zero": abs(float(native["carbon_free_growth"])) <= 1e-8,
        "selected_carbon_free_zero": abs(float(selected["carbon_free_growth"])) <= 1e-8,
        "native_precursors_all_producible": native_prec["precursors_tested"] == native_prec["precursors_producible"],
        "selected_precursors_all_producible": selected_prec["precursors_tested"] == selected_prec["precursors_producible"],
        "mother_unchanged": mother_before == mother_after,
        "bipolaris_reference_unchanged": bipolaris_before == bipolaris_after,
    }

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "mother_model": str(MODEL),
        "mother_sha256_before": mother_before,
        "mother_sha256_after": mother_after,
        "mother_unchanged": mother_before == mother_after,
        "bipolaris_reference": str(S15.S06.BIPOLARIS_FINAL),
        "bipolaris_reference_sha256_before": bipolaris_before,
        "bipolaris_reference_sha256_after": bipolaris_after,
        "bipolaris_reference_unchanged": bipolaris_before == bipolaris_after,
        "selected_correction": {
            "final_trigger_removals": FINAL_TRIGGER_REMOVALS,
            "retained_from_trigger_set": RETAINED_FROM_TRIGGER_SET,
            "harmonized_biomass_id": HARMONIZED_BIOMASS_ID,
        },
        "conditions": condition_rows,
        "precursor_summary": precursor_summary,
        "pass_gate": pass_gate,
        "all_pass_gate_checks_true": all(pass_gate.values()),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 21 compact summary")
    print(f"mother_unchanged: {mother_before == mother_after}")
    print(f"bipolaris_reference_unchanged: {bipolaris_before == bipolaris_after}")
    print(f"final_trigger_removals: {';'.join(FINAL_TRIGGER_REMOVALS)}")
    print(f"retained_from_trigger_set: {';'.join(RETAINED_FROM_TRIGGER_SET)}")
    print(f"native_memote_pass: {native['memote_pass']}")
    print(f"native_unconserved_count: {native['unconserved_count']}")
    print(f"native_unconserved_metabolites: {native['unconserved_metabolites']}")
    print(f"native_growth: {native['growth']}")
    print(f"native_carbon_free_growth: {native['carbon_free_growth']}")
    print(f"selected_memote_pass: {selected['memote_pass']}")
    print(f"selected_unconserved_count: {selected['unconserved_count']}")
    print(f"selected_unconserved_metabolites: {selected['unconserved_metabolites']}")
    print(f"selected_growth: {selected['growth']}")
    print(f"selected_carbon_free_growth: {selected['carbon_free_growth']}")
    print(
        "native_biomass_precursors_producible: "
        f"{native_prec['precursors_producible']}/{native_prec['precursors_tested']}"
    )
    print(
        "selected_biomass_precursors_producible: "
        f"{selected_prec['precursors_producible']}/{selected_prec['precursors_tested']}"
    )
    print(f"reaction_count: {selected['reaction_count']}")
    print(f"metabolite_count: {selected['metabolite_count']}")
    print(f"gene_count: {selected['gene_count']}")
    print(f"all_pass_gate_checks_true: {all(pass_gate.values())}")
    print(f"outputs: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
