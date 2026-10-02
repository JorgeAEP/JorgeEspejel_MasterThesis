#!/usr/bin/env python3
"""Script 15: test the Fusarium analogue of the Bipolaris 23->7 step.

Copy-only diagnostic. The mother SBML is never modified and no SBML is saved.

From Script 14:

    raw Fusarium                         32 unconserved metabolites
    + UF01339_C correction               22 unconserved metabolites
    + remove Bipolaris Batch-E overlap    7 unconserved metabolites

The deletion-only condition had zero growth because it removed the native
artificial biomass objective. This script repeats the same stoichiometric
background but adds the already validated Fusarium 52-component direct biomass
objective from Script 06.

This is the proper Fusarium equivalent of the successful Bipolaris background:

    UF01339 correction
    + direct 52-component biomass
    + remove artificial biomass/macromolecule/polymer context

Outputs are tables only.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import contextlib
from datetime import datetime, timezone
from pathlib import Path

import cobra
from cobra import Reaction
from memote.support import consistency as con
from memote.support.consistency import con_helpers


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "04_models/Fusarium_graminearum_carvefungi.xml"
SCRIPT_06 = ROOT / "03_scripts/06_test_bipolaris_52_component_biomass_and_shared_medium_copy_only.py"
SCRIPT_08 = ROOT / "03_scripts/08_find_grouped_native_medium_rescue_for_shared_medium_growth.py"
OUT = ROOT / "06_reports/15_direct_52_biomass_after_UF01339_and_BatchE_removal_copy_only"

HARMONIZED_BIOMASS_ID = "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF"
ZERO = 1e-9

BIPOLARIS_BATCH_E_CONTEXT = [
    "BIOMASS",
    "UF01078_C",
    "UF01092_C",
    "UF01093_C",
    "UF01094_C",
    "UF01408_C",
    "UF01408_x",
    "UF01729_C",
    "UF01730_O",
    "UF01737_O",
    "UF01822_C",
    "UF01822_x",
    "UF02405_C",
    "UF02406_C",
    "UF02407_C",
    "UF02409_C",
    "UF02414_C",
]

SELECTED_FUSARIUM_SUPPLEMENT_GROUPS = [
    "amino_acids",
    "lipids_fatty_acids_sterols",
    "vitamins_cofactors",
]
SELECTED_FUSARIUM_SINGLE_SUPPLEMENTS = ["UF03044_E"]  # guaiacol


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S06 = load_module(SCRIPT_06, "script06_for_script15")
S08 = load_module(SCRIPT_08, "script08_for_script15")


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


def quiet_memote_state(model: cobra.Model) -> dict[str, object]:
    # MEMOTE caches biomass/internal-reaction detection by model object.
    # This script creates and annotates a direct biomass reaction in memory, so
    # clear the helper caches before each consistency call to avoid stale
    # biomass-detection state.
    for helper in (con_helpers.find_biomass_reaction, con_helpers.get_internals):
        cache_clear = getattr(helper, "cache_clear", None)
        if callable(cache_clear):
            cache_clear()
    with open(os.devnull, "w", encoding="utf-8") as sink:
        with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
            passed = bool(con.check_stoichiometric_consistency(model))
            unconserved = [] if passed else sorted(m.id for m in con.find_unconserved_metabolites(model))
    return {
        "memote_pass": passed,
        "unconserved_count": len(unconserved),
        "unconserved_metabolites": ";".join(unconserved),
    }


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


def exchange_has_carbon(reaction: cobra.Reaction) -> bool:
    return any(met.elements.get("C", 0) for met in reaction.metabolites)


def carbon_free_medium(model: cobra.Model, medium: dict[str, float]) -> dict[str, float]:
    result: dict[str, float] = {}
    for reaction_id, value in medium.items():
        if not model.reactions.has_id(reaction_id):
            continue
        reaction = model.reactions.get_by_id(reaction_id)
        if exchange_has_carbon(reaction):
            continue
        result[reaction_id] = float(value)
    return result


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
    }


def apply_uf01339_correction(model: cobra.Model) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if "UF01339_C" not in model.reactions:
        return rows
    rows.append(reaction_snapshot(model, "UF01339_C", "before_UF01339_C_correction"))
    model.reactions.get_by_id("UF01339_C").add_metabolites(
        {
            model.metabolites.get_by_id("h2o_c"): 1.0,
            model.metabolites.get_by_id("h_c"): -1.0,
            model.metabolites.get_by_id("ppi_c"): -1.0,
        },
        combine=True,
    )
    rows.append(reaction_snapshot(model, "UF01339_C", "after_UF01339_C_correction"))
    return rows


def remove_batch_e_overlap(model: cobra.Model) -> list[str]:
    present = [rid for rid in BIPOLARIS_BATCH_E_CONTEXT if rid in model.reactions]
    if present:
        model.remove_reactions(present, remove_orphans=False)
    return present


def build_shared_plus_37_selected_medium(model: cobra.Model) -> tuple[dict[str, float], list[dict[str, object]]]:
    """Reconstruct the Script 12 selected medium for this in-memory model."""
    equivalents = S06.read_shared_medium_equivalents()
    all_shared_medium = {
        str(row["carvefungi_exchange"]): float(row["full_flux"])
        for row in equivalents
        if str(row["carvefungi_exchange"])
    }
    missing_interfaces = sorted(set(all_shared_medium).difference(rid for rid in all_shared_medium if model.reactions.has_id(rid)))
    interface_rows = S06.add_all_missing_temporary_interfaces(model, missing_interfaces)

    shared_plus_interfaces = {
        reaction_id: flux
        for reaction_id, flux in all_shared_medium.items()
        if model.reactions.has_id(reaction_id)
    }

    native_medium = {
        reaction_id: float(flux)
        for reaction_id, flux in cobra.io.read_sbml_model(str(MODEL)).medium.items()
        if model.reactions.has_id(reaction_id)
    }
    native_only = {
        reaction_id: flux
        for reaction_id, flux in native_medium.items()
        if reaction_id not in shared_plus_interfaces
    }
    selected_additions: dict[str, float] = {}
    supplement_rows: list[dict[str, object]] = []
    for reaction_id, flux in sorted(native_only.items()):
        meta = S08.exchange_metadata(model, reaction_id)
        group = S08.classify_exchange(meta)
        selected = (
            group in SELECTED_FUSARIUM_SUPPLEMENT_GROUPS
            or reaction_id in SELECTED_FUSARIUM_SINGLE_SUPPLEMENTS
        )
        if not selected:
            continue
        selected_additions[reaction_id] = flux
        supplement_rows.append({
            **meta,
            "native_flux": flux,
            "selection_reason": group if group in SELECTED_FUSARIUM_SUPPLEMENT_GROUPS else "single_guaiacol_rescue_trigger",
        })

    medium = {**shared_plus_interfaces, **selected_additions}
    for row in interface_rows:
        supplement_rows.append({
            "reaction_id": row["exchange_id"],
            "name": "temporary shared interface",
            "selection_reason": "shared_medium_interface",
            "native_flux": "",
            "added_reactions": row["added_reactions"],
            "reason": row["reason"],
        })
    return medium, supplement_rows


def biomass_precursor_producibility(
    model: cobra.Model,
    medium: dict[str, float],
    biomass_id: str,
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
    state = quiet_memote_state(model)
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
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mother_before = sha256(MODEL)
    bipolaris_before = sha256(S06.BIPOLARIS_FINAL)

    mother = cobra.io.read_sbml_model(str(MODEL))
    native_medium = dict(mother.medium)

    bipolaris = cobra.io.read_sbml_model(str(S06.BIPOLARIS_FINAL))
    bipolaris_biomass = S06.find_bipolaris_harmonized_biomass(bipolaris)

    baseline_raw = mother.copy()
    baseline_uf01339 = mother.copy()
    uf01339_audit = apply_uf01339_correction(baseline_uf01339)

    harmonized = mother.copy()
    direct_audit: list[dict[str, object]] = []
    direct_audit.extend(apply_uf01339_correction(harmonized))
    removed = remove_batch_e_overlap(harmonized)
    biomass_info = S06.add_harmonized_biomass(harmonized, bipolaris_biomass)
    biomass_id = str(biomass_info["harmonized_biomass_id"])
    # Script 06's biomass helper was originally for growth diagnostics. For
    # MEMOTE-style consistency diagnostics, explicitly mark the direct objective
    # as biomass so it is excluded from the internal stoichiometric-consistency
    # matrix, as MEMOTE does for biomass reactions.
    harmonized.reactions.get_by_id(biomass_id).annotation["sbo"] = "SBO:0000629"
    harmonized.reactions.get_by_id(biomass_id).name = (
        "biomass Fusarium direct 52-component harmonized CarveFungi biomass"
    )
    direct_audit.append(reaction_snapshot(harmonized, biomass_id, "after_direct_52_biomass_addition"))

    selected_medium, supplement_rows = build_shared_plus_37_selected_medium(harmonized.copy())

    condition_rows = [
        evaluate_condition(baseline_raw, "raw_mother_native_medium", native_medium),
        evaluate_condition(baseline_uf01339, "UF01339_C_only_native_medium", native_medium),
        evaluate_condition(harmonized, "UF01339_plus_direct_52_biomass_remove_Batch_E_overlap_native_medium", native_medium),
        evaluate_condition(harmonized, "UF01339_plus_direct_52_biomass_remove_Batch_E_overlap_selected_shared_plus_37_medium", selected_medium),
    ]

    native_precursors = biomass_precursor_producibility(harmonized, native_medium, biomass_id)
    selected_precursors = biomass_precursor_producibility(harmonized, selected_medium, biomass_id)
    for row in native_precursors:
        row["medium_condition"] = "native_medium"
    for row in selected_precursors:
        row["medium_condition"] = "selected_shared_plus_37_medium"
    precursor_rows = native_precursors + selected_precursors

    removed_rows = []
    for rid in removed:
        removed_rows.append(reaction_snapshot(mother, rid, "removed_Batch_E_overlap_reaction_before_removal"))

    write_csv(OUT / "condition_results.csv", condition_rows)
    write_csv(OUT / "removed_Batch_E_overlap_reactions.csv", removed_rows)
    write_csv(OUT / "reaction_edit_audit.csv", direct_audit + uf01339_audit)
    write_csv(OUT / "selected_shared_plus_37_medium_rows.csv", supplement_rows)
    write_csv(OUT / "biomass_precursor_producibility.csv", precursor_rows)

    precursor_summary = []
    for medium_condition, rows in [
        ("native_medium", native_precursors),
        ("selected_shared_plus_37_medium", selected_precursors),
    ]:
        precursor_summary.append({
            "medium_condition": medium_condition,
            "precursors_tested": len(rows),
            "precursors_producible": sum(bool(row["producible"]) for row in rows),
        })
    write_csv(OUT / "biomass_precursor_producibility_summary.csv", precursor_summary)

    mother_after = sha256(MODEL)
    bipolaris_after = sha256(S06.BIPOLARIS_FINAL)
    final_native = condition_rows[2]
    final_selected = condition_rows[3]
    native_prec_summary = precursor_summary[0]
    selected_prec_summary = precursor_summary[1]

    pass_gate = {
        "native_growth_positive": float(final_native["growth"]) > ZERO,
        "native_carbon_free_zero": abs(float(final_native["carbon_free_growth"])) <= 1e-8,
        "selected_growth_positive": float(final_selected["growth"]) > ZERO,
        "selected_carbon_free_zero": abs(float(final_selected["carbon_free_growth"])) <= 1e-8,
        "native_precursors_all_producible": native_prec_summary["precursors_tested"] == native_prec_summary["precursors_producible"],
        "selected_precursors_all_producible": selected_prec_summary["precursors_tested"] == selected_prec_summary["precursors_producible"],
        "direct_background_native_unconserved_is_7": int(final_native["unconserved_count"]) == 7,
        "direct_background_selected_unconserved_is_7": int(final_selected["unconserved_count"]) == 7,
        "mother_unchanged": mother_before == mother_after,
        "bipolaris_reference_unchanged": bipolaris_before == bipolaris_after,
    }

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "mother_model": str(MODEL),
        "mother_sha256_before": mother_before,
        "mother_sha256_after": mother_after,
        "mother_unchanged": mother_before == mother_after,
        "bipolaris_reference": str(S06.BIPOLARIS_FINAL),
        "bipolaris_reference_sha256_before": bipolaris_before,
        "bipolaris_reference_sha256_after": bipolaris_after,
        "bipolaris_reference_unchanged": bipolaris_before == bipolaris_after,
        "removed_Batch_E_overlap_count": len(removed),
        "removed_Batch_E_overlap_reactions": removed,
        "harmonized_biomass": biomass_info,
        "conditions": condition_rows,
        "precursor_summary": precursor_summary,
        "pass_gate": pass_gate,
        "all_pass_gate_checks_true": all(pass_gate.values()),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 15 compact summary")
    print(f"mother_unchanged: {mother_before == mother_after}")
    print(f"bipolaris_reference_unchanged: {bipolaris_before == bipolaris_after}")
    print(f"removed_Batch_E_overlap_count: {len(removed)}")
    print(f"removed_Batch_E_overlap_reactions: {';'.join(removed)}")
    print(f"raw_unconserved_count: {condition_rows[0]['unconserved_count']}")
    print(f"UF01339_only_unconserved_count: {condition_rows[1]['unconserved_count']}")
    print(f"direct_background_unconserved_count: {final_native['unconserved_count']}")
    print(f"direct_background_unconserved_metabolites: {final_native['unconserved_metabolites']}")
    print(f"direct_background_native_growth: {final_native['growth']}")
    print(f"direct_background_native_carbon_free_growth: {final_native['carbon_free_growth']}")
    print(f"direct_background_selected_medium_growth: {final_selected['growth']}")
    print(f"direct_background_selected_medium_carbon_free_growth: {final_selected['carbon_free_growth']}")
    print(
        "native_biomass_precursors_producible: "
        f"{native_prec_summary['precursors_producible']}/{native_prec_summary['precursors_tested']}"
    )
    print(
        "selected_medium_biomass_precursors_producible: "
        f"{selected_prec_summary['precursors_producible']}/{selected_prec_summary['precursors_tested']}"
    )
    print(f"all_pass_gate_checks_true: {all(pass_gate.values())}")
    print(f"outputs: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
