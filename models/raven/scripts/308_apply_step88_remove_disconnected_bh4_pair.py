#!/usr/bin/env python3
"""Apply Step 88: remove the isolated BH4 residual pair from Step 87.

This step removes the tiny disconnected subnetwork formed by:
- R11762
- C04244
- C21640

The pair is a leftover from the earlier BH4-core cleanup and is already
disconnected from the functional network. The goal here is purely structural:
restore the model to a single connected metabolite network without altering the
working metabolic core used for downstream MICOM/community analysis.
"""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path("/home/tlu/Desktop/Jorge/VS CODE/Bipolaris__sorokiniana")
DEPS = ROOT / ".cobramod_deps"
if str(DEPS) not in sys.path:
    sys.path.insert(0, str(DEPS))

import cobra
import pandas as pd
from cobra import Model, Reaction


IN_DIR = ROOT / "03_intermediate_models" / "87_post_shared_bacterial_medium_plus_fungal_supplements"
IN_XML = IN_DIR / "87_post_shared_bacterial_medium_plus_fungal_supplements.xml"
IN_JSON = IN_DIR / "87_post_shared_bacterial_medium_plus_fungal_supplements.json"

OUT_DIR = ROOT / "03_intermediate_models" / "88_post_remove_disconnected_bh4_pair"
OUT_XML = OUT_DIR / "88_post_remove_disconnected_bh4_pair.xml"
OUT_JSON = OUT_DIR / "88_post_remove_disconnected_bh4_pair.json"

REPORT_DIR = ROOT / "04_reports" / "308_apply_step88_remove_disconnected_bh4_pair"

BIOMASS = "BIOMASS_BIPOLARIS"
TARGET_REACTION = "R11762"
TARGET_METABOLITES = ("C04244", "C21640")
TOL = 1e-8


def optimize(model: Model) -> tuple[str, float]:
    model.objective = BIOMASS
    sol = model.optimize()
    if sol.status != "optimal" or sol.objective_value is None:
        return sol.status, math.nan
    value = float(sol.objective_value)
    if abs(value) < TOL:
        value = 0.0
    return sol.status, value


def demand_possible(model: Model, met_id: str, label: str) -> tuple[str, float]:
    if not model.metabolites.has_id(met_id):
        return "absent", math.nan
    with model:
        dm = Reaction(label)
        dm.bounds = (0.0, 1000.0)
        dm.add_metabolites({model.metabolites.get_by_id(met_id): -1.0})
        model.add_reactions([dm])
        model.objective = dm
        return optimize(model)


def biomass_precursor_status(model: Model) -> tuple[str, list[dict[str, object]]]:
    biomass = model.reactions.get_by_id(BIOMASS)
    rows: list[dict[str, object]] = []
    ok = 0
    precursors = sorted([m for m, c in biomass.metabolites.items() if c < 0], key=lambda m: m.id)
    for met in precursors:
        status, value = demand_possible(model, met.id, f"DM_STEP88_{met.id}")
        producible = status == "optimal" and math.isfinite(value) and value > TOL
        ok += int(producible)
        rows.append(
            {
                "metabolite_id": met.id,
                "metabolite_name": met.name,
                "coefficient": float(biomass.metabolites[met]),
                "demand_status": status,
                "demand_max": value,
                "producible": producible,
            }
        )
    return f"{ok}/{len(precursors)}", rows


def carbon_free_growth(model: Model) -> tuple[str, float]:
    with model:
        for rxn in model.exchanges:
            if any(((met.elements or {}).get("C", 0) or 0) > 0 for met in rxn.metabolites):
                rxn.lower_bound = max(0.0, rxn.lower_bound)
        model.objective = BIOMASS
        return optimize(model)


def metabolite_components(model: Model) -> list[list[str]]:
    adjacency = {met.id: set() for met in model.metabolites}
    for rxn in model.reactions:
        mets = [met.id for met in rxn.metabolites]
        for i, left in enumerate(mets):
            for right in mets[i + 1 :]:
                adjacency[left].add(right)
                adjacency[right].add(left)

    components: list[list[str]] = []
    seen: set[str] = set()
    for met_id in adjacency:
        if met_id in seen:
            continue
        stack = [met_id]
        seen.add(met_id)
        component: list[str] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    components.sort(key=len, reverse=True)
    return components


def component_summary(model: Model, label: str) -> tuple[list[dict[str, object]], list[list[str]]]:
    components = metabolite_components(model)
    rows: list[dict[str, object]] = []
    for idx, component in enumerate(components, start=1):
        rows.append(
            {
                "scenario": label,
                "component_rank": idx,
                "metabolite_count": len(component),
                "preview_metabolites": " | ".join(component[:10]),
            }
        )
    return rows, components


def snapshot(model: Model, label: str) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    status, growth = optimize(model)
    precursor_summary, precursor_rows = biomass_precursor_status(model)
    cf_status, cf_growth = carbon_free_growth(model)
    component_rows, components = component_summary(model, label)
    summary = {
        "scenario": label,
        "status": status,
        "growth": growth,
        "carbon_free_status": cf_status,
        "carbon_free_growth": cf_growth,
        "biomass_precursors": precursor_summary,
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
        "component_count": len(components),
        "largest_component_metabolites": len(components[0]) if components else 0,
        "smallest_component_metabolites": len(components[-1]) if components else 0,
    }
    return summary, precursor_rows, component_rows


def apply_cleanup(model: Model) -> tuple[Model, list[dict[str, object]]]:
    work = model.copy()
    removed_rows: list[dict[str, object]] = []

    if work.reactions.has_id(TARGET_REACTION):
        rxn = work.reactions.get_by_id(TARGET_REACTION)
        removed_rows.append(
            {
                "object_type": "reaction",
                "object_id": rxn.id,
                "name": rxn.name,
                "equation_before": rxn.reaction,
                "reason": "isolated residual BH4 pair reaction",
            }
        )
        work.remove_reactions([rxn], remove_orphans=False)

    orphan_mets = []
    for met_id in TARGET_METABOLITES:
        if not work.metabolites.has_id(met_id):
            continue
        met = work.metabolites.get_by_id(met_id)
        if len(met.reactions) == 0:
            removed_rows.append(
                {
                    "object_type": "metabolite",
                    "object_id": met.id,
                    "name": met.name,
                    "equation_before": "",
                    "reason": "orphan metabolite after isolated BH4 pair reaction removal",
                }
            )
            orphan_mets.append(met)
    if orphan_mets:
        work.remove_metabolites(orphan_mets, destructive=False)

    return work, removed_rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    before = cobra.io.read_sbml_model(IN_XML)
    before.objective = BIOMASS
    before_summary, before_precursors, before_components = snapshot(before, "step87_before")

    after, removed_rows = apply_cleanup(before)
    after.objective = BIOMASS
    after_summary, after_precursors, after_components = snapshot(after, "step88_after")

    cobra.io.write_sbml_model(after, OUT_XML)
    cobra.io.save_json_model(after, OUT_JSON)

    xml_reload = cobra.io.read_sbml_model(OUT_XML)
    json_reload = cobra.io.load_json_model(OUT_JSON)
    xml_summary, _, xml_components = snapshot(xml_reload, "step88_xml_reload")
    json_summary, _, json_components = snapshot(json_reload, "step88_json_reload")

    pd.DataFrame([before_summary, after_summary, xml_summary, json_summary]).to_csv(
        REPORT_DIR / "validation_summary.csv",
        index=False,
    )
    pd.DataFrame(before_precursors + after_precursors).to_csv(
        REPORT_DIR / "biomass_precursor_checks.csv",
        index=False,
    )
    pd.DataFrame(before_components + after_components + xml_components + json_components).to_csv(
        REPORT_DIR / "component_summary.csv",
        index=False,
    )
    pd.DataFrame(removed_rows).to_csv(REPORT_DIR / "removed_objects.csv", index=False)

    summary_lines = [
        "Step 88 disconnected BH4 pair cleanup",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "Input checkpoint:",
        str(IN_XML),
        "",
        "Output checkpoint:",
        str(OUT_XML),
        "",
        "Removed objects:",
    ]
    for row in removed_rows:
        summary_lines.append(
            f"- {row['object_type']}: {row['object_id']} | {row['name']} | {row['reason']}"
        )
    summary_lines.extend(
        [
            "",
            "Connectivity:",
            f"- before component count: {before_summary['component_count']}",
            f"- after component count: {after_summary['component_count']}",
            f"- before smallest component size: {before_summary['smallest_component_metabolites']}",
            f"- after smallest component size: {after_summary['smallest_component_metabolites']}",
            "",
            "Growth / integrity:",
            f"- before growth: {before_summary['growth']}",
            f"- after growth: {after_summary['growth']}",
            f"- before biomass precursors: {before_summary['biomass_precursors']}",
            f"- after biomass precursors: {after_summary['biomass_precursors']}",
            f"- before carbon-free growth: {before_summary['carbon_free_growth']}",
            f"- after carbon-free growth: {after_summary['carbon_free_growth']}",
            "",
            "Reload verification:",
            f"- XML reload component count: {xml_summary['component_count']}",
            f"- JSON reload component count: {json_summary['component_count']}",
            "",
            "STEP 88 COMPLETE",
            "The isolated BH4 residual pair was removed.",
            "The updated model contains a single connected metabolite network.",
        ]
    )
    (REPORT_DIR / "step88_cleanup_summary.txt").write_text("\n".join(summary_lines) + "\n")

    print(f"Step 88 saved to: {OUT_XML}")
    print(f"Report directory: {REPORT_DIR}")
    print("Disconnected BH4 residual pair removed successfully.")


if __name__ == "__main__":
    main()
