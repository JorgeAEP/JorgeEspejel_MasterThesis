#!/usr/bin/env python3
"""Script 67: exhaustive b5 minimal-inconsistency and hitting-set search.

Diagnostic only. The mother SBML is never modified and no SBML is saved.

Script 66 localized the residual MEMOTE failure to small interactions inside
the 15-reaction cytochrome-b5 subsystem. No single b5 reaction recreated the
failure, but three pairs did, all involving UF00430_CR. Every leave-one-out b5
condition still failed, implying additional redundant inconsistent supports.

This script therefore treats the 15 b5 reactions as a finite mathematical
subsystem:

1. Build the current viable seven-proton background on model.copy().
2. Apply old_C explicit-H neutralization on that copy.
3. Remove all 15 b5 reactions to confirm the clean passing control.
4. Enumerate b5 add-back subsets with monotonic pruning and collect all
   inclusion-minimal inconsistent subsets.
5. Solve the minimum hitting-set problem over those minimal bad subsets.
6. Test candidate b5 deletion sets in the full H-neutralized b5-containing
   background, ranking MEMOTE pass and growth retention.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import itertools
import json
import os
from pathlib import Path
from typing import Iterable

import cobra
import pandas as pd
from cobra import Metabolite, Reaction
from memote.support import consistency as con
from memote.support import consistency_helpers as con_helpers


ROOT = Path("/home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme")
MODEL = ROOT / "04_models/Bipolaris_sorokiniana_carvefungi.xml"
SCRIPT25_AUDIT = ROOT / "06_reports/25_heavy_element_unbalanced_reaction_audit/post_UF01339_internal_reaction_balance_audit.csv"
SCRIPT62_INDIVIDUAL_INVENTORY = ROOT / "06_reports/62_old_C_individual_addback_from_clean_background/old_C_addback_reaction_inventory.csv"
RELAXED_BIOMASS_MAPPING = ROOT / "06_reports/relaxed_curated_biomass_variant/relaxed_biomass_mapping.csv"
SCRIPT39_TESTS = ROOT / "06_reports/39_audit_batch_E_biomass_macromolecule_polymer_pools/submodule_and_cumulative_ablation_tests.csv"
MEDIUM_MAP = ROOT / "06_reports/micom_carvefungi_emma_full_phosphate/reports/full_carvefungi_medium_mapping.csv"
OUT = ROOT / "06_reports/67_b5_minimal_inconsistency_hitting_set"

UF01339 = ["UF01339_C", "UF01339_m", "UF01339_x"]
DIRECT_BIOMASS_ID = "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF"
E_CONTEXT_CONDITION = "cumulative_through_E08_cell_wall_carbohydrate_pools"
RESIDUAL_PROTONS = {"h_c", "h_e", "h_m", "h_n", "h_r", "h_x", "hi_m"}
ZERO = 1e-9
DEFAULT_MAX_SUBSET_SIZE = 3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def set_solver(model: cobra.Model) -> None:
    if "gurobi" in cobra.util.solver.solvers:
        model.solver = "gurobi"
    elif "glpk" in cobra.util.solver.solvers:
        model.solver = "glpk"
    try:
        model.solver.configuration.verbosity = 0
    except Exception:
        pass


def met_id(base: str, rid: str) -> str:
    return f"{base}_{rid.rsplit('_', 1)[1].lower()}"


def apply_uf01339_repair(model: cobra.Model) -> None:
    for rid in UF01339:
        if rid not in model.reactions:
            raise SystemExit(f"Expected reaction missing: {rid}")
        rxn = model.reactions.get_by_id(rid)
        desired = {
            model.metabolites.get_by_id(met_id("datp", rid)): -1.0,
            model.metabolites.get_by_id(met_id("h2o", rid)): -1.0,
            model.metabolites.get_by_id(met_id("damp", rid)): 1.0,
            model.metabolites.get_by_id(met_id("h", rid)): 1.0,
            model.metabolites.get_by_id(met_id("ppi", rid)): 1.0,
        }
        rxn.add_metabolites({met: -coeff for met, coeff in rxn.metabolites.items()})
        rxn.add_metabolites(desired)


def load_e_context_reactions() -> list[str]:
    table = pd.read_csv(SCRIPT39_TESTS)
    match = table.loc[table["condition"].astype(str) == E_CONTEXT_CONDITION]
    if match.empty:
        raise SystemExit(f"Missing Script 39 condition: {E_CONTEXT_CONDITION}")
    return sorted(x for x in str(match.iloc[0]["removed_reactions"]).split(";") if x and x != "nan")


def remove_reactions_if_present(model: cobra.Model, rids: Iterable[str]) -> list[str]:
    present = [rid for rid in sorted(set(rids)) if rid in model.reactions]
    if present:
        model.remove_reactions(present, remove_orphans=False)
    return present


def add_direct_biomass(model: cobra.Model) -> None:
    table = pd.read_csv(RELAXED_BIOMASS_MAPPING)
    included = table.loc[table["status"].astype(str) == "included"]
    if included.empty:
        raise SystemExit("No included direct-biomass rows found")
    missing = [str(mid) for mid in included["carve_id"] if str(mid) not in model.metabolites]
    if missing:
        raise SystemExit(f"Missing direct-biomass metabolites: {missing}")
    rxn = Reaction(DIRECT_BIOMASS_ID)
    rxn.name = "Direct near-curated Bipolaris biomass for Script 67 tests"
    rxn.bounds = (0.0, 1000.0)
    rxn.annotation = {"sbo": "SBO:0000629", "comment": "Diagnostic copy-only biomass."}
    rxn.add_metabolites({
        model.metabolites.get_by_id(str(row["carve_id"])): float(row["coefficient"])
        for _, row in included.iterrows()
    })
    model.add_reactions([rxn])
    model.objective = rxn


def build_harmonized_copy(mother: cobra.Model) -> cobra.Model:
    work = mother.copy()
    set_solver(work)
    apply_uf01339_repair(work)
    remove_reactions_if_present(work, load_e_context_reactions())
    add_direct_biomass(work)
    return work


def load_medium(model: cobra.Model) -> dict[str, float]:
    medium: dict[str, float] = {}
    if MEDIUM_MAP.exists():
        table = pd.read_csv(MEDIUM_MAP)
        for _, row in table.iterrows():
            rid = str(row.get("carve_model_exchange", ""))
            if rid in model.reactions:
                medium[rid] = max(medium.get(rid, 0.0), float(row.get("flux", 0.0)))
    for rid in ["UF02549_E", "UF03474_E", "UF03376_E"]:
        if rid in model.reactions:
            medium[rid] = max(medium.get(rid, 0.0), 1000.0)
    return medium


def carbon_free_medium(model: cobra.Model, medium: dict[str, float]) -> dict[str, float]:
    result = dict(medium)
    for rid in result:
        if rid not in model.reactions:
            continue
        mets = list(model.reactions.get_by_id(rid).metabolites)
        if mets and mets[0].elements.get("C", 0) > 0:
            result[rid] = 0.0
    return result


def optimize_value(model: cobra.Model, medium: dict[str, float]) -> tuple[str, float | None]:
    set_solver(model)
    try:
        with model:
            model.medium = {rid: value for rid, value in medium.items() if rid in model.reactions}
            with open(os.devnull, "w", encoding="utf-8") as sink:
                with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
                    solution = model.optimize()
            return solution.status, float(solution.objective_value or 0.0)
    except Exception as exc:
        return f"error: {exc!r}", None


def memote_state(model: cobra.Model) -> dict:
    set_solver(model)
    with open(os.devnull, "w", encoding="utf-8") as sink:
        with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
            passed = bool(con.check_stoichiometric_consistency(model))
            metabolites = [] if passed else sorted(m.id for m in con.find_unconserved_metabolites(model))
    return {
        "memote_pass": passed,
        "unconserved_count": len(metabolites),
        "unconserved_metabolites": metabolites,
    }


def memote_internal_ids(model: cobra.Model) -> set[str]:
    return {rxn.id for rxn in con_helpers.get_internals(model)}


def to_bool_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def load_old_C(base: cobra.Model) -> set[str]:
    audit = pd.read_csv(SCRIPT25_AUDIT)
    rids = audit["reaction_id"].astype(str)
    current_internal = memote_internal_ids(base)
    return set(rids[to_bool_series(audit["H_charge_only_signal"])]) & current_internal


def parse_mass_balance(text: str) -> dict[str, float]:
    if pd.isna(text) or not str(text).strip():
        return {}
    value = ast.literal_eval(str(text))
    return {str(k): float(v) for k, v in value.items()}


def load_old_C_inventory(old_C: set[str]) -> pd.DataFrame:
    if not SCRIPT62_INDIVIDUAL_INVENTORY.exists():
        raise SystemExit(f"Missing Script 62 individual inventory: {SCRIPT62_INDIVIDUAL_INVENTORY}")
    table = pd.read_csv(SCRIPT62_INDIVIDUAL_INVENTORY)
    table = table.loc[table["reaction_id"].astype(str).isin(old_C)].copy()
    table["mass_balance_dict"] = table["mass_balance"].apply(parse_mass_balance)
    return table


def get_met(model: cobra.Model, mid: str) -> Metabolite | None:
    if mid not in model.metabolites:
        return None
    return model.metabolites.get_by_id(mid)


def choose_proton(model: cobra.Model, rxn: Reaction) -> Metabolite | None:
    candidates = [met for met in rxn.metabolites if met.id in RESIDUAL_PROTONS]
    if candidates:
        candidates = sorted(candidates, key=lambda met: (met.id == "hi_m", met.id))
        return candidates[0]
    compartments = sorted(rxn.compartments)
    if len(compartments) == 1:
        return get_met(model, f"h_{compartments[0]}")
    return None


def neutralize_reaction_H(model: cobra.Model, rid: str) -> dict:
    if rid not in model.reactions:
        return {"reaction_id": rid, "status": "missing", "before": "", "after": "", "note": ""}
    rxn = model.reactions.get_by_id(rid)
    before = rxn.check_mass_balance()
    h_imbalance = float(before.get("H", 0.0))
    if abs(h_imbalance) <= ZERO:
        return {
            "reaction_id": rid,
            "status": "no_H_imbalance",
            "before": json.dumps(before, sort_keys=True),
            "after": json.dumps(before, sort_keys=True),
            "note": "No H term in local mass balance.",
        }
    proton = choose_proton(model, rxn)
    if proton is None:
        return {
            "reaction_id": rid,
            "status": "skipped_no_proton_available",
            "before": json.dumps(before, sort_keys=True),
            "after": json.dumps(before, sort_keys=True),
            "note": "No suitable proton metabolite found.",
        }
    rxn.add_metabolites({proton: -h_imbalance})
    after = rxn.check_mass_balance()
    return {
        "reaction_id": rid,
        "status": "changed",
        "before": json.dumps(before, sort_keys=True),
        "after": json.dumps(after, sort_keys=True),
        "note": f"Adjusted {proton.id} coefficient by {-h_imbalance:g}.",
        "after_equation": rxn.reaction,
    }


def apply_old_C_explicit_H_neutralization(model: cobra.Model, inventory: pd.DataFrame, old_C: set[str]) -> list[dict]:
    explicit_H = set(inventory.loc[
        inventory["mass_balance_dict"].apply(lambda mb: abs(float(mb.get("H", 0.0))) > ZERO),
        "reaction_id",
    ].astype(str)) & old_C
    edits = []
    for rid in sorted(explicit_H):
        edits.append(neutralize_reaction_H(model, rid))
    return edits


def b5_reaction_ids(inventory: pd.DataFrame) -> set[str]:
    marker = inventory["missing_formula_metabolites"].fillna("").astype(str)
    return set(inventory.loc[marker.str.contains("ficytb5|focytb5"), "reaction_id"].astype(str))


def evaluate_addback(
    h_neutralized_full: cobra.Model,
    b5_rids: set[str],
    add_back: frozenset[str],
    baseline_growth: float | None,
    baseline_cf: float | None,
    with_growth: bool = False,
) -> dict:
    work = h_neutralized_full.copy()
    removed = remove_reactions_if_present(work, b5_rids - set(add_back))
    state = memote_state(work)
    row = {
        "n_b5_add_back": len(add_back),
        "b5_add_back_reactions": ";".join(sorted(add_back)),
        "n_b5_removed": len(removed),
        "b5_removed_reactions": ";".join(removed),
        **state,
        "unconserved_metabolites": ";".join(state["unconserved_metabolites"]),
    }
    if with_growth:
        medium = load_medium(work)
        growth_status, growth = optimize_value(work, medium)
        cf_status, cf_growth = optimize_value(work, carbon_free_medium(work, medium))
        row.update({
            "growth_status": growth_status,
            "growth": growth,
            "growth_delta_vs_original_baseline": growth - baseline_growth if growth is not None and baseline_growth is not None else None,
            "growth_ratio_vs_original_baseline": growth / baseline_growth if growth is not None and baseline_growth and abs(baseline_growth) > ZERO else None,
            "carbon_free_status": cf_status,
            "carbon_free_growth": cf_growth,
            "carbon_free_delta_vs_original_baseline": cf_growth - baseline_cf if cf_growth is not None and baseline_cf is not None else None,
        })
    return row


def enumerate_minimal_bad_subsets(
    h_neutralized_full: cobra.Model,
    b5_rids: set[str],
    baseline_growth: float | None,
    baseline_cf: float | None,
    max_subset_size: int,
) -> tuple[list[frozenset[str]], pd.DataFrame]:
    ordered = sorted(b5_rids)
    minimal_bad: list[frozenset[str]] = []
    rows: list[dict] = []
    memote_tests_done = 0

    for size in range(min(len(ordered), max_subset_size) + 1):
        print(f"Enumerating b5 add-back subsets of size {size} ...", flush=True)
        for combo in itertools.combinations(ordered, size):
            subset = frozenset(combo)
            containing_known_bad = [bad for bad in minimal_bad if bad.issubset(subset)]
            if containing_known_bad:
                rows.append({
                    "subset_size": size,
                    "b5_add_back_reactions": ";".join(sorted(subset)),
                    "tested_memote": False,
                    "skipped_reason": "contains_known_minimal_bad_subset",
                    "witness_minimal_bad_subset": ";".join(sorted(containing_known_bad[0])),
                    "memote_pass": False,
                    "unconserved_count": None,
                    "unconserved_metabolites": "",
                    "is_inclusion_minimal_bad": False,
                })
                continue

            result = evaluate_addback(
                h_neutralized_full,
                b5_rids,
                subset,
                baseline_growth,
                baseline_cf,
                with_growth=False,
            )
            memote_tests_done += 1
            if memote_tests_done % 500 == 0:
                print(
                    f"  MEMOTE subset tests performed: {memote_tests_done}; "
                    f"minimal bad subsets found: {len(minimal_bad)}; "
                    f"rows recorded: {len(rows)}",
                    flush=True,
                )
            is_bad = not bool(result["memote_pass"])
            if is_bad:
                minimal_bad.append(subset)
                print(
                    "  Found minimal inconsistent subset: "
                    f"size={size}; reactions={';'.join(sorted(subset))}",
                    flush=True,
                )
            rows.append({
                "subset_size": size,
                "b5_add_back_reactions": ";".join(sorted(subset)),
                "tested_memote": True,
                "skipped_reason": "",
                "witness_minimal_bad_subset": "",
                "is_inclusion_minimal_bad": is_bad,
                **result,
            })

    return minimal_bad, pd.DataFrame(rows)


def configured_max_subset_size(b5_count: int) -> int:
    raw = os.environ.get("SCRIPT67_MAX_SUBSET_SIZE", str(DEFAULT_MAX_SUBSET_SIZE)).strip()
    if raw.lower() in {"all", "exhaustive", "full"}:
        return b5_count
    try:
        value = int(raw)
    except ValueError as exc:
        raise SystemExit(
            "SCRIPT67_MAX_SUBSET_SIZE must be an integer or one of: all, exhaustive, full"
        ) from exc
    if value < 0:
        raise SystemExit("SCRIPT67_MAX_SUBSET_SIZE must be non-negative")
    return min(value, b5_count)


def hits_all(candidate: frozenset[str], bad_subsets: list[frozenset[str]]) -> bool:
    return all(candidate & bad for bad in bad_subsets)


def find_minimum_hitting_sets(b5_rids: set[str], bad_subsets: list[frozenset[str]]) -> list[frozenset[str]]:
    ordered = sorted(b5_rids)
    for size in range(len(ordered) + 1):
        hits = [frozenset(combo) for combo in itertools.combinations(ordered, size) if hits_all(frozenset(combo), bad_subsets)]
        if hits:
            return hits
    return []


def reaction_audit(model: cobra.Model, rids: Iterable[str]) -> pd.DataFrame:
    rows = []
    for rid in sorted(rids):
        if rid not in model.reactions:
            rows.append({"reaction_id": rid, "status": "missing"})
            continue
        rxn = model.reactions.get_by_id(rid)
        rows.append({
            "reaction_id": rid,
            "status": "present",
            "reaction_name": rxn.name,
            "reaction_equation": rxn.reaction,
            "bounds": str(rxn.bounds),
            "compartments": ";".join(sorted(rxn.compartments)),
            "residual_protons": ";".join(sorted(m.id for m in rxn.metabolites if m.id in RESIDUAL_PROTONS)),
            "metabolites": "; ".join(
                f"{met.id}:{coeff:g}" for met, coeff in sorted(rxn.metabolites.items(), key=lambda item: item[0].id)
            ),
        })
    return pd.DataFrame(rows)


def test_deletion_hitting_sets(
    h_neutralized_full: cobra.Model,
    b5_rids: set[str],
    hitting_sets: list[frozenset[str]],
    baseline_growth: float | None,
    baseline_cf: float | None,
) -> pd.DataFrame:
    rows = []
    for index, deletion_set in enumerate(hitting_sets, start=1):
        kept = frozenset(b5_rids - set(deletion_set))
        result = evaluate_addback(
            h_neutralized_full,
            b5_rids,
            kept,
            baseline_growth,
            baseline_cf,
            with_growth=True,
        )
        rows.append({
            "candidate": f"hitting_set_{index:03d}",
            "n_b5_deleted": len(deletion_set),
            "b5_deleted_reactions": ";".join(sorted(deletion_set)),
            "n_b5_kept": len(kept),
            "b5_kept_reactions": ";".join(sorted(kept)),
            **result,
            "unconserved_metabolites": result["unconserved_metabolites"],
            "interpretation": "COPY_ONLY_B5_MINIMUM_HITTING_SET_DELETION_TEST",
        })
    results = pd.DataFrame(rows)
    if not results.empty:
        results = results.sort_values(
            ["memote_pass", "unconserved_count", "growth", "n_b5_deleted", "candidate"],
            ascending=[False, True, False, True, True],
        )
    return results


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    before_checksum = sha256(MODEL)
    mother = cobra.io.read_sbml_model(str(MODEL))
    base = build_harmonized_copy(mother)
    old_C = load_old_C(base)
    inventory = load_old_C_inventory(old_C)
    b5_rids = b5_reaction_ids(inventory)

    original_state = memote_state(base)
    medium = load_medium(base)
    baseline_status, baseline_growth = optimize_value(base, medium)
    baseline_cf_status, baseline_cf = optimize_value(base, carbon_free_medium(base, medium))

    h_neutralized_full = base.copy()
    edits = apply_old_C_explicit_H_neutralization(h_neutralized_full, inventory, old_C)
    h_neutralized_full_state = memote_state(h_neutralized_full)

    clean = h_neutralized_full.copy()
    clean_removed = remove_reactions_if_present(clean, b5_rids)
    clean_state = memote_state(clean)

    max_subset_size = configured_max_subset_size(len(b5_rids))
    print(
        "Script 67 subset search limit: "
        f"size <= {max_subset_size} "
        "(set SCRIPT67_MAX_SUBSET_SIZE=all for full exhaustive search)",
        flush=True,
    )
    minimal_bad, subset_tests = enumerate_minimal_bad_subsets(
        h_neutralized_full,
        b5_rids,
        baseline_growth,
        baseline_cf,
        max_subset_size,
    )
    minimal_bad_rows = []
    for index, subset in enumerate(sorted(minimal_bad, key=lambda s: (len(s), sorted(s))), start=1):
        minimal_bad_rows.append({
            "minimal_bad_subset": f"minimal_bad_{index:03d}",
            "size": len(subset),
            "b5_add_back_reactions": ";".join(sorted(subset)),
        })
    minimal_bad_table = pd.DataFrame(minimal_bad_rows)

    hitting_sets = find_minimum_hitting_sets(b5_rids, minimal_bad)
    hitting_rows = []
    for index, hit in enumerate(sorted(hitting_sets, key=lambda s: (len(s), sorted(s))), start=1):
        hitting_rows.append({
            "hitting_set": f"hitting_set_{index:03d}",
            "size": len(hit),
            "b5_deleted_reactions": ";".join(sorted(hit)),
        })
    hitting_table = pd.DataFrame(hitting_rows)
    deletion_tests = test_deletion_hitting_sets(
        h_neutralized_full,
        b5_rids,
        sorted(hitting_sets, key=lambda s: (len(s), sorted(s))),
        baseline_growth,
        baseline_cf,
    )

    b5_audit = reaction_audit(h_neutralized_full, b5_rids)
    edit_table = pd.DataFrame(edits)

    subset_tests.to_csv(OUT / "all_b5_subset_tests_with_pruning.csv", index=False)
    minimal_bad_table.to_csv(OUT / "minimal_inconsistent_b5_subsets.csv", index=False)
    hitting_table.to_csv(OUT / "minimum_b5_hitting_sets.csv", index=False)
    deletion_tests.to_csv(OUT / "minimum_hitting_set_deletion_tests.csv", index=False)
    b5_audit.to_csv(OUT / "b5_reaction_audit.csv", index=False)
    edit_table.to_csv(OUT / "old_C_explicit_H_neutralization_edit_audit.csv", index=False)
    inventory.to_csv(OUT / "old_C_inventory_used.csv", index=False)

    with pd.ExcelWriter(OUT / "b5_minimal_inconsistency_hitting_set.xlsx", engine="xlsxwriter") as writer:
        subset_tests.to_excel(writer, sheet_name="all_subset_tests", index=False)
        minimal_bad_table.to_excel(writer, sheet_name="minimal_bad_subsets", index=False)
        hitting_table.to_excel(writer, sheet_name="minimum_hitting_sets", index=False)
        deletion_tests.to_excel(writer, sheet_name="deletion_tests", index=False)
        b5_audit.to_excel(writer, sheet_name="b5_reaction_audit", index=False)
        edit_table.to_excel(writer, sheet_name="H_edit_audit", index=False)
        wrap = writer.book.add_format({"text_wrap": True, "valign": "top"})
        for sheet in writer.sheets.values():
            sheet.freeze_panes(1, 0)
            sheet.set_column(0, 10, 32)
            sheet.set_column(11, 80, 44, wrap)

    after_checksum = sha256(MODEL)
    top_deletion_tests = deletion_tests.head(20).to_dict(orient="records") if not deletion_tests.empty else []
    summary = {
        "mother_model": str(MODEL),
        "mother_checksum_before": before_checksum,
        "mother_checksum_after": after_checksum,
        "mother_model_unchanged": before_checksum == after_checksum,
        "diagnostic_only_no_sbml_saved": True,
        "background": "UF01339_plus_direct_biomass_remove_Batch_E_context_17",
        "original_baseline": original_state,
        "original_baseline_growth_status": baseline_status,
        "original_baseline_growth": baseline_growth,
        "original_baseline_carbon_free_status": baseline_cf_status,
        "original_baseline_carbon_free_growth": baseline_cf,
        "old_C_size": len(old_C),
        "old_C_explicit_H_edits_changed": sum(1 for edit in edits if edit["status"] == "changed"),
        "b5_reaction_count": len(b5_rids),
        "b5_reactions": sorted(b5_rids),
        "h_neutralized_full_state": h_neutralized_full_state,
        "clean_background": {
            "b5_removed_reactions": clean_removed,
            **clean_state,
        },
        "subset_search": {
            "total_possible_subsets": 2 ** len(b5_rids),
            "max_subset_size_tested": max_subset_size,
            "is_full_exhaustive_search": max_subset_size == len(b5_rids),
            "rows_recorded": len(subset_tests),
            "memote_tests_performed": int(subset_tests["tested_memote"].sum()),
            "subsets_pruned": int((~subset_tests["tested_memote"]).sum()),
            "minimal_inconsistent_subset_count": len(minimal_bad),
            "minimal_inconsistent_subsets": minimal_bad_table.to_dict(orient="records"),
        },
        "minimum_hitting_set_size": int(hitting_table["size"].min()) if not hitting_table.empty else None,
        "minimum_hitting_set_count": len(hitting_sets),
        "minimum_hitting_sets": hitting_table.to_dict(orient="records"),
        "deletion_tests_top_ranked": top_deletion_tests,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    deletion_display_cols = [
        "candidate", "n_b5_deleted", "b5_deleted_reactions", "memote_pass",
        "unconserved_count", "unconserved_metabolites", "growth",
        "growth_ratio_vs_original_baseline", "carbon_free_growth",
    ]
    deletion_display = deletion_tests[deletion_display_cols] if not deletion_tests.empty else pd.DataFrame(columns=deletion_display_cols)
    md = "\n".join([
        "# Script 67 b5 minimal-inconsistency hitting-set search", "",
        "Diagnostic copy-only exhaustive b5 subset and hitting-set tests. No SBML was saved.", "",
        "## Minimal inconsistent b5 subsets", "",
        f"Subset search limit: size <= {max_subset_size}. "
        "Set `SCRIPT67_MAX_SUBSET_SIZE=all` for a full exhaustive search.", "",
        "```csv",
        minimal_bad_table.to_csv(index=False).strip(),
        "```", "",
        "## Minimum hitting sets", "",
        "```csv",
        hitting_table.to_csv(index=False).strip(),
        "```", "",
        "## Deletion tests for minimum hitting sets", "",
        "```csv",
        deletion_display.to_csv(index=False).strip(),
        "```", "",
        "## Integrity", "",
        f"- mother checksum before: `{before_checksum}`",
        f"- mother checksum after: `{after_checksum}`",
        f"- mother model unchanged: `{before_checksum == after_checksum}`", "",
    ])
    (OUT / "summary.md").write_text(md, encoding="utf-8")

    print("Script 67 compact summary")
    print(f"original_baseline_unconserved_count: {original_state['unconserved_count']}")
    print(f"original_baseline_unconserved_metabolites: {';'.join(original_state['unconserved_metabolites'])}")
    print(f"original_baseline_growth: {baseline_growth}")
    print(f"original_baseline_carbon_free_growth: {baseline_cf}")
    print(f"old_C_size: {len(old_C)}")
    print(f"old_C_explicit_H_edits_changed: {sum(1 for edit in edits if edit['status'] == 'changed')}")
    print(f"b5_reaction_count: {len(b5_rids)}")
    print(f"H_neutralized_full_unconserved_count: {h_neutralized_full_state['unconserved_count']}")
    print(f"clean_background_pass: {clean_state['memote_pass']}")
    print(f"clean_background_unconserved_count: {clean_state['unconserved_count']}")
    print(f"total_possible_subsets: {2 ** len(b5_rids)}")
    print(f"max_subset_size_tested: {max_subset_size}")
    print(f"is_full_exhaustive_search: {max_subset_size == len(b5_rids)}")
    print(f"memote_tests_performed: {int(subset_tests['tested_memote'].sum())}")
    print(f"subsets_pruned: {int((~subset_tests['tested_memote']).sum())}")
    print(f"minimal_inconsistent_subset_count: {len(minimal_bad)}")
    if not minimal_bad_table.empty:
        print("minimal_inconsistent_subsets:")
        for _, row in minimal_bad_table.iterrows():
            print(f"  {row['minimal_bad_subset']} | size={row['size']} | add_back={row['b5_add_back_reactions']}")
    print(f"minimum_hitting_set_size: {summary['minimum_hitting_set_size']}")
    print(f"minimum_hitting_set_count: {len(hitting_sets)}")
    if not hitting_table.empty:
        print("minimum_hitting_sets:")
        for _, row in hitting_table.iterrows():
            print(f"  {row['hitting_set']} | size={row['size']} | delete={row['b5_deleted_reactions']}")
    passing = deletion_tests.loc[deletion_tests["memote_pass"]] if not deletion_tests.empty else pd.DataFrame()
    print(f"minimum_hitting_set_deletion_tests_passing_memote: {len(passing)}")
    if not passing.empty:
        best = passing.sort_values(["growth", "n_b5_deleted"], ascending=[False, True]).iloc[0]
        print(f"best_passing_deletion_set: {best['candidate']} | delete={best['b5_deleted_reactions']} | growth={best['growth']} | carbon_free_growth={best['carbon_free_growth']}")
    print(f"mother_model_unchanged: {before_checksum == after_checksum}")
    print(f"outputs: {OUT}")


if __name__ == "__main__":
    main()
