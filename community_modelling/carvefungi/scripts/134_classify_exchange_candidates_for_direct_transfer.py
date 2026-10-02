#!/usr/bin/env python3
"""Classify re-audited exchange candidates before biological interpretation.

This post-processing step uses the identifier-audit outputs from Script 133.
It separates background compounds supplied by the medium or representing
currency/boundary chemistry from non-medium chemical matches. It also records
whether donor and recipient used the same exported exchange identifier, the
same external metabolite identifier, or only a curated chemical crosswalk.

It does not claim transfer: a crosswalk match between different model interfaces
is weaker evidence than a common modelled external pool.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
AUDIT = PACKAGE / "reports/133_fungal_bacterial_exchange_identifier_audit"
MASTER = PACKAGE / "reports/96_fungal_species_syncom_factorial_panel/summaries/master_run_index.csv"
CROSSWALK = PACKAGE / "reports/88_dual_namespace_carvefungi_media/dual_namespace_media_audit_long.csv"
DEFAULT_OUT = PACKAGE / "reports/134_exchange_candidate_classification"

# These compounds are not informative as direct cross-feeding candidates in
# this setting, regardless of an opposite-direction flux pattern.
BOUNDARY_KEYS = {"c00001", "c00007", "c00011", "c00080"}  # H2O, O2, CO2, H+


def base_id(value: object) -> str:
    value = str(value).strip().removeprefix("R_").removeprefix("EX_")
    return re.sub(r"_(?:e|E|m|M)$", "", value).lower()


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, value: str) -> str:
        self.parent.setdefault(value, value)
        if self.parent[value] != value:
            self.parent[value] = self.find(self.parent[value])
        return self.parent[value]

    def union(self, first: str, second: str) -> None:
        first, second = self.find(first), self.find(second)
        if first != second:
            self.parent[second] = first


def load_crosswalk() -> UnionFind:
    union_find = UnionFind()
    table = pd.read_csv(CROSSWALK)
    for _, row in table.iterrows():
        first, second = base_id(row.get("mapped_from", "")), base_id(row.get("carvefungi_metabolite", ""))
        if first and second and first != "nan" and second != "nan":
            union_find.union(first, second)
    return union_find


def supplied_medium_keys(union_find: UnionFind, scenario: str) -> pd.DataFrame:
    master = pd.read_csv(MASTER)
    master = master[(master["biomass_mode"] == "shared_52") & (master["scenario"] == scenario)]
    rows: list[dict[str, object]] = []
    for _, run in master.iterrows():
        medium = pd.read_csv(run["medium_csv"])
        reaction_column = "reaction" if "reaction" in medium.columns else medium.columns[0]
        for _, entry in medium.iterrows():
            raw = str(entry[reaction_column])
            rows.append({"run_id": run["run_id"], "chemical_key": f"id:{union_find.find(base_id(raw))}", "medium_reaction": raw})
    return pd.DataFrame(rows).drop_duplicates(["run_id", "chemical_key"])


def community_medium_residuals(master: pd.DataFrame, chemical_keys: set[str]) -> pd.DataFrame:
    """Read net community-medium fluxes for the candidate chemical keys."""
    rows: list[dict[str, object]] = []
    for _, run in master.iterrows():
        path = (
            PACKAGE / "reports/96_fungal_species_syncom_factorial_panel/raw_micom_runs"
            / str(run["manifest"])
            / str(run["scenario"])
            / f"{run['manifest']}_{run['scenario']}_emma_style_fluxes.csv"
        )
        fluxes = pd.read_csv(path)
        organism_column = fluxes.columns[0]
        medium = fluxes[fluxes[organism_column] == "medium"].iloc[0]
        for chemical_key in chemical_keys:
            column = f"EX_{chemical_key.removeprefix('id:')}_m"
            value = float(medium[column]) if column in medium.index and pd.notna(medium[column]) else 0.0
            rows.append({"run_id": run["run_id"], "chemical_key": chemical_key, "medium_net_flux": value})
    return pd.DataFrame(rows)


def interface_evidence(row: pd.Series) -> str:
    if row["raw_exchange_donor"] == row["raw_exchange_recipient"]:
        return "same_exported_exchange_identifier"
    if row["external_metabolite_id_donor"] == row["external_metabolite_id_recipient"]:
        return "same_external_metabolite_identifier"
    return "chemical_crosswalk_only"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario", default="proline_relaxed_1")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    union_find = load_crosswalk()
    medium = supplied_medium_keys(union_find, args.scenario)
    pairs = pd.read_csv(AUDIT / "04_candidate_donor_recipient_pairs_reaudited.csv")
    active = pd.read_csv(AUDIT / "01_active_exchange_fluxes_mapped.csv")

    pairs = pairs.merge(medium[["run_id", "chemical_key"]].assign(medium_associated=True), on=["run_id", "chemical_key"], how="left")
    pairs["medium_associated"] = pairs["medium_associated"].fillna(False)
    pairs["base_key"] = pairs["chemical_key"].str.removeprefix("id:")
    pairs["interface_evidence"] = pairs.apply(interface_evidence, axis=1)
    pairs["interpretation_class"] = "nonmedium_candidate"
    pairs.loc[pairs["medium_associated"], "interpretation_class"] = "medium_associated"
    pairs.loc[pairs["base_key"].isin(BOUNDARY_KEYS), "interpretation_class"] = "currency_or_boundary"
    pairs.to_csv(output / "01_all_candidate_pairs_classified.csv", index=False)

    # The relevant balance comparison sums all bacterial secretion and fungal
    # uptake for each run and chemical key, rather than comparing one selected
    # donor-recipient pair in isolation.
    active["base_key"] = active["chemical_key"].str.removeprefix("id:")
    active = active.merge(medium[["run_id", "chemical_key"]].assign(medium_associated=True), on=["run_id", "chemical_key"], how="left")
    active["medium_associated"] = active["medium_associated"].fillna(False)
    balance_rows: list[dict[str, object]] = []
    for (run_id, chemical_key), group in active.groupby(["run_id", "chemical_key"], dropna=False):
        bacterial_secretion = group[(group["role"] == "bacteria") & (group["direction"] == "secretion")]["abs_flux"].sum()
        fungal_uptake = group[(group["role"] == "fungus") & (group["direction"] == "uptake")]["abs_flux"].sum()
        if bacterial_secretion <= 0 or fungal_uptake <= 0:
            continue
        balance_rows.append(
            {
                "run_id": run_id,
                "chemical_key": chemical_key,
                "medium_associated": bool(group["medium_associated"].iloc[0]),
                "bacterial_total_secretion": bacterial_secretion,
                "fungal_total_uptake": fungal_uptake,
                "smaller_to_larger_flux_ratio": min(bacterial_secretion, fungal_uptake) / max(bacterial_secretion, fungal_uptake),
            }
        )
    balance = pd.DataFrame(balance_rows)
    master = pd.read_csv(MASTER)
    master = master[(master["biomass_mode"] == "shared_52") & (master["scenario"] == args.scenario)]
    residuals = community_medium_residuals(master, set(balance["chemical_key"]))
    balance = balance.merge(residuals, on=["run_id", "chemical_key"], how="left")
    balance["medium_residual_fraction"] = balance["medium_net_flux"].abs() / balance[["bacterial_total_secretion", "fungal_total_uptake"]].max(axis=1)
    balance.to_csv(output / "02_bacterial_secretion_fungal_uptake_balance.csv", index=False)

    direct = pairs[(pairs["donor_role"] == "bacteria") & (pairs["recipient_role"] == "fungus")].copy()
    direct_summary = (
        direct.groupby(["interpretation_class", "interface_evidence"], dropna=False)
        .agg(candidate_pairs=("pair_flux", "size"), chemical_keys=("chemical_key", "nunique"), runs=("run_id", "nunique"))
        .reset_index()
    )
    direct_summary.to_csv(output / "03_bacteria_to_fungus_evidence_summary.csv", index=False)

    informative = direct[direct["interpretation_class"] == "nonmedium_candidate"].copy()
    informative.to_csv(output / "04_nonmedium_bacteria_to_fungus_candidates.csv", index=False)
    validated = informative.merge(balance, on=["run_id", "chemical_key"], how="left")
    validated = validated[
        (validated["interface_evidence"] == "same_external_metabolite_identifier")
        & (validated["medium_residual_fraction"] <= 1e-4)
    ].copy()
    validated["bacterial_secretion_covers_fungal_uptake"] = (
        validated["bacterial_total_secretion"] >= validated["fungal_total_uptake"]
    )
    validated["similar_flux_scale"] = validated["smaller_to_larger_flux_ratio"] >= 0.5
    validated.to_csv(output / "05_shared_pool_compatible_nonmedium_candidates.csv", index=False)
    (output / "summary.json").write_text(json.dumps({
        "scenario": args.scenario,
        "all_candidate_pairs": int(len(pairs)),
        "bacteria_to_fungus_pairs": int(len(direct)),
        "nonmedium_bacteria_to_fungus_pairs": int(len(informative)),
        "shared_pool_compatible_nonmedium_pairs": int(len(validated)),
        "note": "Non-medium candidates remain model-derived hypotheses. The shared-pool-compatible subset has the same external metabolite identifier and negligible net medium residual, but it is not direct biological proof of transfer.",
    }, indent=2), encoding="utf-8")
    print((output / "summary.json").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
