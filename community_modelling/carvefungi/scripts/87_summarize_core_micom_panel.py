#!/usr/bin/env python3
"""Script 87: summarize the core MICOM panel for RAVEN vs CarveFungi.

Run this after Script 86. It converts the raw MICOM outputs into comparison
tables for analyses 2--12:

2  fungal monoculture baseline
3  bacterial community baseline
4  fungal+SynCom community growth
5  relative fungal-growth effect
6  bacterial member growth changes after fungus addition
7  resource-competition / substrate-utilization candidates
8  metabolic niche overlap
9  putative bacteria -> fungus exchanges
10 putative fungus -> bacteria exchanges
11 HvSC1 vs HvSC2 exchange-profile comparison
12 RAVEN vs CarveFungi exchange-profile comparison

The exchange-transfer tables are hypothesis-generating only: direction is
inferred from one MICOM optimum by sign-opposed exchange fluxes and should be
confirmed later with FVA/blocking/drop-out tests.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

import pandas as pd


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
REPORT_ROOT = PACKAGE / "reports/86_core_micom_panel"
OUT = PACKAGE / "reports/87_core_micom_panel_summary"

TOL = 1e-9

CORE_MANIFESTS = [
    "fungus_alone__raven_independent",
    "fungus_alone__carvefungi",
    "hvsc1_bacteria_only",
    "hvsc2_bacteria_only",
    "hvsc1_plus_bipolaris__raven_independent",
    "hvsc2_plus_bipolaris__raven_independent",
    "hvsc1_plus_bipolaris__carvefungi",
    "hvsc2_plus_bipolaris__carvefungi",
]


def classify_manifest(manifest: str) -> dict[str, str]:
    if "raven_independent" in manifest:
        track = "raven_independent"
    elif "carvefungi" in manifest:
        track = "carvefungi"
    else:
        track = ""

    if manifest.startswith("hvsc1"):
        syncom = "HvSC1"
    elif manifest.startswith("hvsc2"):
        syncom = "HvSC2"
    else:
        syncom = ""

    if manifest.startswith("fungus_alone"):
        run_type = "fungus_alone"
    elif manifest.endswith("bacteria_only"):
        run_type = "bacteria_only"
    elif "plus_bipolaris" in manifest:
        run_type = "fungus_plus_syncom"
    else:
        run_type = "unknown"

    return {"track": track, "syncom": syncom, "run_type": run_type}


def result_paths(manifest: str, scenario: str) -> dict[str, Path]:
    base = REPORT_ROOT / manifest / scenario
    prefix = f"{manifest}_{scenario}_emma_style"
    return {
        "summary": base / f"{prefix}_community_growth_summary.csv",
        "growth": base / f"{prefix}_individualtaxon_growth.csv",
        "fluxes": base / f"{prefix}_fluxes.csv",
        "run_summary": base / "run_summary.json",
    }


def discover_runs() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if not REPORT_ROOT.exists():
        raise FileNotFoundError(
            f"{REPORT_ROOT} does not exist. Run Script 86 before Script 87."
        )

    for manifest_dir in sorted(p for p in REPORT_ROOT.iterdir() if p.is_dir() and not p.name.startswith("_")):
        manifest = manifest_dir.name
        for scenario_dir in sorted(p for p in manifest_dir.iterdir() if p.is_dir()):
            scenario = scenario_dir.name
            paths = result_paths(manifest, scenario)
            if paths["summary"].exists() and paths["growth"].exists() and paths["fluxes"].exists():
                row = {"manifest": manifest, "scenario": scenario}
                row.update(classify_manifest(manifest))
                row.update({k: str(v) for k, v in paths.items()})
                rows.append(row)
    return rows


def read_summary(path: str | Path) -> dict[str, object]:
    table = pd.read_csv(path)
    if table.empty:
        return {"community_growth": math.nan, "objective_value": math.nan, "status": "empty"}
    row = table.iloc[0].to_dict()
    return {
        "community_growth": row.get("community_growth", math.nan),
        "objective_value": row.get("objective_value", math.nan),
        "status": row.get("status", ""),
        "tradeoff": row.get("tradeoff", math.nan),
    }


def read_growth(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    if "compartments" not in table.columns:
        raise ValueError(f"Expected column 'compartments' in {path}")
    return table


def is_fungus_id(value: object) -> bool:
    return "bipolaris" in str(value).lower()


def is_medium_id(value: object) -> bool:
    return str(value).lower() == "medium"


def numeric(value: object) -> float:
    try:
        if pd.isna(value):
            return 0.0
        return float(value)
    except Exception:
        return 0.0


def growth_for_fungus(growth: pd.DataFrame) -> float:
    subset = growth.loc[growth["compartments"].map(is_fungus_id)]
    if subset.empty:
        return math.nan
    return float(subset["growth_rate"].iloc[0])


def exchange_columns(columns: Iterable[str]) -> list[str]:
    return [c for c in columns if c.startswith("EX_")]


def read_fluxes_long(path: str | Path, tol: float = TOL) -> pd.DataFrame:
    wide = pd.read_csv(path)
    if "compartment" not in wide.columns:
        raise ValueError(f"Expected column 'compartment' in {path}")
    cols = exchange_columns(wide.columns)
    if not cols:
        return pd.DataFrame(columns=["compartment", "reaction", "flux"])

    long = wide.melt(
        id_vars=["compartment"],
        value_vars=cols,
        var_name="reaction",
        value_name="flux",
    )
    long["flux"] = pd.to_numeric(long["flux"], errors="coerce").fillna(0.0)
    long = long.loc[long["flux"].abs() > tol].copy()
    if long.empty:
        return long
    long["role"] = long["compartment"].map(
        lambda c: "medium" if is_medium_id(c) else ("fungus" if is_fungus_id(c) else "bacterium")
    )
    long["direction_inferred"] = long["flux"].map(lambda x: "secretion_or_export" if x > 0 else "uptake_or_import")
    return long


def safe_ratio(numerator: float, denominator: float) -> float:
    if denominator is None or pd.isna(denominator) or abs(float(denominator)) <= TOL:
        return math.nan
    return float(numerator) / float(denominator)


def write_csv(table: pd.DataFrame, name: str) -> Path:
    path = OUT / name
    table.to_csv(path, index=False)
    return path


def total_growth_by_role(growth: pd.DataFrame) -> dict[str, float]:
    """Return simple and abundance-weighted growth totals by role."""
    table = growth.copy()
    table["role"] = table["compartments"].map(
        lambda c: "fungus" if is_fungus_id(c) else ("medium" if is_medium_id(c) else "bacterium")
    )
    table["growth_rate"] = pd.to_numeric(table["growth_rate"], errors="coerce").fillna(0.0)
    if "abundance" in table.columns:
        table["abundance"] = pd.to_numeric(table["abundance"], errors="coerce").fillna(0.0)
    else:
        table["abundance"] = 0.0

    bacteria = table.loc[table["role"] == "bacterium"].copy()
    fungus = table.loc[table["role"] == "fungus"].copy()
    return {
        "fungal_growth": float(fungus["growth_rate"].sum()) if not fungus.empty else math.nan,
        "bacterial_growth_sum_unweighted": float(bacteria["growth_rate"].sum()) if not bacteria.empty else 0.0,
        "bacterial_growth_sum_abundance_weighted": float((bacteria["growth_rate"] * bacteria["abundance"]).sum()) if not bacteria.empty else 0.0,
        "bacterial_member_count": int(bacteria.shape[0]),
    }


def producer_consumer_table(exchange_long: pd.DataFrame) -> pd.DataFrame:
    """Summarize bacterial producer/consumer identities for every exchange."""
    if exchange_long.empty:
        return pd.DataFrame(
            columns=[
                "reaction",
                "bacterial_total_flux",
                "bacterial_producer_strains",
                "bacterial_consumer_strains",
                "bacterial_producer_count",
                "bacterial_consumer_count",
            ]
        )

    bacteria = exchange_long.loc[exchange_long["role"] == "bacterium"].copy()
    if bacteria.empty:
        return pd.DataFrame(
            columns=[
                "reaction",
                "bacterial_total_flux",
                "bacterial_producer_strains",
                "bacterial_consumer_strains",
                "bacterial_producer_count",
                "bacterial_consumer_count",
            ]
        )

    rows: list[dict[str, object]] = []
    for reaction, group in bacteria.groupby("reaction"):
        producers = sorted(map(str, group.loc[group["flux"] > TOL, "compartment"].unique()))
        consumers = sorted(map(str, group.loc[group["flux"] < -TOL, "compartment"].unique()))
        rows.append(
            {
                "reaction": reaction,
                "bacterial_total_flux": float(group["flux"].sum()),
                "bacterial_producer_strains": ";".join(producers),
                "bacterial_consumer_strains": ";".join(consumers),
                "bacterial_producer_count": len(producers),
                "bacterial_consumer_count": len(consumers),
            }
        )
    return pd.DataFrame(rows).sort_values("reaction", ignore_index=True)


def write_standardized_run_exports(
    manifest: str,
    scenario: str,
    meta: dict[str, str],
    summary: dict[str, object],
    growth: pd.DataFrame,
    exchange_long: pd.DataFrame,
) -> dict[str, object]:
    """Write a complete per-run export bundle so MICOM need not be rerun."""
    export_dir = OUT / "standardized_per_run_exports" / manifest / scenario
    export_dir.mkdir(parents=True, exist_ok=True)

    totals = total_growth_by_role(growth)
    run_summary = {
        "manifest": manifest,
        "scenario": scenario,
        **meta,
        **summary,
        **totals,
    }
    pd.DataFrame([run_summary]).to_csv(export_dir / "standardized_growth_summary.csv", index=False)

    member = growth.copy()
    member.insert(0, "manifest", manifest)
    member.insert(1, "scenario", scenario)
    member.insert(2, "track", meta["track"])
    member.insert(3, "syncom", meta["syncom"])
    member.insert(4, "run_type", meta["run_type"])
    member["role"] = member["compartments"].map(
        lambda c: "fungus" if is_fungus_id(c) else ("medium" if is_medium_id(c) else "bacterium")
    )
    member.to_csv(export_dir / "individual_growth_rates.csv", index=False)

    ex = exchange_long.copy()
    if not ex.empty:
        ex.insert(0, "manifest", manifest)
        ex.insert(1, "scenario", scenario)
        ex.insert(2, "track", meta["track"])
        ex.insert(3, "syncom", meta["syncom"])
        ex.insert(4, "run_type", meta["run_type"])
    ex.to_csv(export_dir / "all_exchange_fluxes_long.csv", index=False)

    medium_ex = ex.loc[ex["role"] == "medium"].copy() if not ex.empty else pd.DataFrame()
    fungal_ex = ex.loc[ex["role"] == "fungus"].copy() if not ex.empty else pd.DataFrame()
    bacterial_ex = ex.loc[ex["role"] == "bacterium"].copy() if not ex.empty else pd.DataFrame()

    medium_ex.to_csv(export_dir / "medium_exchange_fluxes.csv", index=False)
    fungal_ex.to_csv(export_dir / "fungal_uptake_secretion_fluxes.csv", index=False)
    bacterial_ex.to_csv(export_dir / "bacterial_exchange_fluxes_by_strain.csv", index=False)

    producer_consumer = producer_consumer_table(exchange_long)
    producer_consumer.to_csv(export_dir / "bacterial_producer_consumer_identities.csv", index=False)

    # A joined table with fungal flux plus bacterial producers/consumers is the
    # most convenient file for cross-feeding interpretation.
    if not fungal_ex.empty:
        fungal_by_rxn = fungal_ex.groupby("reaction", as_index=False)["flux"].sum()
        fungal_by_rxn = fungal_by_rxn.rename(columns={"flux": "fungal_flux"})
        joined = fungal_by_rxn.merge(producer_consumer, on="reaction", how="left")
        joined["fungal_direction_inferred"] = joined["fungal_flux"].map(
            lambda x: "secretion_or_export" if x > 0 else "uptake_or_import"
        )
        joined["putative_bacteria_to_fungus"] = (
            (joined["fungal_flux"] < -TOL)
            & (joined["bacterial_producer_count"].fillna(0).astype(int) > 0)
        )
        joined["putative_fungus_to_bacteria"] = (
            (joined["fungal_flux"] > TOL)
            & (joined["bacterial_consumer_count"].fillna(0).astype(int) > 0)
        )
    else:
        joined = pd.DataFrame()
    joined.to_csv(export_dir / "fungal_bacterial_exchange_interface.csv", index=False)

    run_summary["standardized_export_dir"] = str(export_dir)
    run_summary["exchange_flux_rows"] = int(ex.shape[0])
    run_summary["medium_exchange_rows"] = int(medium_ex.shape[0])
    run_summary["fungal_exchange_rows"] = int(fungal_ex.shape[0])
    run_summary["bacterial_exchange_rows"] = int(bacterial_ex.shape[0])
    run_summary["producer_consumer_rows"] = int(producer_consumer.shape[0])
    return run_summary


def exchange_set(long: pd.DataFrame, role: str | None = None) -> set[str]:
    if long.empty:
        return set()
    table = long if role is None else long.loc[long["role"] == role]
    return set(table["reaction"].astype(str))


def summarize_set_overlap(a: set[str], b: set[str]) -> dict[str, object]:
    shared = a & b
    union = a | b
    denom_min = min(len(a), len(b))
    return {
        "shared_count": len(shared),
        "first_unique_count": len(a - b),
        "second_unique_count": len(b - a),
        "union_count": len(union),
        "jaccard": len(shared) / len(union) if union else math.nan,
        "overlap_coefficient": len(shared) / denom_min if denom_min else math.nan,
        "shared_reactions": ";".join(sorted(shared)),
        "first_unique_reactions": ";".join(sorted(a - b)),
        "second_unique_reactions": ";".join(sorted(b - a)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    runs = discover_runs()
    if not runs:
        raise RuntimeError(f"No complete MICOM runs found under {REPORT_ROOT}")

    run_table = pd.DataFrame(runs)
    run_table.to_csv(OUT / "00_discovered_runs.csv", index=False)

    growth_summary_rows: list[dict[str, object]] = []
    member_growth_rows: list[dict[str, object]] = []
    exchange_rows: list[pd.DataFrame] = []
    standardized_run_rows: list[dict[str, object]] = []

    # Keep per-run tables in memory for cross-run comparisons.
    run_growth: dict[tuple[str, str], pd.DataFrame] = {}
    run_exchange: dict[tuple[str, str], pd.DataFrame] = {}

    for run in runs:
        manifest = str(run["manifest"])
        scenario = str(run["scenario"])
        meta = classify_manifest(manifest)
        summary = read_summary(str(run["summary"]))
        growth = read_growth(str(run["growth"]))
        flux_long = read_fluxes_long(str(run["fluxes"]))

        run_growth[(manifest, scenario)] = growth
        run_exchange[(manifest, scenario)] = flux_long

        standardized_run_rows.append(
            write_standardized_run_exports(
                manifest=manifest,
                scenario=scenario,
                meta=meta,
                summary=summary,
                growth=growth,
                exchange_long=flux_long,
            )
        )

        growth_summary_rows.append(
            {
                "manifest": manifest,
                "scenario": scenario,
                **meta,
                **summary,
                **total_growth_by_role(growth),
                "member_count": int(growth.shape[0]),
            }
        )

        member = growth.copy()
        member.insert(0, "manifest", manifest)
        member.insert(1, "scenario", scenario)
        member.insert(2, "track", meta["track"])
        member.insert(3, "syncom", meta["syncom"])
        member.insert(4, "run_type", meta["run_type"])
        member["role"] = member["compartments"].map(
            lambda c: "fungus" if is_fungus_id(c) else ("medium" if is_medium_id(c) else "bacterium")
        )
        member_growth_rows.append(member.to_dict(orient="records"))

        if not flux_long.empty:
            flux_long = flux_long.copy()
            flux_long.insert(0, "manifest", manifest)
            flux_long.insert(1, "scenario", scenario)
            flux_long.insert(2, "track", meta["track"])
            flux_long.insert(3, "syncom", meta["syncom"])
            flux_long.insert(4, "run_type", meta["run_type"])
            exchange_rows.append(flux_long)

    growth_summary = pd.DataFrame(growth_summary_rows)
    member_growth = pd.DataFrame([row for chunk in member_growth_rows for row in chunk])
    exchange_long = (
        pd.concat(exchange_rows, ignore_index=True)
        if exchange_rows
        else pd.DataFrame(columns=["manifest", "scenario", "track", "syncom", "run_type", "compartment", "reaction", "flux", "role"])
    )

    write_csv(growth_summary, "01_all_run_growth_summary.csv")
    write_csv(member_growth, "01_all_member_growth_long.csv")
    write_csv(exchange_long, "07_exchange_flux_long.csv")
    standardized_runs = pd.DataFrame(standardized_run_rows)
    write_csv(standardized_runs, "00_standardized_per_run_export_index.csv")

    fungal_baseline = growth_summary.loc[growth_summary["run_type"] == "fungus_alone"].copy()
    write_csv(fungal_baseline, "02_fungal_monoculture_baseline.csv")

    fungal_baseline_ex = exchange_long.loc[exchange_long["run_type"] == "fungus_alone"].copy()
    write_csv(fungal_baseline_ex, "02_fungal_monoculture_exchange_fluxes.csv")

    bacterial_baseline = growth_summary.loc[growth_summary["run_type"] == "bacteria_only"].copy()
    write_csv(bacterial_baseline, "03_bacterial_community_baseline.csv")

    bacterial_member_baseline = member_growth.loc[member_growth["run_type"] == "bacteria_only"].copy()
    write_csv(bacterial_member_baseline, "03_bacterial_baseline_member_growth.csv")

    bacterial_baseline_ex = exchange_long.loc[exchange_long["run_type"] == "bacteria_only"].copy()
    write_csv(bacterial_baseline_ex, "03_bacterial_baseline_exchange_fluxes.csv")

    community_growth = growth_summary.loc[growth_summary["run_type"] == "fungus_plus_syncom"].copy()
    write_csv(community_growth, "04_community_growth_simulation.csv")

    community_member_growth = member_growth.loc[member_growth["run_type"] == "fungus_plus_syncom"].copy()
    write_csv(community_member_growth, "04_community_member_growth.csv")

    # 5. Relative fungal-growth effect.
    alone_lookup = {
        (str(row["track"]), str(row["scenario"])): numeric(row["fungal_growth"])
        for _, row in fungal_baseline.iterrows()
    }
    rel_rows: list[dict[str, object]] = []
    for _, row in community_growth.iterrows():
        alone = alone_lookup.get((str(row["track"]), str(row["scenario"])), math.nan)
        community = numeric(row["fungal_growth"])
        ratio = safe_ratio(community, alone)
        rel_rows.append(
            {
                "scenario": row["scenario"],
                "track": row["track"],
                "syncom": row["syncom"],
                "fungal_growth_alone": alone,
                "fungal_growth_in_community": community,
                "community_growth": row["community_growth"],
                "ratio_community_vs_alone": ratio,
                "percent_change": (ratio - 1.0) * 100.0 if not pd.isna(ratio) else math.nan,
            }
        )
    rel = pd.DataFrame(rel_rows)
    write_csv(rel, "05_relative_fungal_growth_effect.csv")

    # 6. Bacterial growth changes after fungus addition.
    bacteria_only_lookup: dict[tuple[str, str, str], float] = {}
    for _, row in bacterial_member_baseline.iterrows():
        bacteria_only_lookup[(str(row["syncom"]), str(row["scenario"]), str(row["compartments"]))] = numeric(row["growth_rate"])

    delta_rows: list[dict[str, object]] = []
    for _, row in community_member_growth.iterrows():
        if str(row["role"]) != "bacterium":
            continue
        base = bacteria_only_lookup.get((str(row["syncom"]), str(row["scenario"]), str(row["compartments"])), math.nan)
        with_fungus = numeric(row["growth_rate"])
        ratio = safe_ratio(with_fungus, base)
        delta_rows.append(
            {
                "scenario": row["scenario"],
                "track": row["track"],
                "syncom": row["syncom"],
                "bacterium": row["compartments"],
                "growth_bacteria_only": base,
                "growth_with_fungus": with_fungus,
                "delta": with_fungus - base if not pd.isna(base) else math.nan,
                "ratio_with_fungus_vs_bacteria_only": ratio,
                "percent_change": (ratio - 1.0) * 100.0 if not pd.isna(ratio) else math.nan,
            }
        )
    member_deltas = pd.DataFrame(delta_rows)
    write_csv(member_deltas, "06_bacterial_member_growth_changes_after_fungus.csv")

    # 7--10. Exchange interaction heuristics inside each fungus+SynCom run.
    competition_rows: list[dict[str, object]] = []
    niche_rows: list[dict[str, object]] = []
    bac_to_fun_rows: list[dict[str, object]] = []
    fun_to_bac_rows: list[dict[str, object]] = []

    for _, row in community_growth.iterrows():
        manifest = str(row["manifest"])
        scenario = str(row["scenario"])
        track = str(row["track"])
        syncom = str(row["syncom"])
        ex = run_exchange.get((manifest, scenario), pd.DataFrame())
        if ex.empty:
            continue

        fungal = ex.loc[ex["role"] == "fungus"]
        bacteria = ex.loc[ex["role"] == "bacterium"]
        fset = exchange_set(ex, "fungus")
        bset = exchange_set(ex, "bacterium")
        overlap = summarize_set_overlap(fset, bset)
        niche_rows.append(
            {
                "scenario": scenario,
                "track": track,
                "syncom": syncom,
                "fungal_exchange_count": len(fset),
                "bacterial_exchange_count": len(bset),
                **overlap,
            }
        )

        fungal_flux = fungal.groupby("reaction")["flux"].sum()
        bacterial_sum = bacteria.groupby("reaction")["flux"].sum()
        bacterial_min = bacteria.groupby("reaction")["flux"].min()
        bacterial_max = bacteria.groupby("reaction")["flux"].max()
        bacterial_producers = bacteria.loc[bacteria["flux"] > TOL].groupby("reaction")["compartment"].apply(lambda x: ";".join(sorted(map(str, x))))
        bacterial_consumers = bacteria.loc[bacteria["flux"] < -TOL].groupby("reaction")["compartment"].apply(lambda x: ";".join(sorted(map(str, x))))

        shared_reactions = sorted(set(fungal_flux.index) & set(bacterial_sum.index))
        for reaction in shared_reactions:
            ff = float(fungal_flux.get(reaction, 0.0))
            bsum = float(bacterial_sum.get(reaction, 0.0))
            bmin = float(bacterial_min.get(reaction, 0.0))
            bmax = float(bacterial_max.get(reaction, 0.0))

            if ff < -TOL and bmin < -TOL:
                competition_rows.append(
                    {
                        "scenario": scenario,
                        "track": track,
                        "syncom": syncom,
                        "reaction": reaction,
                        "fungal_flux": ff,
                        "bacterial_total_flux": bsum,
                        "bacterial_min_member_flux": bmin,
                        "bacterial_max_member_flux": bmax,
                        "interpretation": "putative_shared_uptake_competition",
                    }
                )

            if ff < -TOL and bmax > TOL:
                bac_to_fun_rows.append(
                    {
                        "scenario": scenario,
                        "track": track,
                        "syncom": syncom,
                        "reaction": reaction,
                        "fungal_uptake_flux": ff,
                        "bacterial_total_flux": bsum,
                        "producer_strains": bacterial_producers.get(reaction, ""),
                        "interpretation": "putative_bacteria_to_fungus_transfer",
                    }
                )

            if ff > TOL and bmin < -TOL:
                fun_to_bac_rows.append(
                    {
                        "scenario": scenario,
                        "track": track,
                        "syncom": syncom,
                        "reaction": reaction,
                        "fungal_secretion_flux": ff,
                        "bacterial_total_flux": bsum,
                        "receiver_strains": bacterial_consumers.get(reaction, ""),
                        "interpretation": "putative_fungus_to_bacteria_transfer",
                    }
                )

    competition = pd.DataFrame(competition_rows).sort_values(
        ["scenario", "track", "syncom", "reaction"], ignore_index=True
    ) if competition_rows else pd.DataFrame()
    niche = pd.DataFrame(niche_rows)
    bac_to_fun = pd.DataFrame(bac_to_fun_rows)
    fun_to_bac = pd.DataFrame(fun_to_bac_rows)

    write_csv(competition, "07_resource_competition_candidates.csv")
    write_csv(niche, "08_metabolic_niche_overlap.csv")
    write_csv(bac_to_fun, "09_candidate_bacteria_to_fungus_exchanges.csv")
    write_csv(fun_to_bac, "10_candidate_fungus_to_bacteria_exchanges.csv")

    # 11. HvSC1 vs HvSC2 exchange-profile comparisons within a fungal model.
    hv_rows: list[dict[str, object]] = []
    for scenario in sorted(community_growth["scenario"].dropna().unique()):
        for track in sorted(community_growth["track"].dropna().unique()):
            key1 = (f"hvsc1_plus_bipolaris__{track}", scenario)
            key2 = (f"hvsc2_plus_bipolaris__{track}", scenario)
            if key1 not in run_exchange or key2 not in run_exchange:
                continue
            set1 = exchange_set(run_exchange[key1])
            set2 = exchange_set(run_exchange[key2])
            overlap = summarize_set_overlap(set1, set2)
            hv_rows.append(
                {
                    "scenario": scenario,
                    "track": track,
                    "first": "HvSC1",
                    "second": "HvSC2",
                    **overlap,
                }
            )
    hv_compare = pd.DataFrame(hv_rows)
    write_csv(hv_compare, "11_hvsc1_vs_hvsc2_exchange_profile_comparison.csv")

    # 12. RAVEN vs CarveFungi exchange-profile comparisons within a SynCom.
    model_rows: list[dict[str, object]] = []
    for scenario in sorted(community_growth["scenario"].dropna().unique()):
        for syncom in sorted(community_growth["syncom"].dropna().unique()):
            prefix = syncom.lower()
            key1 = (f"{prefix}_plus_bipolaris__raven_independent", scenario)
            key2 = (f"{prefix}_plus_bipolaris__carvefungi", scenario)
            if key1 not in run_exchange or key2 not in run_exchange:
                continue
            set1 = exchange_set(run_exchange[key1])
            set2 = exchange_set(run_exchange[key2])
            overlap = summarize_set_overlap(set1, set2)
            model_rows.append(
                {
                    "scenario": scenario,
                    "syncom": syncom,
                    "first": "RAVEN",
                    "second": "CarveFungi",
                    **overlap,
                }
            )
    model_compare = pd.DataFrame(model_rows)
    write_csv(model_compare, "12_raven_vs_carvefungi_exchange_profile_comparison.csv")

    # Candidate lists for Scripts 88+ (FVA, blocking, supplementation, drop-out, pairwise).
    candidate_frames: list[pd.DataFrame] = []
    if not bac_to_fun.empty:
        tmp = bac_to_fun.copy()
        tmp["candidate_type"] = "bacteria_to_fungus_exchange_for_FVA_or_blocking"
        tmp["priority_score"] = tmp["fungal_uptake_flux"].abs()
        candidate_frames.append(tmp)
    if not fun_to_bac.empty:
        tmp = fun_to_bac.copy()
        tmp["candidate_type"] = "fungus_to_bacteria_exchange_for_FVA_or_blocking"
        tmp["priority_score"] = tmp["fungal_secretion_flux"].abs()
        candidate_frames.append(tmp)
    if not competition.empty:
        tmp = competition.copy()
        tmp["candidate_type"] = "shared_resource_for_competition_supplementation"
        tmp["priority_score"] = tmp["fungal_flux"].abs() + tmp["bacterial_total_flux"].abs()
        candidate_frames.append(tmp)
    if not member_deltas.empty:
        tmp = member_deltas.copy()
        tmp["candidate_type"] = "bacterial_strain_for_dropout_or_pairwise"
        tmp["reaction"] = ""
        tmp["priority_score"] = tmp["delta"].abs()
        candidate_frames.append(tmp)

    candidates = (
        pd.concat(candidate_frames, ignore_index=True, sort=False)
        .sort_values("priority_score", ascending=False, ignore_index=True)
        if candidate_frames
        else pd.DataFrame()
    )
    write_csv(candidates, "13_to_18_candidate_followup_list.csv")

    summary = {
        "package": str(PACKAGE),
        "input_report_root": str(REPORT_ROOT),
        "output_report_root": str(OUT),
        "complete_runs_discovered": int(len(runs)),
        "scenarios": sorted(run_table["scenario"].unique()),
        "manifests": sorted(run_table["manifest"].unique()),
        "fungal_monoculture_rows": int(fungal_baseline.shape[0]),
        "bacterial_baseline_rows": int(bacterial_baseline.shape[0]),
        "fungus_plus_syncom_rows": int(community_growth.shape[0]),
        "exchange_flux_rows": int(exchange_long.shape[0]),
        "standardized_per_run_exports": int(standardized_runs.shape[0]),
        "resource_competition_candidates": int(competition.shape[0]),
        "bacteria_to_fungus_candidates": int(bac_to_fun.shape[0]),
        "fungus_to_bacteria_candidates": int(fun_to_bac.shape[0]),
        "followup_candidate_rows": int(candidates.shape[0]),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = [
        "# Core MICOM panel summary",
        "",
        "This report summarizes the batch outputs from Script 86 for the RAVEN and CarveFungi *Bipolaris sorokiniana* tracks.",
        "",
        "## What this summary covers",
        "",
        "- Analysis 2: fungal monoculture baseline.",
        "- Analysis 3: bacterial SynCom baseline.",
        "- Analysis 4: fungal-SynCom community growth.",
        "- Analysis 5: relative fungal-growth effect.",
        "- Analysis 6: bacterial member growth changes after fungus addition.",
        "- Analysis 7: shared resource-competition candidates.",
        "- Analysis 8: metabolic niche-overlap metrics.",
        "- Analyses 9-10: putative directional exchange candidates.",
        "- Analyses 11-12: HvSC1/HvSC2 and RAVEN/CarveFungi exchange-profile comparisons.",
        "",
        "## Important limitation",
        "",
        "The exchange-transfer calls are inferred from one MICOM optimal flux distribution. They are useful for prioritization, but should be confirmed with FVA, blocking, nutrient supplementation, species drop-out, and pairwise simulations.",
        "",
        "## Output tables",
        "",
        "- `00_standardized_per_run_export_index.csv`",
        "- `02_fungal_monoculture_baseline.csv`",
        "- `03_bacterial_community_baseline.csv`",
        "- `04_community_growth_simulation.csv`",
        "- `05_relative_fungal_growth_effect.csv`",
        "- `06_bacterial_member_growth_changes_after_fungus.csv`",
        "- `07_resource_competition_candidates.csv`",
        "- `08_metabolic_niche_overlap.csv`",
        "- `09_candidate_bacteria_to_fungus_exchanges.csv`",
        "- `10_candidate_fungus_to_bacteria_exchanges.csv`",
        "- `11_hvsc1_vs_hvsc2_exchange_profile_comparison.csv`",
        "- `12_raven_vs_carvefungi_exchange_profile_comparison.csv`",
        "- `13_to_18_candidate_followup_list.csv`",
        "",
        "## Per-run standardized exports",
        "",
        "Each completed MICOM run also gets a self-contained export folder under:",
        "",
        "`standardized_per_run_exports/<manifest>/<scenario>/`",
        "",
        "Each folder contains:",
        "",
        "- `standardized_growth_summary.csv`: fungal growth, total bacterial growth, community growth, objective/status.",
        "- `individual_growth_rates.csv`: growth rate of every fungal/bacterial member.",
        "- `all_exchange_fluxes_long.csv`: all nonzero exchange-like fluxes in long format.",
        "- `medium_exchange_fluxes.csv`: medium uptake/secretion fluxes.",
        "- `fungal_uptake_secretion_fluxes.csv`: fungal exchange fluxes.",
        "- `bacterial_exchange_fluxes_by_strain.csv`: bacterial exchange fluxes by strain.",
        "- `bacterial_producer_consumer_identities.csv`: producer/consumer strain identities per exchange.",
        "- `fungal_bacterial_exchange_interface.csv`: fungal flux joined to bacterial producers/consumers.",
        "",
        "## Compact counts",
        "",
        f"- Complete runs discovered: {summary['complete_runs_discovered']}",
        f"- Standardized per-run exports: {summary['standardized_per_run_exports']}",
        f"- Exchange flux rows: {summary['exchange_flux_rows']}",
        f"- Resource-competition candidates: {summary['resource_competition_candidates']}",
        f"- Putative bacteria-to-fungus candidates: {summary['bacteria_to_fungus_candidates']}",
        f"- Putative fungus-to-bacteria candidates: {summary['fungus_to_bacteria_candidates']}",
        f"- Candidate follow-up rows for analyses 13-18: {summary['followup_candidate_rows']}",
        "",
    ]
    (OUT / "summary.md").write_text("\n".join(md), encoding="utf-8")

    print("Script 87 compact summary")
    print(f"complete_runs_discovered: {summary['complete_runs_discovered']}")
    print(f"scenarios: {';'.join(summary['scenarios'])}")
    print(f"manifests: {';'.join(summary['manifests'])}")
    print(f"fungus_plus_syncom_rows: {summary['fungus_plus_syncom_rows']}")
    print(f"exchange_flux_rows: {summary['exchange_flux_rows']}")
    print(f"followup_candidate_rows: {summary['followup_candidate_rows']}")
    print(f"outputs: {OUT}")


if __name__ == "__main__":
    main()
