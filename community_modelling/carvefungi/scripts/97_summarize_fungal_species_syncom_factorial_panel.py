#!/usr/bin/env python3
"""Summarize the fungal-species × SynCom factorial MICOM panel.

Reads outputs from Script 96 and creates analysis-ready long-format tables:

* master_run_index.csv
* growth_summary_all_runs.csv
* fungal_growth_effects.csv
* bacterial_growth_effects.csv
* exchange_flux_long_nonzero.csv
* candidate_opposite_direction_exchanges.csv

The exchange tables are intentionally "candidate" tables. They identify
opposite-direction flux patterns between fungi and bacteria, but they do not
replace later FVA/blocking tests.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
REPORT_ROOT = PACKAGE / "reports/96_fungal_species_syncom_factorial_panel"
RAW_ROOT = REPORT_ROOT / "raw_micom_runs"
SUMMARY_DIR = REPORT_ROOT / "summaries"


def infer_metadata(manifest: str, scenario: str) -> dict[str, object]:
    if "specific" in manifest:
        biomass_mode = "species_specific"
    elif manifest.endswith("_52") or "_52" in manifest:
        biomass_mode = "shared_52"
    else:
        biomass_mode = "not_applicable"

    if manifest.startswith("hvsc1"):
        syncom = "HvSC1"
    elif manifest.startswith("hvsc2"):
        syncom = "HvSC2"
    else:
        syncom = "none"

    if "bipolaris" in manifest and "fusarium" in manifest:
        fungus_context = "Bs+Fg"
    elif "fungi_pair" in manifest:
        fungus_context = "Bs+Fg"
    elif "bipolaris" in manifest:
        fungus_context = "Bs"
    elif "fusarium" in manifest:
        fungus_context = "Fg"
    else:
        fungus_context = "none"

    if "bacteria_only" in manifest:
        community_type = "bacteria_only"
    elif manifest.startswith("fungus_alone"):
        community_type = "fungus_alone"
    elif manifest.startswith("fungi_pair"):
        community_type = "fungi_pair"
    elif "plus_fungi_pair" in manifest:
        community_type = "dual_fungus_syncom"
    elif "plus_bipolaris" in manifest or "plus_fusarium" in manifest:
        community_type = "single_fungus_syncom"
    else:
        community_type = "unknown"

    return {
        "run_id": f"{manifest}__{scenario}",
        "manifest": manifest,
        "scenario": scenario,
        "biomass_mode": biomass_mode,
        "syncom": syncom,
        "fungus_context": fungus_context,
        "community_type": community_type,
    }


def organism_role(compartment: str) -> str:
    text = str(compartment)
    if text == "medium":
        return "medium"
    if "Bipolaris" in text or "Fusarium" in text:
        return "fungus"
    return "bacterium"


def fungus_label(compartment: str) -> str:
    text = str(compartment)
    if "Bipolaris" in text:
        return "Bs"
    if "Fusarium" in text:
        return "Fg"
    return "none"


def read_job_rows() -> pd.DataFrame:
    status_path = REPORT_ROOT / "job_status.csv"
    if not status_path.exists():
        raise FileNotFoundError(f"{status_path}. Run Script 96 first.")
    status = pd.read_csv(status_path)
    rows = []
    for _, row in status.iterrows():
        meta = infer_metadata(str(row["manifest"]), str(row["scenario"]))
        meta.update(
            {
                "returncode": int(row["returncode"]),
                "status": str(row["status"]),
                "medium_csv": str(row["medium_csv"]),
                "summary_path": str(row["summary_path"]),
                "stdout_log": str(row.get("stdout_log", "")),
                "stderr_log": str(row.get("stderr_log", "")),
            }
        )
        rows.append(meta)
    return pd.DataFrame(rows)


def result_paths(manifest: str, scenario: str) -> dict[str, Path]:
    base = RAW_ROOT / manifest / scenario / f"{manifest}_{scenario}_emma_style"
    return {
        "community": base.with_name(base.name + "_community_growth_summary.csv"),
        "growth": base.with_name(base.name + "_individualtaxon_growth.csv"),
        "fluxes": base.with_name(base.name + "_fluxes.csv"),
    }


def build_growth_summary(master: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for _, meta in master.iterrows():
        paths = result_paths(meta["manifest"], meta["scenario"])
        if not paths["growth"].exists():
            continue
        growth = pd.read_csv(paths["growth"])
        first = growth.columns[0]
        growth = growth.rename(columns={first: "organism"})
        growth["role"] = growth["organism"].map(organism_role)
        growth["fungus"] = growth["organism"].map(fungus_label)

        community_growth = None
        community_status = None
        if paths["community"].exists():
            community = pd.read_csv(paths["community"])
            if not community.empty:
                community_growth = community.loc[0, "community_growth"]
                community_status = community.loc[0, "status"]

        for key in [
            "run_id",
            "manifest",
            "scenario",
            "biomass_mode",
            "syncom",
            "fungus_context",
            "community_type",
        ]:
            growth[key] = meta[key]
        growth["community_growth"] = community_growth
        growth["community_status"] = community_status
        frames.append(growth)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_fungal_growth_effects(growth: pd.DataFrame) -> pd.DataFrame:
    if growth.empty:
        return pd.DataFrame()
    fungi = growth[growth["role"] == "fungus"].copy()
    baseline_rows = []

    alone = fungi[fungi["community_type"] == "fungus_alone"]
    for _, row in alone.iterrows():
        baseline_rows.append(
            {
                "baseline_type": "fungus_alone",
                "scenario": row["scenario"],
                "biomass_mode": row["biomass_mode"],
                "fungus": row["fungus"],
                "baseline_growth": row["growth_rate"],
            }
        )

    pair = fungi[fungi["community_type"] == "fungi_pair"]
    for _, row in pair.iterrows():
        baseline_rows.append(
            {
                "baseline_type": "fungi_pair",
                "scenario": row["scenario"],
                "biomass_mode": row["biomass_mode"],
                "fungus": row["fungus"],
                "baseline_growth": row["growth_rate"],
            }
        )

    baselines = pd.DataFrame(baseline_rows)
    effect_rows = []
    for _, row in fungi.iterrows():
        if row["community_type"] == "single_fungus_syncom":
            baseline_type = "fungus_alone"
        elif row["community_type"] == "dual_fungus_syncom":
            baseline_type = "fungi_pair"
        else:
            continue
        match = baselines[
            (baselines["baseline_type"] == baseline_type)
            & (baselines["scenario"] == row["scenario"])
            & (baselines["biomass_mode"] == row["biomass_mode"])
            & (baselines["fungus"] == row["fungus"])
        ]
        baseline_growth = float(match.iloc[0]["baseline_growth"]) if not match.empty else float("nan")
        growth_value = float(row["growth_rate"])
        ratio = growth_value / baseline_growth if baseline_growth and baseline_growth == baseline_growth else float("nan")
        effect_rows.append(
            {
                "run_id": row["run_id"],
                "scenario": row["scenario"],
                "biomass_mode": row["biomass_mode"],
                "community_type": row["community_type"],
                "syncom": row["syncom"],
                "fungus": row["fungus"],
                "baseline_type": baseline_type,
                "baseline_growth": baseline_growth,
                "community_fungal_growth": growth_value,
                "delta_growth": growth_value - baseline_growth if baseline_growth == baseline_growth else float("nan"),
                "ratio_vs_baseline": ratio,
                "percent_change_vs_baseline": (ratio - 1.0) * 100 if ratio == ratio else float("nan"),
            }
        )
    return pd.DataFrame(effect_rows)


def build_bacterial_growth_effects(growth: pd.DataFrame) -> pd.DataFrame:
    if growth.empty:
        return pd.DataFrame()
    bacteria = growth[growth["role"] == "bacterium"].copy()
    baseline = bacteria[bacteria["community_type"] == "bacteria_only"][
        ["scenario", "syncom", "organism", "growth_rate"]
    ].rename(columns={"growth_rate": "bacteria_only_growth"})
    comparisons = bacteria[bacteria["community_type"].isin(["single_fungus_syncom", "dual_fungus_syncom"])].copy()
    joined = comparisons.merge(baseline, on=["scenario", "syncom", "organism"], how="left")
    joined["with_fungus_growth"] = joined["growth_rate"]
    joined["delta_growth"] = joined["with_fungus_growth"] - joined["bacteria_only_growth"]
    joined["ratio_vs_bacteria_only"] = joined["with_fungus_growth"] / joined["bacteria_only_growth"]
    joined["percent_change_vs_bacteria_only"] = (joined["ratio_vs_bacteria_only"] - 1.0) * 100
    keep = [
        "run_id",
        "scenario",
        "biomass_mode",
        "community_type",
        "syncom",
        "fungus_context",
        "organism",
        "bacteria_only_growth",
        "with_fungus_growth",
        "delta_growth",
        "ratio_vs_bacteria_only",
        "percent_change_vs_bacteria_only",
    ]
    return joined[keep]


def build_exchange_flux_tables(master: pd.DataFrame, threshold: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    flux_frames = []
    for _, meta in master.iterrows():
        paths = result_paths(meta["manifest"], meta["scenario"])
        if not paths["fluxes"].exists():
            continue
        fluxes = pd.read_csv(paths["fluxes"])
        first = fluxes.columns[0]
        fluxes = fluxes.rename(columns={first: "organism"})
        exchange_cols = [col for col in fluxes.columns if str(col).startswith("EX_")]
        if not exchange_cols:
            continue
        long = fluxes[["organism"] + exchange_cols].melt(
            id_vars=["organism"],
            var_name="exchange",
            value_name="flux",
        )
        long = long[long["flux"].abs() > threshold].copy()
        if long.empty:
            continue
        long["role"] = long["organism"].map(organism_role)
        long["fungus"] = long["organism"].map(fungus_label)
        long["flux_direction_sign"] = long["flux"].map(lambda x: "positive" if x > 0 else "negative")
        for key in [
            "run_id",
            "manifest",
            "scenario",
            "biomass_mode",
            "syncom",
            "fungus_context",
            "community_type",
        ]:
            long[key] = meta[key]
        flux_frames.append(long)

    if not flux_frames:
        return pd.DataFrame(), pd.DataFrame()
    flux_long = pd.concat(flux_frames, ignore_index=True)

    candidate_rows = []
    for (run_id, exchange), group in flux_long.groupby(["run_id", "exchange"]):
        fungi = group[group["role"] == "fungus"]
        bacteria = group[group["role"] == "bacterium"]
        if fungi.empty or bacteria.empty:
            continue
        for _, frow in fungi.iterrows():
            opposite = bacteria[bacteria["flux"] * frow["flux"] < 0]
            for _, brow in opposite.iterrows():
                direction = "bacteria_to_fungus" if brow["flux"] > 0 and frow["flux"] < 0 else "fungus_to_bacteria"
                candidate_rows.append(
                    {
                        "run_id": run_id,
                        "scenario": frow["scenario"],
                        "biomass_mode": frow["biomass_mode"],
                        "community_type": frow["community_type"],
                        "syncom": frow["syncom"],
                        "fungus": frow["fungus"],
                        "exchange": exchange,
                        "candidate_direction": direction,
                        "fungal_organism": frow["organism"],
                        "fungal_flux": frow["flux"],
                        "bacterial_organism": brow["organism"],
                        "bacterial_flux": brow["flux"],
                    }
                )
    return flux_long, pd.DataFrame(candidate_rows)


def write_readme(master: pd.DataFrame) -> None:
    readme = REPORT_ROOT / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Fungal-species × SynCom factorial MICOM panel",
                "",
                "This folder contains the main community modelling panel for the CarveMe fungal models.",
                "",
                "Main ecological medium:",
                "",
                "- `original_proline_0p1`: original shared barley/community medium with `EX_proL_m = 0.1`.",
                "",
                "Targeted control:",
                "",
                "- `proline_relaxed_1`: Fusarium-specific monoculture control with `EX_proL_m = 1.0`.",
                "",
                "The main design compares:",
                "",
                "- Bipolaris alone",
                "- Fusarium alone",
                "- Bipolaris + Fusarium",
                "- HvSC1 alone",
                "- HvSC2 alone",
                "- each single fungus + each SynCom",
                "- both fungi together + each SynCom",
                "",
                "Each fungal context is tested with the shared 52-component biomass and the species-specific biomass.",
                "",
                f"Runs summarized: {master.shape[0]}",
                "",
                "Important summary tables are in `summaries/`.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-flux-summary", action="store_true", help="Only summarize growth tables.")
    parser.add_argument("--flux-threshold", type=float, default=1e-9)
    args = parser.parse_args()

    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    master = read_job_rows()
    master.to_csv(SUMMARY_DIR / "master_run_index.csv", index=False)

    growth = build_growth_summary(master)
    growth.to_csv(SUMMARY_DIR / "growth_summary_all_runs.csv", index=False)

    fungal_effects = build_fungal_growth_effects(growth)
    fungal_effects.to_csv(SUMMARY_DIR / "fungal_growth_effects.csv", index=False)

    bacterial_effects = build_bacterial_growth_effects(growth)
    bacterial_effects.to_csv(SUMMARY_DIR / "bacterial_growth_effects.csv", index=False)

    exchange_rows = 0
    candidate_rows = 0
    if not args.skip_flux_summary:
        flux_long, candidates = build_exchange_flux_tables(master, args.flux_threshold)
        flux_long.to_csv(SUMMARY_DIR / "exchange_flux_long_nonzero.csv", index=False)
        candidates.to_csv(SUMMARY_DIR / "candidate_opposite_direction_exchanges.csv", index=False)
        exchange_rows = int(flux_long.shape[0])
        candidate_rows = int(candidates.shape[0])

    summary = {
        "report_root": str(REPORT_ROOT),
        "summary_dir": str(SUMMARY_DIR),
        "runs": int(master.shape[0]),
        "growth_rows": int(growth.shape[0]),
        "fungal_growth_effect_rows": int(fungal_effects.shape[0]),
        "bacterial_growth_effect_rows": int(bacterial_effects.shape[0]),
        "exchange_flux_rows": exchange_rows,
        "candidate_opposite_direction_exchange_rows": candidate_rows,
    }
    (SUMMARY_DIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_readme(master)

    print("Script 97 compact summary")
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
