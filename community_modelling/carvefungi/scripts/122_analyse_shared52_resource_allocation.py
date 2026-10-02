#!/usr/bin/env python3
"""Analyse member-level medium resource use in shared-52 MICOM solutions.

Uses the complete Script 96 long exchange export. Medium rows are kept as a
community-input reference, while fungal and bacterial uptake are calculated
from organism rows directly. Positive flux is secretion and negative flux is
uptake. This script does not rerun MICOM or infer intracellular uptake.
"""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import cobra


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
SOURCE = PACKAGE / "reports/96_fungal_species_syncom_factorial_panel/summaries/exchange_flux_long_nonzero.csv"
REPORT_ROOT = PACKAGE / "reports/96_fungal_species_syncom_factorial_panel"
OUT = PACKAGE / "reports/122_shared52_resource_allocation"
FUNGAL_MODEL_PATHS = [
    PACKAGE / "fungal_models/biomass_variant_panel/Bipolaris_sorokiniana_carvefungi_52biomass.xml",
    PACKAGE / "fungal_models/biomass_variant_panel/Fusarium_graminearum_carvefungi_52biomass.xml",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--flux-threshold",
        type=float,
        default=1e-4,
        help="Minimum member uptake flux retained for overlap analysis (default: 1e-4).",
    )
    args = parser.parse_args()
    if args.flux_threshold < 0:
        raise ValueError("--flux-threshold must be non-negative")

    OUT.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(REPORT_ROOT / "summaries/master_run_index.csv")
    master = master[(master["biomass_mode"] == "shared_52") & (master["scenario"] == "proline_relaxed_1")]
    fungal_exchange_ids = set()
    fungal_exchange_metabolites = {}
    metabolite_names = {}
    for model_path in FUNGAL_MODEL_PATHS:
        model = cobra.io.read_sbml_model(str(model_path))
        for metabolite in model.metabolites:
            base_id = metabolite.id
            if base_id.endswith("_e") or base_id.endswith("_m"):
                base_id = base_id[:-2]
            if metabolite.name and base_id not in metabolite_names:
                metabolite_names[base_id] = metabolite.name
        for reaction in model.exchanges:
            fungal_exchange_ids.add(reaction.id)
            metabolites = list(reaction.metabolites)
            if metabolites:
                metabolite_id = metabolites[0].id
                if metabolite_id.endswith("_e"):
                    metabolite_id = metabolite_id[:-2]
                fungal_exchange_metabolites[reaction.id] = metabolite_id
    frames = []
    for _, meta in master.iterrows():
        flux_path = REPORT_ROOT / "raw_micom_runs" / str(meta["manifest"]) / str(meta["scenario"]) / f"{meta['manifest']}_{meta['scenario']}_emma_style_fluxes.csv"
        if not flux_path.exists():
            continue
        fluxes = pd.read_csv(flux_path)
        first = fluxes.columns[0]
        fluxes = fluxes.rename(columns={first: "organism"})
        exchange_cols = [
            column for column in fluxes.columns
            if str(column).startswith("EX_") or str(column) in fungal_exchange_ids
        ]
        if not exchange_cols:
            continue
        long = fluxes[["organism"] + exchange_cols].melt(
            id_vars=["organism"], var_name="metabolite", value_name="flux"
        )
        for key in ["run_id", "manifest", "scenario", "biomass_mode", "syncom", "fungus_context", "community_type"]:
            long[key] = meta[key]
        frames.append(long)
    if not frames:
        raise FileNotFoundError("No shared-52 raw MICOM flux files were found.")
    raw = pd.concat(frames, ignore_index=True)
    raw["role"] = raw["organism"].map(
        lambda value: "medium" if str(value) == "medium" else ("fungus" if ("Bipolaris" in str(value) or "Fusarium" in str(value)) else "bacteria")
    )
    raw["fungus"] = raw["organism"].map(
        lambda value: "Bs" if "Bipolaris" in str(value) else ("Fg" if "Fusarium" in str(value) else "none")
    )
    def harmonize_exchange(exchange: str) -> str:
        exchange = str(exchange)
        if exchange in fungal_exchange_metabolites:
            return fungal_exchange_metabolites[exchange]
        if exchange.startswith("EX_"):
            return exchange[3:].removesuffix("_m").removesuffix("_e")
        return exchange

    raw["metabolite"] = raw["metabolite"].map(harmonize_exchange)
    raw["metabolite_name"] = raw["metabolite"].map(metabolite_names).fillna("")
    raw["flux"] = pd.to_numeric(raw["flux"], errors="coerce")
    raw["uptake"] = (-raw["flux"]).clip(lower=0.0)

    organism = raw[raw["role"].isin(["fungus", "bacteria"])].copy()
    uptake = (
        organism.groupby(["run_id", "manifest", "syncom", "fungus_context", "community_type", "role", "fungus", "organism", "metabolite", "metabolite_name"], dropna=False)
        .agg(uptake_flux=("uptake", "sum"))
        .reset_index()
    )
    # Apply the cutoff after harmonising and aggregating exchange rows. This
    # avoids counting tiny solver residuals as resource use.
    uptake = uptake[uptake["uptake_flux"] >= args.flux_threshold].copy()
    uptake.to_csv(OUT / "01_member_uptake_long.csv", index=False)

    role_uptake = (
        uptake.groupby(["run_id", "manifest", "syncom", "fungus_context", "community_type", "role", "metabolite", "metabolite_name"], dropna=False)
        .agg(uptake_flux=("uptake_flux", "sum"))
        .reset_index()
    )
    role_uptake.to_csv(OUT / "02_role_uptake_long.csv", index=False)

    medium = raw[raw["role"] == "medium"].copy()
    # In the MICOM export, negative medium flux supplies the shared pool.
    medium["community_medium_supply"] = (-medium["flux"]).clip(lower=0.0)
    medium = medium[medium["community_medium_supply"] >= args.flux_threshold]
    medium = (
        medium.groupby(["run_id", "manifest", "syncom", "fungus_context", "community_type", "metabolite", "metabolite_name"], dropna=False)
        .agg(community_medium_supply=("community_medium_supply", "sum"))
        .reset_index()
    )
    medium.to_csv(OUT / "03_community_medium_supply_long.csv", index=False)

    contexts = role_uptake[["run_id", "syncom", "fungus_context", "community_type"]].drop_duplicates()
    overlap_rows = []
    for _, context in contexts.iterrows():
        subset = role_uptake[role_uptake["run_id"] == context["run_id"]]
        fungal = subset[subset["role"] == "fungus"].groupby("metabolite")["uptake_flux"].sum()
        bacterial = subset[subset["role"] == "bacteria"].groupby("metabolite")["uptake_flux"].sum()
        metabolites = sorted(set(fungal.index) | set(bacterial.index))
        f = fungal.reindex(metabolites, fill_value=0.0).to_numpy(dtype=float)
        b = bacterial.reindex(metabolites, fill_value=0.0).to_numpy(dtype=float)
        shared = int(((f > 0) & (b > 0)).sum())
        denom = float(np.linalg.norm(f) * np.linalg.norm(b))
        cosine = float(np.dot(f, b) / denom) if denom else 0.0
        overlap_rows.append({
            "run_id": context["run_id"],
            "syncom": context["syncom"],
            "fungus_context": context["fungus_context"],
            "community_type": context["community_type"],
            "fungal_uptake_metabolites": int((f > 0).sum()),
            "bacterial_uptake_metabolites": int((b > 0).sum()),
            "shared_uptake_metabolites": shared,
            "fungal_bacterial_weighted_overlap_cosine": cosine,
            "fungal_total_uptake": float(f.sum()),
            "bacterial_total_uptake": float(b.sum()),
        })
    pd.DataFrame(overlap_rows).to_csv(OUT / "04_fungal_bacterial_resource_overlap.csv", index=False)

    fungal_fungal_rows = []
    for run_id, run_data in uptake[uptake["role"] == "fungus"].groupby("run_id"):
        profiles = run_data.groupby(["fungus", "metabolite"])["uptake_flux"].sum().unstack(fill_value=0.0)
        if not {"Bs", "Fg"}.issubset(profiles.index):
            continue
        bs = profiles.loc["Bs"]
        fg = profiles.loc["Fg"]
        metabolites = sorted(set(bs.index) | set(fg.index))
        bs_values = bs.reindex(metabolites, fill_value=0.0).to_numpy(dtype=float)
        fg_values = fg.reindex(metabolites, fill_value=0.0).to_numpy(dtype=float)
        shared_mask = (bs_values > 0) & (fg_values > 0)
        denom = float(np.linalg.norm(bs_values) * np.linalg.norm(fg_values))
        metadata = run_data.iloc[0]
        fungal_fungal_rows.append({
            "run_id": run_id,
            "syncom": metadata["syncom"],
            "fungus_context": metadata["fungus_context"],
            "community_type": metadata["community_type"],
            "bipolaris_uptake_metabolites": int((bs_values > 0).sum()),
            "fusarium_uptake_metabolites": int((fg_values > 0).sum()),
            "shared_uptake_metabolites": int(shared_mask.sum()),
            "fungal_fungal_weighted_overlap_cosine": float(np.dot(bs_values, fg_values) / denom) if denom else 0.0,
            "bipolaris_total_uptake": float(bs_values.sum()),
            "fusarium_total_uptake": float(fg_values.sum()),
        })
    pd.DataFrame(fungal_fungal_rows).to_csv(OUT / "13_bipolaris_fusarium_resource_overlap.csv", index=False)

    member_overlap_rows = []
    member_overlap_metabolite_rows = []
    for run_id, run_data in uptake.groupby("run_id"):
        fungi = run_data[run_data["role"] == "fungus"]
        bacteria = run_data[run_data["role"] == "bacteria"]
        for fungal_member, fungal_data in fungi.groupby("organism"):
            fungal_profile = fungal_data.groupby("metabolite")["uptake_flux"].sum()
            for bacterial_member, bacterial_data in bacteria.groupby("organism"):
                bacterial_profile = bacterial_data.groupby("metabolite")["uptake_flux"].sum()
                metabolites = sorted(set(fungal_profile.index) | set(bacterial_profile.index))
                f = fungal_profile.reindex(metabolites, fill_value=0.0).to_numpy(dtype=float)
                b = bacterial_profile.reindex(metabolites, fill_value=0.0).to_numpy(dtype=float)
                shared_mask = (f > 0) & (b > 0)
                denom = float(np.linalg.norm(f) * np.linalg.norm(b))
                for metabolite, fungal_flux, bacterial_flux in zip(metabolites, f, b):
                    if fungal_flux > 0 and bacterial_flux > 0:
                        member_overlap_metabolite_rows.append({
                            "run_id": run_id,
                            "syncom": fungal_data["syncom"].iloc[0],
                            "fungus_context": fungal_data["fungus_context"].iloc[0],
                            "community_type": fungal_data["community_type"].iloc[0],
                            "fungal_member": fungal_member,
                            "bacterial_member": bacterial_member,
                            "metabolite": metabolite,
                            "metabolite_name": fungal_data.loc[fungal_data["metabolite"] == metabolite, "metabolite_name"].iloc[0] if not fungal_data.loc[fungal_data["metabolite"] == metabolite, "metabolite_name"].empty else "",
                            "fungal_uptake_flux": float(fungal_flux),
                            "bacterial_uptake_flux": float(bacterial_flux),
                            "shared_min_flux": float(min(fungal_flux, bacterial_flux)),
                        })
                member_overlap_rows.append({
                    "run_id": run_id,
                    "syncom": fungal_data["syncom"].iloc[0],
                    "fungus_context": fungal_data["fungus_context"].iloc[0],
                    "community_type": fungal_data["community_type"].iloc[0],
                    "fungal_member": fungal_member,
                    "bacterial_member": bacterial_member,
                    "shared_metabolites": int(shared_mask.sum()),
                    "weighted_overlap_min_flux": float(np.minimum(f[shared_mask], b[shared_mask]).sum()),
                    "weighted_overlap_cosine": float(np.dot(f, b) / denom) if denom else 0.0,
                    "overlap_fraction_of_fungal_uptake": float(np.minimum(f[shared_mask], b[shared_mask]).sum() / f.sum()) if f.sum() else 0.0,
                    "overlap_fraction_of_bacterial_uptake": float(np.minimum(f[shared_mask], b[shared_mask]).sum() / b.sum()) if b.sum() else 0.0,
                    "fungal_total_uptake": float(f.sum()),
                    "bacterial_total_uptake": float(b.sum()),
                })
    member_overlap = pd.DataFrame(member_overlap_rows).sort_values(
        ["run_id", "weighted_overlap_min_flux"], ascending=[True, False]
    )
    member_overlap.to_csv(OUT / "06_bacterial_fungal_overlap_by_member.csv", index=False)
    shared_resources = pd.DataFrame(member_overlap_metabolite_rows)
    shared_resources.to_csv(OUT / "07_shared_resource_use_by_member_metabolite.csv", index=False)

    # Ranked tables are the compact thesis-facing outputs: they retain the
    # named bacterial member and the metabolites contributing to its overlap.
    top_members = (
        member_overlap.sort_values(
            ["run_id", "weighted_overlap_min_flux", "shared_metabolites"],
            ascending=[True, False, False],
        )
        .groupby("run_id", sort=False)
        .head(20)
        .reset_index(drop=True)
    )
    top_members.to_csv(OUT / "08_top_bacterial_fungal_overlap_members.csv", index=False)

    if shared_resources.empty:
        top_resources = shared_resources.copy()
    else:
        # The fungal uptake value is a context-level quantity and must not be
        # summed once for every bacterial member sharing that metabolite.
        fungal_resource_totals = (
            uptake[uptake["role"] == "fungus"]
            .groupby(
                ["run_id", "syncom", "fungus_context", "community_type", "metabolite", "metabolite_name"],
                as_index=False,
            )
            .agg(fungal_uptake_flux=("uptake_flux", "sum"))
        )
        top_resources = (
            shared_resources.groupby(
                ["run_id", "syncom", "fungus_context", "community_type", "metabolite"],
                as_index=False,
            )
            .agg(
                bacterial_members=("bacterial_member", "nunique"),
                shared_min_flux=("shared_min_flux", "sum"),
                bacterial_uptake_flux=("bacterial_uptake_flux", "sum"),
            )
            .merge(
                fungal_resource_totals,
                on=["run_id", "syncom", "fungus_context", "community_type", "metabolite"],
                how="left",
            )
            .sort_values(
                ["run_id", "shared_min_flux", "bacterial_members"],
                ascending=[True, False, False],
            )
            .groupby("run_id", sort=False)
            .head(20)
            .reset_index(drop=True)
        )
    top_resources.to_csv(OUT / "09_top_shared_resources.csv", index=False)

    if shared_resources.empty:
        resource_allocation = pd.DataFrame()
    else:
        resource_allocation = (
            shared_resources.groupby(
                ["run_id", "syncom", "fungus_context", "community_type", "metabolite", "metabolite_name"],
                as_index=False,
            )
            .agg(
                bacterial_members=("bacterial_member", "nunique"),
                shared_min_flux=("shared_min_flux", "sum"),
                bacterial_uptake_flux=("bacterial_uptake_flux", "sum"),
            )
            .merge(
                uptake[uptake["role"] == "fungus"].groupby(
                    ["run_id", "syncom", "fungus_context", "community_type", "metabolite", "metabolite_name"],
                    as_index=False,
                ).agg(fungal_uptake_flux=("uptake_flux", "sum")),
                on=["run_id", "syncom", "fungus_context", "community_type", "metabolite", "metabolite_name"],
                how="left",
            )
        )
        resource_allocation["fungal_share_of_shared_uptake"] = (
            resource_allocation["fungal_uptake_flux"] /
            (resource_allocation["fungal_uptake_flux"] + resource_allocation["bacterial_uptake_flux"])
        )
        resource_allocation["bacterial_share_of_shared_uptake"] = 1.0 - resource_allocation["fungal_share_of_shared_uptake"]
        resource_allocation = resource_allocation.sort_values(
            ["run_id", "shared_min_flux"], ascending=[True, False]
        )
    resource_allocation.to_csv(OUT / "12_shared_resource_allocation.csv", index=False)

    if shared_resources.empty:
        competition_matrix = pd.DataFrame()
    else:
        competition_matrix = (
            shared_resources.groupby(
                ["run_id", "fungal_member", "bacterial_member"], as_index=False
            )
            .agg(
                shared_metabolites=("metabolite", "nunique"),
                weighted_overlap_min_flux=("shared_min_flux", "sum"),
            )
            .pivot_table(
                index=["fungal_member", "bacterial_member"],
                columns="run_id",
                values="weighted_overlap_min_flux",
                fill_value=0.0,
            )
            .reset_index()
        )
    competition_matrix.to_csv(OUT / "10_member_overlap_matrix.csv", index=False)

    member_consumption = uptake.copy()
    member_totals = member_consumption.groupby(
        ["run_id", "organism"], as_index=False
    )["uptake_flux"].sum().rename(columns={"uptake_flux": "member_total_uptake"})
    member_consumption = member_consumption.merge(
        member_totals, on=["run_id", "organism"], how="left"
    )
    member_consumption["uptake_fraction_of_member"] = (
        member_consumption["uptake_flux"] / member_consumption["member_total_uptake"]
    )
    member_consumption["metabolite_rank"] = member_consumption.groupby(
        ["run_id", "organism"]
    )["uptake_flux"].rank(method="first", ascending=False).astype(int)
    member_consumption.sort_values(
        ["run_id", "role", "organism", "metabolite_rank"]
    ).to_csv(OUT / "11_member_consumption_ranked.csv", index=False)

    matrix = role_uptake.pivot_table(index="metabolite", columns="run_id", values="uptake_flux", aggfunc="sum", fill_value=0.0).reset_index()
    matrix.to_csv(OUT / "05_role_uptake_matrix.csv", index=False)

    metadata = {
        "source": str(SOURCE),
        "rows_after_filter": int(len(raw)),
        "member_uptake_rows": int(len(uptake)),
        "contexts": int(contexts.shape[0]),
        "member_overlap_rows": int(len(member_overlap)),
        "member_overlap_metabolite_rows": int(len(member_overlap_metabolite_rows)),
        "flux_threshold": args.flux_threshold,
        "top_member_rows": int(len(top_members)),
        "top_resource_rows": int(len(top_resources)),
        "member_consumption_rows": int(len(member_consumption)),
        "shared_resource_allocation_rows": int(len(resource_allocation)),
        "fungus_fungus_overlap_rows": int(len(fungal_fungal_rows)),
        "note": "Uptake is calculated from member-level exchange rows; intracellular consumption is not inferred. Tiny uptake values below the threshold are excluded.",
    }
    (OUT / "summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
