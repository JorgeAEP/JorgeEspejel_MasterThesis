#!/usr/bin/env python3
"""Screen harmonized member-level MICOM exchange fluxes for interactions.

The screen uses the existing shared-52, proline-relaxed MICOM exchange export.
It excludes medium rows from donor/recipient pairs, harmonizes exchange IDs to
the metabolite_key field, and treats positive flux as secretion and negative
flux as uptake. Results are candidate interactions, not proof of transfer.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
SOURCE = PACKAGE / "reports/96_fungal_species_syncom_factorial_panel/summaries/exchange_flux_long_nonzero.csv"
BACTERIAL_SOURCE = PACKAGE / "reports/103_exchange_interpretation_tables/05_bacterial_exchange_fluxes_by_run.csv"
OUT = PACKAGE / "reports/121_member_level_interaction_screen"


def harmonize_metabolite(row: pd.Series) -> str:
    """Use the already harmonized metabolite_key, with an exchange fallback."""
    value = str(row.get("metabolite_key", ""))
    if value and value != "nan":
        return value
    value = str(row.get("exchange", ""))
    return value.removeprefix("EX_").removesuffix("_m").removesuffix("_e")


def add_direction(data: pd.DataFrame, threshold: float) -> pd.DataFrame:
    data = data.copy()
    data["flux"] = pd.to_numeric(data["flux"], errors="coerce")
    data["metabolite_harmonized"] = data.apply(harmonize_metabolite, axis=1)
    data["direction"] = data["flux"].map(lambda value: "secretion" if value > 0 else "uptake")
    data["abs_flux"] = data["flux"].abs()
    return data[data["abs_flux"] >= threshold].copy()


def make_pairs(data: pd.DataFrame, donor_role: str, recipient_role: str) -> pd.DataFrame:
    donors = data[(data["role"] == donor_role) & (data["direction"] == "secretion")]
    recipients = data[(data["role"] == recipient_role) & (data["direction"] == "uptake")]
    keys = ["run_id", "metabolite_harmonized"]
    donor_cols = keys + ["organism", "fungus_context", "community_type", "syncom", "flux", "abs_flux"]
    recipient_cols = keys + ["organism", "fungus_context", "community_type", "syncom", "flux", "abs_flux"]
    merged = donors[donor_cols].merge(recipients[recipient_cols], on=keys, suffixes=("_donor", "_recipient"))
    if merged.empty:
        return pd.DataFrame(columns=["run_id", "metabolite_harmonized", "donor_role", "recipient_role"])
    merged["donor_role"] = donor_role
    merged["recipient_role"] = recipient_role
    # A member cannot be counted as its own donor and recipient. Self-pairs
    # commonly arise when one organism has reversible or duplicated exchange
    # records for the same harmonized metabolite.
    if donor_role == recipient_role:
        merged = merged[merged["organism_donor"] != merged["organism_recipient"]].copy()
    merged["donor_flux"] = merged["abs_flux_donor"]
    merged["recipient_flux"] = merged["abs_flux_recipient"]
    merged["pair_flux"] = merged[["donor_flux", "recipient_flux"]].min(axis=1)
    merged["pool_balance_difference"] = merged["donor_flux"] - merged["recipient_flux"]
    return merged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", default=str(SOURCE))
    parser.add_argument("--output-dir", default=str(OUT))
    parser.add_argument("--threshold", type=float, default=1e-6)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(args.source)
    # Script 96 already contains complete member-level rows. The fallback
    # supports older exports where bacterial rows were stored separately.
    if "bacterium" not in set(raw["role"].dropna().unique()) and BACTERIAL_SOURCE.exists():
        raw_bacterial = pd.read_csv(BACTERIAL_SOURCE)
        raw = pd.concat([raw, raw_bacterial], ignore_index=True)
    raw["role"] = raw["role"].replace({"bacterium": "bacteria"})
    raw = raw[(raw["biomass_mode"] == "shared_52") & (raw["scenario"] == "proline_relaxed_1")]
    raw = raw[raw["role"].isin(["bacteria", "fungus"])].copy()
    data = add_direction(raw, args.threshold)

    pairs = pd.concat([
        make_pairs(data, "bacteria", "fungus"),
        make_pairs(data, "fungus", "bacteria"),
        make_pairs(data, "bacteria", "bacteria"),
        make_pairs(data, "fungus", "fungus"),
    ], ignore_index=True)
    pairs.to_csv(output_dir / "01_candidate_donor_recipient_pairs.csv", index=False)

    unique_columns = [
        "run_id", "metabolite_harmonized", "organism_donor", "organism_recipient",
        "donor_role", "recipient_role",
    ]
    if pairs.empty:
        unique_pairs = pairs.copy()
    else:
        unique_pairs = (
            pairs.groupby(unique_columns, dropna=False)
            .agg(
                donor_flux=("donor_flux", "max"),
                recipient_flux=("recipient_flux", "max"),
                pair_flux=("pair_flux", "max"),
                pool_balance_difference=("pool_balance_difference", "mean"),
            )
            .reset_index()
        )
    unique_pairs.to_csv(output_dir / "02_unique_donor_recipient_metabolite_pairs.csv", index=False)

    summary_columns = [
        "run_id", "fungus_context_donor", "community_type_donor", "syncom_donor",
        "donor_role", "recipient_role", "candidate_pairs", "unique_metabolites",
        "summed_pair_flux", "mean_pair_flux",
    ]
    if unique_pairs.empty:
        summary = pd.DataFrame(columns=summary_columns)
    else:
        summary = (
            pairs.groupby(["run_id", "fungus_context_donor", "community_type_donor", "syncom_donor", "donor_role", "recipient_role"], dropna=False)
            .agg(candidate_pairs=("pair_flux", "size"), unique_metabolites=("metabolite_harmonized", "nunique"), summed_pair_flux=("pair_flux", "sum"), mean_pair_flux=("pair_flux", "mean"))
            .reset_index()
        )
    summary.to_csv(output_dir / "02_interaction_summary_by_condition.csv", index=False)

    category_summary = (
        unique_pairs.groupby(["donor_role", "recipient_role"], dropna=False)
        .agg(
            unique_donor_recipient_metabolite_pairs=("pair_flux", "size"),
            unique_metabolites=("metabolite_harmonized", "nunique"),
            summed_pair_flux=("pair_flux", "sum"),
            mean_pair_flux=("pair_flux", "mean"),
        )
        .reset_index()
        if not unique_pairs.empty
        else pd.DataFrame(columns=["donor_role", "recipient_role", "unique_donor_recipient_metabolite_pairs", "unique_metabolites", "summed_pair_flux", "mean_pair_flux"])
    )
    categories = pd.DataFrame([
        {"donor_role": "bacteria", "recipient_role": "bacteria"},
        {"donor_role": "bacteria", "recipient_role": "fungus"},
        {"donor_role": "fungus", "recipient_role": "bacteria"},
        {"donor_role": "fungus", "recipient_role": "fungus"},
    ])
    category_summary = categories.merge(category_summary, on=["donor_role", "recipient_role"], how="left").fillna(0)
    category_summary.to_csv(output_dir / "04_interaction_category_summary.csv", index=False)

    global_columns = ["organism_donor", "organism_recipient", "metabolite_harmonized", "donor_role", "recipient_role"]
    global_unique = unique_pairs.drop_duplicates(global_columns) if not unique_pairs.empty else unique_pairs.copy()
    global_unique.to_csv(output_dir / "05_global_unique_donor_recipient_metabolite_pairs.csv", index=False)

    member_counts = data.groupby(["run_id", "role"], dropna=False).agg(
        member_exchange_rows=("metabolite_harmonized", "size"),
        unique_metabolites=("metabolite_harmonized", "nunique"),
        secretion_rows=("direction", lambda values: int((values == "secretion").sum())),
        uptake_rows=("direction", lambda values: int((values == "uptake").sum())),
        total_absolute_flux=("abs_flux", "sum"),
    ).reset_index()
    member_counts.to_csv(output_dir / "03_member_exchange_diagnostics.csv", index=False)

    metadata = {
        "source": str(Path(args.source).resolve()),
        "bacterial_source": str(BACTERIAL_SOURCE.resolve()),
        "threshold": args.threshold,
        "rows_after_shared52_filter": int(len(raw)),
        "rows_after_threshold": int(len(data)),
        "candidate_pair_rows": int(len(pairs)),
        "unique_candidate_pair_rows": int(len(unique_pairs)),
        "global_unique_candidate_pair_rows": int(len(global_unique)),
        "roles_included": ["bacteria", "fungus"],
        "medium_rows_excluded": True,
        "direction_rule": "positive flux = secretion; negative flux = uptake",
        "interpretation": "Candidate donor-recipient pairs require same run and harmonized metabolite identifier; they are not proof of transfer.",
    }
    (output_dir / "summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
