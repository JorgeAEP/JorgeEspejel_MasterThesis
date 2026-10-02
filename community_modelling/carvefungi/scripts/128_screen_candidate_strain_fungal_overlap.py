#!/usr/bin/env python3
"""Screen literature-prioritized SynCom strains for fungal resource overlap.

This is a post-processing analysis of Script 122 outputs. It does not rerun
MICOM. Overlap means co-consumption of the same harmonised metabolite in the
same exported community solution; it is not proof of competition, transfer or
fungal suppression.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
SOURCE = PACKAGE / "reports/122_shared52_resource_allocation"
TAXONOMY = PACKAGE.parent / "Bipolaris__sorokiniana/08_reference_tables/syncom_metadata/SynComs_taxonomy.csv"
OUT = PACKAGE / "reports/128_candidate_strain_fungal_overlap"

DEFAULT_STRAINS = ["397", "1338", "2774", "163", "793", "790", "997", "892"]


def parse_strain_ids(values: list[str] | None) -> list[str]:
    raw = values or [",".join(DEFAULT_STRAINS)]
    ids: list[str] = []
    for value in raw:
        ids.extend(part.strip() for part in value.split(",") if part.strip())
    invalid = [value for value in ids if not value.isdigit()]
    if invalid:
        raise ValueError(f"Strain IDs must be numeric: {', '.join(invalid)}")
    return list(dict.fromkeys(ids))


def strain_id_from_member(value: object) -> str:
    match = re.match(r"^(\d+)_", str(value))
    return match.group(1) if match else ""


def load_taxonomy() -> pd.DataFrame:
    taxonomy = pd.read_csv(TAXONOMY, dtype=str).fillna("")
    taxonomy.columns = [str(column).strip() for column in taxonomy.columns]
    taxonomy["Strain ID"] = taxonomy["Strain ID"].astype(str).str.strip()
    keep = [
        "Strain ID",
        "Syncom",
        "Genus",
        "Species",
        "Phylum",
        "Class",
        "Order",
        "Family",
    ]
    return taxonomy[[column for column in keep if column in taxonomy.columns]].rename(
        columns={"Strain ID": "strain_id", "Syncom": "taxonomy_syncom"}
    )


def add_rank(overlap: pd.DataFrame) -> pd.DataFrame:
    """Rank candidates against all bacterial members in each fungal context."""
    out = overlap.copy()
    rank_keys = ["run_id", "fungal_member"]
    out["overlap_rank"] = (
        out.groupby(rank_keys)["weighted_overlap_min_flux"]
        .rank(method="min", ascending=False)
        .astype("Int64")
    )
    out["members_in_context"] = out.groupby(rank_keys)["bacterial_member"].transform("nunique")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strain-id",
        action="append",
        help="Candidate strain ID; repeat or provide comma-separated IDs.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top shared metabolites retained per candidate/context (default: 5).",
    )
    args = parser.parse_args()
    if args.top_n < 1:
        raise ValueError("--top-n must be positive")

    strain_ids = parse_strain_ids(args.strain_id)
    overlap_path = SOURCE / "06_bacterial_fungal_overlap_by_member.csv"
    shared_path = SOURCE / "07_shared_resource_use_by_member_metabolite.csv"
    if not overlap_path.exists() or not shared_path.exists():
        raise FileNotFoundError(
            "Script 122 outputs were not found. Run scripts/122_analyse_shared52_resource_allocation.py first."
        )

    OUT.mkdir(parents=True, exist_ok=True)
    overlap = pd.read_csv(overlap_path)
    shared = pd.read_csv(shared_path)
    overlap["strain_id"] = overlap["bacterial_member"].map(strain_id_from_member)
    shared["strain_id"] = shared["bacterial_member"].map(strain_id_from_member)

    taxonomy = load_taxonomy()
    ranked_overlap = add_rank(overlap)
    selected_overlap = ranked_overlap[ranked_overlap["strain_id"].isin(strain_ids)].copy()
    selected_shared = shared[shared["strain_id"].isin(strain_ids)].copy()

    selected_overlap = selected_overlap.merge(taxonomy, on="strain_id", how="left")
    selected_shared = selected_shared.merge(taxonomy, on="strain_id", how="left")

    overlap_columns = [
        "strain_id",
        "taxonomy_syncom",
        "Genus",
        "Species",
        "run_id",
        "syncom",
        "fungus_context",
        "community_type",
        "fungal_member",
        "bacterial_member",
        "shared_metabolites",
        "weighted_overlap_min_flux",
        "weighted_overlap_cosine",
        "overlap_fraction_of_fungal_uptake",
        "overlap_fraction_of_bacterial_uptake",
        "fungal_total_uptake",
        "bacterial_total_uptake",
        "overlap_rank",
        "members_in_context",
    ]
    selected_overlap[overlap_columns].sort_values(
        ["strain_id", "fungus_context", "syncom", "fungal_member"]
    ).to_csv(OUT / "01_candidate_overlap_summary.csv", index=False)

    shared_columns = [
        "strain_id",
        "taxonomy_syncom",
        "Genus",
        "Species",
        "run_id",
        "syncom",
        "fungus_context",
        "community_type",
        "fungal_member",
        "bacterial_member",
        "metabolite",
        "metabolite_name",
        "fungal_uptake_flux",
        "bacterial_uptake_flux",
        "shared_min_flux",
    ]
    selected_shared = selected_shared.sort_values(
        ["strain_id", "run_id", "fungal_member", "shared_min_flux"],
        ascending=[True, True, True, False],
    )
    selected_shared["metabolite_rank_within_candidate"] = selected_shared.groupby(
        ["strain_id", "run_id", "fungal_member"]
    )["shared_min_flux"].rank(method="first", ascending=False).astype("Int64")
    selected_shared = selected_shared[
        selected_shared["metabolite_rank_within_candidate"] <= args.top_n
    ]
    selected_shared[shared_columns + ["metabolite_rank_within_candidate"]].to_csv(
        OUT / "02_candidate_shared_metabolites.csv", index=False
    )

    missing = [strain_id for strain_id in strain_ids if strain_id not in set(selected_overlap["strain_id"])]
    summary = {
        "source": str(SOURCE),
        "threshold": 1e-4,
        "candidate_strains": strain_ids,
        "overlap_rows": int(len(selected_overlap)),
        "shared_metabolite_rows_retained": int(len(selected_shared)),
        "candidate_strains_without_overlap_rows": missing,
        "interpretation": (
            "Overlap is harmonised co-consumption in the same MICOM solution. "
            "It is not evidence of causal competition, metabolite transfer or fungal suppression."
        ),
    }
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    report_lines = [
        "# Candidate-strain fungal resource-overlap screen",
        "",
        "This report reuses Script 122 member-level uptake results. It does not rerun MICOM.",
        "Overlap means that a bacterial member and a fungal member both retained uptake",
        "for the same harmonised metabolite in the same exported community solution.",
        "It is not proof of competition, transfer, suppression or causality.",
        "",
        "## Candidate interpretation",
        "",
    ]
    if selected_overlap.empty:
        report_lines.append("No selected candidate strain occurred in the Script 122 overlap table.")
    else:
        compact = (
            selected_overlap.groupby(["strain_id", "Genus", "Species"], dropna=False)
            .agg(
                overlap_rows=("strain_id", "size"),
                maximum_overlap=("weighted_overlap_min_flux", "max"),
                maximum_fungal_fraction=("overlap_fraction_of_fungal_uptake", "max"),
                best_rank=("overlap_rank", "min"),
            )
            .reset_index()
            .sort_values("maximum_overlap", ascending=False)
        )
        for _, row in compact.iterrows():
            report_lines.append(
                f"- Strain {row['strain_id']} ({row['Genus']} {row['Species']}): "
                f"maximum overlap score {row['maximum_overlap']:.4g}, "
                f"maximum fungal-uptake fraction {100 * row['maximum_fungal_fraction']:.3g}%, "
                f"best rank {int(row['best_rank'])}."
            )
    report_lines.extend(
        [
            "",
            "## Suggested next step",
            "",
            "Prioritise a perturbation only when a candidate combines substantial overlap",
            "with a clear biological rationale. A dropout or add-in MICOM experiment is",
            "still needed to test whether the candidate changes community or fungal growth.",
        ]
    )
    (OUT / "candidate_overlap_report.md").write_text("\n".join(report_lines) + "\n")

    print(json.dumps(summary, indent=2))
    print(f"outputs: {OUT}")


if __name__ == "__main__":
    main()
