#!/usr/bin/env python3
"""Plot corrected fungal--bacterial metabolic niche overlap."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
SOURCE = PACKAGE / "reports/135_fungal_syncom_interface_metrics/02_shared_fungal_bacterial_uptake_metabolites.csv"
AUDIT = PACKAGE / "reports/133_fungal_bacterial_exchange_identifier_audit/01_active_exchange_fluxes_mapped.csv"
OUT = PACKAGE / "reports/137_fungal_bacterial_niche_overlap"
THESIS_FIGURES = PACKAGE.parent / "Writing_Master_Thesis/figures/results"

CONTEXT_ORDER = [
    ("HvSC1", "Bs", "B. sorokiniana", "Bs"),
    ("HvSC1", "Bs+Fg", "B. sorokiniana", "Bs + Fg: Bs"),
    ("HvSC1", "Bs+Fg", "F. graminearum", "Bs + Fg: Fg"),
    ("HvSC1", "Fg", "F. graminearum", "Fg"),
    ("HvSC2", "Bs", "B. sorokiniana", "Bs"),
    ("HvSC2", "Bs+Fg", "B. sorokiniana", "Bs + Fg: Bs"),
    ("HvSC2", "Bs+Fg", "F. graminearum", "Bs + Fg: Fg"),
    ("HvSC2", "Fg", "F. graminearum", "Fg"),
]


def load_matrix() -> tuple[pd.DataFrame, list[str]]:
    detail = pd.read_csv(SOURCE)
    audit = pd.read_csv(AUDIT, usecols=["chemical_key", "external_metabolite_name"])
    names = audit.dropna(subset=["external_metabolite_name"]).drop_duplicates("chemical_key")
    detail = detail.merge(names, on="chemical_key", how="left")
    detail["external_metabolite_name"] = detail["external_metabolite_name"].fillna(detail["chemical_key"])
    boundary_terms = ("water", "h2o", "dioxygen", "carbon dioxide", "proton", "h+")
    boundary_mask = detail["external_metabolite_name"].str.lower().apply(
        lambda value: any(term in value for term in boundary_terms)
    )
    detail = detail[~boundary_mask].copy()
    detail["profile"] = detail.apply(
        lambda row: f"{row['syncom']}|{row['fungus_context']}|{row['fungal_recipient']}", axis=1
    )
    matrix = detail.pivot_table(
        index=["chemical_key", "external_metabolite_name"],
        columns="profile",
        values="shared_min_uptake",
        aggfunc="sum",
        fill_value=0.0,
    )
    profile_keys = [f"{syncom}|{context}|{recipient}" for syncom, context, recipient, label in CONTEXT_ORDER]
    matrix = matrix.reindex(columns=profile_keys, fill_value=0.0)
    matrix["maximum"] = matrix.max(axis=1)
    matrix = matrix.sort_values("maximum", ascending=False).head(25).drop(columns="maximum")
    return matrix, profile_keys


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    THESIS_FIGURES.mkdir(parents=True, exist_ok=True)
    matrix, profile_keys = load_matrix()
    matrix.columns = [f"{syncom}\n{label}" for syncom, context, recipient, label in CONTEXT_ORDER]
    matrix.index = [str(name).capitalize() for _, name in matrix.index]
    matrix.to_csv(OUT / "01_fungal_bacterial_niche_overlap_matrix.csv")

    values = np.log10(1.0 + matrix.to_numpy(float))
    figure, axis = plt.subplots(figsize=(11.2, 8.0), dpi=180)
    image = axis.imshow(values, aspect="auto", cmap="YlGnBu")
    axis.set_xticks(np.arange(len(matrix.columns)), matrix.columns, rotation=40, ha="right")
    axis.set_yticks(np.arange(len(matrix.index)), matrix.index)
    axis.tick_params(axis="both", labelsize=7.5, length=0)
    axis.set_xlabel("Fungal recipient and bacterial SynCom context", fontweight="bold", labelpad=12)
    axis.set_ylabel("Metabolite jointly consumed by fungus and SynCom", fontweight="bold", labelpad=12)
    axis.set_title("Fungal--bacterial metabolic niche overlap", fontweight="bold", pad=12)
    axis.set_xticks(np.arange(-0.5, len(matrix.columns), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=0.8)
    axis.tick_params(which="minor", bottom=False, left=False)
    for spine in axis.spines.values():
        spine.set_visible(False)
    colorbar = figure.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
    colorbar.set_label("log$_{10}$(1 + pairwise minimum uptake flux)", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)
    figure.subplots_adjust(left=0.28, right=0.94, bottom=0.20, top=0.92)
    for path in (OUT / "fungal_bacterial_niche_overlap.png", THESIS_FIGURES / "fungal_bacterial_niche_overlap.png"):
        figure.savefig(path, dpi=400, bbox_inches="tight", facecolor="white")
    figure.savefig(OUT / "fungal_bacterial_niche_overlap.pdf", bbox_inches="tight", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
