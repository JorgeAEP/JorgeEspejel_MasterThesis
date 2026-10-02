#!/usr/bin/env python3
"""Plot non-medium bacteria-to-fungus exchange candidates for the thesis.

The plot aggregates donor-level records by candidate metabolite and fungal
recipient. Bubble area shows the smaller of total bacterial secretion and
fungal uptake, avoiding inflation when several bacteria secrete the same
metabolite. A black outline marks candidates with similar total secretion and
uptake scales.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
SOURCE = PACKAGE / "reports/134_exchange_candidate_classification/05_shared_pool_compatible_nonmedium_candidates.csv"
OUT = PACKAGE / "reports/136_fungal_bacterial_exchange_candidate_figure"
THESIS_FIGURES = PACKAGE.parent / "Writing_Master_Thesis/figures/results"

SYNCOMS = ("HvSC1", "HvSC2")
CONTEXT_ORDER = ("Bs", "Bs+Fg", "Fg")
CONTEXT_LABELS = {
    ("Bs", "B. sorokiniana"): "Bs",
    ("Bs+Fg", "B. sorokiniana"): "Bs in pair",
    ("Bs+Fg", "F. graminearum"): "Fg in pair",
    ("Fg", "F. graminearum"): "Fg",
}
METABOLITE_LABELS = {
    "id:2mbald": "2-Methylbutyraldehyde",
    "id:acald": "Acetaldehyde",
    "id:cytd": "Cytidine",
    "id:gal": "D-Galactose",
    "id:gcald": "Glycolaldehyde",
    "id:gua": "Guanine",
    "id:hqn": "Hydroquinone",
    "id:hxan": "Hypoxanthine",
    "id:inost": "Myo-inositol",
    "id:ins": "Inosine",
    "id:ump": "UMP",
    "id:ura": "Uracil",
    "id:uri": "Uridine",
}
COLORS = {"B. sorokiniana": "#188977", "F. graminearum": "#CC6A4B"}


def fungal_label(organism: str) -> str:
    if "Bipolaris" in organism:
        return "B. sorokiniana"
    if "Fusarium" in organism:
        return "F. graminearum"
    raise ValueError(f"Unrecognised fungal recipient: {organism}")


def load_candidates() -> pd.DataFrame:
    raw = pd.read_csv(SOURCE)
    raw["fungal_recipient"] = raw["organism_recipient"].map(fungal_label)
    grouped = (
        raw.groupby(
            ["syncom_donor", "fungus_context_donor", "fungal_recipient", "chemical_key"],
            as_index=False,
        )
        .agg(
            bacterial_total_secretion=("bacterial_total_secretion", "first"),
            fungal_total_uptake=("fungal_total_uptake", "first"),
            candidate_donor_count=("organism_donor", "nunique"),
            has_similar_flux_scale=("similar_flux_scale", "any"),
        )
        .rename(columns={"syncom_donor": "syncom", "fungus_context_donor": "fungus_context"})
    )
    grouped["candidate_flux"] = grouped[["bacterial_total_secretion", "fungal_total_uptake"]].min(axis=1)
    grouped["context_label"] = grouped.apply(
        lambda row: CONTEXT_LABELS[(row["fungus_context"], row["fungal_recipient"])], axis=1
    )
    grouped["metabolite"] = grouped["chemical_key"].map(METABOLITE_LABELS)
    return grouped


def bubble_size(values: pd.Series) -> np.ndarray:
    # The candidate-flux range spans several orders of magnitude.
    return 28 + 42 * (np.log10(values.to_numpy() + 1e-4) + 4)


def draw_figure(data: pd.DataFrame) -> plt.Figure:
    metabolites = (
        data.groupby(["chemical_key", "metabolite"])["candidate_flux"].max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    y_positions = {key: index for index, (key, _) in enumerate(metabolites)}
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 6.6), sharey=True, dpi=180)

    for axis, syncom in zip(axes, SYNCOMS):
        panel = data[data["syncom"] == syncom].copy()
        contexts = [
            label
            for context in CONTEXT_ORDER
            for label in ("Bs", "Bs in pair", "Fg in pair", "Fg")
            if any(
                (panel["fungus_context"] == context)
                & (panel["context_label"] == label)
            )
        ]
        # The comprehension above preserves biological ordering and omits labels
        # that cannot occur for a given fungal context.
        contexts = list(dict.fromkeys(contexts))
        x_positions = {label: index for index, label in enumerate(contexts)}

        for _, row in panel.iterrows():
            edgecolor = "#222222" if row["has_similar_flux_scale"] else "white"
            axis.scatter(
                x_positions[row["context_label"]],
                y_positions[row["chemical_key"]],
                s=bubble_size(pd.Series([row["candidate_flux"]]))[0],
                color=COLORS[row["fungal_recipient"]],
                edgecolor=edgecolor,
                linewidth=0.8,
                alpha=0.92,
                zorder=3,
            )

        axis.set_title(syncom, fontweight="bold", y=1.03)
        axis.set_xticks(range(len(contexts)), contexts, rotation=35, ha="right")
        axis.set_xlim(-0.6, len(contexts) - 0.4)
        axis.set_ylim(len(metabolites) - 0.5, -0.5)
        axis.grid(axis="both", color="#DFE4E2", linewidth=0.75)
        axis.set_axisbelow(True)
        axis.spines[["top", "right", "left", "bottom"]].set_visible(False)
        axis.tick_params(axis="both", length=0, labelsize=8)

    axes[0].set_yticks(range(len(metabolites)), [label for _, label in metabolites])
    axes[0].set_ylabel("Non-medium candidate metabolite", fontweight="bold")
    axes[1].tick_params(axis="y", labelleft=False)
    figure.text(
        0.5,
        0.055,
        "Fungal recipient and community context",
        ha="center",
        va="center",
        fontweight="bold",
    )
    size_values = [0.1, 10, 100]
    handles = [
        Line2D([0], [0], marker="o", color="w", label=r"$B.\ sorokiniana$ recipient", markerfacecolor=COLORS["B. sorokiniana"], markersize=7),
        Line2D([0], [0], marker="o", color="w", label=r"$F.\ graminearum$ recipient", markerfacecolor=COLORS["F. graminearum"], markersize=7),
        Line2D([0], [0], marker="o", color="#222222", label="Similar secretion and uptake", markerfacecolor="white", markersize=7, linewidth=0),
    ]
    handles.extend(
        Line2D(
            [0], [0], marker="o", color="w", label=f"Flux {value:g}",
            markerfacecolor="#8A8A8A", markeredgecolor="white",
            markersize=np.sqrt(bubble_size(pd.Series([value]))[0]),
        )
        for value in size_values
    )
    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.015),
        frameon=False,
        fontsize=7.6,
        columnspacing=1.6,
        handletextpad=0.55,
    )
    figure.subplots_adjust(left=0.23, right=0.99, bottom=0.22, top=0.92, wspace=0.06)
    return figure


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    THESIS_FIGURES.mkdir(parents=True, exist_ok=True)
    data = load_candidates()
    data.to_csv(OUT / "01_exchange_candidate_fluxes_by_context.csv", index=False)
    figure = draw_figure(data)
    for path in (
        OUT / "fungal_bacterial_exchange_candidates.png",
        THESIS_FIGURES / "fungal_bacterial_exchange_candidates.png",
    ):
        figure.savefig(path, dpi=400, bbox_inches="tight", facecolor="white")
    figure.savefig(OUT / "fungal_bacterial_exchange_candidates.pdf", bbox_inches="tight", facecolor="white")
    plt.close(figure)
    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "source": str(SOURCE),
                "candidate_metabolites": int(data["chemical_key"].nunique()),
                "fungal_recipient_contexts": int(data[["syncom", "fungus_context", "fungal_recipient"]].drop_duplicates().shape[0]),
                "interpretation": "Bubble area is the smaller of total bacterial secretion and fungal uptake for non-medium, shared-pool-compatible candidate metabolites. It does not prove direct transfer.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
