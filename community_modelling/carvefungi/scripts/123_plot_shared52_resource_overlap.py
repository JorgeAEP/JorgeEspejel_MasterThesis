#!/usr/bin/env python3
"""Create a publication-style shared-52 resource-overlap figure."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
REPORT = PACKAGE / "reports/122_shared52_resource_allocation"
THESIS_FIGURES = PACKAGE.parent / "Writing_Master_Thesis/figures/results"
OUT = REPORT / "figures"

DUAL_RUNS = {
    "Bs + Fg": "fungi_pair__carvefungi_52__proline_relaxed_1",
    "Bs + Fg + HvSC1": "hvsc1_plus_fungi_pair__carvefungi_52__proline_relaxed_1",
    "Bs + Fg + HvSC2": "hvsc2_plus_fungi_pair__carvefungi_52__proline_relaxed_1",
}
SYNCOM_CONTEXTS = [
    ("Bs + HvSC1", "HvSC1_Bs"),
    ("Fg + HvSC1", "HvSC1_Fg"),
    ("Bs + Fg + HvSC1", "HvSC1_Bs+Fg"),
    ("Bs + HvSC2", "HvSC2_Bs"),
    ("Fg + HvSC2", "HvSC2_Fg"),
    ("Bs + Fg + HvSC2", "HvSC2_Bs+Fg"),
]


def compact_name(name: str) -> str:
    replacements = {
        "hydrogenphosphate": "Phosphate",
        "dioxygen": "Oxygen",
        "citrate(3-)": "Citrate",
        "water": "Water",
        "alpha-D-glucose": "Glucose",
        "S-adenosyl-L-methionine zwitterion": "SAM",
        "L-argininium": "Arginine",
        "hydroquinone": "Hydroquinone",
        "L-serine zwitterion": "Serine",
        "D-fructose": "Fructose",
    }
    return replacements.get(name, name)


def load_fungal_overlap() -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(REPORT / "13_bipolaris_fusarium_resource_overlap.csv")
    summary = summary.set_index("run_id").loc[list(DUAL_RUNS.values())].reset_index()
    summary["context"] = list(DUAL_RUNS)

    uptake = pd.read_csv(REPORT / "11_member_consumption_ranked.csv")
    uptake = uptake[
        uptake["run_id"].isin(DUAL_RUNS.values()) & (uptake["role"] == "fungus")
    ].copy()
    uptake["context"] = uptake["run_id"].map({v: k for k, v in DUAL_RUNS.items()})

    rows = []
    for context, group in uptake.groupby("context"):
        profiles = group.pivot_table(
            index=["metabolite", "metabolite_name"],
            columns="fungus",
            values="uptake_flux",
            aggfunc="sum",
            fill_value=0.0,
        )
        if not {"Bs", "Fg"}.issubset(profiles.columns):
            continue
        for (metabolite, name), row in profiles.iterrows():
            bs = float(row.get("Bs", 0.0))
            fg = float(row.get("Fg", 0.0))
            if bs > 0 and fg > 0:
                rows.append(
                    {
                        "context": context,
                        "metabolite": metabolite,
                        "name": compact_name(name),
                        "overlap": min(bs, fg),
                    }
                )
    detail = pd.DataFrame(rows)
    matrix = detail.pivot_table(
        index=["metabolite", "name"],
        columns="context",
        values="overlap",
        fill_value=0.0,
    )
    matrix = matrix.reindex(columns=list(DUAL_RUNS), fill_value=0.0)
    matrix["maximum"] = matrix.max(axis=1)
    # Keep the strongest overlap patterns for a readable main-text figure.
    matrix = (
        matrix.sort_values("maximum", ascending=False)
        .head(25)
        .drop(columns="maximum")
    )
    return summary, matrix


def load_bacterial_overlap() -> pd.DataFrame:
    detail = pd.read_csv(REPORT / "07_shared_resource_use_by_member_metabolite.csv")
    detail["context"] = detail["syncom"] + "_" + detail["fungus_context"]
    detail = (
        detail.groupby(["context", "metabolite", "metabolite_name"], as_index=False)
        .agg(overlap=("shared_min_flux", "sum"))
    )
    selected = ["acald", "gua", "ura", "orot", "ade", "thm"]
    matrix = detail[detail["metabolite"].isin(selected)].pivot_table(
        index="context", columns="metabolite", values="overlap", fill_value=0.0
    )
    matrix = matrix.reindex(
        index=[key for _, key in SYNCOM_CONTEXTS], columns=selected, fill_value=0.0
    )
    matrix.columns = [
        "Acetaldehyde",
        "Guanine",
        "Uracil",
        "Orotate",
        "Adenine",
        "Thiamine",
    ]
    matrix.index = [label for label, _ in SYNCOM_CONTEXTS]
    return matrix


def heatmap(ax, values, row_labels, col_labels, cmap, colorbar_label, log=True):
    plot_values = np.log10(1.0 + values) if log else values
    image = ax.imshow(plot_values, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(col_labels)), labels=col_labels, rotation=40, ha="right")
    ax.set_yticks(np.arange(len(row_labels)), labels=row_labels)
    ax.tick_params(axis="both", labelsize=8, length=0)
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = ax.figure.colorbar(image, ax=ax, fraction=0.045, pad=0.03)
    colorbar.ax.tick_params(labelsize=7)
    colorbar.set_label(colorbar_label, fontsize=8)
    return image


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    THESIS_FIGURES.mkdir(parents=True, exist_ok=True)
    summary, fungal_matrix = load_fungal_overlap()
    bacterial_matrix = load_bacterial_overlap()

    summary.to_csv(OUT / "01_fungal_fungal_summary.csv", index=False)
    fungal_matrix.to_csv(OUT / "02_fungal_fungal_heatmap.csv")
    bacterial_matrix.to_csv(OUT / "03_fungal_bacterial_heatmap.csv")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
        }
    )
    fungal_values = fungal_matrix.to_numpy()
    fungal_height = max(8.0, 2.4 + 0.23 * len(fungal_matrix))
    fungal_figure, fungal_axis = plt.subplots(figsize=(7.4, fungal_height), dpi=180)
    fungal_figure.subplots_adjust(left=0.28, right=0.88, bottom=0.12, top=0.94)
    heatmap(
        fungal_axis,
        fungal_values,
        [name for _, name in fungal_matrix.index],
        ["Bs + Fg", "Bs + Fg\n+ HvSC1", "Bs + Fg\n+ HvSC2"],
        "PuRd",
        r"log$_{10}$(1 + minimum uptake flux)",
    )
    fungal_axis.set_title("Fungal--fungal resource overlap", fontweight="bold", pad=38)
    shared_annotations = [
        (34, 0.580),
        (41, 0.626),
        (49, 0.633),
    ]
    for column, (shared_count, cosine) in enumerate(shared_annotations):
        fungal_axis.text(
            column,
            1.01,
            f"{shared_count} metabolites\ncosine = {cosine:.3f}",
            transform=fungal_axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=7.5,
            linespacing=1.05,
        )
    for label in fungal_axis.get_yticklabels():
        label.set_fontsize(6.5)

    bacterial_values = bacterial_matrix.to_numpy()
    bacterial_figure, bacterial_axis = plt.subplots(figsize=(8.0, 5.8), dpi=180)
    bacterial_figure.subplots_adjust(left=0.25, right=0.88, bottom=0.22, top=0.88)
    heatmap(
        bacterial_axis,
        bacterial_values,
        list(bacterial_matrix.index),
        list(bacterial_matrix.columns),
        "PuRd",
        r"log$_{10}$(1 + summed overlap flux)",
    )
    bacterial_axis.set_title("Fungal--bacterial resource overlap", fontweight="bold", pad=12)

    fungal_paths = [
        OUT / "fungal_fungal_resource_overlap.png",
        THESIS_FIGURES / "fungal_fungal_resource_overlap.png",
    ]
    bacterial_paths = [
        OUT / "fungal_bacterial_resource_overlap.png",
        THESIS_FIGURES / "fungal_bacterial_resource_overlap.png",
    ]
    for path in fungal_paths:
        fungal_figure.savefig(path, dpi=400, bbox_inches="tight", facecolor="white")
    for path in bacterial_paths:
        bacterial_figure.savefig(path, dpi=400, bbox_inches="tight", facecolor="white")
    fungal_figure.savefig(OUT / "fungal_fungal_resource_overlap.pdf", bbox_inches="tight")
    bacterial_figure.savefig(OUT / "fungal_bacterial_resource_overlap.pdf", bbox_inches="tight")
    plt.close(fungal_figure)
    plt.close(bacterial_figure)
    print(
        {
            "fungal_figure": str(THESIS_FIGURES / "fungal_fungal_resource_overlap.png"),
            "bacterial_figure": str(THESIS_FIGURES / "fungal_bacterial_resource_overlap.png"),
            "fungal_heatmap_rows": int(len(fungal_matrix)),
            "bacterial_heatmap_rows": int(len(bacterial_matrix)),
            "bacterial_heatmap_columns": int(len(bacterial_matrix.columns)),
        }
    )


if __name__ == "__main__":
    main()
