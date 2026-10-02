#!/usr/bin/env python3
"""Audit fungal--bacterial exchange matching in the shared-biomass MICOM panel.

The earlier donor--recipient screen inspected only ``EX_*`` columns. This is
insufficient for the derived CarveFungi models, whose native exchange reactions
are predominantly named ``UF...``. This audit reads the model exchange reactions
for every member in each MICOM run, maps raw flux columns back to those reactions,
and derives a chemical key from the exchanged external metabolite.

Existing curated medium crosswalks are used to join KEGG compound identifiers to
their CarveFungi/biGG-style equivalents. Potential name/formula aliases not
covered by that crosswalk are exported for manual review and are *not* counted as
donor--recipient matches. Consequently, the candidate screen is conservative.

The script performs no optimisation and does not modify source simulations.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import cobra
import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
REPORT_ROOT = PACKAGE / "reports/96_fungal_species_syncom_factorial_panel"
MASTER = REPORT_ROOT / "summaries/master_run_index.csv"
CROSSWALK = PACKAGE / "reports/88_dual_namespace_carvefungi_media/dual_namespace_media_audit_long.csv"
DEFAULT_OUT = PACKAGE / "reports/133_fungal_bacterial_exchange_identifier_audit"


def strip_exchange_prefix_and_compartment(value: str) -> str:
    value = str(value).strip()
    value = value.removeprefix("R_").removeprefix("EX_")
    return re.sub(r"_(?:e|E|m|M)$", "", value)


def normalise_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        self.parent.setdefault(item, item)
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, first: str, second: str) -> None:
        first_root, second_root = self.find(first), self.find(second)
        if first_root != second_root:
            self.parent[second_root] = first_root


def crosswalk_clusters() -> tuple[UnionFind, set[str]]:
    """Return identifier clusters from the already curated medium crosswalk."""
    union_find = UnionFind()
    mapped_tokens: set[str] = set()
    if not CROSSWALK.exists():
        return union_find, mapped_tokens

    crosswalk = pd.read_csv(CROSSWALK)
    for _, row in crosswalk.iterrows():
        source = strip_exchange_prefix_and_compartment(row.get("mapped_from", ""))
        target = strip_exchange_prefix_and_compartment(row.get("carvefungi_metabolite", ""))
        if not source or not target or source.lower() == "nan" or target.lower() == "nan":
            continue
        source_token, target_token = f"id:{source.lower()}", f"id:{target.lower()}"
        union_find.union(source_token, target_token)
        mapped_tokens.update([source_token, target_token])
    return union_find, mapped_tokens


def model_exchange_index(model_path: Path) -> dict[str, dict[str, str]]:
    """Index model exchange reactions by their raw MICOM-export aliases."""
    model = cobra.io.read_sbml_model(str(model_path))
    aliases: dict[str, dict[str, str]] = {}
    for reaction in model.exchanges:
        metabolites = list(reaction.metabolites)
        external = next(
            (met for met in metabolites if str(met.id).endswith(("_e", "_E"))),
            metabolites[0] if metabolites else None,
        )
        if external is None:
            continue
        record = {
            "reaction_id": reaction.id,
            "external_metabolite_id": external.id,
            "external_metabolite_name": external.name or "",
            "formula": external.formula or "",
            "external_base_id": strip_exchange_prefix_and_compartment(external.id),
        }
        for alias in {reaction.id, reaction.id.removeprefix("R_")}:
            aliases[alias] = record
    return aliases


def flux_path(manifest: str, scenario: str) -> Path:
    stem = f"{manifest}_{scenario}_emma_style_fluxes.csv"
    return REPORT_ROOT / "raw_micom_runs" / manifest / scenario / stem


def collect_active_exchanges(master: pd.DataFrame, threshold: float, scenario: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    cache: dict[Path, dict[str, dict[str, str]]] = {}
    mapped_rows: list[dict[str, object]] = []
    unrecognised_rows: list[dict[str, object]] = []

    selected = master[(master["biomass_mode"] == "shared_52") & (master["scenario"] == scenario)].copy()
    for _, run in selected.iterrows():
        manifest = str(run["manifest"])
        members = pd.read_csv(PACKAGE / "manifests" / f"{manifest}.csv")
        member_info = {
            str(row["id"]): {"role": "fungus" if str(row["kind"]) == "fungus" else "bacteria", "file": Path(row["file"])}
            for _, row in members.iterrows()
        }
        indexes: dict[str, dict[str, dict[str, str]]] = {}
        for organism, info in member_info.items():
            path = info["file"]
            if path not in cache:
                cache[path] = model_exchange_index(path)
            indexes[organism] = cache[path]

        path = flux_path(manifest, scenario)
        fluxes = pd.read_csv(path)
        organism_column = fluxes.columns[0]
        for _, row in fluxes.iterrows():
            organism = str(row[organism_column])
            if organism == "medium" or organism not in member_info:
                continue
            role = member_info[organism]["role"]
            index = indexes[organism]
            for raw_exchange, value in row.iloc[1:].items():
                if raw_exchange not in index or pd.isna(value) or abs(float(value)) < threshold:
                    continue
                reaction = index[raw_exchange]
                base_id = reaction["external_base_id"]
                mapped_rows.append(
                    {
                        "run_id": run["run_id"],
                        "manifest": manifest,
                        "syncom": run["syncom"],
                        "fungus_context": run["fungus_context"],
                        "community_type": run["community_type"],
                        "organism": organism,
                        "role": role,
                        "raw_exchange": raw_exchange,
                        "reaction_id": reaction["reaction_id"],
                        "external_metabolite_id": reaction["external_metabolite_id"],
                        "external_metabolite_name": reaction["external_metabolite_name"],
                        "formula": reaction["formula"],
                        "external_base_id": base_id,
                        "flux": float(value),
                        "abs_flux": abs(float(value)),
                        "direction": "secretion" if float(value) > 0 else "uptake",
                    }
                )

            # Only external-looking exchange columns not found in the member's
            # actual exchange set are documented. Internal UF reactions use the
            # same prefix, but are deliberately outside this audit's scope.
            for raw_exchange, value in row.iloc[1:].items():
                external_like = str(raw_exchange).startswith("EX_") or str(raw_exchange).endswith("_E")
                if not external_like or raw_exchange in index or pd.isna(value) or abs(float(value)) < threshold:
                    continue
                unrecognised_rows.append(
                    {"run_id": run["run_id"], "organism": organism, "role": role, "raw_column": raw_exchange, "flux": float(value)}
                )
    return pd.DataFrame(mapped_rows), pd.DataFrame(unrecognised_rows)


def add_chemical_keys(data: pd.DataFrame, union_find: UnionFind, mapped_tokens: set[str]) -> pd.DataFrame:
    data = data.copy()
    data["identifier_token"] = data["external_base_id"].map(lambda value: f"id:{str(value).lower()}")
    data["chemical_key"] = data["identifier_token"].map(union_find.find)
    data["mapping_basis"] = data["identifier_token"].map(
        lambda token: "curated_crosswalk" if token in mapped_tokens else "exact_external_identifier"
    )
    return data


def candidate_pairs(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for donor_role, recipient_role in [("bacteria", "fungus"), ("fungus", "bacteria"), ("fungus", "fungus"), ("bacteria", "bacteria")]:
        donors = data[(data["role"] == donor_role) & (data["direction"] == "secretion")]
        recipients = data[(data["role"] == recipient_role) & (data["direction"] == "uptake")]
        merged = donors.merge(recipients, on=["run_id", "chemical_key"], suffixes=("_donor", "_recipient"))
        if donor_role == recipient_role:
            merged = merged[merged["organism_donor"] != merged["organism_recipient"]]
        if merged.empty:
            continue
        merged["donor_role"] = donor_role
        merged["recipient_role"] = recipient_role
        merged["pair_flux"] = merged[["abs_flux_donor", "abs_flux_recipient"]].min(axis=1)
        rows.append(merged)
    if not rows:
        return pd.DataFrame()
    pairs = pd.concat(rows, ignore_index=True)
    return pairs.sort_values(["donor_role", "recipient_role", "run_id", "pair_flux"], ascending=[True, True, True, False])


def unresolved_name_formula_aliases(data: pd.DataFrame) -> pd.DataFrame:
    """Flag plausible aliases for review without treating them as matches."""
    fungi = data[data["role"] == "fungus"].copy()
    bacteria = data[data["role"] == "bacteria"].copy()
    for frame in (fungi, bacteria):
        frame["name_key"] = frame["external_metabolite_name"].map(normalise_name)
        frame["formula_key"] = frame["formula"].fillna("").astype(str)
    # Joining first on run, name and formula avoids a very large all-by-all
    # bacteria--fungus merge while retaining only chemically plausible aliases.
    fungi = fungi[(fungi["name_key"] != "") & (fungi["formula_key"] != "")]
    bacteria = bacteria[(bacteria["name_key"] != "") & (bacteria["formula_key"] != "")]
    fungi = fungi.drop_duplicates(
        ["run_id", "organism", "name_key", "formula_key", "direction", "chemical_key"]
    )
    bacteria = bacteria.drop_duplicates(
        ["run_id", "organism", "name_key", "formula_key", "direction", "chemical_key"]
    )
    merged = bacteria.merge(
        fungi,
        on=["run_id", "name_key", "formula_key"],
        suffixes=("_bacteria", "_fungus"),
    )
    plausible = merged[
        (merged["direction_bacteria"] != merged["direction_fungus"])
        & (merged["chemical_key_bacteria"] != merged["chemical_key_fungus"])
    ].copy()
    keep = [
        "run_id", "organism_bacteria", "raw_exchange_bacteria", "external_metabolite_id_bacteria", "external_metabolite_name_bacteria", "flux_bacteria",
        "organism_fungus", "raw_exchange_fungus", "external_metabolite_id_fungus", "external_metabolite_name_fungus", "flux_fungus",
    ]
    return plausible[keep].drop_duplicates() if not plausible.empty else pd.DataFrame(columns=keep)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario", default="proline_relaxed_1", help="Existing panel scenario to audit.")
    parser.add_argument("--flux-threshold", type=float, default=1e-4)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(MASTER)
    union_find, mapped_tokens = crosswalk_clusters()
    active, unrecognised = collect_active_exchanges(master, args.flux_threshold, args.scenario)
    if active.empty:
        raise RuntimeError("No active model exchange reactions were found for the selected panel.")
    active = add_chemical_keys(active, union_find, mapped_tokens)
    active.to_csv(output / "01_active_exchange_fluxes_mapped.csv", index=False)
    unrecognised.to_csv(output / "02_unrecognised_exchange_like_flux_columns.csv", index=False)

    coverage = (
        active.groupby("role", dropna=False)
        .agg(
            active_flux_rows=("raw_exchange", "size"),
            organisms=("organism", "nunique"),
            reactions=("reaction_id", "nunique"),
            chemical_keys=("chemical_key", "nunique"),
            curated_crosswalk_rows=("mapping_basis", lambda values: int((values == "curated_crosswalk").sum())),
        )
        .reset_index()
    )
    coverage.to_csv(output / "03_mapping_coverage_by_role.csv", index=False)

    pairs = candidate_pairs(active)
    pairs.to_csv(output / "04_candidate_donor_recipient_pairs_reaudited.csv", index=False)
    if pairs.empty:
        category = pd.DataFrame(columns=["donor_role", "recipient_role", "candidate_pairs", "unique_metabolites"])
    else:
        category = (
            pairs.groupby(["donor_role", "recipient_role"], dropna=False)
            .agg(candidate_pairs=("pair_flux", "size"), unique_metabolites=("chemical_key", "nunique"))
            .reset_index()
        )
    categories = pd.DataFrame(
        [("bacteria", "bacteria"), ("bacteria", "fungus"), ("fungus", "bacteria"), ("fungus", "fungus")],
        columns=["donor_role", "recipient_role"],
    )
    category = categories.merge(category, on=["donor_role", "recipient_role"], how="left").fillna(0)
    category.to_csv(output / "05_candidate_categories_reaudited.csv", index=False)

    if pairs.empty:
        metabolite_summary = pd.DataFrame()
    else:
        metabolite_summary = (
            pairs.groupby(
                ["donor_role", "recipient_role", "chemical_key", "external_metabolite_name_donor"],
                dropna=False,
            )
            .agg(
                candidate_pairs=("pair_flux", "size"),
                runs=("run_id", "nunique"),
                donors=("organism_donor", "nunique"),
                recipients=("organism_recipient", "nunique"),
                maximum_pair_flux=("pair_flux", "max"),
            )
            .reset_index()
            .sort_values(["donor_role", "recipient_role", "maximum_pair_flux"], ascending=[True, True, False])
        )
    metabolite_summary.to_csv(output / "06_candidate_metabolite_summary_reaudited.csv", index=False)

    bacteria_to_fungus = pairs[(pairs["donor_role"] == "bacteria") & (pairs["recipient_role"] == "fungus")].copy()
    recipient_rows: list[dict[str, object]] = []
    for (syncom, recipient), group in bacteria_to_fungus.groupby(["syncom_donor", "organism_recipient"], dropna=False):
        globally_unique = group.drop_duplicates(["organism_donor", "organism_recipient", "chemical_key"])
        recipient_rows.append(
            {
                "syncom": syncom,
                "fungal_recipient": recipient,
                "run_specific_candidate_pairs": len(group),
                "globally_unique_donor_recipient_metabolite_combinations": len(globally_unique),
                "unique_chemical_keys": group["chemical_key"].nunique(),
                "unique_bacterial_donors": group["organism_donor"].nunique(),
                "community_runs": group["run_id"].nunique(),
            }
        )
    pd.DataFrame(recipient_rows).to_csv(output / "08_bacteria_to_fungus_by_syncom_and_fungal_recipient.csv", index=False)

    by_run = (
        bacteria_to_fungus.groupby(["run_id", "syncom_donor", "fungus_context_donor", "organism_recipient"], dropna=False)
        .agg(
            candidate_pairs=("pair_flux", "size"),
            unique_chemical_keys=("chemical_key", "nunique"),
            unique_bacterial_donors=("organism_donor", "nunique"),
        )
        .reset_index()
    )
    by_run.to_csv(output / "09_bacteria_to_fungus_by_community_context.csv", index=False)

    aliases = unresolved_name_formula_aliases(active)
    aliases.to_csv(output / "07_unresolved_name_formula_aliases_for_review.csv", index=False)

    metadata = {
        "source_panel": str(REPORT_ROOT),
        "scenario": args.scenario,
        "flux_threshold": args.flux_threshold,
        "active_exchange_flux_rows": int(len(active)),
        "active_fungal_exchange_reactions": int(active.loc[active["role"] == "fungus", "reaction_id"].nunique()),
        "active_bacterial_exchange_reactions": int(active.loc[active["role"] == "bacteria", "reaction_id"].nunique()),
        "unrecognised_exchange_like_flux_columns": int(len(unrecognised)),
        "unresolved_name_formula_aliases": int(len(aliases)),
        "interpretation": "Only exact external-metabolite identifiers and curated crosswalk relationships are matched. Name/formula aliases require manual review and are excluded from candidate exchange counts.",
    }
    (output / "summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (output / "README.md").write_text(
        "# Fungal--bacterial exchange identifier audit\n\n"
        "This audit checks whether candidate exchange screening can include every active fungal and bacterial exchange reaction. "
        "It uses model reaction stoichiometry and curated crosswalks, without rerunning MICOM.\n\n"
        "Candidate pairs remain flux-pattern hypotheses, not proof of transfer.\n",
        encoding="utf-8",
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
