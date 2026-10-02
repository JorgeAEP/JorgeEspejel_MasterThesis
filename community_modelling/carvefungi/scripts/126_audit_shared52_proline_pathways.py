#!/usr/bin/env python3
"""Audit proline reactions and fluxes in the existing shared-52 MICOM runs.

The Script 96 flux exports contain full taxon-by-reaction flux matrices, not
only exchange reactions. This script uses those exports to connect:

* the proline coefficient in each biomass reaction;
* fungal biomass flux and reported member growth;
* intracellular proline-related reaction fluxes;
* proline transport and external uptake; and
* the biomass-associated proline requirement.

It does not infer intracellular synthesis from exchange summaries and does not
modify models or rerun MICOM.
"""

from __future__ import annotations

import json
from pathlib import Path

import cobra
import pandas as pd


SCRIPT = Path(__file__).resolve()
PACKAGE = SCRIPT.parents[1]
REPORTS = PACKAGE / "reports"
RAW = REPORTS / "96_fungal_species_syncom_factorial_panel/raw_micom_runs"
OUT = REPORTS / "126_shared52_proline_pathway_audit"
SCENARIO = "proline_relaxed_1"

MODELS = {
    "Bs": {
        "organism_token": "Bipolaris_sorokiniana",
        "path": PACKAGE / "fungal_models/biomass_variant_panel/Bipolaris_sorokiniana_carvefungi_52biomass.xml",
        "biomass": "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF",
        "standalone_capacity": 952.2298245614036,
    },
    "Fg": {
        "organism_token": "Fusarium_graminearum",
        "path": PACKAGE / "fungal_models/biomass_variant_panel/Fusarium_graminearum_carvefungi_52biomass.xml",
        "biomass": "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF",
        "standalone_capacity": 1.0,
    },
}

CONTEXTS = {
    "fungus_alone__bipolaris_carvefungi_52": "Bs",
    "hvsc1_plus_bipolaris_carvefungi_52": "Bs",
    "hvsc2_plus_bipolaris_carvefungi_52": "Bs",
    "fungus_alone__fusarium_carvefungi_52": "Fg",
    "hvsc1_plus_fusarium_carvefungi_52": "Fg",
    "hvsc2_plus_fusarium_carvefungi_52": "Fg",
    "fungi_pair__carvefungi_52": "both",
    "hvsc1_plus_fungi_pair__carvefungi_52": "both",
    "hvsc2_plus_fungi_pair__carvefungi_52": "both",
}


def model_audit(label: str, config: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    model = cobra.io.read_sbml_model(str(config["path"]))
    biomass = model.reactions.get_by_id(str(config["biomass"]))
    proline_terms = {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in biomass.metabolites.items()
        if "prol" in metabolite.id.lower()
    }
    rows = []
    for metabolite_id, coefficient in proline_terms.items():
        rows.append({
            "fungus": label,
            "biomass_reaction": biomass.id,
            "biomass_proline_metabolite": metabolite_id,
            "biomass_proline_coefficient": coefficient,
            "standalone_maximum_demand_capacity": config["standalone_capacity"],
        })

    reaction_rows = []
    for reaction in model.reactions:
        proline_metabolites = [
            metabolite for metabolite in reaction.metabolites
            if "prol" in metabolite.id.lower()
        ]
        if not proline_metabolites:
            continue
        if reaction.id == biomass.id:
            role = "biomass"
        elif reaction.id.endswith("_E") or reaction.id.startswith("EX_"):
            role = "external_exchange"
        elif "TCE" in reaction.id or "TCM" in reaction.id or "transport" in reaction.name.lower():
            role = "transport_or_export"
        elif any(reaction.metabolites[metabolite] > 0 for metabolite in proline_metabolites):
            role = "candidate_intracellular_production"
        else:
            role = "proline_consumption_or_processing"
        reaction_rows.append({
            "fungus": label,
            "reaction": reaction.id,
            "reaction_name": reaction.name,
            "role": role,
            "lower_bound": reaction.lower_bound,
            "upper_bound": reaction.upper_bound,
            "equation": reaction.reaction,
            "proline_metabolites": ";".join(metabolite.id for metabolite in proline_metabolites),
        })
    return rows, reaction_rows


def flux_path(manifest: str) -> Path:
    return RAW / manifest / SCENARIO / f"{manifest}_{SCENARIO}_emma_style_fluxes.csv"


def growth_path(manifest: str) -> Path:
    return RAW / manifest / SCENARIO / f"{manifest}_{SCENARIO}_emma_style_individualtaxon_growth.csv"


def extract_context(manifest: str, context: str, models: dict[str, cobra.Model]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    path = flux_path(manifest)
    growth_file = growth_path(manifest)
    if not path.exists() or not growth_file.exists():
        raise FileNotFoundError(f"Missing shared-52 export for {manifest}: {path}")
    fluxes = pd.read_csv(path)
    growth = pd.read_csv(growth_file)
    organism_column = fluxes.columns[0]
    growth_column = growth.columns[0]
    summary_rows = []
    reaction_rows = []
    labels = [context] if context != "both" else ["Bs", "Fg"]
    for label in labels:
        config = MODELS[label]
        token = str(config["organism_token"])
        matching = fluxes[fluxes[organism_column].astype(str).str.contains(token, regex=False)]
        if matching.empty:
            raise ValueError(f"No {label} taxon row found in {path}")
        flux_row = matching.iloc[0]
        growth_matching = growth[growth[growth_column].astype(str).str.contains(token, regex=False)]
        member_growth = float(growth_matching.iloc[0]["growth_rate"]) if not growth_matching.empty else None
        biomass = models[label].reactions.get_by_id(str(config["biomass"]))
        biomass_flux = float(flux_row.get(biomass.id, 0.0))
        coefficient = float(next(
            coefficient for metabolite, coefficient in biomass.metabolites.items()
            if "prol" in metabolite.id.lower()
        ))
        biomass_requirement = biomass_flux * abs(coefficient)

        model_reactions = []
        for reaction in models[label].reactions:
            if any("prol" in metabolite.id.lower() for metabolite in reaction.metabolites):
                value = flux_row.get(reaction.id, 0.0)
                value = None if pd.isna(value) else float(value)
                model_reactions.append((reaction, value))
                reaction_rows.append({
                    "manifest": manifest,
                    "context": context,
                    "fungus": label,
                    "reaction": reaction.id,
                    "reaction_name": reaction.name,
                    "flux": value,
                    "role": classify_reaction(reaction, biomass.id),
                })

        synthesis = [value for reaction, value in model_reactions if classify_reaction(reaction, biomass.id) == "candidate_intracellular_production" and value is not None]
        transport = [value for reaction, value in model_reactions if classify_reaction(reaction, biomass.id) == "transport_or_export" and value is not None]
        exchange = [value for reaction, value in model_reactions if classify_reaction(reaction, biomass.id) == "external_exchange" and value is not None]
        uptake = sum(max(0.0, -value) for value in exchange)
        summary_rows.append({
            "manifest": manifest,
            "context": context,
            "fungus": label,
            "member_growth": member_growth,
            "biomass_reaction_flux": biomass_flux,
            "biomass_proline_coefficient": coefficient,
            "biomass_associated_proline_requirement": biomass_requirement,
            "external_proline_uptake": uptake,
            "candidate_intracellular_production_flux_sum": sum(max(0.0, value) for value in synthesis),
            "candidate_intracellular_production_reaction_count": len(synthesis),
            "proline_transport_or_export_flux_sum": sum(transport),
            "standalone_demand_capacity": config["standalone_capacity"],
            "requirement_to_standalone_demand_capacity": biomass_requirement / config["standalone_capacity"],
            "interpretation_note": "Reaction-level flux audit; pathway limitation requires a dedicated constrained rerun.",
        })
    return summary_rows, reaction_rows


def classify_reaction(reaction: cobra.Reaction, biomass_id: str) -> str:
    if reaction.id == biomass_id:
        return "biomass"
    if reaction.id.endswith("_E") or reaction.id.startswith("EX_"):
        return "external_exchange"
    if "TCE" in reaction.id or "TCM" in reaction.id or "transport" in reaction.name.lower():
        return "transport_or_export"
    proline_metabolites = [m for m in reaction.metabolites if "prol" in m.id.lower()]
    if any(reaction.metabolites[m] > 0 for m in proline_metabolites):
        return "candidate_intracellular_production"
    return "proline_consumption_or_processing"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audits = []
    reaction_audits = []
    models = {}
    for label, config in MODELS.items():
        models[label] = cobra.io.read_sbml_model(str(config["path"]))
        audit, reactions = model_audit(label, config)
        audits.extend(audit)
        reaction_audits.extend(reactions)

    summary_rows = []
    flux_rows = []
    for manifest, context in CONTEXTS.items():
        current_summary, current_flux = extract_context(manifest, context, models)
        summary_rows.extend(current_summary)
        flux_rows.extend(current_flux)

    pd.DataFrame(audits).to_csv(OUT / "01_model_proline_biomass_audit.csv", index=False)
    pd.DataFrame(reaction_audits).to_csv(OUT / "02_model_proline_reaction_audit.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(OUT / "03_shared52_proline_flux_summary.csv", index=False)
    pd.DataFrame(flux_rows).to_csv(OUT / "04_shared52_proline_reaction_fluxes_long.csv", index=False)
    metadata = {
        "output_dir": str(OUT),
        "scenario": SCENARIO,
        "contexts": len(CONTEXTS),
        "source": str(RAW),
        "note": "Existing Script 96 full taxon-by-reaction MICOM exports were audited. Intracellular synthesis is reported only from explicit reaction fluxes; no exchange-only inference is used.",
    }
    (OUT / "summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({
        "output_dir": str(OUT),
        "model_audit_rows": len(audits),
        "reaction_audit_rows": len(reaction_audits),
        "community_summary_rows": len(summary_rows),
        "community_reaction_rows": len(flux_rows),
    }, indent=2))


if __name__ == "__main__":
    main()
