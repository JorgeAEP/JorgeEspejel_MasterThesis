#!/usr/bin/env python3
"""
Original shared 52-component biomass baseline
=============================================

Purpose
-------
This file preserves the EXACT 52-component direct biomass composition currently
encoded in both of the uploaded fungal models:

- Bipolaris sorokiniana:
  R_BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF

- Fusarium graminearum:
  R_BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF

The two reactions were checked directly in the SBML files and have identical
reactant identities and identical coefficients.

Recommended reaction ID
-----------------------
BIOMASS_ORIGINAL_52_COMPONENT_BASELINE

What this biomass represents
----------------------------
This is a direct biomass drain composed of 52 consumed metabolites. It has no
explicit biomass pseudometabolite product and no explicit ATP -> ADP + Pi growth-
associated energy term in the reaction itself.

This script is intended to preserve the baseline in a portable and auditable
form so that it can be compared against:

1. BIOMASS_BIPOLARIS_ALTERNARIA_HARMONIZED
2. BIOMASS_FUSARIUM_FUJIKUROI_HARMONIZED

Important provenance note
-------------------------
This script records the exact stoichiometry present in the current models.
It does NOT assign a literature origin to the 52-component formulation.
The literature provenance should be documented separately from the model
stoichiometry once the exact source is confirmed.

Why keep this baseline
----------------------
Even if the new taxonomically closer biomass candidates are ultimately preferred,
the original 52-component objective should be retained for sensitivity analysis.
It provides a common reference because both fungal models currently use the same
52-component material demand.

This makes it possible to ask:

- Does changing biomass composition alter predicted growth?
- Does a taxonomically closer biomass introduce new pathway bottlenecks?
- Which conclusions are robust to biomass formulation?
- Are differences between Bipolaris and Fusarium driven by the networks, rather
  than by using different objective functions?

Component count
---------------
52 consumed metabolites
0 reaction products
52 total metabolites represented in the reaction

COBRApy / SBML ID note
----------------------
The SBML file stores metabolite IDs with the prefix "M_" and reaction IDs with
"R_". COBRApy commonly strips these SBML prefixes on import. The helper below
therefore accepts either representation automatically.

Usage
-----
From another script:

    import cobra
    from biomass_original_52_component_baseline import add_biomass

    model = cobra.io.read_sbml_model("your_model.xml")
    add_biomass(model, set_objective=True)
    print(model.optimize())

CLI:
    python biomass_original_52_component_baseline.py input.xml output.xml --set-objective
"""

BIOMASS_ID = "BIOMASS_ORIGINAL_52_COMPONENT_BASELINE"
BIOMASS_NAME = "Original shared 52-component fungal biomass baseline"

SOURCE_REACTIONS = {
    "Bipolaris_sorokiniana": "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF",
    "Fusarium_graminearum": "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF",
}

# Signed stoichiometry:
# negative = consumed by biomass
# positive = produced
#
# These coefficients were extracted directly from the current uploaded SBML files.
STOICHIOMETRY = {
    "M_12dgr_c": -0.0186362284446717,
    "M_1aglyc_c": -0.0294408196006602,
    "M_1agpc_c": -0.0064217979613458,
    "M_1gr3petha_c": -0.0102236713718235,
    "M_9Z12Z15Zlnlncg_c": -0.0007281194240842,
    "M_9Z12Zlnlc_c": -0.0145179514974092,
    "M_9Zocda_c": -0.021651412664999,
    "M_alaL_c": -1.40316755110735,
    "M_amp_c": -0.0278122418792401,
    "M_argL_c": -0.023934424463854,
    "M_asnL_c": -0.020987291324898,
    "M_aspL_c": -0.0916085176736937,
    "M_clpn_m": -0.003625691059926,
    "M_cmp_c": -0.0212935474239356,
    "M_coa_c": -0.0094588631306413,
    "M_cysL_c": -0.0164219605782342,
    "M_damp_c": -0.0046210726467883,
    "M_dgmp_c": -0.0046725257559158,
    "M_dtmp_c": -0.0043816348730713,
    "M_epist_c": -0.0002155109978423,
    "M_ergst_c": -0.008828451209328601,
    "M_glnL_c": -0.160239928635564,
    "M_gluL_c": -0.953163180812749,
    "M_gly_c": -0.07832529261919779,
    "M_gmp_c": -0.0310692730533455,
    "M_hisL_c": -0.0038333158420551,
    "M_hxda_c": -0.014955598297248,
    "M_ileL_c": -0.0087751029887381,
    "M_lanost_c": -0.0001764456612167,
    "M_leuL_c": -0.0098074680462367,
    "M_lysL_c": -0.035316501473541,
    "M_metL_c": -0.0013356222931388,
    "M_pc_c": -0.011711645223229,
    "M_pe_c": -0.011159882072025,
    "M_pheL_c": -0.0021827808278501,
    "M_proL_c": -0.0866064142634991,
    "M_ps_c": -0.0054635336115283,
    "M_ptd1ino_c": -0.0049797500548819,
    "M_ribflv_c": -0.0192898415602559,
    "M_serL_c": -0.06439454106010579,
    "M_tda_c": -0.0004438612800492,
    "M_thf_c": -0.0162988865570555,
    "M_thmpp_c": -0.0170697253608154,
    "M_thrL_c": -0.0415943771796608,
    "M_tre_c": -0.009012711226120001,
    "M_triglyc_c": -0.0213981747750414,
    "M_trpL_c": -0.0004705016777462,
    "M_tyrL_c": -0.0044744465431786,
    "M_uacgam_c": -0.442829187390135,
    "M_udpg_c": -1.43781930636124,
    "M_ump_c": -0.0284192826668287,
    "M_valL_c": -0.07188710001310469,
}


def _candidate_ids(sbml_met_id):
    """Return SBML-style and COBRApy-style candidate IDs."""
    if sbml_met_id.startswith("M_"):
        return [sbml_met_id, sbml_met_id[2:]]
    return [sbml_met_id, "M_" + sbml_met_id]


def get_metabolite(model, sbml_met_id):
    """Find a metabolite whether COBRApy kept or removed the SBML M_ prefix."""
    for candidate in _candidate_ids(sbml_met_id):
        try:
            return model.metabolites.get_by_id(candidate)
        except KeyError:
            pass
    raise KeyError(f"Required metabolite not found: {sbml_met_id}")


def validate_metabolites(model):
    """Return SBML-style IDs that cannot be resolved in the supplied model."""
    missing = []
    for mid in STOICHIOMETRY:
        try:
            get_metabolite(model, mid)
        except KeyError:
            missing.append(mid)
    return missing


def add_biomass(model, set_objective=False, replace_existing=False):
    """
    Add a copy of the original 52-component baseline to a COBRApy model.

    The model's existing biomass reaction is NOT removed.

    Parameters
    ----------
    model : cobra.Model
        Loaded COBRApy model.
    set_objective : bool
        If True, set this added baseline copy as the objective.
    replace_existing : bool
        If a reaction with BIOMASS_ID already exists, remove that copy first.

    Returns
    -------
    cobra.Reaction
        Newly added baseline reaction.
    """
    import cobra

    if BIOMASS_ID in model.reactions:
        if not replace_existing:
            raise ValueError(
                f"{BIOMASS_ID} already exists. "
                "Use replace_existing=True only if you intentionally want to replace the copy."
            )
        model.remove_reactions([model.reactions.get_by_id(BIOMASS_ID)])

    missing = validate_metabolites(model)
    if missing:
        raise KeyError("Missing required metabolites: " + ", ".join(missing))

    rxn = cobra.Reaction(BIOMASS_ID)
    rxn.name = BIOMASS_NAME
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0

    rxn.add_metabolites({
        get_metabolite(model, mid): coeff
        for mid, coeff in STOICHIOMETRY.items()
    })

    rxn.notes = {
        "role": "Baseline copy for biomass sensitivity analysis",
        "source_model_reactions":
            "Bipolaris BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF and "
            "Fusarium BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF",
        "component_count": "52 consumed metabolites",
        "energetic_note":
            "Direct biomass drain; no explicit ATP hydrolysis/GAM term in this reaction.",
        "provenance_note":
            "Exact stoichiometry copied from current models. Literature origin intentionally "
            "not asserted here until independently confirmed.",
        "comparison_use":
            "Retain as a common baseline when testing the Alternaria-derived Bipolaris "
            "and F. fujikuroi-derived Fusarium biomass candidates."
    }

    model.add_reactions([rxn])

    if set_objective:
        model.objective = rxn

    return rxn


def print_report():
    print(f"Reaction: {BIOMASS_ID}")
    print(f"Consumed components: {sum(1 for x in STOICHIOMETRY.values() if x < 0)}")
    print(f"Produced metabolites: {sum(1 for x in STOICHIOMETRY.values() if x > 0)}")
    print("Baseline is identical in the uploaded Bipolaris and Fusarium models.")


if __name__ == "__main__":
    import argparse
    import cobra

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_sbml")
    parser.add_argument("output_sbml")
    parser.add_argument("--set-objective", action="store_true")
    parser.add_argument("--replace-existing", action="store_true")
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(args.input_sbml)
    add_biomass(
        model,
        set_objective=args.set_objective,
        replace_existing=args.replace_existing,
    )
    cobra.io.write_sbml_model(model, args.output_sbml)
    print_report()
    print(f"Saved: {args.output_sbml}")
