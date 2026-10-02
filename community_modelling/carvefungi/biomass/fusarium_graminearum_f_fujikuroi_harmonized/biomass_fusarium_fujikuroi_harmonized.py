#!/usr/bin/env python3
"""
Harmonized Fusarium fujikuroi-derived biomass candidate for Fusarium graminearum
===============================================================================

Purpose
-------
This file defines a candidate biomass reaction for the CarveFungi-derived
Fusarium graminearum model. It is intended for comparative analysis and does
NOT delete or overwrite the existing 52-component biomass.

Recommended reaction ID
-----------------------
BIOMASS_FUSARIUM_FUJIKUROI_HARMONIZED

Reference donor
---------------
Fusarium fujikuroi model iCY1235:
Li YW et al. (2023), Journal of Agricultural and Food Chemistry
71(48):18890-18897.
DOI: 10.1021/acs.jafc.3c05309

Why this donor is biologically useful
-------------------------------------
Fusarium fujikuroi is substantially closer to Fusarium graminearum than the
previous Aspergillus fumigatus-derived biomass baseline.

    Fusarium graminearum:
        Kingdom Fungi
        Phylum Ascomycota
        Class Sordariomycetes
        Order Hypocreales
        Family Nectriaceae
        Genus Fusarium

    Fusarium fujikuroi:
        Kingdom Fungi
        Phylum Ascomycota
        Class Sordariomycetes
        Order Hypocreales
        Family Nectriaceae
        Genus Fusarium

The donor and target therefore share genus, family, order and class. By
contrast, Aspergillus fumigatus is in Eurotiomycetes / Eurotiales /
Aspergillaceae. iCY1235 is therefore a much closer phylogenetic proxy for
F. graminearum biomass composition.

Harmonization performed
-----------------------
The published iCY1235 biomass was translated to the metabolite resolution
present in the F. graminearum CarveFungi model.

1. Directly retained:
   - dATP, dCTP, dGTP, dTTP
   - ATP, CTP, GTP, UTP
   - all 20 amino acids
   - trehalose, UDP-N-acetylglucosamine, UDP-galactose
   - 10-formyl-THF, SAM, biotin, CoA, 5,10-methylene-THF
   - NAD, NADH, NADP, NADPH
   - NH4+, heme b, putrescine
   - riboflavin, TPP, PLP, THF, Fe2+
   - growth-associated ATP hydrolysis and polymerization by-products

2. Lipid-resolution harmonization:
   The donor uses explicit PA and PE molecular species. The CarveFungi model
   uses generic phosphatidate and phosphatidylethanolamine pools.

       PA14:0 + PA16:0 + PA16:1
       0.000068 + 0.000778 + 0.000200 = 0.001046
       -> M_paSC_c

       PE14:0 + PE16:0 + PE16:1
       0.000349 + 0.003997 + 0.001026 = 0.005372
       -> M_pe_c

   This preserves the summed donor molar demand while matching the resolution
   of the target reconstruction.

3. Omitted because no defensible target-model representation exists:
   - 5-methyl-THF
   - cytosolic FAD demand
   - glycogen
   - spermidine
   - sulfate
   - Ca2+, Cl-, Co2+, Cu2+, Fe3+, K+, Mg2+, Mn2+, Zn2+

The omitted ions/cofactors are MODEL-RESOLUTION omissions, not biological
claims that F. graminearum does not require them. Unsupported reactions or
artificial uptake were deliberately not added merely to make biomass feasible.

Sulfate was initially retained from the donor biomass because cytosolic sulfate
(`M_so4_c` / `so4_c`) exists in the target model. A later feasibility check
showed that the corrected F. graminearum CarveFungi model lacks a supported
extracellular sulfate exchange/interface (`so4_e` or equivalent), so the
retained sulfate term blocked the otherwise feasible biomass. It was therefore
omitted rather than adding a blind artificial sulfate source.

Glycogen and spermidine were also initially retained from the donor biomass.
Under the final shared barley/community medium, they are not producible in the
corrected F. graminearum CarveFungi model. They are therefore omitted in this
final shared-medium-compatible biomass candidate rather than adding unsupported
medium supplements or artificial sources.

Normalization note
------------------
The retained coefficients preserve the published iCY1235 molar coefficients.
The reaction is NOT mass-renormalized after omission of unsupported trace
components. It is appropriate for feasibility/sensitivity analysis. If strict
absolute gDW interpretation is required, perform a later mass-normalization
after the final retained composition is fixed.

Important source detail
-----------------------
The coefficients below follow the actual iCY1235 model reaction table rather
than relying only on the supplementary PDF table. This matters because the
model table contains small differences from the PDF in some lipid/DNA entries.

Existing comparison biomass
----------------------------
Keep the current biomass as a baseline:
R_BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF

Usage
-----
From another Python script:

    import cobra
    from biomass_fusarium_fujikuroi_harmonized import add_biomass

    model = cobra.io.read_sbml_model("your_model.xml")
    add_biomass(model, set_objective=True)
    print(model.optimize())

CLI:
    python biomass_fusarium_fujikuroi_harmonized.py input.xml output.xml --set-objective
"""

BIOMASS_ID = "BIOMASS_FUSARIUM_FUJIKUROI_HARMONIZED"
BIOMASS_NAME = "Fusarium fujikuroi-derived harmonized biomass for Fusarium graminearum"

REFERENCE = {
    "donor_species": "Fusarium fujikuroi",
    "donor_model": "iCY1235",
    "citation": "Li et al. (2023) J Agric Food Chem 71(48):18890-18897",
    "doi": "10.1021/acs.jafc.3c05309",
    "target_species": "Fusarium graminearum",
    "baseline_reaction": "BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF",
}

LIPID_AGGREGATION = {
    "M_paSC_c": {
        "donor_components": {
            "pa140": 0.000068,
            "pa160": 0.000778,
            "pa161": 0.000200,
        },
        "summed_coefficient": 0.001046,
        "reason": "CarveFungi represents a generic phosphatidate pool."
    },
    "M_pe_c": {
        "donor_components": {
            "pe140": 0.000349,
            "pe160": 0.003997,
            "pe161": 0.001026,
        },
        "summed_coefficient": 0.005372,
        "reason": "CarveFungi represents a generic phosphatidylethanolamine pool."
    },
}

OMITTED_COMPONENTS = {
    "5-methyl-THF": {
        "donor_coefficient": 0.004431,
        "reason": "No defensible equivalent metabolite represented in the current target model."
    },
    "FAD_c": {
        "donor_coefficient": 0.002592,
        "reason": "Target model does not represent a supported cytosolic FAD biomass pool."
    },
    "sulfate": {
        "donor_coefficient": 0.004126,
        "reason": (
            "M_so4_c/so4_c exists, but the corrected Fusarium CarveFungi model has no "
            "supported extracellular sulfate exchange/interface; omitted rather than "
            "adding a blind artificial sulfate source."
        ),
    },
    "glycogen": {
        "donor_coefficient": 0.154187,
        "reason": (
            "Retained donor glycogen term is not producible under the final shared "
            "barley/community medium in the corrected Fusarium CarveFungi model; "
            "omitted rather than expanding the medium or adding an artificial source."
        ),
    },
    "spermidine": {
        "donor_coefficient": 0.013699,
        "reason": (
            "Retained donor spermidine term is not producible under the final shared "
            "barley/community medium in the corrected Fusarium CarveFungi model; "
            "omitted rather than expanding the medium or adding an artificial source."
        ),
    },
    "Ca2+": {"donor_coefficient": 0.004952, "reason": "Corresponding cytosolic biomass species absent."},
    "Cl-": {"donor_coefficient": 0.004952, "reason": "Corresponding cytosolic biomass species absent."},
    "Co2+": {"donor_coefficient": 0.000024, "reason": "Corresponding cytosolic biomass species absent."},
    "Cu2+": {"donor_coefficient": 0.000674, "reason": "Corresponding cytosolic biomass species absent."},
    "Fe3+": {"donor_coefficient": 0.007428, "reason": "Corresponding cytosolic biomass species absent."},
    "K+": {"donor_coefficient": 0.185690, "reason": "Corresponding cytosolic biomass species absent."},
    "Mg2+": {"donor_coefficient": 0.008253, "reason": "Corresponding cytosolic biomass species absent."},
    "Mn2+": {"donor_coefficient": 0.000658, "reason": "Corresponding cytosolic biomass species absent."},
    "Zn2+": {"donor_coefficient": 0.000324, "reason": "Corresponding cytosolic biomass species absent."},
}

# Signed stoichiometry:
# negative = biomass reactant / consumed
# positive = product
STOICHIOMETRY = {
    # DNA
    "M_datp_c": -0.016297,
    "M_dctp_c": -0.042326,
    "M_dgtp_c": -0.042326,
    "M_dttp_c": -0.016297,

    # RNA + growth-associated ATP demand
    # iCY1235 global ATP coefficient = 75.7 GAM + 0.11045 RNA ATP
    "M_atp_c": -75.81045,
    "M_ctp_c": -0.13226,
    "M_gtp_c": -0.10035,
    "M_utp_c": -0.17554,

    # Protein
    "M_alaL_c": -0.60978,
    "M_argL_c": -0.10195,
    "M_asnL_c": -0.18343,
    "M_aspL_c": -0.18599,
    "M_cysL_c": -0.08069,
    "M_glnL_c": -0.14180,
    "M_gluL_c": -0.25085,
    "M_gly_c": -0.55317,
    "M_hisL_c": -0.07991,
    "M_ileL_c": -0.24831,
    "M_leuL_c": -0.13981,
    "M_lysL_c": -0.07080,
    "M_metL_c": -0.12059,
    "M_pheL_c": -0.13521,
    "M_proL_c": -0.37332,
    "M_serL_c": -0.21339,
    "M_thrL_c": -0.26691,
    "M_trpL_c": -0.03275,
    "M_tyrL_c": -0.06539,
    "M_valL_c": -0.22152,

    # Carbohydrate / cell-associated pools
    "M_tre_c": -0.087643,
    "M_uacgam_c": -0.083227,
    "M_udpgal_c": -0.166453,

    # Harmonized generic lipid pools
    "M_paSC_c": -0.001046,
    "M_pe_c": -0.005372,

    # Cofactors and small molecules
    "M_10fthf_c": -0.004309,
    "M_amet_c": -0.005085,
    "M_btn_c": -0.008349,
    "M_coa_c": -0.002660,
    "M_mlthf_c": -0.004460,
    "M_nad_c": -0.003066,
    "M_nadh_c": -0.003062,
    "M_nadp_c": -0.002744,
    "M_nadph_c": -0.002740,
    "M_nh4_c": -0.012379,
    "M_pheme_c": -0.003306,
    "M_ptrc_c": -0.022527,
    "M_ribflv_c": -0.005397,
    "M_thmpp_c": -0.004810,
    "M_pydx5p_c": -0.008287,
    "M_thf_c": -0.004581,
    "M_fe2_c": -0.006388,

    # Water
    "M_h2o_c": -71.62443,

    # Products
    "M_adp_c": 75.700000,
    "M_h_c": 76.585527,
    "M_pi_c": 75.700000,
    "M_ppi_c": 0.635847,
    "M_udp_c": 0.249680,
}


def validate_metabolites(model):
    """Return missing metabolite IDs required by this candidate biomass."""
    return [mid for mid in STOICHIOMETRY if mid not in model.metabolites]


def add_biomass(model, set_objective=False, replace_existing=False):
    """
    Add the harmonized biomass reaction to a COBRApy model.

    Parameters
    ----------
    model : cobra.Model
        Loaded COBRApy model.
    set_objective : bool
        If True, make this candidate the model objective.
    replace_existing : bool
        If True and a reaction with BIOMASS_ID already exists, remove it first.

    Returns
    -------
    cobra.Reaction
        Added biomass reaction.
    """
    import cobra

    if BIOMASS_ID in model.reactions:
        if not replace_existing:
            raise ValueError(
                f"{BIOMASS_ID} already exists. Use replace_existing=True to overwrite it."
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
        model.metabolites.get_by_id(mid): coeff
        for mid, coeff in STOICHIOMETRY.items()
    })

    rxn.notes = {
        "donor_species": REFERENCE["donor_species"],
        "donor_model": REFERENCE["donor_model"],
        "citation": REFERENCE["citation"],
        "doi": REFERENCE["doi"],
        "target_species": REFERENCE["target_species"],
        "harmonization_summary":
            "iCY1235 biomass adapted to Fusarium graminearum CarveFungi resolution; "
            "PA/PE molecular species collapsed to generic pools; unsupported trace ions, "
            "5-methyl-THF and cytosolic FAD demand omitted.",
        "biological_rationale":
            "F. fujikuroi and F. graminearum share the genus Fusarium, family Nectriaceae, "
            "order Hypocreales and class Sordariomycetes, making iCY1235 a closer biomass "
            "proxy than the previous Aspergillus-derived baseline.",
        "normalization_warning":
            "Published retained molar coefficients are preserved; reaction was not renormalized after omissions."
    }

    model.add_reactions([rxn])

    if set_objective:
        model.objective = rxn

    return rxn


def print_report():
    print(f"Reaction: {BIOMASS_ID}")
    print(f"Number of represented metabolites: {len(STOICHIOMETRY)}")
    print("\nAggregated lipid pools:")
    for mid, value in LIPID_AGGREGATION.items():
        print(f"  - {mid}: {value['summed_coefficient']}")
        print(f"    {value['reason']}")
    print("\nOmitted donor components:")
    for key, value in OMITTED_COMPONENTS.items():
        print(f"  - {key}: {value['reason']}")


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
    print(f"\nSaved: {args.output_sbml}")
