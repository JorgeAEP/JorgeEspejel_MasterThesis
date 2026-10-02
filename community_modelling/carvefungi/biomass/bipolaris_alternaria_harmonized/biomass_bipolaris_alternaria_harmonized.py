#!/usr/bin/env python3
"""
Harmonized Alternaria-derived biomass candidate for Bipolaris sorokiniana
========================================================================

Purpose
-------
This file defines a candidate biomass reaction for the CarveFungi-derived
Bipolaris sorokiniana model. It is intended for comparative analysis and does
NOT delete or overwrite the existing biomass reaction.

Recommended reaction ID
-----------------------
BIOMASS_BIPOLARIS_ALTERNARIA_HARMONIZED

Reference donor
---------------
Alternaria sp. MG1 model iYL1539:
Lu Y. et al. (2019), Microbial Cell Factories 18:13.
DOI: 10.1186/s12934-019-1063-7

Why this donor is biologically useful
-------------------------------------
Bipolaris and Alternaria are much closer taxonomically than the previously
used Aspergillus fumigatus-derived biomass baseline.

    Bipolaris sorokiniana:
        Kingdom Fungi
        Phylum Ascomycota
        Class Dothideomycetes
        Order Pleosporales
        Family Pleosporaceae

    Alternaria sp.:
        Kingdom Fungi
        Phylum Ascomycota
        Class Dothideomycetes
        Order Pleosporales
        Family Pleosporaceae

Thus, the donor and target share class, order, and family. In contrast,
Aspergillus fumigatus belongs to Eurotiomycetes / Eurotiales / Aspergillaceae.
This does not make the Alternaria biomass experimentally species-specific for
B. sorokiniana, but it makes it a more phylogenetically relevant proxy.

Harmonization performed
-----------------------
The published Alternaria biomass was translated to metabolite IDs present in
the B. sorokiniana CarveFungi model.

1. Directly retained:
   - all 20 amino acids
   - AMP, GMP, CMP, UMP
   - dAMP, dGMP, dTMP
   - triglyceride
   - ergosterol
   - phosphatidylethanolamine
   - phosphatidylserine
   - phosphatidylinositol
   - phosphatidylcholine
   - chitin
   - beta-1,3-glucan
   - ATP hydrolysis / growth-associated energetic cost

2. Resolution/compartment adaptations:
   - donor phosphatidate -> generic CarveFungi phosphatidate M_paSC_c
   - donor phosphatidylglycerol -> mitochondrial M_pg_m
   - donor cardiolipin -> mitochondrial M_clpn_m

   The mitochondrial placement of PG/cardiolipin is preferable to inventing
   unsupported cytosolic production because cardiolipin is characteristically
   associated with mitochondrial membranes in eukaryotes.

3. Omitted:
   - dCMP: metabolite is absent from the current Bipolaris CarveFungi model.
   - mannitol: M_mnl_c exists, but was not producible in the tested model state.

These omissions are MODEL-RESOLUTION decisions. They must not be interpreted
as evidence that B. sorokiniana biologically lacks dCMP or mannitol.

Normalization note
------------------
The retained coefficients preserve the published donor molar coefficients.
They are NOT renormalized after omission of dCMP and mannitol. Therefore this
candidate is best used first for feasibility, sensitivity analysis, and
relative comparison. If absolute growth units are to be interpreted strictly
as gDW, a later biomass mass-normalization step should be performed using the
final retained composition.

Existing comparison biomass
----------------------------
The current model biomass should be kept as a baseline:
R_BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF

Usage
-----
From another Python script:

    import cobra
    from biomass_bipolaris_alternaria_harmonized import add_biomass

    model = cobra.io.read_sbml_model("your_model.xml")
    add_biomass(model, set_objective=True)
    print(model.optimize())

CLI:
    python biomass_bipolaris_alternaria_harmonized.py input.xml output.xml --set-objective
"""

BIOMASS_ID = "BIOMASS_BIPOLARIS_ALTERNARIA_HARMONIZED"
BIOMASS_NAME = "Alternaria-derived harmonized biomass for Bipolaris sorokiniana"

REFERENCE = {
    "donor_species": "Alternaria sp. MG1",
    "donor_model": "iYL1539",
    "citation": "Lu et al. (2019) Microbial Cell Factories 18:13",
    "doi": "10.1186/s12934-019-1063-7",
    "target_species": "Bipolaris sorokiniana",
    "baseline_reaction": "BIOMASS_BIPOLARIS_DIRECT_HARMONIZED_CF",
}

ADAPTATIONS = {
    "phosphatidate": {
        "donor_term": "pa[c]",
        "target_metabolite": "M_paSC_c",
        "reason": "Generic phosphatidate pool available in CarveFungi."
    },
    "phosphatidylglycerol": {
        "donor_term": "pg[c]",
        "target_metabolite": "M_pg_m",
        "reason": "Only supported PG representation is mitochondrial."
    },
    "cardiolipin": {
        "donor_term": "cl[c]",
        "target_metabolite": "M_clpn_m",
        "reason": "Mitochondrial cardiolipin is represented and biologically preferable."
    },
}

OMITTED_COMPONENTS = {
    "dCMP": {
        "donor_coefficient": 0.0048,
        "reason": "M_dcmp_c is absent from the current model."
    },
    "D-mannitol": {
        "donor_coefficient": 0.8304,
        "reason": "M_mnl_c exists but was not producible in the tested model state."
    },
}

# Signed stoichiometry:
# negative = biomass reactant / consumed
# positive = product
STOICHIOMETRY = {
    # Energy
    "M_h2o_c": -71.4,
    "M_atp_c": -71.4,

    # Amino acids
    "M_alaL_c": -0.0985,
    "M_argL_c": -0.0423,
    "M_asnL_c": -0.1842,
    "M_aspL_c": -0.0940,
    "M_cysL_c": -0.0599,
    "M_gluL_c": -0.1118,
    "M_glnL_c": -0.1527,
    "M_gly_c": -0.0846,
    "M_hisL_c": -0.0293,
    "M_ileL_c": -0.0454,
    "M_leuL_c": -0.0738,
    "M_lysL_c": -0.0683,
    "M_metL_c": -0.1190,
    "M_pheL_c": -0.0372,
    "M_proL_c": -0.0385,
    "M_serL_c": -0.0672,
    "M_thrL_c": -0.0627,
    "M_trpL_c": -0.0705,
    "M_tyrL_c": -0.0334,
    "M_valL_c": -0.0465,

    # DNA/RNA
    "M_damp_c": -0.0046,
    "M_dgmp_c": -0.0048,
    # dCMP intentionally omitted
    "M_dtmp_c": -0.0046,
    "M_amp_c": -0.0234,
    "M_gmp_c": -0.0237,
    "M_cmp_c": -0.0222,
    "M_ump_c": -0.0225,

    # Lipids
    "M_triglyc_c": -0.117159,
    "M_ergst_c": -0.020289,
    "M_pe_c": -0.000506,
    "M_ps_c": -0.002524,
    "M_ptd1ino_c": -0.002333,
    "M_pg_m": -0.009747,
    "M_pc_c": -0.010624,
    "M_paSC_c": -0.000804,
    "M_clpn_m": -0.001676,

    # Cell wall
    "M_chitin_c": -0.7195,
    "M_13BDglcn_c": -1.0734,

    # Mannitol intentionally omitted

    # Energy products
    "M_h_c": 71.4,
    "M_pi_c": 71.4,
    "M_adp_c": 71.4,
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
            "Alternaria biomass adapted to Bipolaris CarveFungi metabolite resolution; "
            "mitochondrial PG/cardiolipin used; dCMP and non-producible mannitol omitted.",
        "biological_rationale":
            "Alternaria and Bipolaris are both Pleosporaceae/Pleosporales/Dothideomycetes, "
            "making this donor phylogenetically closer than the previous Aspergillus-derived baseline.",
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
    print("\nOmitted donor components:")
    for key, value in OMITTED_COMPONENTS.items():
        print(f"  - {key}: {value['reason']}")
    print("\nAdapted donor components:")
    for key, value in ADAPTATIONS.items():
        print(f"  - {key}: {value['donor_term']} -> {value['target_metabolite']}")
        print(f"    {value['reason']}")


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
