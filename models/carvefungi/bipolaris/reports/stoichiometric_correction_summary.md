# Stoichiometric correction summary for the Bipolaris sorokiniana CarveFungi-derived GEM

This document summarizes the final stoichiometric curation path used to resolve the MEMOTE stoichiometric-consistency failure in the CarveFungi-derived *Bipolaris sorokiniana* genome-scale metabolic model.

All diagnostic work was performed on `model.copy()`. The mother model was never modified:

```text
04_models/Bipolaris_sorokiniana_carvefungi.xml
```

The mother-model SHA256 checksum remained unchanged throughout the investigation.

## Starting point

The original model failed MEMOTE stoichiometric consistency with:

```text
35 unconserved metabolites
```

The curation reduced this in three main correction stages:

```text
35 -> 23 -> 7 -> 0
```

Each stage was validated separately, with growth and carbon-free growth monitored to avoid introducing biologically implausible free biomass production.

## Step 1: UF01339 nucleotide hydrolysis correction

The first validated correction involved the UF01339 reaction family in multiple compartments:

```text
UF01339_C
UF01339_m
UF01339_x
```

The corrected reaction form was:

```text
dATP + H2O -> dAMP + H+ + PPi
```

This correction addressed an incorrect nucleotide hydrolysis stoichiometry and reduced the MEMOTE unconserved metabolite count:

```text
35 -> 23
```

This step was retained because it represented a direct stoichiometric correction to a chemically interpretable reaction family.

## Step 2: Direct biomass and Batch E representation harmonization

The second major correction was not a single biochemical reaction edit. It was a representation-level harmonization of the biomass/macromolecule/polymer layer.

The final viable background used:

```text
validated UF01339 correction
direct biomass representation
removal of the Batch E artificial biomass/macromolecule/polymer context
```

The retained background was:

```text
UF01339_plus_direct_biomass_remove_Batch_E_context_17
```

This reduced the remaining MEMOTE unconserved metabolites:

```text
23 -> 7
```

while preserving the key model-level controls:

```text
growth ≈ 11.1463500607
carbon-free growth = 0
```

The remaining seven unconserved metabolites after this step were:

```text
h_c;h_e;h_m;h_n;h_r;h_x;hi_m
```

These were interpreted as residual proton-bookkeeping and representation-convention issues rather than biomass composition failure.

## Step 3: Localized cytochrome-b5/desaturase representation correction

The final seven unconserved metabolites were localized through a series of copy-only diagnostics. Broad deletion controls, generic metabolite splitting, respiratory alternatives, and source-driven family harmonization did not resolve the final inconsistency.

The decisive localization identified a small cytochrome-b5/desaturase-associated subsystem. The final selected correction was the minimal among-tested two-reaction stoichiometric representation correction that restored MEMOTE stoichiometric consistency while preserving growth and carbon-free controls.

The selected correction modifies:

```text
UF00430_CR
UF00435_C
```

### UF00430_CR correction

Before:

```text
4fzym_r + 2 focytb5_c + h_r + o2_r
-> 2 ficytb5_c + h2o_r + zym_int1_r
```

After:

```text
4fzym_r + 2 focytb5_c + 2 h_r + o2_r
-> 2 ficytb5_c + h2o_r + zym_int1_r
```

In coefficient form:

```text
h_r: -1 -> -2
```

### UF00435_C correction

Before:

```text
2 focytb5_c + 2 h_c + o2_c + stcoa_c
-> 9Zocdcoa_c + 2 ficytb5_c + 2 h2o_c
```

After:

```text
2 focytb5_c + o2_c + stcoa_c
-> 9Zocdcoa_c + 2 ficytb5_c + h2o_c
```

In coefficient form:

```text
h_c:   -2 -> 0
h2o_c: +2 -> +1
```

This final correction reduced the MEMOTE unconserved metabolite count:

```text
7 -> 0
```

## Final validation

The selected correction was validated in Script 73 as a copy-only final gate:

```text
03_scripts/73_final_validate_selected_two_reaction_stoichiometric_correction_copy_only.py
```

Validation output:

```text
corrected MEMOTE stoichiometric consistency: pass
corrected unconserved metabolite count: 0
corrected unconserved metabolites: none
corrected growth: 11.146350061412337
growth ratio vs baseline: 1.000000000065361
corrected carbon-free growth: 0
biomass reaction unchanged: True
objective unchanged: True
biomass precursors producible: 52/52
reaction count baseline/corrected: 2057/2057
metabolite count baseline/corrected: 1564/1564
gene count baseline/corrected: 551/551
mother model unchanged: True
all pass-gate checks true: True
```

The final validation confirmed that the correction resolved the remaining stoichiometric inconsistency without altering model size, biomass reaction, objective, growth capacity, carbon-free growth behavior, or the mother model.

## Interpretation

The stoichiometric curation can be summarized as:

```text
1. Correct nucleotide hydrolysis stoichiometry in UF01339:
   35 -> 23

2. Harmonize direct biomass representation and remove the artificial Batch E layer:
   23 -> 7

3. Apply a localized cytochrome-b5/desaturase representation correction:
   7 -> 0
```

The final two-reaction correction is best described as a localized stoichiometric representation harmonization selected after systematic exclusion of broader alternatives. It restored MEMOTE stoichiometric consistency while preserving the validated model-level constraints.

## Important safeguard

No modified SBML was saved during this validation workflow. The final correction is ready to be applied to a derived curated SBML only after explicit approval. The original mother model must remain unchanged.
