# Fusarium graminearum CarveFungi stoichiometric correction summary

This report summarizes the copy-only stoichiometric curation path used for the CarveFungi-derived *Fusarium graminearum* genome-scale metabolic model.

The mother model was never modified during diagnostics:

```text
04_models/Fusarium_graminearum_carvefungi.xml
SHA256: a62db5dde09f71282e7b13c84de440f5049d5cc47b94ef3ec94af9b6566b33da
```

All tests were performed on `model.copy()`.

## Starting point

The original real MEMOTE report showed stoichiometric inconsistency:

```text
MEMOTE stoichiometric consistency: FAIL
unconserved metabolites: 32
```

The unconserved metabolites were mainly currency/bookkeeping species:

```text
H / H2O / O2 / H2O2
Pi / PPi / polyPi
hi_m
```

This pattern was very similar to the early CarveFungi-derived *Bipolaris sorokiniana* model, so the Fusarium curation reused the successful Bipolaris strategy: first correct the nucleotide hydrolysis error, then harmonize the artificial biomass/pool representation, then localize the final proton bookkeeping layer.

## Step 1: UF01339_C nucleotide hydrolysis correction

The original Fusarium `UF01339_C` reaction was:

```text
datp_c + 2 h2o_c -> damp_c + 2 h_c + 2 ppi_c
```

The corrected form tested was:

```text
datp_c + h2o_c -> damp_c + h_c + ppi_c
```

This correction reduced the unconserved metabolite count:

```text
32 -> 22
```

It also removed the Pi/PPi/polyPi part of the inconsistency while preserving growth and carbon-free controls.

## Step 2: direct 52-component biomass and Batch-E-overlap removal

The next correction followed the successful Bipolaris representation-level harmonization. The native artificial biomass/macromolecule/polymer layer was replaced by a direct 52-component biomass objective:

```text
BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF
```

The direct biomass uses the same 52-component CarveFungi-compatible biomass composition used for the curated Bipolaris CarveFungi model.

The following 15 artificial biomass/pool-context reactions were removed:

```text
BIOMASS
UF01078_C
UF01092_C
UF01093_C
UF01094_C
UF01408_C
UF01729_C
UF01730_O
UF01737_O
UF01822_C
UF02405_C
UF02406_C
UF02407_C
UF02409_C
UF02414_C
```

This reduced the remaining unconserved metabolites:

```text
22 -> 7
```

The remaining seven were:

```text
h_c; h_e; h_m; h_n; h_r; h_x; hi_m
```

The model retained growth:

```text
native Fusarium medium growth: 69.54976856798156
selected shared+37 medium growth: 1.1546488888888882
carbon-free growth: 0.0
biomass precursors producible: 52/52
```

## Step 3: localization of the final seven proton metabolites

The final seven were localized using the same approach that worked for Bipolaris:

```text
viable 7-proton background
remove focused reaction classes/subsystems
if clean, add pieces back
identify the smallest growth-compatible clean representation
```

Script 16 showed that the broad 60-reaction light H/charge imbalance layer contained the final inconsistency.

Script 18 reduced the causal class to an eight-reaction interaction between cytochrome/proton bookkeeping and glycerolipid acyltransferase bookkeeping:

```text
UF00283_CM
UF00284_CM
UF00794_c
UF00806_M
UF00874_C
UF00969_C
UF00970_C
UF01036_C
```

Removing all eight restored stoichiometric consistency but eliminated growth.

Script 19 then added subsets back from the clean background. The minimal clean and growth-restoring add-back was:

```text
UF00969_C
UF00970_C
```

Script 20 tested whether source-backed reformulations could preserve more of the eight-reaction set. The only growth-positive MEMOTE-clean conditions retained:

```text
UF00969_C
UF00970_C
```

or:

```text
UF00794_c
UF00969_C
UF00970_C
```

The source-backed variants of `UF00874_C` and `UF01036_C` still recreated the final seven-proton inconsistency. Therefore the final selected representation retains `UF00794_c`, `UF00969_C`, and `UF00970_C`, and removes the remaining five trigger/interface reactions.

## Final selected correction

The final selected correction is:

```text
1. Correct UF01339_C.

2. Remove the 15 Batch-E-overlap artificial biomass/pool reactions.

3. Add direct 52-component biomass:
   BIOMASS_FUSARIUM_DIRECT_HARMONIZED_CF

4. Remove final five trigger/interface reactions:
   UF00283_CM
   UF00284_CM
   UF00806_M
   UF00874_C
   UF01036_C

5. Retain from the eight-reaction trigger set:
   UF00794_c
   UF00969_C
   UF00970_C
```

The five final removals are best interpreted as representation-level incompatibilities in this CarveFungi stoichiometric background, not as proof that all associated biochemical activities are absent from *F. graminearum*. In particular, `UF00283_CM`, `UF00284_CM`, and `UF00806_M` resemble plausible cytochrome/proton-transfer chemistry, but their compartment/proton-carrier representation was incompatible with the selected consistent model background.

## Final copy-only validation

Script:

```text
03_scripts/21_final_validate_selected_fusarium_stoichiometric_correction_copy_only.py
```

Output folder:

```text
06_reports/21_final_selected_fusarium_stoichiometric_correction_validation
```

Final validation:

```text
MEMOTE support stoichiometric consistency: PASS
unconserved metabolites: 0

native medium growth: 69.54976856798156
native carbon-free growth: 0.0

selected shared+37 medium growth: 1.1546488888888882
selected shared+37 carbon-free growth: 0.0

biomass precursors producible on native medium: 52/52
biomass precursors producible on selected medium: 52/52

reaction count: 1535
metabolite count: 1167
gene count: 547

mother model unchanged: True
all validation gates passed: True
```

## Key output files

```text
summary.json
final_validation_condition_results.csv
selected_correction_reaction_and_medium_audit.csv
biomass_precursor_producibility.csv
biomass_precursor_producibility_summary.csv
```

## Interpretation

The selected curation resolves the formal stoichiometric inconsistency while preserving growth on both the native Fusarium medium and the selected shared+37 comparison medium. It also prevents carbon-free growth and preserves production of all 52 direct biomass precursors.

This gives a validated corrected Fusarium CarveFungi background suitable for saving as a derived model and then checking with real MEMOTE. The original mother model should remain archived unchanged.
