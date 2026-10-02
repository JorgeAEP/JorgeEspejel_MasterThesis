# Final shared barley/community medium

This folder contains the **single final shared medium** to use for the fungal and community comparisons involving:

- *Bipolaris sorokiniana* CarveFungi
- *Fusarium graminearum* CarveFungi
- the RAVEN/independent *Bipolaris* reconstruction, where applicable
- HvSC1 and HvSC2 bacterial communities

## Final medium file

```text
final_shared_barley_community_medium_dual_namespace.csv
```

Rows: 125 medium reactions.

This is the latest shared barley/community medium. It corresponds to the previous expanded medium:

```text
../shared_barley_medium_versions/04_shared_barley_medium_with_fusarium_expansion_dual_namespace.csv
```

with one final bound adjustment:

```text
EX_proL_m: 0.1 -> 1.0
```

This proline adjustment was introduced after biomass-sensitivity testing showed
that the finalized *Fusarium fujikuroi*-derived biomass was specifically
proline-limited under the shared medium. Increasing proline to 1.0 restored a
more reasonable Fusarium growth rate without increasing *Bipolaris* growth and
without adding any new medium component.

## Why this is the final shared medium

The earlier shared medium supported the *Bipolaris* comparison, but the *Fusarium graminearum* CarveFungi model with the shared 52-component biomass did not grow on that strict medium.

The final medium keeps the original barley/community medium and adds the selected 37 Fusarium-compatible supplement validated in the Fusarium curation workflow:

```text
shared-only Fusarium growth: 0.0
shared + selected 37 supplement growth: 1.1546488888888882
carbon-free growth: 0.0
biomass precursors: 52/52
```

Therefore, this final medium is the common medium to use when comparing *Bipolaris* and *Fusarium* under the same environmental scenario.

The final proline relaxation should be interpreted as a bound adjustment on an
already-present barley/community medium component, not as the addition of a new
nutrient source.

## Important interpretation

This is a **simulation medium file**, not a claim that every organism can consume every listed component.

Each model can only use a medium row if the community/model interface contains a compatible exchange route. Keeping the same final medium for all models means the external scenario is standardized; differences in growth then reflect model/network/interface differences rather than different manually supplied media.

## Provenance files

The older medium versions remain in:

```text
../shared_barley_medium_versions/
```

This folder also includes:

```text
source_medium_version_audit.csv
```

which records the historical medium files used to build the final version.
