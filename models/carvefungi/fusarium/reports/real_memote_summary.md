# Real MEMOTE summary: Fusarium CarveFungi after stoichiometry correction

Date: 2026-09-23

## Model tested

Derived model:

`/home/tlu/Desktop/Jorge/VS CODE/Fusarium_graminearum_carveme/04_models/Fusarium_graminearum_carvefungi_after_stoichiometry_correction.xml`

This is a derived corrected model. The mother model was not overwritten.

## Main result

The real MEMOTE run confirms that the selected correction resolves the global stoichiometric inconsistency.

| MEMOTE test | Result | Details |
|---|---:|---|
| `test_stoichiometric_consistency` | passed | Model stoichiometry is consistent |
| `test_unconserved_metabolites` | passed | 0 unconserved metabolites |
| `test_inconsistent_min_stoichiometry` | passed | No reported minimal unconservable sets |

## Full MEMOTE run status

The full MEMOTE suite still reports many expected non-stoichiometric failures, mainly annotation, SBO, biomass pseudo-reaction, and per-reaction mass/charge-balance categories.

Terminal summary:

```text
178 failed, 110 passed, 12 skipped
```

Important remaining failed categories include:

- reaction/metabolite/gene annotation coverage;
- SBO annotations;
- per-reaction mass balance: 48 reactions;
- per-reaction charge balance: 76 reactions;
- biomass-related MEMOTE checks for several pseudo/biomass-like reactions.

These do not contradict the main curation target: global stoichiometric consistency now passes and MEMOTE reports zero unconserved metabolites.

## Output files

- HTML report: `real_memote_fusarium_after_stoichiometry_correction.html`
- JSON results: `real_memote_fusarium_after_stoichiometry_correction_results.json.gz`
- MEMOTE log: `real_memote_fusarium_after_stoichiometry_correction.log`
- Runner script: `run_real_memote_fusarium_after_stoichiometry_correction.sh`
