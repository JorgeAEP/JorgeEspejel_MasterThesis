# Script 74 after-stoichiometry-correction model creation

This script writes a derived curated model for real MEMOTE testing. The mother model is not modified.

## Derived files

- SBML: `/home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme/04_models/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.xml`
- SBML SHA256: `09e2caeb2a7feeb5cf491b32942f02ca183e47027daf109b17934e6a42066fc1`
- JSON: `/home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme/04_models/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.json`
- JSON SHA256: `d12ff551e0c6559901684b570d596d7639e42f4c0274247b9f29d182d990490e`

## Internal validation before real MEMOTE

```csv
condition,memote_pass,unconserved_count,unconserved_metabolites,growth,growth_ratio_vs_baseline,carbon_free_growth,reaction_count,metabolite_count,gene_count
current_viable_harmonized_baseline,False,7,h_c;h_e;h_m;h_n;h_r;h_x;hi_m,11.146350060683758,0.9999999999999963,0.0,2057,1564,551
after_stoichiometry_correction_model,True,0,,11.146350061412337,1.000000000065361,0.0,2057,1564,551
```

## Pass gate

```json
{
  "saved_derived_xml_exists": true,
  "saved_derived_json_exists": true,
  "internal_memote_unconserved_zero": true,
  "growth_retained_99_999_percent": true,
  "carbon_free_growth_zero": true,
  "all_biomass_precursors_producible": true,
  "biomass_reaction_unchanged": true,
  "objective_unchanged": true,
  "mother_model_unchanged": true
}
```

## Real MEMOTE command

```bash
# Full MEMOTE test collection as JSON:
02_environment/.venv-carveme/bin/memote run \
  --ignore-git \
  --solver gurobi \
  --filename /home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme/06_reports/74_after_stoichiometry_correction_model_creation/real_memote_after_stoichiometry_correction_results.json.gz \
  /home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme/04_models/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.xml

# Full MEMOTE HTML snapshot report:
02_environment/.venv-carveme/bin/memote report snapshot \
  --solver gurobi \
  --filename /home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme/06_reports/74_after_stoichiometry_correction_model_creation/real_memote_after_stoichiometry_correction.html \
  /home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme/04_models/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.xml
```

## Integrity

- mother checksum before: `55adadc17bdf116a62821880a6342cba5804742866280aa3cdb8ce8b700b5541`
- mother checksum after: `55adadc17bdf116a62821880a6342cba5804742866280aa3cdb8ce8b700b5541`
- mother model unchanged: `True`
