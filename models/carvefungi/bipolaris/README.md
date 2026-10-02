# *Bipolaris sorokiniana* CarveFungi model

The native model was reconstructed from the archived protein FASTA and
CarveFungi-compatible functional/localisation inputs. The construction script
is `scripts/02_run_bipolaris_carvefungi_workflow.py`.

`curated/Bipolaris_sorokiniana_carvefungi_after_stoichiometry_correction.xml`
is the final corrected derived model. The final model-creation script is
`scripts/74_create_after_stoichiometry_correction_model_for_real_memote.py`.
It depends on the archived Scripts 67 and 73, which identify and validate the
selected correction set.

The final correction included the direct biomass representation and the
selected stoichiometric adjustments to `UF00430_CR` and `UF00435_C`.
`reports/selected_reaction_before_after_audit.csv` is the concise
reaction-level audit; the validation and checksum files document the saved
model.
