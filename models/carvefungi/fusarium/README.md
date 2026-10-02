# *Fusarium graminearum* CarveFungi model

The native model was reconstructed from the archived protein FASTA and
CarveFungi-compatible functional/localisation inputs. The reconstruction
workflow is documented in `scripts/03_run_fusarium_carvefungi_workflow.py`.

`curated/Fusarium_graminearum_carvefungi_after_stoichiometry_correction.xml`
is the final corrected derived model. The model-creation script is
`scripts/22_create_corrected_fusarium_model_and_prepare_real_memote.py`; it
uses the archived Script 21 and its required Script 15 helper.

The selected final correction is documented in
`reports/selected_correction_reaction_and_medium_audit.csv`. The validation
and MEMOTE summary records are retained in `reports/`.
