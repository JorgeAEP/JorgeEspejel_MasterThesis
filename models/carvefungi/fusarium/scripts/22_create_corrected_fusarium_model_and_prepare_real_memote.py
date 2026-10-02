#!/usr/bin/env python3
"""Create the final corrected Fusarium CarveFungi model and real-MEMOTE runner.

This script does not modify the mother model. It builds the selected corrected
model from Script 21 in memory, writes a derived SBML/JSON model, validates the
main copy-only gates again, and prepares a quoted bash runner for real MEMOTE.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import cobra


ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "04_models/Fusarium_graminearum_carvefungi.xml"
SCRIPT_21 = ROOT / "03_scripts/21_final_validate_selected_fusarium_stoichiometric_correction_copy_only.py"
OUT = ROOT / "06_reports/22_corrected_fusarium_model_creation_and_real_memote"
DERIVED_XML = ROOT / "04_models/Fusarium_graminearum_carvefungi_after_stoichiometry_correction.xml"
DERIVED_JSON = ROOT / "04_models/Fusarium_graminearum_carvefungi_after_stoichiometry_correction.json"
MEMOTE = ROOT.parent / "Bipolaris_sorokiniana_carveme/02_environment/.venv-carveme/bin/memote"
ZERO = 1e-9


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S21 = load_module(SCRIPT_21, "script21_for_script22")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mother_before = sha256(MODEL)

    corrected, native_medium, selected_medium, _audit_rows = S21.build_selected_corrected_copy()

    native_eval = S21.evaluate_condition(corrected, "corrected_native_medium", native_medium)
    selected_eval = S21.evaluate_condition(corrected, "corrected_selected_shared_plus_37_medium", selected_medium)

    native_precursors = S21.biomass_precursor_producibility(
        corrected, native_medium, S21.HARMONIZED_BIOMASS_ID, "native_medium"
    )
    selected_precursors = S21.biomass_precursor_producibility(
        corrected, selected_medium, S21.HARMONIZED_BIOMASS_ID, "selected_shared_plus_37_medium"
    )

    cobra.io.write_sbml_model(corrected, str(DERIVED_XML))
    cobra.io.save_json_model(corrected, str(DERIVED_JSON))

    mother_after = sha256(MODEL)
    derived_xml_sha = sha256(DERIVED_XML)
    derived_json_sha = sha256(DERIVED_JSON)

    real_memote_json = OUT / "real_memote_fusarium_after_stoichiometry_correction_results.json.gz"
    real_memote_html = OUT / "real_memote_fusarium_after_stoichiometry_correction.html"
    real_memote_log = OUT / "real_memote_fusarium_after_stoichiometry_correction.log"
    runner = OUT / "run_real_memote_fusarium_after_stoichiometry_correction.sh"

    runner.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                f"MODEL={json.dumps(str(DERIVED_XML))}",
                f"MEMOTE={json.dumps(str(MEMOTE))}",
                f"JSON_OUT={json.dumps(str(real_memote_json))}",
                f"HTML_OUT={json.dumps(str(real_memote_html))}",
                f"LOG_OUT={json.dumps(str(real_memote_log))}",
                'echo "Running real MEMOTE on: $MODEL" | tee "$LOG_OUT"',
                '"$MEMOTE" run --ignore-git --solver gurobi --filename "$JSON_OUT" "$MODEL" 2>&1 | tee -a "$LOG_OUT"',
                '"$MEMOTE" report snapshot --solver gurobi --filename "$HTML_OUT" "$MODEL" 2>&1 | tee -a "$LOG_OUT"',
                'echo "MEMOTE JSON: $JSON_OUT" | tee -a "$LOG_OUT"',
                'echo "MEMOTE HTML: $HTML_OUT" | tee -a "$LOG_OUT"',
                "",
            ]
        ),
        encoding="utf-8",
    )
    runner.chmod(0o755)

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "mother_model": str(MODEL),
        "mother_sha256_before": mother_before,
        "mother_sha256_after": mother_after,
        "mother_unchanged": mother_before == mother_after,
        "derived_xml": str(DERIVED_XML),
        "derived_xml_sha256": derived_xml_sha,
        "derived_json": str(DERIVED_JSON),
        "derived_json_sha256": derived_json_sha,
        "native_eval": native_eval,
        "selected_eval": selected_eval,
        "native_biomass_precursors_producible": f"{sum(bool(r['producible']) for r in native_precursors)}/{len(native_precursors)}",
        "selected_biomass_precursors_producible": f"{sum(bool(r['producible']) for r in selected_precursors)}/{len(selected_precursors)}",
        "real_memote_runner": str(runner),
        "real_memote_json": str(real_memote_json),
        "real_memote_html": str(real_memote_html),
        "real_memote_log": str(real_memote_log),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 22 compact summary")
    print(f"derived_xml: {DERIVED_XML}")
    print(f"derived_xml_sha256: {derived_xml_sha}")
    print(f"derived_json: {DERIVED_JSON}")
    print(f"derived_json_sha256: {derived_json_sha}")
    print(f"native_memote_pass: {native_eval['memote_pass']}")
    print(f"native_unconserved_count: {native_eval['unconserved_count']}")
    print(f"native_growth: {native_eval['growth']}")
    print(f"native_carbon_free_growth: {native_eval['carbon_free_growth']}")
    print(f"selected_memote_pass: {selected_eval['memote_pass']}")
    print(f"selected_unconserved_count: {selected_eval['unconserved_count']}")
    print(f"selected_growth: {selected_eval['growth']}")
    print(f"selected_carbon_free_growth: {selected_eval['carbon_free_growth']}")
    print(f"native_biomass_precursors_producible: {summary['native_biomass_precursors_producible']}")
    print(f"selected_biomass_precursors_producible: {summary['selected_biomass_precursors_producible']}")
    print(f"mother_unchanged: {mother_before == mother_after}")
    print(f"outputs: {OUT}")
    print("real_memote_runner:")
    print(runner)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
