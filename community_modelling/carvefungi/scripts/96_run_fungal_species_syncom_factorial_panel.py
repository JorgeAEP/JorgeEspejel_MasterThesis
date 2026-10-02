#!/usr/bin/env python3
"""Run the main fungal-species × SynCom MICOM factorial panel.

Core design:

* two fungal pathogens: Bipolaris sorokiniana (Bs), Fusarium graminearum (Fg)
* two bacterial SynComs: HvSC1, HvSC2
* two fungal biomass modes: shared 52-component biomass, species-specific biomass
* main medium: original shared barley/community medium with EX_proL_m = 0.1

The proline-relaxed medium (EX_proL_m = 1.0) is treated as a targeted
Fusarium monoculture support control, not as the default ecological medium.

Run Script 94 first to generate model variants, manifests and media.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
ROOT = PACKAGE.parent

RUNNER = ROOT / "Bipolaris__sorokiniana/02_py_scripts/306_run_phaseA_kegg_medium_emma_style.py"
PYTHON = ROOT / "Bipolaris__sorokiniana/.venv-report/bin/python"

REPORT_ROOT = PACKAGE / "reports/96_fungal_species_syncom_factorial_panel"
RAW_ROOT = REPORT_ROOT / "raw_micom_runs"

MAIN_MANIFESTS = [
    # Bacterial controls.
    "hvsc1_bacteria_only",
    "hvsc2_bacteria_only",
    # Fungal monoculture controls.
    "fungus_alone__bipolaris_carvefungi_52",
    "fungus_alone__bipolaris_carvefungi_specific",
    "fungus_alone__fusarium_carvefungi_52",
    "fungus_alone__fusarium_carvefungi_specific",
    # Dual-fungus controls.
    "fungi_pair__carvefungi_52",
    "fungi_pair__carvefungi_specific",
    # Single fungus + SynCom.
    "hvsc1_plus_bipolaris_carvefungi_52",
    "hvsc2_plus_bipolaris_carvefungi_52",
    "hvsc1_plus_bipolaris_carvefungi_specific",
    "hvsc2_plus_bipolaris_carvefungi_specific",
    "hvsc1_plus_fusarium_carvefungi_52",
    "hvsc2_plus_fusarium_carvefungi_52",
    "hvsc1_plus_fusarium_carvefungi_specific",
    "hvsc2_plus_fusarium_carvefungi_specific",
    # Dual fungi + SynCom.
    "hvsc1_plus_fungi_pair__carvefungi_52",
    "hvsc2_plus_fungi_pair__carvefungi_52",
    "hvsc1_plus_fungi_pair__carvefungi_specific",
    "hvsc2_plus_fungi_pair__carvefungi_specific",
]

PROLINE_RELAXED_CONTROL_MANIFESTS = [
    "fungus_alone__fusarium_carvefungi_specific",
]

SCENARIOS = {
    "original_proline_0p1": PACKAGE / "media/biomass_variant_panel/final_shared_proline_0p1.csv",
    "proline_relaxed_1": PACKAGE / "media/biomass_variant_panel/final_shared_proline_1.csv",
}


def output_summary_path(manifest: str, scenario: str) -> Path:
    return (
        RAW_ROOT
        / manifest
        / scenario
        / f"{manifest}_{scenario}_emma_style_community_growth_summary.csv"
    )


def run_job(
    manifest: str,
    scenario: str,
    medium_csv: Path,
    solver: str,
    fraction: float,
    debug: bool,
    force: bool,
) -> dict[str, object]:
    summary_path = output_summary_path(manifest, scenario)
    if summary_path.exists() and not force:
        return {
            "manifest": manifest,
            "scenario": scenario,
            "medium_csv": str(medium_csv),
            "returncode": 0,
            "status": "skipped_existing",
            "summary_path": str(summary_path),
            "stdout_log": "",
            "stderr_log": "",
        }

    cmd = [
        str(PYTHON),
        str(RUNNER),
        "--package-root",
        str(PACKAGE),
        "--report-root",
        str(RAW_ROOT),
        "--manifest-name",
        manifest,
        "--medium-name",
        "full",
        "--medium-csv",
        str(medium_csv),
        "--medium-label",
        scenario,
        "--solver",
        solver,
        "--fraction",
        str(fraction),
    ]
    if debug:
        cmd.append("--debug")

    log_dir = REPORT_ROOT / "_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_log = log_dir / f"{manifest}__{scenario}.stdout.log"
    stderr_log = log_dir / f"{manifest}__{scenario}.stderr.log"

    with stdout_log.open("w", encoding="utf-8") as out, stderr_log.open("w", encoding="utf-8") as err:
        completed = subprocess.run(cmd, stdout=out, stderr=err, text=True)

    return {
        "manifest": manifest,
        "scenario": scenario,
        "medium_csv": str(medium_csv),
        "returncode": completed.returncode,
        "status": "ok" if completed.returncode == 0 else "failed",
        "summary_path": str(summary_path),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "cmd": " ".join(cmd),
    }


def build_job_plan(include_relaxed_control: bool, manifests: list[str] | None, scenarios: list[str] | None) -> list[tuple[str, str, Path]]:
    if manifests:
        selected = manifests
        selected_scenarios = scenarios or ["original_proline_0p1"]
        return [(manifest, scenario, SCENARIOS[scenario]) for scenario in selected_scenarios for manifest in selected]

    jobs = [(manifest, "original_proline_0p1", SCENARIOS["original_proline_0p1"]) for manifest in MAIN_MANIFESTS]
    if include_relaxed_control:
        jobs.extend(
            [
                (manifest, "proline_relaxed_1", SCENARIOS["proline_relaxed_1"])
                for manifest in PROLINE_RELAXED_CONTROL_MANIFESTS
            ]
        )
    return jobs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", default="cplex")
    parser.add_argument("--fraction", type=float, default=1.0)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--force", action="store_true", help="Rerun jobs even if output summaries already exist.")
    parser.add_argument(
        "--skip-proline-relaxed-control",
        action="store_true",
        help="Do not run the Fusarium-specific proline-relaxed monoculture control.",
    )
    parser.add_argument("--manifest", action="append", choices=MAIN_MANIFESTS + PROLINE_RELAXED_CONTROL_MANIFESTS)
    parser.add_argument("--scenario", action="append", choices=list(SCENARIOS))
    args = parser.parse_args()

    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_ROOT.mkdir(parents=True, exist_ok=True)
    if not PYTHON.exists():
        raise FileNotFoundError(PYTHON)
    if not RUNNER.exists():
        raise FileNotFoundError(RUNNER)

    jobs = build_job_plan(
        include_relaxed_control=not args.skip_proline_relaxed_control,
        manifests=args.manifest,
        scenarios=args.scenario,
    )

    for manifest, scenario, medium_csv in jobs:
        manifest_path = PACKAGE / "manifests" / f"{manifest}.csv"
        if not manifest_path.exists():
            raise FileNotFoundError(f"{manifest_path}. Run Script 94 first.")
        if not medium_csv.exists():
            raise FileNotFoundError(f"{medium_csv}. Run Script 94 first.")

    print("Script 96 fungal-species × SynCom factorial MICOM panel")
    print(f"started: {datetime.now().isoformat(timespec='seconds')}")
    print(f"package: {PACKAGE}")
    print(f"report_root: {REPORT_ROOT}")
    print(f"raw_micom_runs: {RAW_ROOT}")
    print(f"job_count_planned: {len(jobs)}")

    rows = []
    for count, (manifest, scenario, medium_csv) in enumerate(jobs, start=1):
        print(f"[{count}/{len(jobs)}] {manifest} :: {scenario}", flush=True)
        row = run_job(
            manifest=manifest,
            scenario=scenario,
            medium_csv=medium_csv,
            solver=args.solver,
            fraction=args.fraction,
            debug=args.debug,
            force=args.force,
        )
        rows.append(row)
        print(f"  {row['status']} ({row['returncode']})", flush=True)
        pd.DataFrame(rows).to_csv(REPORT_ROOT / "job_status_live.csv", index=False)

    status = pd.DataFrame(rows)
    status.to_csv(REPORT_ROOT / "job_status.csv", index=False)
    summary = {
        "package": str(PACKAGE),
        "report_root": str(REPORT_ROOT),
        "raw_micom_runs": str(RAW_ROOT),
        "solver": args.solver,
        "fraction": args.fraction,
        "job_count": int(len(rows)),
        "failed_jobs": int((status["returncode"] != 0).sum()),
        "all_jobs_ok": bool((status["returncode"] == 0).all()),
        "finished": datetime.now().isoformat(timespec="seconds"),
    }
    (REPORT_ROOT / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 96 compact summary")
    print(f"job_count: {summary['job_count']}")
    print(f"failed_jobs: {summary['failed_jobs']}")
    print(f"all_jobs_ok: {summary['all_jobs_ok']}")
    print(f"outputs: {REPORT_ROOT}")

    if not summary["all_jobs_ok"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
