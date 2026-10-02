#!/usr/bin/env python3
"""Script 86: run the core MICOM panel for RAVEN vs CarveFungi.

This is the overnight batch runner for the main comparison analyses.

Default panel:
    scenarios:
        full
        phosphate_limited
        combined_source_limited
    manifests:
        fungus_alone__raven_independent
        fungus_alone__carvefungi
        hvsc1_bacteria_only
        hvsc2_bacteria_only
        hvsc1_plus_bipolaris__raven_independent
        hvsc2_plus_bipolaris__raven_independent
        hvsc1_plus_bipolaris__carvefungi
        hvsc2_plus_bipolaris__carvefungi

This covers core analyses 2--12. Follow-up analyses 13--18 should be run after
Script 87 identifies important candidate exchanges/strains/nutrients.
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

REPORT_ROOT = PACKAGE / "reports/86_core_micom_panel"

CORE_MANIFESTS = [
    "fungus_alone__raven_independent",
    "fungus_alone__carvefungi",
    "hvsc1_bacteria_only",
    "hvsc2_bacteria_only",
    "hvsc1_plus_bipolaris__raven_independent",
    "hvsc2_plus_bipolaris__raven_independent",
    "hvsc1_plus_bipolaris__carvefungi",
    "hvsc2_plus_bipolaris__carvefungi",
]

DEFAULT_SCENARIOS = {
    "full": PACKAGE / "media/sensitivity_panel/full.csv",
    "phosphate_limited": PACKAGE / "media/sensitivity_panel/phosphate_limited.csv",
    "combined_source_limited": PACKAGE / "media/combined_source_limited_panel/combined_source_limited.csv",
}

ALL_SENSITIVITY_SCENARIOS = {
    "full": PACKAGE / "media/sensitivity_panel/full.csv",
    "minus_folate": PACKAGE / "media/sensitivity_panel/minus_folate.csv",
    "minus_riboflavin": PACKAGE / "media/sensitivity_panel/minus_riboflavin.csv",
    "minus_thiamine": PACKAGE / "media/sensitivity_panel/minus_thiamine.csv",
    "carbon_limited": PACKAGE / "media/sensitivity_panel/carbon_limited.csv",
    "nitrogen_limited": PACKAGE / "media/sensitivity_panel/nitrogen_limited.csv",
    "phosphate_limited": PACKAGE / "media/sensitivity_panel/phosphate_limited.csv",
    "sulfur_limited": PACKAGE / "media/sensitivity_panel/sulfur_limited.csv",
    "combined_source_limited": PACKAGE / "media/combined_source_limited_panel/combined_source_limited.csv",
}


def output_summary_path(manifest: str, scenario: str) -> Path:
    return (
        REPORT_ROOT
        / manifest
        / scenario
        / f"{manifest}_{scenario}_emma_style_community_growth_summary.csv"
    )


def run_job(manifest: str, scenario: str, medium_csv: Path, solver: str, fraction: float, debug: bool, force: bool) -> dict[str, object]:
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
        str(REPORT_ROOT),
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", default="cplex")
    parser.add_argument("--fraction", type=float, default=1.0)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--force", action="store_true", help="Rerun jobs even if output summaries already exist.")
    parser.add_argument(
        "--all-sensitivity",
        action="store_true",
        help="Run all sensitivity-panel scenarios instead of only full/phosphate/combined.",
    )
    parser.add_argument("--manifest", action="append", choices=CORE_MANIFESTS, help="Restrict to one or more manifests.")
    parser.add_argument("--scenario", action="append", help="Restrict to one or more scenario labels.")
    args = parser.parse_args()

    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    if not PYTHON.exists():
        raise FileNotFoundError(PYTHON)
    if not RUNNER.exists():
        raise FileNotFoundError(RUNNER)

    scenario_map = ALL_SENSITIVITY_SCENARIOS if args.all_sensitivity else DEFAULT_SCENARIOS
    if args.scenario:
        unknown = sorted(set(args.scenario).difference(scenario_map))
        if unknown:
            raise ValueError(f"Unknown scenarios: {unknown}. Available: {sorted(scenario_map)}")
        scenario_map = {name: scenario_map[name] for name in args.scenario}

    manifests = args.manifest or CORE_MANIFESTS
    for name, path in scenario_map.items():
        if not path.exists():
            raise FileNotFoundError(path)
    for manifest in manifests:
        manifest_path = PACKAGE / "manifests" / f"{manifest}.csv"
        if not manifest_path.exists():
            raise FileNotFoundError(manifest_path)

    print("Script 86 core MICOM panel")
    print(f"started: {datetime.now().isoformat(timespec='seconds')}")
    print(f"package: {PACKAGE}")
    print(f"report_root: {REPORT_ROOT}")
    print(f"manifests: {';'.join(manifests)}")
    print(f"scenarios: {';'.join(scenario_map)}")

    job_rows = []
    total = len(manifests) * len(scenario_map)
    counter = 0
    for scenario, medium_csv in scenario_map.items():
        for manifest in manifests:
            counter += 1
            print(f"[{counter}/{total}] {manifest} :: {scenario}", flush=True)
            row = run_job(
                manifest=manifest,
                scenario=scenario,
                medium_csv=medium_csv,
                solver=args.solver,
                fraction=args.fraction,
                debug=args.debug,
                force=args.force,
            )
            job_rows.append(row)
            print(f"  {row['status']} ({row['returncode']})", flush=True)
            pd.DataFrame(job_rows).to_csv(REPORT_ROOT / "job_status_live.csv", index=False)

    job_status = pd.DataFrame(job_rows)
    job_status.to_csv(REPORT_ROOT / "job_status.csv", index=False)
    summary = {
        "package": str(PACKAGE),
        "report_root": str(REPORT_ROOT),
        "solver": args.solver,
        "fraction": args.fraction,
        "manifests": manifests,
        "scenarios": list(scenario_map),
        "job_count": int(len(job_rows)),
        "failed_jobs": int((job_status["returncode"] != 0).sum()),
        "all_jobs_ok": bool((job_status["returncode"] == 0).all()),
        "finished": datetime.now().isoformat(timespec="seconds"),
    }
    (REPORT_ROOT / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 86 compact summary")
    print(f"job_count: {summary['job_count']}")
    print(f"failed_jobs: {summary['failed_jobs']}")
    print(f"all_jobs_ok: {summary['all_jobs_ok']}")
    print(f"outputs: {REPORT_ROOT}")

    if not summary["all_jobs_ok"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
