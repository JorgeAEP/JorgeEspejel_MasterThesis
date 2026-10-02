#!/usr/bin/env python3
"""Run the dual-fungi biomass-variant MICOM panel.

This panel is designed for the next comparison layer after the single-fungus
CarveFungi work:

* Bipolaris + Fusarium together
* HvSC1 or HvSC2 added separately, not together
* current shared 52-component biomass vs species-specific biomass
* final shared barley/community medium with proline relaxed to 1.0
* proline-limited control with EX_proL_m reset to 0.1

Script 94 must be run first because it creates the derived model variants,
media files and manifests used here.
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

REPORT_ROOT = PACKAGE / "reports/95_dual_fungi_biomass_variant_micom_panel"

DEFAULT_MANIFESTS = [
    # Controls.
    "hvsc1_bacteria_only",
    "hvsc2_bacteria_only",
    "fungus_alone__bipolaris_carvefungi_52",
    "fungus_alone__bipolaris_carvefungi_specific",
    "fungus_alone__fusarium_carvefungi_52",
    "fungus_alone__fusarium_carvefungi_specific",
    "fungi_pair__carvefungi_52",
    "fungi_pair__carvefungi_specific",
    # Single fungus + bacterial SynCom comparisons.
    "hvsc1_plus_bipolaris_carvefungi_52",
    "hvsc2_plus_bipolaris_carvefungi_52",
    "hvsc1_plus_bipolaris_carvefungi_specific",
    "hvsc2_plus_bipolaris_carvefungi_specific",
    "hvsc1_plus_fusarium_carvefungi_52",
    "hvsc2_plus_fusarium_carvefungi_52",
    "hvsc1_plus_fusarium_carvefungi_specific",
    "hvsc2_plus_fusarium_carvefungi_specific",
    # Main dual-fungi + bacterial SynCom comparisons.
    "hvsc1_plus_fungi_pair__carvefungi_52",
    "hvsc2_plus_fungi_pair__carvefungi_52",
    "hvsc1_plus_fungi_pair__carvefungi_specific",
    "hvsc2_plus_fungi_pair__carvefungi_specific",
]

MAIN_COMPARISON_MANIFESTS = [
    "fungi_pair__carvefungi_52",
    "fungi_pair__carvefungi_specific",
    "hvsc1_plus_bipolaris_carvefungi_52",
    "hvsc2_plus_bipolaris_carvefungi_52",
    "hvsc1_plus_bipolaris_carvefungi_specific",
    "hvsc2_plus_bipolaris_carvefungi_specific",
    "hvsc1_plus_fusarium_carvefungi_52",
    "hvsc2_plus_fusarium_carvefungi_52",
    "hvsc1_plus_fusarium_carvefungi_specific",
    "hvsc2_plus_fusarium_carvefungi_specific",
    "hvsc1_plus_fungi_pair__carvefungi_52",
    "hvsc2_plus_fungi_pair__carvefungi_52",
    "hvsc1_plus_fungi_pair__carvefungi_specific",
    "hvsc2_plus_fungi_pair__carvefungi_specific",
]

SCENARIOS = {
    "final_shared_proline_1": PACKAGE / "media/biomass_variant_panel/final_shared_proline_1.csv",
    "final_shared_proline_0p1": PACKAGE / "media/biomass_variant_panel/final_shared_proline_0p1.csv",
}


def output_summary_path(manifest: str, scenario: str) -> Path:
    return (
        REPORT_ROOT
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
        "--main-only",
        action="store_true",
        help="Run only fungi-pair and fungi-pair-plus-SynCom manifests, skipping fungus-alone and bacteria-only controls.",
    )
    parser.add_argument("--manifest", action="append", choices=DEFAULT_MANIFESTS, help="Restrict to one or more manifests.")
    parser.add_argument("--scenario", action="append", choices=list(SCENARIOS), help="Restrict to one or more scenarios.")
    args = parser.parse_args()

    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    if not PYTHON.exists():
        raise FileNotFoundError(PYTHON)
    if not RUNNER.exists():
        raise FileNotFoundError(RUNNER)

    if args.manifest:
        manifests = args.manifest
    elif args.main_only:
        manifests = MAIN_COMPARISON_MANIFESTS
    else:
        manifests = DEFAULT_MANIFESTS

    scenarios = {name: SCENARIOS[name] for name in (args.scenario or SCENARIOS)}

    for manifest in manifests:
        manifest_path = PACKAGE / "manifests" / f"{manifest}.csv"
        if not manifest_path.exists():
            raise FileNotFoundError(f"{manifest_path}. Run Script 94 first.")
    for scenario, path in scenarios.items():
        if not path.exists():
            raise FileNotFoundError(f"{path}. Run Script 94 first.")

    print("Script 95 dual-fungi biomass-variant MICOM panel")
    print(f"started: {datetime.now().isoformat(timespec='seconds')}")
    print(f"package: {PACKAGE}")
    print(f"report_root: {REPORT_ROOT}")
    print(f"manifests: {';'.join(manifests)}")
    print(f"scenarios: {';'.join(scenarios)}")

    rows = []
    total = len(manifests) * len(scenarios)
    count = 0
    for scenario, medium_csv in scenarios.items():
        for manifest in manifests:
            count += 1
            print(f"[{count}/{total}] {manifest} :: {scenario}", flush=True)
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
        "solver": args.solver,
        "fraction": args.fraction,
        "manifests": manifests,
        "scenarios": list(scenarios),
        "job_count": int(len(rows)),
        "failed_jobs": int((status["returncode"] != 0).sum()),
        "all_jobs_ok": bool((status["returncode"] == 0).all()),
        "finished": datetime.now().isoformat(timespec="seconds"),
    }
    (REPORT_ROOT / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 95 compact summary")
    print(f"job_count: {summary['job_count']}")
    print(f"failed_jobs: {summary['failed_jobs']}")
    print(f"all_jobs_ok: {summary['all_jobs_ok']}")
    print(f"outputs: {REPORT_ROOT}")

    if not summary["all_jobs_ok"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
